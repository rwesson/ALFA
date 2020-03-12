!Copyright (C) 2013- Roger Wesson
!Free under the terms of the GNU General Public License v3

! this is the general fitting routine, which first subtracts a continuum, then fits sky lines if requested, then fits the emission lines
! it is called from within parallelised wrappers for 2D and 3D data, and simply called for 1D data.

! subtract the continuum

if (messages .and. subtractcontinuum) print *,gettime(),"fitting continuum"
call fit_continuum(realspec,spectrumlength, continuum, continuumwindow, subtractcontinuum)

! now do the fitting
! first get guesses for the redshift and resolution

if (.not. resolution_estimated) then
  ! estimate resolution assuming nyquist sampling
  resolutionguess=2*realspec(2)%wavelength/(realspec(3)%wavelength-realspec(1)%wavelength)
  if (messages) print *,gettime(),"estimated spectrograph resolution assuming Nyquist sampling: ",resolutionguess
else
  if (messages) print *,gettime(),"estimated spectrograph resolution from user input: ",resolutionguess
endif

if (rtol1 .eq. 0.d0) then
  rtol1=0.9*resolutionguess ! user didn't specify a value, default behaviour is to allow resolution to vary between 0.1 and 1.9x the initial guess on the first pass
endif

if (messages) print "(X,A,A,F8.1,A,F7.1)",gettime(),"initial guesses for velocity and resolution: ",c*(redshiftguess-1),"km/s, R=",resolutionguess

! first, subtract sky spectrum if requested. do in chunks of 400 units. no overlap necessary because velocity is zero

allocate(skyspectrum(spectrumlength))
skyspectrum%wavelength=realspec%wavelength
skyspectrum%flux=0.d0

if (subtractsky) then
  if (messages) print *,gettime(),"fitting sky emission"

  !get an array for all the sky lines in the range
  call selectlines(skylines_catalogue,realspec(1)%wavelength, realspec(size(realspec))%wavelength, skylines, nlines)
  linearraypos=1

  !if there are any sky lines to fit, then go though in chunks of 400 units
  if (nlines .gt. 0) then
    do i=1,spectrumlength,400

! overlap=nint(2*vtol2/(1-realspec(i)%wavelength/realspec(i+1)%wavelength)) ! this needs fixing, i+1 can be out of bounds
      overlap=20

! avoid refitting final section. if spectrumlength%400<overlap, this would happen

      if ((spectrumlength-i)<overlap) exit

      if (i .eq. 1) then
        startpos=1
        startwlen=realspec(1)%wavelength
      else
        startpos=i-overlap
        startwlen=realspec(i)%wavelength
      endif

      if (i+400+overlap-1 .gt. spectrumlength) then
        endpos=spectrumlength
        endwlen=realspec(spectrumlength)%wavelength
      else
        endpos=i+400+overlap-1
        endwlen=realspec(i+400)%wavelength
      endif

      allocate(spectrumchunk(endpos-startpos+1))
      spectrumchunk = realspec(startpos:endpos)

      call selectlines(skylines_catalogue, startwlen, endwlen, skylines_section, nlines)

      if (nlines .gt. 0) then
        if (messages) print "(' ',A,A,F7.1,A,F7.1,A,I3,A)",gettime(),"fitting sky emission from ",spectrumchunk(1)%wavelength," to ",spectrumchunk(size(spectrumchunk))%wavelength," with ",nlines," lines"
        call fit(spectrumchunk, 1.0, resolutionguess, skylines_section, 3.e-6, rtol2, generations, popsize, pressure) !velocity guess=0km/s, tolerance~1km/s
    !use redshift and resolution from this chunk as initial values for next chunk
!    redshiftguess=skylines_section(1)%redshift
!    resolutionguess=skylines_section(1)%resolution
    !copy results back
        skylines(linearraypos:linearraypos+nlines-1)=skylines_section
        linearraypos=linearraypos+nlines
      endif
      deallocate(spectrumchunk)
    enddo

  ! make full sky spectrum and subtract at end

    call makespectrum(skylines,skyspectrum)
    realspec%flux = realspec%flux - skyspectrum%flux

  else

    print *,gettime(),"no sky lines in wavelength range covered by spectrum"

  endif ! nlines .gt. 0
endif ! subtractsky

! select the lines from the catalogue

call selectlines(stronglines_catalogue, minval(realspec%wavelength), maxval(realspec%wavelength), fittedlines, nlines)

if (nlines .eq. 0) then
  if (messages) print *,gettime(),"Warning: no reference lines detected, using default guesses for velocity and resolution"
  redshiftguess_overall=1.d0
else
  if (messages) print *,gettime(),"estimating resolution and velocity using ",nlines," lines"
  !create an array containing just the regions around lines of interest
  !otherwise far more data is being processed than necessary
  !for each line, take 50 data points nearest to wavelength
  !should be way more than enough and still far more efficient than fitting the whole spectrum
  !todo: calculate the number of points from vtol1

  allocate(stronglines(50*nlines))
  do i=1,nlines
    linelocation=minloc(abs(stronglines_catalogue(i)%wavelength*redshiftguess-realspec%wavelength),1)
    if (linelocation-24 .gt. 0 .and. linelocation+25 .lt. size(realspec)) then
      stronglines(50*(i-1)+1:50*i) = realspec(linelocation-24:linelocation+25)
    endif
  enddo

  !now fit the strong lines

  call fit(stronglines, redshiftguess, resolutionguess, fittedlines, vtol1, rtol1, generations, popsize, pressure)

  if (messages) print *,gettime(),"estimated velocity and resolution: ",c*(fittedlines(1)%redshift-1),fittedlines(1)%resolution
  redshiftguess_overall = fittedlines(1)%redshift ! when fitting chunks, use this redshift to get lines in the right range from the catalogue. if velocity from each chunk is used, then there's a chance that a line could be missed or double counted due to variations in the calculated velocity between chunks.
  redshiftguess=fittedlines(1)%redshift
  resolutionguess=fittedlines(1)%resolution

  deallocate(stronglines)

endif

! then again in chunks with tighter tolerance

linearraypos=1

!get total number of lines and an array to put them all in

call selectlines(deeplines_catalogue, realspec(1)%wavelength/redshiftguess_overall, realspec(size(realspec))%wavelength/redshiftguess_overall, fittedlines, totallines)

if (totallines .eq. 0) then
  print *,gettime(),"Error: no known emission lines in this spectrum."
  print *,gettime(),"       Are your wavelength units correct?  Default catalogues use Angstroms"
  call exit(1)
endif

if (messages) print *, gettime(),"fitting full spectrum with ",totallines," lines"

! process the spectrum so that it only contains line regions

maskedspectrum=realspec
maskedspectrum%flux=0
do i=1,totallines
  where (abs(fittedlines(i)%wavelength-maskedspectrum%wavelength)<6)
    maskedspectrum=realspec
  endwhere
enddo

!todo: replace later use of realspec with maskedspectrum, but write out realspec at the end
realspec%flux=maskedspectrum%flux

!now go through spectrum in chunks of 440 units.  Each one overlaps by 20 units with the previous and succeeding chunk, to avoid the code attempting to fit part of a line profile
!at beginning and end, padding is only to the right and left respectively

do i=1,spectrumlength,400

! overlap=nint(2*vtol2/(1-realspec(i)%wavelength/realspec(i+1)%wavelength)) ! this needs fixing, i+1 can be out of bounds
  overlap=20

! avoid refitting final section. if spectrumlength%400<overlap, this would happen

  if ((spectrumlength-i)<overlap) exit

  if (i .eq. 1) then
    startpos=1
    startwlen=realspec(1)%wavelength/redshiftguess_overall
  else
    startpos=i-overlap
    startwlen=realspec(i)%wavelength/redshiftguess_overall
  endif

  if (i+400+overlap-1 .gt. spectrumlength) then
    endpos=spectrumlength
    endwlen=realspec(spectrumlength)%wavelength/redshiftguess_overall
  else
    endpos=i+400+overlap-1
    endwlen=realspec(i+400)%wavelength/redshiftguess_overall
  endif

  allocate(spectrumchunk(endpos-startpos+1))
  spectrumchunk = realspec(startpos:endpos)

  call selectlines(deeplines_catalogue, startwlen, endwlen, fittedlines_section, nlines)

  if (nlines .gt. 0) then
    if (messages) print "(' ',A,A,F7.1,A,F7.1,A,I3,A)",gettime(),"fitting from ",spectrumchunk(1)%wavelength," to ",spectrumchunk(size(spectrumchunk))%wavelength," with ",nlines," lines"
    call fit(spectrumchunk, redshiftguess, resolutionguess, fittedlines_section, vtol2, rtol2, generations, popsize, pressure)
    !use redshift and resolution from this chunk as initial values for next chunk
!    redshiftguess=fittedlines_section(1)%redshift
!    resolutionguess=fittedlines_section(1)%resolution
    !copy results back
    fittedlines(linearraypos:linearraypos+nlines-1)=fittedlines_section
    linearraypos=linearraypos+nlines
    !report any errors
    if (maxval(fittedlines_section%wavelength*fittedlines_section%redshift) .gt. maxval(spectrumchunk%wavelength) .or. minval(fittedlines_section%wavelength*fittedlines_section%redshift) .lt. minval(spectrumchunk%wavelength)) then
      if (messages) print *,"              Warning: some lines ended up outside the fitting region."
      if (messages) print *,"              Try reducing the value of vtol2, which is currently set to ",c*vtol2
    endif
  endif
  deallocate(spectrumchunk)
enddo

!make the fitted spectrum

call makespectrum(fittedlines, fittedspectrum)

!account for blends - for each line, determine if the subsequent line is separated by less than the resolution
!if it is, then flag that line with the lineid of the first member of the blend

if (messages) print *
if (messages) print *,gettime(),"flagging blends"

fittedlines%blended = 0

do i=1,totallines-1
  if (abs(fittedlines(i)%wavelength-fittedlines(i+1)%wavelength) .lt. 0.75*fittedlines(i)%wavelength/fittedlines(i)%resolution) then
    if (fittedlines(i)%blended .gt. 0) then
      fittedlines(i+1)%blended = fittedlines(i)%blended
    else
      fittedlines(i+1)%blended = i
    endif
  endif
enddo

! now, go through the list and add all blended line fluxes to the first line in the blend
! uncertainty at a given wavelength is independent of the line flux, so the uncertainty calculated
! at the position of the blend will apply to the total flux.

blendpeak = 0.d0

do i=1,totallines-1
  if (fittedlines(i+1)%blended .gt. 0 .and. blendpeak .eq. 0.d0) then
    !line is first in blend. add line flux to temporary variable, store line number
    blendpeak = blendpeak + fittedlines(i)%peak
  elseif (fittedlines(i+1)%blended .gt. 0 .and. blendpeak .gt. 0) then
    !line is part of multiple blend. add line flux to temporary variable, leave line number alone
    blendpeak = blendpeak + fittedlines(i)%peak
    fittedlines(i)%peak = 0.d0
!  elseif (fittedlines(i)%blended .eq. 0 .and. blendpeak .eq. 0.d0)
    !line is isolated, nothing to be done
  elseif (fittedlines(i+1)%blended .eq. 0 .and. blendpeak .gt. 0) then
    !line was last component of blend. total flux needs to be attributed to first line in blend
    blendpeak = blendpeak + fittedlines(i)%peak
    fittedlines(i)%peak = 0.d0
    fittedlines(fittedlines(i)%blended)%peak = blendpeak
    blendpeak = 0.d0
  endif
enddo

! calculate the uncertainties

if (messages) print *,gettime(),"estimating uncertainties"
call get_uncertainties(fittedspectrum, realspec, fittedlines)

! write out the fitted spectrum

open(100+tid,file=trim(outputdirectory)//trim(outputbasename)//"_fit")

write (100+tid,*) "#alfa version ",VERSION
write (100+tid,*) "#fit generated using: ",trim(commandline)
write (100+tid,*) "#""wavelength""  ""input spectrum ""  ""fitted spectrum""  ""cont-subbed orig"" ""continuum""  ""sky lines""  ""residuals""  ""uncertainty"""
do i=1,spectrumlength
  write(100+tid,"(F9.2, 7(ES12.3))") fittedspectrum(i)%wavelength,realspec(i)%flux + continuum(i)%flux, fittedspectrum(i)%flux + continuum(i)%flux + skyspectrum(i)%flux, realspec(i)%flux, continuum(i)%flux, skyspectrum(i)%flux, realspec(i)%flux - fittedspectrum(i)%flux, realspec(i)%uncertainty
enddo

close(100+tid)

! normalise if H beta is present and user did not specify a normalisation

hbetaflux=0.d0
normalisation=1.d0

do i=1,totallines
  if (abs(fittedlines(i)%wavelength - 4861.33) .lt. 0.005) then
    hbetaflux = gaussianflux(fittedlines(i)%peak,(fittedlines(i)%wavelength/fittedlines(i)%resolution))
    exit
  endif
enddo

if (.not. normalise) then

  normalisation = 1.d0

  if (hbetaflux .gt. 0.d0) then
    normalisation = 100./hbetaflux
    if (messages) print "(' ',A,A,ES9.3,A)",gettime(),"H beta detected with flux ",hbetaflux," - normalising to 100.0"
  endif

  if (normalisation .eq. 0.d0) then
    if (messages) print *,gettime(),"no H beta detected, no normalisation applied"
    normalisation = 1.d0
  endif

else

  if (normalisation .eq. 1.d0) then
    if (messages) print *,gettime(),"no normalisation applied, measured fluxes will be reported"
  else
    if (messages) print *,gettime(),"normalising to H beta = 100.0 assuming flux of ",normalisation
    normalisation = 100./normalisation
  endif

endif

fittedlines%peak = fittedlines%peak * normalisation
continuum%flux = continuum%flux * normalisation !for continuum jumps to be scaled
realspec%uncertainty = realspec%uncertainty * normalisation !for continuum jumps to be scaled

! now write out the line list.

if (maxval(fittedlines%peak) .lt. 0.1 .or. maxval(fittedlines%peak) .gt. 1.e7) then
  fluxformat="ES12.3"
else
  fluxformat="F12.3"
endif

if (messages) print *,gettime(),"writing output files ",trim(outputdirectory),trim(outputbasename),"_lines.tex and ",trim(outputdirectory),trim(outputbasename),"_fit"

if (messages) open(100+tid,file=trim(outputdirectory)//trim(outputbasename)//"_lines.tex")

open(200+tid,file=trim(outputdirectory)//trim(outputbasename)//"_lines")
write (200+tid,*) "#wlen (obs)  (rest)      flux    uncertainty     peak        FWHM"

if (messages) then
  write(100+tid,*) "\\ \hline"
  write(100+tid,*) "Observed wavelength & Rest wavelength & Flux & Uncertainty & Peak & FWHM & Ion & Multiplet & Lower term & Upper term & g$_1$ & g$_2$ \\"
  write(100+tid,*) "\hline"
endif

!check whether to write out continuum fluxes

writeb1=0
writeb2=0
writep1=0
writep2=0

if (minval(continuum%wavelength)-3630.0 .lt.0..and. maxval(continuum%wavelength)-3630. .gt.0.) writeb1=minloc(abs(fittedlines%wavelength-3630.0),1)
if (minval(continuum%wavelength)-3700.0 .lt.0..and. maxval(continuum%wavelength)-3700. .gt.0.) writeb2=minloc(abs(fittedlines%wavelength-3700.0),1)
if (minval(continuum%wavelength)-8100.0 .lt.0..and. maxval(continuum%wavelength)-8100. .gt.0.) writep1=minloc(abs(fittedlines%wavelength-8100.0),1)
if (minval(continuum%wavelength)-8400.0 .lt.0..and. maxval(continuum%wavelength)-8400. .gt.0.) writep2=minloc(abs(fittedlines%wavelength-8400.0),1)

do i=1,totallines

  if (i.eq.writeb1) then
    if (messages) write (100+tid,"(F9.2,' &           & ',"//fluxformat//",' & ',"//fluxformat//",' & Balmer cont.\\')") 3630.0, continuum(minloc(abs(continuum%wavelength-3630.)))%flux, realspec(minloc(abs(continuum%wavelength-3630.)))%uncertainty
    write (200+tid,"(2(F9.2),"//fluxformat//","//fluxformat//",'         0.0         0.0')") 3630.0, 3630.0, continuum(minloc(abs(continuum%wavelength-3630.)))%flux, realspec(minloc(abs(continuum%wavelength-3630.)))%uncertainty
  endif

  if (i.eq.writeb2) then
    if (messages) write (100+tid,"(F9.2,' &           & ',"//fluxformat//",' & ',"//fluxformat//",' & Paschen cont.\\')") 3700.0, continuum(minloc(abs(continuum%wavelength-3700.)))%flux, realspec(minloc(abs(continuum%wavelength-3700.)))%uncertainty
    write (200+tid,"(2(F9.2),"//fluxformat//","//fluxformat//",'         0.0         0.0')") 3700.0, 3700.0, continuum(minloc(abs(continuum%wavelength-3700.)))%flux, realspec(minloc(abs(continuum%wavelength-3700.)))%uncertainty
  endif

  if (i.eq.writep1) then
    if (messages) write (100+tid,"(F9.2,' &           & ',"//fluxformat//",' & ',"//fluxformat//",' & Paschen cont.\\')") 8100.0, continuum(minloc(abs(continuum%wavelength-8100.)))%flux, realspec(minloc(abs(continuum%wavelength-8100.)))%uncertainty
    write (200+tid,"(2(F9.2),"//fluxformat//","//fluxformat//",'         0.0         0.0')") 8100.0, 8100.0, continuum(minloc(abs(continuum%wavelength-8100.)))%flux, realspec(minloc(abs(continuum%wavelength-8100.)))%uncertainty
  endif

  if (i.eq.writep2) then
    if (messages) write (100+tid,"(F9.2,' &           & ',"//fluxformat//",' & ',"//fluxformat//",' & Brackett cont.\\')") 8400.0, continuum(minloc(abs(continuum%wavelength-8400.)))%flux, realspec(minloc(abs(continuum%wavelength-8400.)))%uncertainty
    write (200+tid,"(2(F9.2),"//fluxformat//","//fluxformat//",'         0.0         0.0')") 8400.0, 8400.0, continuum(minloc(abs(continuum%wavelength-8400.)))%flux, realspec(minloc(abs(continuum%wavelength-8400.)))%uncertainty
  endif

  if (fittedlines(i)%blended .eq. 0 .and. fittedlines(i)%uncertainty .gt. detectionlimit) then
    if (messages) write (100+tid,"(F9.2,' & ',F9.2,' & ',"//fluxformat//",' & ',"//fluxformat//",' & ',"//fluxformat//",' & ',"//fluxformat//",A85)") fittedlines(i)%wavelength*fittedlines(i)%redshift,fittedlines(i)%wavelength,gaussianflux(fittedlines(i)%peak,(fittedlines(i)%wavelength/fittedlines(i)%resolution)), gaussianflux(fittedlines(i)%peak,(fittedlines(i)%wavelength/fittedlines(i)%resolution))/fittedlines(i)%uncertainty, fittedlines(i)%peak/normalisation, fittedlines(i)%wavelength/fittedlines(i)%resolution * 2.35482, fittedlines(i)%linedata
    write (200+tid,"(2(F9.2),4("//fluxformat//"))") fittedlines(i)%wavelength*fittedlines(i)%redshift, fittedlines(i)%wavelength, gaussianflux(fittedlines(i)%peak,(fittedlines(i)%wavelength/fittedlines(i)%resolution)), gaussianflux(fittedlines(i)%peak,(fittedlines(i)%wavelength/fittedlines(i)%resolution))/fittedlines(i)%uncertainty,fittedlines(i)%peak/normalisation, fittedlines(i)%wavelength/fittedlines(i)%resolution * 2.35482
  elseif (fittedlines(i)%blended .ne. 0) then
    if (fittedlines(fittedlines(i)%blended)%uncertainty .gt. detectionlimit) then
      if (messages) write (100+tid,"(F9.2,' & ',F9.2,' &            * &            * &            * &            *',A85)") fittedlines(i)%wavelength*fittedlines(i)%redshift,fittedlines(i)%wavelength,fittedlines(i)%linedata
      write (200+tid,"(F9.2,F9.2,'           *           *           *           *')") fittedlines(i)%wavelength*fittedlines(i)%redshift,fittedlines(i)%wavelength
    endif
! write out 3 sigma upper limit for non-detections if upperlimits flag is set
  elseif (fittedlines(i)%uncertainty .le. detectionlimit .and. upperlimits .eqv. .true.) then
    if (messages) write (100+tid,"(F9.2,' & ',F9.2,' & ',"//fluxformat//",' & upper limit ',A85)") fittedlines(i)%wavelength*fittedlines(i)%redshift,fittedlines(i)%wavelength, 3.*gaussianflux(fittedlines(i)%peak,(fittedlines(i)%wavelength/fittedlines(i)%resolution))/fittedlines(i)%uncertainty, fittedlines(i)%linedata
    write (200+tid,"(2(F9.2),"//fluxformat//",' upper limit')") fittedlines(i)%wavelength*fittedlines(i)%redshift, fittedlines(i)%wavelength, 3.*gaussianflux(fittedlines(i)%peak,(fittedlines(i)%wavelength/fittedlines(i)%resolution))/fittedlines(i)%uncertainty
  endif
enddo

!write out measured Hbeta flux to latex table, if normalisation was applied and if output is required

if (hbetaflux .gt. 0.d0 .and. normalisation .ne. 1.d0 .and. messages) then
  write (100+tid,*) "\hline"
  write (100+tid,"(A,ES8.2)") "\multicolumn{10}{l}{Measured flux of H$\beta$: ",hbetaflux,"} \\"
endif

if (messages) write (100+tid,*) "\hline"

!done, close files

close(200+tid)
if (messages) close(100+tid) 
