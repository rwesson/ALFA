module mod_output
use mod_functions
use mod_globals

contains

subroutine write_plaintext

  open(100+threadnumber,file=trim(outputdirectory)//trim(outputbasename)//"_fit")

  write (100+threadnumber,*) "#alfa version ",VERSION
  write (100+threadnumber,*) "#fit generated using: ",trim(commandline)
  write (100+threadnumber,*) "#""wavelength""  ""input spectrum ""  ""fitted spectrum""  ""cont-subbed orig"" ""continuum""  ""sky lines""  ""residuals""  ""uncertainty"""
  do i=1,spectrumlength
    write(100+threadnumber,"(F9.2, 7(ES12.3))") fittedspectrum(i)%wavelength,realspec(i)%flux + continuum(i)%flux, fittedspectrum(i)%flux + continuum(i)%flux + skyspectrum(i)%flux, realspec(i)%flux, continuum(i)%flux, skyspectrum(i)%flux, realspec(i)%flux - fittedspectrum(i)%flux, maskedspectrum(i)%uncertainty
  enddo

  close(100+threadnumber)

! now write out the line list.

  if (maxval(fittedlines%peak) .lt. 0.1 .or. maxval(fittedlines%peak) .gt. 1.e7) then
    fluxformat="ES12.3"
  else
    fluxformat="F12.3"
  endif

  if (messages) print *,gettime(),"writing output files ",trim(outputdirectory),trim(outputbasename),"_lines and ",trim(outputdirectory),trim(outputbasename),"_fit"

  open(200+threadnumber,file=trim(outputdirectory)//trim(outputbasename)//"_lines")
  write (200+threadnumber,*) "#wlen (obs)  (rest)      flux    uncertainty     peak        FWHM"

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
      write (200+threadnumber,"(2(F9.2),"//fluxformat//","//fluxformat//",'         0.0         0.0')") 3630.0, 3630.0, continuum(minloc(abs(continuum%wavelength-3630.)))%flux, realspec(minloc(abs(continuum%wavelength-3630.)))%uncertainty
    endif

    if (i.eq.writeb2) then
      write (200+threadnumber,"(2(F9.2),"//fluxformat//","//fluxformat//",'         0.0         0.0')") 3700.0, 3700.0, continuum(minloc(abs(continuum%wavelength-3700.)))%flux, realspec(minloc(abs(continuum%wavelength-3700.)))%uncertainty
    endif

    if (i.eq.writep1) then
      write (200+threadnumber,"(2(F9.2),"//fluxformat//","//fluxformat//",'         0.0         0.0')") 8100.0, 8100.0, continuum(minloc(abs(continuum%wavelength-8100.)))%flux, realspec(minloc(abs(continuum%wavelength-8100.)))%uncertainty
    endif

    if (i.eq.writep2) then
      write (200+threadnumber,"(2(F9.2),"//fluxformat//","//fluxformat//",'         0.0         0.0')") 8400.0, 8400.0, continuum(minloc(abs(continuum%wavelength-8400.)))%flux, realspec(minloc(abs(continuum%wavelength-8400.)))%uncertainty
    endif

    if (fittedlines(i)%blended .eq. 0 .and. fittedlines(i)%uncertainty .gt. detectionlimit) then
      write (200+threadnumber,"(2(F9.2),4("//fluxformat//"))") fittedlines(i)%wavelength*fittedlines(i)%redshift, fittedlines(i)%wavelength, gaussianflux(fittedlines(i)%peak,(fittedlines(i)%wavelength/fittedlines(i)%resolution)), gaussianflux(fittedlines(i)%peak,(fittedlines(i)%wavelength/fittedlines(i)%resolution))/fittedlines(i)%uncertainty,fittedlines(i)%peak/normalisation, fittedlines(i)%wavelength/fittedlines(i)%resolution * 2.35482
    elseif (fittedlines(i)%blended .ne. 0) then
      if (fittedlines(fittedlines(i)%blended)%uncertainty .gt. detectionlimit) then
        write (200+threadnumber,"(F9.2,F9.2,'           *           *           *           *')") fittedlines(i)%wavelength*fittedlines(i)%redshift,fittedlines(i)%wavelength
      endif
! write out 3 sigma upper limit for non-detections if upperlimits flag is set
    elseif (fittedlines(i)%uncertainty .le. detectionlimit .and. upperlimits .eqv. .true.) then
      write (200+threadnumber,"(2(F9.2),"//fluxformat//",' upper limit')") fittedlines(i)%wavelength*fittedlines(i)%redshift, fittedlines(i)%wavelength, 3.*gaussianflux(fittedlines(i)%peak,(fittedlines(i)%wavelength/fittedlines(i)%resolution))/fittedlines(i)%uncertainty
    endif
  enddo

!write out measured Hbeta flux if normalisation was applied

  if (hbetaflux .gt. 0.d0 .and. normalisation .ne. 1.d0) then
    write (200+threadnumber,"(A,ES8.2)") "# measured flux of Hbeta: ",hbetaflux
  endif

!done, close files

  close(200+threadnumber)

end subroutine write_plaintext




subroutine write_latex

  open(100+threadnumber,file=trim(outputdirectory)//trim(outputbasename)//"_fit.tex")

  write (100+threadnumber,*) "%alfa version ",VERSION
  write (100+threadnumber,*) "%fit generated using: ",trim(commandline)
  write (100+threadnumber,*) "#""wavelength""  ""input spectrum ""  ""fitted spectrum""  ""cont-subbed orig"" ""continuum""  ""sky lines""  ""residuals""  ""uncertainty"""
  do i=1,spectrumlength
    write(100+threadnumber,"(F9.2, 7(ES12.3))") fittedspectrum(i)%wavelength,realspec(i)%flux + continuum(i)%flux, fittedspectrum(i)%flux + continuum(i)%flux + skyspectrum(i)%flux, realspec(i)%flux, continuum(i)%flux, skyspectrum(i)%flux, realspec(i)%flux - fittedspectrum(i)%flux, maskedspectrum(i)%uncertainty
  enddo

  close(100+threadnumber)

! now write out the line list.

  if (maxval(fittedlines%peak) .lt. 0.1 .or. maxval(fittedlines%peak) .gt. 1.e7) then
    fluxformat="ES12.3"
  else
    fluxformat="F12.3"
  endif

  if (messages) print *,gettime(),"writing output files ",trim(outputdirectory),trim(outputbasename),"_lines.tex and ",trim(outputdirectory),trim(outputbasename),"_fit.tex"

  open(200+threadnumber,file=trim(outputdirectory)//trim(outputbasename)//"_lines.tex")

  write(200+threadnumber,*) "\\ \hline"
  write(200+threadnumber,*) "Observed wavelength & Rest wavelength & Flux & Uncertainty & Peak & FWHM & Ion & Multiplet & Lower term & Upper term & g$_1$ & g$_2$ \\"
  write(200+threadnumber,*) "\hline"

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
      write (200+threadnumber,"(F9.2,' &           & ',"//fluxformat//",' & ',"//fluxformat//",' & Balmer cont.\\')") 3630.0, continuum(minloc(abs(continuum%wavelength-3630.)))%flux, realspec(minloc(abs(continuum%wavelength-3630.)))%uncertainty
    endif

    if (i.eq.writeb2) then
      write (200+threadnumber,"(F9.2,' &           & ',"//fluxformat//",' & ',"//fluxformat//",' & Paschen cont.\\')") 3700.0, continuum(minloc(abs(continuum%wavelength-3700.)))%flux, realspec(minloc(abs(continuum%wavelength-3700.)))%uncertainty
    endif

    if (i.eq.writep1) then
      write (200+threadnumber,"(F9.2,' &           & ',"//fluxformat//",' & ',"//fluxformat//",' & Paschen cont.\\')") 8100.0, continuum(minloc(abs(continuum%wavelength-8100.)))%flux, realspec(minloc(abs(continuum%wavelength-8100.)))%uncertainty
    endif

    if (i.eq.writep2) then
      write (200+threadnumber,"(F9.2,' &           & ',"//fluxformat//",' & ',"//fluxformat//",' & Brackett cont.\\')") 8400.0, continuum(minloc(abs(continuum%wavelength-8400.)))%flux, realspec(minloc(abs(continuum%wavelength-8400.)))%uncertainty
    endif

    if (fittedlines(i)%blended .eq. 0 .and. fittedlines(i)%uncertainty .gt. detectionlimit) then
      write (200+threadnumber,"(F9.2,' & ',F9.2,' & ',"//fluxformat//",' & ',"//fluxformat//",' & ',"//fluxformat//",' & ',"//fluxformat//",A85)") fittedlines(i)%wavelength*fittedlines(i)%redshift,fittedlines(i)%wavelength,gaussianflux(fittedlines(i)%peak,(fittedlines(i)%wavelength/fittedlines(i)%resolution)), gaussianflux(fittedlines(i)%peak,(fittedlines(i)%wavelength/fittedlines(i)%resolution))/fittedlines(i)%uncertainty, fittedlines(i)%peak/normalisation, fittedlines(i)%wavelength/fittedlines(i)%resolution * 2.35482, fittedlines(i)%linedata
    elseif (fittedlines(i)%blended .ne. 0) then
      if (fittedlines(fittedlines(i)%blended)%uncertainty .gt. detectionlimit) then
        write (200+threadnumber,"(F9.2,' & ',F9.2,' &            * &            * &            * &            *',A85)") fittedlines(i)%wavelength*fittedlines(i)%redshift,fittedlines(i)%wavelength,fittedlines(i)%linedata
      endif
! write out 3 sigma upper limit for non-detections if upperlimits flag is set
    elseif (fittedlines(i)%uncertainty .le. detectionlimit .and. upperlimits .eqv. .true.) then
      write (200+threadnumber,"(F9.2,' & ',F9.2,' & ',"//fluxformat//",' & upper limit ',A85)") fittedlines(i)%wavelength*fittedlines(i)%redshift,fittedlines(i)%wavelength, 3.*gaussianflux(fittedlines(i)%peak,(fittedlines(i)%wavelength/fittedlines(i)%resolution))/fittedlines(i)%uncertainty, fittedlines(i)%linedata
    endif
  enddo

!write out measured Hbeta flux to latex table, if normalisation was applied and if output is required

  if (hbetaflux .gt. 0.d0 .and. normalisation .ne. 1.d0) then
    write (200+threadnumber,*) "\hline"
    write (200+threadnumber,"(A,ES8.2)") "\multicolumn{10}{l}{Measured flux of H$\beta$: ",hbetaflux,"} \\"
  endif

  write (200+threadnumber,*) "\hline"

!done, close files

  close(200+threadnumber)

end subroutine write_latex




subroutine write_csv

  open(100+threadnumber,file=trim(outputdirectory)//trim(outputbasename)//"_fit.csv")

  write (100+threadnumber,*) "#alfa version ",VERSION
  write (100+threadnumber,*) "#fit generated using: ",trim(commandline)
  write (100+threadnumber,*) "#wavelength,input spectrum,fitted spectrum,cont-subbed orig,continuum,sky lines,residuals,uncertainty"
  do i=1,spectrumlength
    write(100+threadnumber,"(F9.2, ',', 6(ES12.3,','), ES12.3)") fittedspectrum(i)%wavelength,realspec(i)%flux + continuum(i)%flux, fittedspectrum(i)%flux + continuum(i)%flux + skyspectrum(i)%flux, realspec(i)%flux, continuum(i)%flux, skyspectrum(i)%flux, realspec(i)%flux - fittedspectrum(i)%flux, maskedspectrum(i)%uncertainty
  enddo

  close(100+threadnumber)

! now write out the line list.

  if (maxval(fittedlines%peak) .lt. 0.1 .or. maxval(fittedlines%peak) .gt. 1.e7) then
    fluxformat="ES12.3"
  else
    fluxformat="F12.3"
  endif

  if (messages) print *,gettime(),"writing output files ",trim(outputdirectory),trim(outputbasename),"_lines.csv and ",trim(outputdirectory),trim(outputbasename),"_fit.csv"

  open(200+threadnumber,file=trim(outputdirectory)//trim(outputbasename)//"_lines.csv")
  write (200+threadnumber,*) "#wlen (obs),(rest),flux,uncertainty,peak,FWHM"

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
      write (200+threadnumber,"(2(F9.2,','),"//fluxformat//",',',"//fluxformat//",',         0.0,         0.0')") 3630.0, 3630.0, continuum(minloc(abs(continuum%wavelength-3630.)))%flux, realspec(minloc(abs(continuum%wavelength-3630.)))%uncertainty
    endif

    if (i.eq.writeb2) then
      write (200+threadnumber,"(2(F9.2,','),"//fluxformat//",',',"//fluxformat//",',         0.0,         0.0')") 3700.0, 3700.0, continuum(minloc(abs(continuum%wavelength-3700.)))%flux, realspec(minloc(abs(continuum%wavelength-3700.)))%uncertainty
    endif

    if (i.eq.writep1) then
      write (200+threadnumber,"(2(F9.2,','),"//fluxformat//",',',"//fluxformat//",',         0.0,         0.0')") 8100.0, 8100.0, continuum(minloc(abs(continuum%wavelength-8100.)))%flux, realspec(minloc(abs(continuum%wavelength-8100.)))%uncertainty
    endif

    if (i.eq.writep2) then
      write (200+threadnumber,"(2(F9.2,','),"//fluxformat//",',',"//fluxformat//",',         0.0,         0.0')") 8400.0, 8400.0, continuum(minloc(abs(continuum%wavelength-8400.)))%flux, realspec(minloc(abs(continuum%wavelength-8400.)))%uncertainty
    endif

    if (fittedlines(i)%blended .eq. 0 .and. fittedlines(i)%uncertainty .gt. detectionlimit) then
      write (200+threadnumber,"(2(F9.2,','),3("//fluxformat//",','),"//fluxformat//")") fittedlines(i)%wavelength*fittedlines(i)%redshift, fittedlines(i)%wavelength, gaussianflux(fittedlines(i)%peak,(fittedlines(i)%wavelength/fittedlines(i)%resolution)), gaussianflux(fittedlines(i)%peak,(fittedlines(i)%wavelength/fittedlines(i)%resolution))/fittedlines(i)%uncertainty,fittedlines(i)%peak/normalisation, fittedlines(i)%wavelength/fittedlines(i)%resolution * 2.35482
    elseif (fittedlines(i)%blended .ne. 0) then
      if (fittedlines(fittedlines(i)%blended)%uncertainty .gt. detectionlimit) then
        write (200+threadnumber,"(F9.2,','F9.2,',',3(A12,','),A12)") fittedlines(i)%wavelength*fittedlines(i)%redshift,fittedlines(i)%wavelength,"*","*","*","*"
      endif
! write out 3 sigma upper limit for non-detections if upperlimits flag is set
    elseif (fittedlines(i)%uncertainty .le. detectionlimit .and. upperlimits .eqv. .true.) then
      write (200+threadnumber,"(2(F9.2,','),"//fluxformat//",' upper limit')") fittedlines(i)%wavelength*fittedlines(i)%redshift, fittedlines(i)%wavelength, 3.*gaussianflux(fittedlines(i)%peak,(fittedlines(i)%wavelength/fittedlines(i)%resolution))/fittedlines(i)%uncertainty
    endif
  enddo

!write out measured Hbeta flux if normalisation was applied

  if (hbetaflux .gt. 0.d0 .and. normalisation .ne. 1.d0) then
    write (200+threadnumber,"(A,ES8.2)") "# measured flux of Hbeta: ",hbetaflux
  endif

!done, close files

  close(200+threadnumber)


end subroutine write_csv

subroutine write_fits
end subroutine write_fits

end module mod_output
