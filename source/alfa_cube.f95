program alfa_cube

use mod_readfiles
use mod_routines
use mod_types
use mod_continuum
use mod_fit
use mod_uncertainties

!alfa wrapper for analysing data cubes
!only tested on MUSE data so far
!todo: better checks for different keywords, non-linear WCS, etc etc
  implicit none

!cfitsio variables

  integer :: status,unit,readwrite,blocksize,naxes(3),nfound, hdutype
  integer :: group,firstpix,cube_i,cube_j,cube_k
  real :: nullval

  real, dimension(:,:,:), allocatable :: cubedata

  logical :: anynull,file_exists
  character :: filename*80
  integer :: alloc_err
  character(len=80) :: outname

  real :: wavelength, dispersion

! alfa variables

  integer :: I, spectrumlength, nlines, linearraypos, totallines, startpos, endpos, copystartpos, copyendpos
  real :: startwlen, endwlen
  character (len=512) :: spectrumfile,stronglinelistfile,deeplinelistfile,skylinelistfile

  type(linelist), dimension(:), allocatable :: referencelinelist, fittedlines, fittedlines_section, skylines, skylines_section
  type(spectrum), dimension(:), allocatable :: realspec, fittedspectrum, spectrumchunk, fittedchunk, skyspectrum, continuum, stronglines

  CHARACTER(len=2048), DIMENSION(:), allocatable :: options
  CHARACTER(len=2048) :: commandline
  integer :: narg

  real :: redshiftguess, resolutionguess, redshiftguess_overall
  real :: vtol1, vtol2, rtol1, rtol2
  real :: blendpeak
  real :: normalisation, hbetaflux
  real :: c
  integer :: linelocation, overlap

  logical :: normalise=.false. !false means spectrum normalised to whatever H beta is detected, true means spectrum normalised to user specified value
  logical :: resolution_estimated=.false. !true means user specified a value, false means estimate from sampling
  logical :: subtractsky=.false. !attempt to fit night sky emission lines

  c=299792.458 !km/s
  !default values in absence of user specificed guess
  redshiftguess=0.0 !km/s
  rtol1=0.d0 !variation allowed in resolution on first pass.  determined later, either from user input, or to be equal to resolution guess.
  rtol2=500. !second pass
  vtol1=0.003 !variation allowed in velocity (expressed as redshift) on first pass. 0.003 = 900 km/s
  vtol2=0.0002 !second pass. 0.0002 = 60 km/s

  stronglinelistfile="linelists/strong_optical"
  deeplinelistfile="linelists/deep_full"
  skylinelistfile="linelists/skylines"

  ! start

  print *,"ALFA, the Automated Line Fitting Algorithm"
  print *,"------------------------------------------"

  print *
  print *,gettime(),": starting code"

  ! random seed

  call init_random_seed()

  ! read command line

  narg = IARGC() !count input arguments
  if (narg .eq. 0) then
    print *,"Usage: alfa [file] [options]"
    print *,"  [file] is an ascii file with columns for wavelength and flux"
    print *,"  [options]:"
    print *,"  -n / --normalise [value]: normalise to Hb=100 assuming that F(Hb)=value"
    print *,"  -vg / --velocity-guess: initial guess for the velocity of the object [km/s]"
    print *,"  -rg / --resolution-guess: initial guess for the resolution [lambda/delta lambda]"
    print *,"  -vtol1 / --velocity-tolerance-1: variation allowed in velocity in first pass (default: 900km/s)"
    print *,"  -vtol2 / --velocity-tolerance-2: variation allowed in velocity in second pass (default: 60km/s)"
    print *,"  -rtol1 / --resolution-tolerance-1: variation allowed in resolution in first pass (default: equal to resolution guess)"
    print *,"  -rtol2 / --resolution-tolerance-2: variation allowed in resolution in second pass (default: 500.)"
    print *,"  -ss / --subtract-sky: attempt to remove night sky emission lines"
    stop
  endif

  call get_command(commandline)
  ALLOCATE (options(Narg))
  if (narg .gt. 1) then
    do i=1,Narg
      call get_command_argument(i,options(i))
    enddo
  endif

  do i=1,narg
    if ((trim(options(i))=="-n" .or. trim(options(i))=="--normalise") .and. (i+1) .le. Narg) then
      read (options(i+1),*) normalisation
      normalise=.true.
    endif
    if ((trim(options(i))=="-vg" .or. trim(options(i))=="--velocity-guess") .and. (i+1) .le. Narg) then
      read (options(i+1),*) redshiftguess
    endif
    if ((trim(options(i))=="-rg" .or. trim(options(i))=="--resolution-guess") .and. (i+1) .le. Narg) then
      read (options(i+1),*) resolutionguess
      resolution_estimated=.true.
    endif
    if ((trim(options(i))=="-vtol1" .or. trim(options(i))=="--velocity-tolerance-1") .and. (i+1) .le. Narg) then
      read (options(i+1),*) vtol1
      vtol1 = vtol1/c
    endif
    if ((trim(options(i))=="-vtol2" .or. trim(options(i))=="--velocity-tolerance-2") .and. (i+1) .le. Narg) then
      read (options(i+1),*) vtol2
      vtol2 = vtol2/c
    endif
    if ((trim(options(i))=="-rtol1" .or. trim(options(i))=="--resolution-tolerance-1") .and. (i+1) .le. Narg) then
      read (options(i+1),*) rtol1
    endif
    if ((trim(options(i))=="-rtol2" .or. trim(options(i))=="--resolution-tolerance-2") .and. (i+1) .le. Narg) then
      read (options(i+1),*) rtol2
    endif
    if ((trim(options(i))=="-ss" .or. trim(options(i))=="--subtract-sky")) then
      subtractsky=.true.
    endif
  enddo

  ! convert from velocity to redshift

  redshiftguess=1.+(redshiftguess/c)

  print *,gettime(),": command line: ",trim(commandline)

  call get_command_argument(1,spectrumfile)
  status=0
  !  Get an unused Logical Unit Number to use to open the FITS file.
  call ftgiou(unit,status)
  !  Open the FITS file
  inquire(file=spectrumfile, exist=file_exists)
  if (.not. file_exists) then
    print *,gettime(),trim(spectrumfile)," does not exist"
    stop
  endif

  readwrite=0
  call ftopen(unit,spectrumfile,readwrite,blocksize,status)

  ! check we have 3 axes
  nfound=0
  call ftgknj(unit,'NAXIS',1,3,naxes,nfound,status)
  do while (nfound .eq. 0)
    call ftmrhd(unit,1,hdutype,status)
    call ftgknj(unit,'NAXIS',1,3,naxes,nfound,status)
  end do

  if (nfound .ne. 3) then
    print *,gettime(),spectrumfile," is not a cube"
    stop
  endif

  group=1
  firstpix=1
  nullval=-999

  if (nfound .eq. 3) then
    allocate(cubedata(naxes(1),naxes(2),naxes(3)), stat=alloc_err)
    if (alloc_err .eq. 0) print *,gettime(), ": reading data cube into memory"
  endif

! find wavelength dispersion

  call ftgkye(unit,"CRVAL3",wavelength,"",status)
  call ftgkye(unit,"CD3_3",dispersion,"",status)

! read data cube into memory

  call ftg3de(unit,group,nullval,naxes(1),naxes(2),naxes(1),naxes(2),naxes(3),cubedata,anynull,status)

  if (status .eq. 0) then
    print *,gettime(), ": successfully read cube into memory"
    print *,gettime(), ": will fit ",naxes(1)*naxes(2)," pixels"
    print *,gettime(), ": details of fitting will be written to the file 'cubeanalysis.log'"
    print *
  else
    print *,gettime(), ": couldn't read cube into memory"
    stop
  endif

  open(4425,file="cubeanalysis.log")

! process cube
  do cube_i=1,naxes(1)
    do cube_j=1,naxes(2)

      write (spectrumfile,"(A5,I3.3,A1,I3.3,A4)") "spec_",cube_i,"_",cube_j,".dat"
      allocate(realspec(naxes(3)))
      spectrumlength=naxes(3)
      realspec%flux=cubedata(cube_i,cube_j,:)
      do cube_k=1,naxes(3)
        realspec(cube_k)%wavelength=wavelength+(cube_k-1)*dispersion
      enddo

!check for valid data
!ultra crude and tailored for NGC 7009 at the moment

      if (maxval(realspec%flux) .lt. 20000.) then
        print *,gettime(), ": skipping pixel ",cube_i,cube_j
        write (4425,*) gettime(), ": skipping pixel ",cube_i,cube_j
        deallocate(realspec)
        cycle
      else
        print *,gettime(), ": fitting pixel ",cube_i,cube_j
        write (4425,*) gettime(), ": fitting pixel ",cube_i,cube_j
      endif

      allocate (fittedspectrum(spectrumlength))
      fittedspectrum%wavelength=realspec%wavelength
      fittedspectrum%flux=0.d0

!now do the fitting
!----start of code taken from alfa.f95
!subtract the continuum

write(4425,*) gettime(),": fitting continuum"
call fit_continuum(realspec,spectrumlength, continuum)

! now do the fitting
! first get guesses for the redshift and resolution

if (.not. resolution_estimated) then
  ! estimate resolution assuming nyquist sampling
  resolutionguess=2*realspec(2)%wavelength/(realspec(3)%wavelength-realspec(1)%wavelength)
  write(4425,*) gettime(),": estimated spectrograph resolution assuming Nyquist sampling: ",resolutionguess
else
  write(4425,*) gettime(),": estimated spectrograph resolution from user input: ",resolutionguess
endif

if (rtol1 .eq. 0.d0) then
  rtol1=0.9*resolutionguess ! user didn't specify a value, default behaviour is to allow resolution to vary between 0 and 2x the initial guess on the first pass
endif

write(4425,"(X,A,A,F8.1,A,F7.1)") gettime(),": initial guesses for velocity and resolution: ",c*(redshiftguess-1),"km/s, R=",resolutionguess

! first, subtract sky spectrum if requested. do in chunks of 400 units. no overlap necessary because velocity is zero

allocate(skyspectrum(spectrumlength))
skyspectrum%wavelength=realspec%wavelength
skyspectrum%flux=0.d0

if (subtractsky) then
  write(4425,*) gettime(),": fitting sky emission"

  !get an array for all the sky lines in the range
  call readlinelist(skylinelistfile, referencelinelist, nlines, skylines, realspec(1)%wavelength, realspec(size(realspec))%wavelength)

  linearraypos=1

  !go though in chunks of 400 units
  do i=1,spectrumlength,400
    if (i+400 .gt. spectrumlength) then
      endpos=spectrumlength
    else
      endpos=i+399
    endif

    allocate(spectrumchunk(endpos-i+1))
    spectrumchunk = realspec(i:endpos)

    !read in sky lines in chunk
    call readlinelist(skylinelistfile, referencelinelist, nlines, skylines_section, realspec(i)%wavelength, realspec(i+399)%wavelength)
    if (nlines .gt. 0) then
      call fit(spectrumchunk, referencelinelist, 1., resolutionguess, skylines_section, 0., rtol2)
      skylines(linearraypos:linearraypos+nlines-1)=skylines_section!(1:nlines)
      linearraypos=linearraypos+nlines
    endif
    deallocate(spectrumchunk)
  enddo
! make full sky spectrum and subtract at end
  call makespectrum(skylines,skyspectrum)
  realspec%flux = realspec%flux - skyspectrum%flux

endif

! read in line catalogue

write(4425,*) gettime(),": reading in line catalogue ",trim(stronglinelistfile)
call readlinelist(stronglinelistfile, referencelinelist, nlines, fittedlines, minval(realspec%wavelength), maxval(realspec%wavelength))

if (nlines .eq. 0) then
  write(4425,*) gettime(),": Warning: no reference lines detected, using default guesses for velocity and resolution"
  redshiftguess_overall=1.d0
else
  write(4425,*) gettime(),": estimating resolution and velocity using ",nlines," lines"
  !create an array containing just the regions around lines of interest
  !otherwise far more data is being processed than necessary
  !for each line, take 50 data points nearest to wavelength
  !should be way more than enough and still far more efficient than fitting the whole spectrum
  !todo: calculate the number of points from vtol1

  allocate(stronglines(50*nlines))
  do i=1,nlines
    linelocation=minloc(abs(referencelinelist(i)%wavelength*redshiftguess-realspec%wavelength),1)
    stronglines(50*(i-1)+1:50*i) = realspec(linelocation-24:linelocation+25)
  enddo

  !now fit the strong lines

  call fit(stronglines, referencelinelist, redshiftguess, resolutionguess, fittedlines, vtol1, rtol1)

  write(4425,*) gettime(),": estimated redshift and resolution: ",c*(fittedlines(1)%redshift-1),fittedlines(1)%resolution
  redshiftguess_overall = fittedlines(1)%redshift ! when fitting chunks, use this redshift to get lines in the right range from the catalogue. if velocity from each chunk is used, then there's a chance that a line could be missed or double counted due to variations in the calculated velocity between chunks.
  redshiftguess=fittedlines(1)%redshift
  resolutionguess=fittedlines(1)%resolution

  deallocate(stronglines)

endif

! then again in chunks with tighter tolerance

linearraypos=1

!call readlinelist with full wavelength range to get total number of lines and an array to put them all in

write(4425,*) gettime(),": reading in line catalogue ",trim(deeplinelistfile)
call readlinelist(deeplinelistfile, referencelinelist, totallines, fittedlines, realspec(1)%wavelength/redshiftguess_overall, realspec(size(realspec))%wavelength/redshiftguess_overall)

write(4425,*)  gettime(), ": fitting full spectrum with ",totallines," lines"

!now go through spectrum in chunks of 440 units.  Each one overlaps by 20 units with the previous and succeeding chunk, to avoid the code attempting to fit part of a line profile
!at beginning and end, padding is only to the right and left respectively

do i=1,spectrumlength,400

  overlap=nint(2*vtol2/(1-realspec(i)%wavelength/realspec(i+1)%wavelength))

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

  call readlinelist(deeplinelistfile, referencelinelist, nlines, fittedlines_section, startwlen, endwlen)

  if (nlines .gt. 0) then
    write(4425,"(' ',A,A,F7.1,A,F7.1,A,I3,A)") gettime(),": fitting from ",spectrumchunk(1)%wavelength," to ",spectrumchunk(size(spectrumchunk))%wavelength," with ",nlines," lines"
    call fit(spectrumchunk, referencelinelist, redshiftguess, resolutionguess, fittedlines_section, vtol2, rtol2)
    !use redshift and resolution from this chunk as initial values for next chunk
    redshiftguess=fittedlines_section(1)%redshift
    resolutionguess=fittedlines_section(1)%resolution
    !copy results back
    fittedlines(linearraypos:linearraypos+nlines-1)=fittedlines_section
    linearraypos=linearraypos+nlines
    !report any errors
    if (maxval(fittedlines_section%wavelength*fittedlines_section%redshift) .gt. maxval(spectrumchunk%wavelength) .or. minval(fittedlines_section%wavelength*fittedlines_section%redshift) .lt. minval(spectrumchunk%wavelength)) then
      write(4425,*) "              Warning: some lines ended up outside the fitting region."
      write(4425,*) "              Try reducing the value of vtol2, which is currently set to ",c*vtol2
    endif
  endif
  deallocate(spectrumchunk)
enddo

!make the fitted spectrum

call makespectrum(fittedlines, fittedspectrum)

!account for blends - for each line, determine if the subsequent line is separated by less than the resolution
!if it is, then flag that line with the lineid of the first member of the blend

write(4425,*) gettime(),": flagging blends"

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

write(4425,*) gettime(),": estimating uncertainties"

call get_uncertainties(fittedspectrum, realspec, fittedlines)

! write out the fitted spectrum

open(100,file=trim(spectrumfile)//"_fit")

write (100,*) """wavelength""  ""fitted spectrum""  ""cont-subbed orig"" ""continuum""  ""sky lines""  ""residuals"""
do i=1,spectrumlength
  write(100,"(F8.2, 5(ES12.3))") fittedspectrum(i)%wavelength,fittedspectrum(i)%flux + continuum(i)%flux + skyspectrum(i)%flux, realspec(i)%flux, continuum(i)%flux, skyspectrum(i)%flux, realspec(i)%flux - fittedspectrum(i)%flux
enddo

close(100)

! normalise if H beta is present and user did not specify a normalisation

hbetaflux=0.d0
do i=1,totallines
  if (abs(fittedlines(i)%wavelength - 4861.33) .lt. 0.005) then
    hbetaflux = gaussianflux(fittedlines(i)%peak,(fittedlines(i)%wavelength/fittedlines(i)%resolution))
  endif
enddo

if (.not. normalise) then

  normalisation = 0.d0

  if (hbetaflux .gt. 0.d0) then
    normalisation = 100./hbetaflux
    write(4425,"(' ',A,A,ES9.3,A)") gettime(),": H beta detected with flux ",hbetaflux," - normalising to 100.0"
  endif

  if (normalisation .eq. 0.d0) then
    write(4425,*) gettime(),": no H beta detected, no normalisation applied"
    normalisation = 1.d0
  endif

else

  if (normalisation .eq. 0.0) then
    write(4425,*) gettime(),": no normalisation applied, measured fluxes will be reported"
    normalisation = 1.d0
  else
    write(4425,*) gettime(),": normalising to H beta = 100.0 assuming flux of ",normalisation
    normalisation = 100./normalisation
  endif

endif
fittedlines%peak = fittedlines%peak * normalisation
continuum%flux = continuum%flux * normalisation !for continuum jumps to be scaled
realspec%uncertainty = realspec%uncertainty * normalisation !for continuum jumps to be scaled

! now write out the line list.

write(4425,*) gettime(),": writing output files ",trim(spectrumfile),"_lines.tex and ",trim(spectrumfile),"_fit"

open(100,file=trim(spectrumfile)//"_lines.tex")
open(101,file=trim(spectrumfile)//"_lines")
write(100,*) "Observed wavelength & Rest wavelength & Flux & Uncertainty & Ion & Multiplet & Lower term & Upper term & g$_1$ & g$_2$ \\"
do i=1,totallines
  if (fittedlines(i)%blended .eq. 0 .and. fittedlines(i)%uncertainty .gt. 3.0) then
    write (100,"(F8.2,' & ',F8.2,' & ',F12.3,' & ',F12.3,A85)") fittedlines(i)%wavelength*fittedlines(i)%redshift,fittedlines(i)%wavelength,gaussianflux(fittedlines(i)%peak,(fittedlines(i)%wavelength/fittedlines(i)%resolution)), gaussianflux(fittedlines(i)%peak,(fittedlines(i)%wavelength/fittedlines(i)%resolution))/fittedlines(i)%uncertainty, fittedlines(i)%linedata
    write (101,"(F8.2,2(F12.3))") fittedlines(i)%wavelength, gaussianflux(fittedlines(i)%peak,(fittedlines(i)%wavelength/fittedlines(i)%resolution)), gaussianflux(fittedlines(i)%peak,(fittedlines(i)%wavelength/fittedlines(i)%resolution))/fittedlines(i)%uncertainty
  elseif (fittedlines(i)%blended .ne. 0) then
    if (fittedlines(fittedlines(i)%blended)%uncertainty .gt. 3.0) then
      write (100,"(F8.2,' & ',F8.2,' &            * &            *',A85)") fittedlines(i)%wavelength*fittedlines(i)%redshift,fittedlines(i)%wavelength,fittedlines(i)%linedata
      write (101,"(F8.2,'           *           *')") fittedlines(i)%wavelength
    endif
  endif
enddo
!write out continuum jump fluxes.  Take continuum at 3630 and 3700A as representative, but write them out as 3645.5 and 3646.5 for input into NEAT

!Balmer
if (minval(abs(continuum%wavelength-3630.)) .lt. 3630./fittedlines(1)%resolution) then
  write (100,"(F8.2,' &          & ',F12.3,' & ',F12.3,' & Balmer jump-\\')") 3645.5, continuum(minloc(abs(continuum%wavelength-3630.)))%flux, realspec(minloc(abs(continuum%wavelength-3630.)))%uncertainty
  write (101,"(F8.2,F12.3,F12.3)") 3645.5, continuum(minloc(abs(continuum%wavelength-3630.)))%flux, realspec(minloc(abs(continuum%wavelength-3630.)))%uncertainty
endif

if (minval(abs(continuum%wavelength-3700.)) .lt. 3700./fittedlines(1)%resolution) then
  write (100,"(F8.2,' &          & ',F12.3,' & ',F12.3,' & Balmer jump+\\')") 3646.5, continuum(minloc(abs(continuum%wavelength-3700.)))%flux, realspec(minloc(abs(continuum%wavelength-3700.          )))%uncertainty
  write (101,"(F8.2,F12.3,F12.3)") 3646.5, continuum(minloc(abs(continuum%wavelength-3700.)))%flux, realspec(minloc(abs(continuum%wavelength-3700.)))%uncertainty
endif

!paschen

if (minval(abs(continuum%wavelength-8100.)) .lt. 8100./fittedlines(1)%resolution) then
  write (100,"(F8.2,' &          & ',F12.3,' & ',F12.3,' & Paschen jump-\\')") 8100.0, continuum(minloc(abs(continuum%wavelength-8100.)))%flux, realspec(minloc(abs(continuum%wavelength-8100.          )))%uncertainty
  write (101,"(F8.2,F12.3,F12.3)") 8100.0, continuum(minloc(abs(continuum%wavelength-8100.)))%flux, realspec(minloc(abs(continuum%wavelength-8100.)))%uncertainty
endif

if (minval(abs(continuum%wavelength-8400.)) .lt. 8400./fittedlines(1)%resolution) then
  write (100,"(F8.2,' &          & ',F12.3,' & ',F12.3,' & Paschen jump+\\')") 8400.0, continuum(minloc(abs(continuum%wavelength-8400.)))%flux, realspec(minloc(abs(continuum%wavelength-8400.          )))%uncertainty
  write (101,"(F8.2,F12.3,F12.3)") 8400.0, continuum(minloc(abs(continuum%wavelength-8400.)))%flux, realspec(minloc(abs(continuum%wavelength-8400.)))%uncertainty
endif

!write out measured Hbeta flux to latex table

if (hbetaflux .gt. 0.d0) then
  write (100,*) "\hline"
  write (100,"(A,ES8.2)") "Measured flux of H$\beta$: ",hbetaflux
  write (100,*) "\hline"
endif

!done, close files

close(101)
close(100)

write(4425,*) gettime(),": all done"
write(4425,*)

!--end of code taken from alfa.f95
!deallocate arrays ready for the next pixel
      deallocate(realspec)
      deallocate(fittedspectrum)
      deallocate(continuum)
      if (allocated(skyspectrum)) deallocate(skyspectrum)
    enddo
  enddo

  deallocate(cubedata)
  close(4425)

!  The FITS file must always be closed before exiting the program. 
!  Any unit numbers allocated with FTGIOU must be freed with FTFIOU.
  call ftclos(unit, status)
  call ftfiou(unit, status)

end program alfa_cube 
