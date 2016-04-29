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
integer :: alloc_err

real :: wavelength, dispersion, cdelt3, cd3_3

! alfa variables

integer :: I, spectrumlength, nlines, linearraypos, totallines, startpos, endpos
real :: startwlen, endwlen, baddata
character (len=512) :: spectrumfile,stronglinelistfile,deeplinelistfile,skylinelistfile,outputdirectory
character(len=12) :: fluxformat !for writing out the line list

type(linelist), dimension(:), allocatable :: skylines_catalogue, stronglines_catalogue, deeplines_catalogue
type(linelist), dimension(:), allocatable :: fittedlines, fittedlines_section, skylines, skylines_section
type(spectrum), dimension(:), allocatable :: realspec, fittedspectrum, spectrumchunk, skyspectrum, continuum, stronglines

CHARACTER(len=2048), DIMENSION(:), allocatable :: options
CHARACTER(len=2048) :: commandline
integer :: narg, nargused

real :: redshiftguess, resolutionguess, redshiftguess_overall
real :: vtol1, vtol2, rtol1, rtol2
real :: blendpeak
real :: normalisation, hbetaflux
real :: c
integer :: linelocation, overlap
integer :: generations, popsize
real :: pressure

logical :: normalise=.false. !false means spectrum normalised to whatever H beta is detected, true means spectrum normalised to user specified value
logical :: resolution_estimated=.false. !true means user specified a value, false means estimate from sampling
logical :: subtractsky=.false. !attempt to fit night sky emission lines

! openmp variables

integer :: tid, nprocessors, omp_get_thread_num, omp_get_num_procs

c=299792.458 !km/s
!default values in absence of user specificed guess
redshiftguess=0.0 !km/s
rtol1=0.d0 !variation allowed in resolution on first pass.  determined later, either from user input, or to be equal to resolution guess.
rtol2=500. !second pass
vtol1=0.003 !variation allowed in velocity (expressed as redshift) on first pass. 0.003 = 900 km/s
vtol2=0.0002 !second pass. 0.0002 = 60 km/s
cube_i=1

stronglinelistfile=trim(PREFIX)//"/share/alfa/strong.cat"
deeplinelistfile=trim(PREFIX)//"/share/alfa/deep.cat"
skylinelistfile=trim(PREFIX)//"/share/alfa/sky.cat"

outputdirectory="./"

popsize=30
pressure=0.3 !pressure * popsize needs to be an integer
generations=500

! start

print *,"ALFACUBE, the Automated Line Fitting Algorithm for data cubes"
print *,"-------------------------------------------------------------"

print *
print *,gettime(),": starting code"

! random seed

call init_random_seed()

! read command line

narg = 0
nargused = 0 !to count options specified
narg = IARGC() !count input arguments

if (narg .eq. 0) then
  print *,"Usage: alfacube [options] [file]"
  print *,"  [file] is a FITS file with three dimensions"
  print *,"  see the man page or online documentation for details of the options"
  stop
endif

include "commandline.f90"

! convert from velocity to redshift

redshiftguess=1.+(redshiftguess/c)

print *,gettime(),": command line: ",trim(commandline)

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
! checks two keywords, and as cfitsio gives the variable a value of zero
! if not found, take the maximum value of the two as the dispersion

call ftgkye(unit,"CRVAL3",wavelength,"",status)
call ftgkye(unit,"CDELT3",cdelt3,"",status) !CD3_3
call ftgkye(unit,"CD3_3",cd3_3,"",status)

dispersion=max(cd3_3,cdelt3)
! reset status
status=0
! read data cube into memory

call ftg3de(unit,group,nullval,naxes(1),naxes(2),naxes(1),naxes(2),naxes(3),cubedata,anynull,status)

if (status .eq. 0) then
  print "(X,A,A,I7,A)",gettime(), ": successfully read ",naxes(1)*naxes(2)," pixels into memory"
else
  print *,gettime(), ": couldn't read cube into memory"
  stop
endif

! process cube
! each pixel is completely independent so it's embarrassingly parallel

nprocessors = omp_get_num_procs()
if (nprocessors .eq. 1) then
  print *,gettime(),": running in serial mode"
else
  print "(X,A,A,I2,A)", gettime(), ": running in parallel mode with ",nprocessors," processors"
endif

!read in catalogues

print *,gettime(),": reading in line catalogues"
call readlinelist(skylinelistfile, skylines_catalogue, nlines,wavelength,wavelength+(naxes(3)-1)*dispersion)
call readlinelist(stronglinelistfile, stronglines_catalogue, nlines,wavelength,wavelength+(naxes(3)-1)*dispersion)
call readlinelist(deeplinelistfile, deeplines_catalogue, nlines,wavelength,wavelength+(naxes(3)-1)*dispersion)

!process cube

print *,gettime(),": processing cube"
!$OMP PARALLEL private(spectrumfile,realspec,fittedspectrum,spectrumlength,continuum,nlines,spectrumchunk,linearraypos,overlap,startpos,startwlen,endpos,endwlen,skylines,skylines_section,stronglines,fittedlines,fittedlines_section,blendpeak,hbetaflux,totallines,skyspectrum,redshiftguess_overall,cube_i,cube_j,tid) firstprivate(redshiftguess,resolutionguess) shared(skylines_catalogue,stronglines_catalogue,deeplines_catalogue, naxes)

!$OMP DO schedule(dynamic)
do cube_i=1,naxes(1)
  do cube_j=1,naxes(2)

    tid=OMP_GET_THREAD_NUM()

    write (spectrumfile,"(A5,I3.3,A1,I3.3,A4)") "spec_",cube_i,"_",cube_j,".dat"
    allocate(realspec(naxes(3)))
    spectrumlength=naxes(3)
    realspec%flux=cubedata(cube_i,cube_j,:)
    do cube_k=1,naxes(3)
      realspec(cube_k)%wavelength=wavelength+(cube_k-1)*dispersion
    enddo

!check for valid data
!ultra crude and tailored for NGC 7009 at the moment
!    baddata=20000.
    baddata=0.
    inquire(file=trim(outputdirectory)//trim(spectrumfile)//"_lines", exist=file_exists)

    if (maxval(realspec%flux) .lt. baddata .or. file_exists) then
      print "(X,A,A,I2,A,I3.3,A,I3.3)",gettime(), "(thread ",tid,") : skipped pixel  ",cube_i,",",cube_j
      deallocate(realspec)
      cycle
    endif

    allocate (fittedspectrum(spectrumlength))
    fittedspectrum%wavelength=realspec%wavelength
    fittedspectrum%flux=0.d0

!now do the fitting
!----start of code taken from alfa.f90
!subtract the continuum

call fit_continuum(realspec,spectrumlength, continuum)

! now do the fitting
! first get guesses for the redshift and resolution

if (.not. resolution_estimated) then
  ! estimate resolution assuming nyquist sampling
  resolutionguess=2*realspec(2)%wavelength/(realspec(3)%wavelength-realspec(1)%wavelength)
else
endif

if (rtol1 .eq. 0.d0) then
  rtol1=0.9*resolutionguess ! user didn't specify a value, default behaviour is to allow resolution to vary between 0.1 and 1.9x the initial guess on the first pass
endif


! first, subtract sky spectrum if requested. do in chunks of 400 units. no overlap necessary because velocity is zero

allocate(skyspectrum(spectrumlength))
skyspectrum%wavelength=realspec%wavelength
skyspectrum%flux=0.d0

if (subtractsky) then

  !get an array for all the sky lines in the range
  call selectlines(skylines_catalogue,realspec(1)%wavelength, realspec(size(realspec))%wavelength, skylines, nlines)
  linearraypos=1

  !go though in chunks of 400 units
  do i=1,spectrumlength,400
    if (i+400 .gt. spectrumlength) then
      endpos=spectrumlength
    else
      endpos=i+400
    endif

    allocate(spectrumchunk(endpos-i+1))
    spectrumchunk = realspec(i:endpos)

    !read in sky lines in chunk
    call selectlines(skylines_catalogue, realspec(i)%wavelength, realspec(endpos)%wavelength, skylines_section, nlines)
    if (nlines .gt. 0) then
      call fit(spectrumchunk, 1., resolutionguess, skylines_section, 0., rtol2, generations, popsize, pressure)
      skylines(linearraypos:linearraypos+nlines-1)=skylines_section!(1:nlines)
      linearraypos=linearraypos+nlines
    endif
    deallocate(spectrumchunk)
  enddo
! make full sky spectrum and subtract at end

  call makespectrum(skylines,skyspectrum)
  realspec%flux = realspec%flux - skyspectrum%flux

endif

! select the lines from the catalogue

call selectlines(stronglines_catalogue, minval(realspec%wavelength), maxval(realspec%wavelength), fittedlines, nlines)

if (nlines .eq. 0) then
  redshiftguess_overall=1.d0
else
  !create an array containing just the regions around lines of interest
  !otherwise far more data is being processed than necessary
  !for each line, take 50 data points nearest to wavelength
  !should be way more than enough and still far more efficient than fitting the whole spectrum
  !todo: calculate the number of points from vtol1

  allocate(stronglines(50*nlines))
  do i=1,nlines
    linelocation=minloc(abs(stronglines_catalogue(i)%wavelength*redshiftguess-realspec%wavelength),1)
    stronglines(50*(i-1)+1:50*i) = realspec(linelocation-24:linelocation+25)
  enddo

  !now fit the strong lines

  call fit(stronglines, redshiftguess, resolutionguess, fittedlines, vtol1, rtol1, generations, popsize, pressure)

  redshiftguess_overall = fittedlines(1)%redshift ! when fitting chunks, use this redshift to get lines in the right range from the catalogue. if velocity from each chunk is used, then there's a chance that a line could be missed or double counted due to variations in the calculated velocity between chunks.
  redshiftguess=fittedlines(1)%redshift
  resolutionguess=fittedlines(1)%resolution

  deallocate(stronglines)

endif

! then again in chunks with tighter tolerance

linearraypos=1

!get total number of lines and an array to put them all in

call selectlines(deeplines_catalogue, realspec(1)%wavelength/redshiftguess_overall, realspec(size(realspec))%wavelength/redshiftguess_overall, fittedlines, totallines)


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

  call selectlines(deeplines_catalogue, startwlen, endwlen, fittedlines_section, nlines)

  if (nlines .gt. 0) then
    call fit(spectrumchunk, redshiftguess, resolutionguess, fittedlines_section, vtol2, rtol2, generations, popsize, pressure)
    !use redshift and resolution from this chunk as initial values for next chunk
    redshiftguess=fittedlines_section(1)%redshift
    resolutionguess=fittedlines_section(1)%resolution
    !copy results back
    fittedlines(linearraypos:linearraypos+nlines-1)=fittedlines_section
    linearraypos=linearraypos+nlines
    !report any errors
    if (maxval(fittedlines_section%wavelength*fittedlines_section%redshift) .gt. maxval(spectrumchunk%wavelength) .or. minval(fittedlines_section%wavelength*fittedlines_section%redshift) .lt. minval(spectrumchunk%wavelength)) then
    endif
  endif
  deallocate(spectrumchunk)
enddo

!make the fitted spectrum

call makespectrum(fittedlines, fittedspectrum)

!account for blends - for each line, determine if the subsequent line is separated by less than the resolution
!if it is, then flag that line with the lineid of the first member of the blend


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


call get_uncertainties(fittedspectrum, realspec, fittedlines)

! write out the fitted spectrum
! use thread number to avoid clashing in parallel mode

open(100+tid,file=trim(outputdirectory)//trim(spectrumfile)//"_fit")
write (100+tid,*) """wavelength""  ""fitted spectrum""  ""cont-subbed orig"" ""continuum""  ""sky lines""  ""residuals"""
do i=1,spectrumlength
  write(100+tid,"(F8.2, 5(ES12.3))") fittedspectrum(i)%wavelength,fittedspectrum(i)%flux + continuum(i)%flux + skyspectrum(i)%flux, realspec(i)%flux, continuum(i)%flux, skyspectrum(i)%flux, realspec(i)%flux - fittedspectrum(i)%flux
enddo

close(100+tid)

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
  endif

  if (normalisation .eq. 0.d0) then
    normalisation = 1.d0
  endif

else

  if (normalisation .eq. 0.0) then
    normalisation = 1.d0
  else
    normalisation = 100./normalisation
  endif

endif
fittedlines%peak = fittedlines%peak * normalisation
continuum%flux = continuum%flux * normalisation !for continuum jumps to be scaled
realspec%uncertainty = realspec%uncertainty * normalisation !for continuum jumps to be scaled

! now write out the line list.

if (maxval(fittedlines%peak) .lt. 0.1) then
  fluxformat="ES12.3"
else
  fluxformat="F12.3"
endif

open(200+tid,file=trim(outputdirectory)//trim(spectrumfile)//"_lines")

do i=1,totallines
  if (fittedlines(i)%blended .eq. 0 .and. fittedlines(i)%uncertainty .gt. 3.0) then
    write (200+tid,"(2(F8.2),2("//fluxformat//"))") fittedlines(i)%wavelength*fittedlines(i)%redshift,fittedlines(i)%wavelength, gaussianflux(fittedlines(i)%peak,(fittedlines(i)%wavelength/fittedlines(i)%resolution)), gaussianflux(fittedlines(i)%peak,(fittedlines(i)%wavelength/fittedlines(i)%resolution))/fittedlines(i)%uncertainty
  elseif (fittedlines(i)%blended .ne. 0) then
    if (fittedlines(fittedlines(i)%blended)%uncertainty .gt. 3.0) then
      write (200+tid,"(F8.2,F8.2,'           *           *')") fittedlines(i)%wavelength*fittedlines(i)%redshift,fittedlines(i)%wavelength
    endif
  endif
enddo
!write out continuum jump fluxes.  Take continuum at 3630 and 3700A as representative, but write them out as 3645.5 and 3646.5 for input into NEAT

!Balmer
if (minval(abs(continuum%wavelength-3630.)) .lt. 3630./fittedlines(1)%resolution) then
  write (200+tid,"(2(F8.2),"//fluxformat//","//fluxformat//")"), 3645.5, 3645.5, continuum(minloc(abs(continuum%wavelength-3630.)))%flux, realspec(minloc(abs(continuum%wavelength-3630.)))%uncertainty
endif

if (minval(abs(continuum%wavelength-3700.)) .lt. 3700./fittedlines(1)%resolution) then
  write (200+tid,"(2(F8.2),"//fluxformat//","//fluxformat//")"), 3646.5, 3646.5, continuum(minloc(abs(continuum%wavelength-3700.)))%flux, realspec(minloc(abs(continuum%wavelength-3700.)))%uncertainty
endif

!paschen

if (minval(abs(continuum%wavelength-8100.)) .lt. 8100./fittedlines(1)%resolution) then
  write (200+tid,"(2(F8.2),"//fluxformat//","//fluxformat//")"), 8100.0, 8100.0, continuum(minloc(abs(continuum%wavelength-8100.)))%flux, realspec(minloc(abs(continuum%wavelength-8100.)))%uncertainty
endif

if (minval(abs(continuum%wavelength-8400.)) .lt. 8400./fittedlines(1)%resolution) then
  write (200+tid,"(2(F8.2),"//fluxformat//","//fluxformat//")"), 8400.0, 8400.0, continuum(minloc(abs(continuum%wavelength-8400.)))%flux, realspec(minloc(abs(continuum%wavelength-8400.)))%uncertainty
endif

!write out measured Hbeta flux to latex table

if (hbetaflux .gt. 0.d0) then
  write (200+tid,"(A,ES8.2)") "#Measured flux of Hbeta: ",hbetaflux
endif

!done, close files

close(200+tid)

!--end of code taken from alfa.f90
!deallocate arrays ready for the next pixel
      deallocate(realspec)
      deallocate(fittedspectrum)
      deallocate(continuum)
      if (allocated(skyspectrum)) deallocate(skyspectrum)

      print "(X,A,A,I2,A,I3.3,A,I3.3)",gettime(), "(thread ",tid,") : finished pixel ",cube_i,",",cube_j

    enddo
  enddo

!$OMP END DO
!$OMP END PARALLEL

  print *,gettime(), ": all processing finished"

  deallocate(cubedata)

!  The FITS file must always be closed before exiting the program. 
!  Any unit numbers allocated with FTGIOU must be freed with FTFIOU.
  call ftclos(unit, status)
  call ftfiou(unit, status)

end program alfa_cube 