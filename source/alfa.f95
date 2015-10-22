program alfa

use mod_readfiles
use mod_routines
use mod_types
use mod_continuum
use mod_fit
use mod_uncertainties

implicit none
integer :: I, spectrumlength, nlines, linearraypos, totallines, startpos, endpos, copystartpos, copyendpos
real :: startwlen, endwlen
character (len=512) :: spectrumfile,linelistfile

type(linelist), dimension(:), allocatable :: referencelinelist, fittedlines, fittedlines_section
type(spectrum), dimension(:), allocatable :: realspec, fittedspectrum, spectrumchunk, fittedchunk
type(spectrum), dimension(:), allocatable :: continuum

CHARACTER(len=2048), DIMENSION(:), allocatable :: options
CHARACTER(len=2048) :: commandline
integer :: narg

real :: redshiftguess, resolutionguess, redshifttolerance, resolutiontolerance
real :: redshiftguess_overall
real :: blendpeak
real :: normalisation, hbetaflux
real :: c

logical :: normalise=.false. !false means spectrum normalised to whatever H beta is detected, true means spectrum normalised to user specified value

c=299792.458 !km/s
!default values in absence of user specificed guess
redshiftguess=0.0 !km/s
resolutionguess=15000.

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
  if ((trim(options(i))=="-vg" .or. trim(options(i))=="--velocity-guess:") .and. (i+1) .le. Narg) then
    read (options(i+1),*) redshiftguess
  endif
  if ((trim(options(i))=="-rg" .or. trim(options(i))=="--resolution-guess") .and. (i+1) .le. Narg) then
    read (options(i+1),*) resolutionguess
  endif
enddo

! convert from velocity to redshift

redshiftguess=1.+(redshiftguess/c)

print *,gettime(),": command line: ",trim(commandline)

call get_command_argument(1,spectrumfile)

! read in spectrum to fit and line list

print *,gettime(),": reading in spectrum ",trim(spectrumfile)
call readspectrum(spectrumfile, realspec, spectrumlength, fittedspectrum)

! then subtract the continuum

print *,gettime(),": fitting continuum"
call fit_continuum(realspec,spectrumlength, continuum)

! now do the fitting
! first get guesses for the redshift and resolution

redshifttolerance=0.003 ! maximum absolute change in 1+v/c allowed from initial guess.  0.003 = 900 km/s
resolutiontolerance=resolutionguess ! maximum absolute change in lambda/delta lambda allowed from initial guess
linelistfile="linelists/strong_optical"

print "(X,A,A,F8.1,A,F7.1)",gettime(),": initial guesses for velocity and resolution: ",c*(redshiftguess-1),"km/s, R=",resolutionguess

print *,gettime(),": reading in line catalogue ",trim(linelistfile)
call readlinelist(linelistfile, referencelinelist, nlines, fittedlines, minval(realspec%wavelength), maxval(realspec%wavelength))

if (nlines .eq. 0) then
  print *,gettime(),": Warning: no reference lines detected, using default guesses for velocity and resolution"
  redshiftguess_overall=1.d0
else
  print *,gettime(),": estimating resolution and velocity using ",nlines," lines"
  call fit(realspec, referencelinelist, redshiftguess, resolutionguess, fittedlines, redshifttolerance, resolutiontolerance)
  print *,gettime(),": estimated redshift and resolution: ",c*(fittedlines(1)%redshift-1),fittedlines(1)%resolution
  redshiftguess_overall = fittedlines(1)%redshift ! when fitting chunks, use this redshift to get lines in the right range from the catalogue. if velocity from each chunk is used, then there's a chance that a line could be missed or double counted due to variations in the calculated velocity between chunks.
  redshiftguess=fittedlines(1)%redshift
  resolutionguess=fittedlines(1)%resolution
endif

! then again in chunks with tighter tolerance

redshifttolerance=0.0002 ! 60 km/s
resolutiontolerance=500.
linelistfile="linelists/deep_full"

linearraypos=1
!call readlinelist with full wavelength range to get total number of lines and an array to put them all in
print *,gettime(),": reading in line catalogue ",trim(linelistfile)
call readlinelist(linelistfile, referencelinelist, totallines, fittedlines, realspec(1)%wavelength/redshiftguess_overall, realspec(size(realspec))%wavelength/redshiftguess_overall)

print *, gettime(), ": fitting full spectrum with ",totallines," lines"

!now go through spectrum in chunks of 420 units.  Each one overlaps by 10 units with the previous and succeeding chunk, to avoid the code attempting to fit part of a line profile
!from the chunk of 420 units, only the middle 400 units is copied into the fit file, so that the overlap regions can't overwrite previously fitted lines with zeroes
!at beginning and end, padding is only to the right and left respectively

do i=1,spectrumlength,400

  if (i .eq. 1) then
    startpos=1
    startwlen=realspec(1)%wavelength/redshiftguess_overall
  else
    startpos=i-10
    startwlen=realspec(i)%wavelength/redshiftguess_overall
  endif

  if (i+409 .gt. spectrumlength) then
    endpos=spectrumlength
    endwlen=realspec(spectrumlength)%wavelength/redshiftguess_overall
  else
    endpos=i+409
    endwlen=realspec(i+400)%wavelength/redshiftguess_overall
  endif

  allocate(spectrumchunk(endpos-startpos+1))
  spectrumchunk = realspec(startpos:endpos)

  call readlinelist(linelistfile, referencelinelist, nlines, fittedlines_section, startwlen, endwlen)

  if (nlines .gt. 0) then
    print "(' ',A,A,F7.1,A,F7.1,A,I3,A)",gettime(),": fitting from ",spectrumchunk(1)%wavelength," to ",spectrumchunk(size(spectrumchunk))%wavelength," with ",nlines," lines"
!    print *,"Best fitting model parameters:       Resolution    Redshift    RMS min      RMS max"
    call fit(spectrumchunk, referencelinelist, redshiftguess, resolutionguess, fittedlines_section, redshifttolerance, resolutiontolerance)
  !use redshift and resolution from this chunk as initial values for next chunk
    redshiftguess=fittedlines_section(1)%redshift
    resolutionguess=fittedlines_section(1)%resolution
  endif

  !copy line fitting results from chunk to main array
  fittedlines(linearraypos:linearraypos+nlines-1)=fittedlines_section
  linearraypos=linearraypos+nlines

  deallocate(spectrumchunk)

enddo

!make the fitted spectrum

do i=1,size(fittedlines)
  where (abs(fittedlines(i)%redshift*fittedlines(i)%wavelength - fittedspectrum%wavelength) .lt. (5*fittedlines(i)%wavelength/fittedlines(i)%resolution))
    fittedspectrum%flux = fittedspectrum%flux + &
    &fittedlines(i)%peak*exp((-(fittedspectrum%wavelength-fittedlines(i)%redshift*fittedlines(i)%wavelength)**2)/(2*(fittedlines(i)%wavelength/fittedlines(i)%resolution)**2))
  end where
enddo

!account for blends - for each line, determine if the subsequent line is separated by less than the resolution
!if it is, then flag that line with the lineid of the first member of the blend

print *
print *,gettime(),": flagging blends"

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

print *,gettime(),": estimating uncertainties"
call get_uncertainties(fittedspectrum, realspec, fittedlines)

! write out the fitted spectrum

open(100,file=trim(spectrumfile)//"_fit")

write (100,*) """wavelength""  ""fitted spectrum""  ""cont-subbed orig"" ""continuum""  ""residuals"""
do i=1,spectrumlength
  write(100,"(F8.2, 4(ES12.3))") fittedspectrum(i)%wavelength,fittedspectrum(i)%flux + continuum(i)%flux, realspec(i)%flux, continuum(i)%flux, realspec(i)%flux - fittedspectrum(i)%flux
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
    print "(' ',A,A,ES9.3,A)",gettime(),": H beta detected with flux ",hbetaflux," - normalising to 100.0"
  endif

  if (normalisation .eq. 0.d0) then
    print *,gettime(),": no H beta detected, no normalisation applied"
    normalisation = 1.d0
  endif

else

  if (normalisation .eq. 0.0) then
    print *,gettime(),": no normalisation applied, measured fluxes will be reported"
    normalisation = 1.d0
  else
    print *,gettime(),": normalising to H beta = 100.0 assuming flux of ",normalisation
    normalisation = 100./normalisation
  endif

endif

fittedlines%peak = fittedlines%peak * normalisation
continuum%flux = continuum%flux * normalisation !for continuum jumps to be scaled
realspec%uncertainty = realspec%uncertainty * normalisation !for continuum jumps to be scaled

! now write out the line list.

print *,gettime(),": writing output files ",trim(spectrumfile),"_lines.tex and ",trim(spectrumfile),"_fit"

open(100,file=trim(spectrumfile)//"_lines.tex")
open(101,file=trim(spectrumfile)//"_lines")
write(100,*) "Observed wavelength & Rest wavelength & Flux & Uncertainty & Ion & Multiplet & Lower term & Upper term & g$_1$ & g$_2$ \\"
do i=1,totallines
print *,fittedlines(i)%wavelength,fittedlines(i)%redshift
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

print *,gettime(),": all done"
print *

end program alfa
