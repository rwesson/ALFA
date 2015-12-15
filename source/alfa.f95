program alfa

use mod_readfiles
use mod_routines
use mod_types
use mod_continuum
use mod_fit
use mod_uncertainties

implicit none
integer :: I, spectrumlength, nlines, linearraypos, totallines, startpos, endpos
real :: startwlen, endwlen
character (len=512) :: spectrumfile,stronglinelistfile,deeplinelistfile,skylinelistfile,outputdirectory

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

logical :: normalise=.false. !false means spectrum normalised to whatever H beta is detected, true means spectrum normalised to user specified value
logical :: resolution_estimated=.false. !true means user specified a value, false means estimate from sampling
logical :: subtractsky=.false. !attempt to fit night sky emission lines
logical :: file_exists

c=299792.458 !km/s
!default values in absence of user specificed guess
redshiftguess=0.0 !km/s
rtol1=0.d0 !variation allowed in resolution on first pass.  determined later, either from user input, or to be equal to resolution guess.
rtol2=500. !second pass
vtol1=0.003 !variation allowed in velocity (expressed as redshift) on first pass. 0.003 = 900 km/s
vtol2=0.0002 !second pass. 0.0002 = 60 km/s

stronglinelistfile="/etc/alfa/strong_optical"
deeplinelistfile="/etc/alfa/deep_full"
skylinelistfile="/etc/alfa/skylines"

outputdirectory="."

! start

print *,"ALFA, the Automated Line Fitting Algorithm"
print *,"------------------------------------------"

print *
print *,gettime(),": starting code"

! random seed

call init_random_seed()

! read command line

narg = IARGC() !count input arguments
nargused = 0 !count options specified
if (narg .eq. 0) then
  print *,"Usage: alfa [options] [file]"
  print *,"  [file] is an ascii file with columns for wavelength and flux"
  print *,"  see the man page or online documentation for details of the options"
  stop
endif

include "commandline.f95"

! convert from velocity to redshift

redshiftguess=1.+(redshiftguess/c)

print *,gettime(),": command line: ",trim(commandline)

! read in spectrum to fit and line list

print *,gettime(),": reading in spectrum ",trim(spectrumfile)
call readspectrum(spectrumfile, realspec, spectrumlength, fittedspectrum)

!read in catalogues

print *,gettime(),": reading in line catalogues ",trim(skylinelistfile),", ",trim(stronglinelistfile),", ",trim(deeplinelistfile)
call readlinelist(skylinelistfile, skylines_catalogue, nlines,minval(realspec%wavelength),maxval(realspec%wavelength))
call readlinelist(stronglinelistfile, stronglines_catalogue, nlines,minval(realspec%wavelength),maxval(realspec%wavelength))
call readlinelist(deeplinelistfile, deeplines_catalogue, nlines,minval(realspec%wavelength),maxval(realspec%wavelength))

! then subtract the continuum

print *,gettime(),": fitting continuum"
call fit_continuum(realspec,spectrumlength, continuum)

! now do the fitting
! first get guesses for the redshift and resolution

if (.not. resolution_estimated) then
  ! estimate resolution assuming nyquist sampling
  resolutionguess=2*realspec(2)%wavelength/(realspec(3)%wavelength-realspec(1)%wavelength)
  print *,gettime(),": estimated spectrograph resolution assuming Nyquist sampling: ",resolutionguess
else
  print *,gettime(),": estimated spectrograph resolution from user input: ",resolutionguess
endif

if (rtol1 .eq. 0.d0) then
  rtol1=0.9*resolutionguess ! user didn't specify a value, default behaviour is to allow resolution to vary between 0.1 and 1.9x the initial guess on the first pass
endif

print "(X,A,A,F8.1,A,F7.1)",gettime(),": initial guesses for velocity and resolution: ",c*(redshiftguess-1),"km/s, R=",resolutionguess

! first, subtract sky spectrum if requested. do in chunks of 400 units. no overlap necessary because velocity is zero

allocate(skyspectrum(spectrumlength))
skyspectrum%wavelength=realspec%wavelength
skyspectrum%flux=0.d0

if (subtractsky) then
  print *,gettime(),": fitting sky emission"

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
      call fit(spectrumchunk, 1., resolutionguess, skylines_section, 0., rtol2)
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
  print *,gettime(),": Warning: no reference lines detected, using default guesses for velocity and resolution"
  redshiftguess_overall=1.d0
else
  print *,gettime(),": estimating resolution and velocity using ",nlines," lines"
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

  call fit(stronglines, redshiftguess, resolutionguess, fittedlines, vtol1, rtol1)

  print *,gettime(),": estimated redshift and resolution: ",c*(fittedlines(1)%redshift-1),fittedlines(1)%resolution
  redshiftguess_overall = fittedlines(1)%redshift ! when fitting chunks, use this redshift to get lines in the right range from the catalogue. if velocity from each chunk is used, then there's a chance that a line could be missed or double counted due to variations in the calculated velocity between chunks.
  redshiftguess=fittedlines(1)%redshift
  resolutionguess=fittedlines(1)%resolution

  deallocate(stronglines)

endif

! then again in chunks with tighter tolerance

linearraypos=1

!get total number of lines and an array to put them all in

call selectlines(deeplines_catalogue, realspec(1)%wavelength/redshiftguess_overall, realspec(size(realspec))%wavelength/redshiftguess_overall, fittedlines, totallines)

print *, gettime(), ": fitting full spectrum with ",totallines," lines"

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
    print "(' ',A,A,F7.1,A,F7.1,A,I3,A)",gettime(),": fitting from ",spectrumchunk(1)%wavelength," to ",spectrumchunk(size(spectrumchunk))%wavelength," with ",nlines," lines"
    call fit(spectrumchunk, redshiftguess, resolutionguess, fittedlines_section, vtol2, rtol2)
    !use redshift and resolution from this chunk as initial values for next chunk
    redshiftguess=fittedlines_section(1)%redshift
    resolutionguess=fittedlines_section(1)%resolution
    !copy results back
    fittedlines(linearraypos:linearraypos+nlines-1)=fittedlines_section
    linearraypos=linearraypos+nlines
    !report any errors
    if (maxval(fittedlines_section%wavelength*fittedlines_section%redshift) .gt. maxval(spectrumchunk%wavelength) .or. minval(fittedlines_section%wavelength*fittedlines_section%redshift) .lt. minval(spectrumchunk%wavelength)) then
      print *,"              Warning: some lines ended up outside the fitting region."
      print *,"              Try reducing the value of vtol2, which is currently set to ",c*vtol2
    endif
  endif
  deallocate(spectrumchunk)
enddo

!make the fitted spectrum

call makespectrum(fittedlines, fittedspectrum)

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

open(100,file=trim(outputdirectory)//trim(spectrumfile)//"_fit")

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

open(100,file=trim(outputdirectory)//trim(spectrumfile)//"_lines.tex")
open(101,file=trim(outputdirectory)//trim(spectrumfile)//"_lines")
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

print *,gettime(),": all done"
print *

end program alfa
