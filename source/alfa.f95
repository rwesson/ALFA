program alfa

use mod_readfiles
use mod_routines
use mod_types
use mod_continuum
use mod_fit
use mod_uncertainties

implicit none
integer :: I, spectrumlength, nlines
character (len=512) :: spectrumfile,linelistfile

type(linelist) :: referencelinelist, fittedlines
type(spectrum), dimension(:), allocatable :: realspec, fittedspectrum
type(spectrum), dimension(:), allocatable :: continuum
character(len=85), dimension(:), allocatable :: linedata

real :: normalisation

CHARACTER*2048, DIMENSION(:), allocatable :: options
CHARACTER*2048 :: commandline
integer :: narg

logical :: normalise

real :: redshiftguess, resolutionguess, tolerance

!temp XXXX

open (999,file="intermediate",status="replace")

!set defaults

normalise = .false.

! start

print *,"ALFA, the Automated Line Fitting Algorithm"
print *,"------------------------------------------"

print *
print *,gettime(),": starting code"

! random seed

call init_random_seed()

! read command line

narg = IARGC() !count input arguments
if (narg .lt. 2) then
  print *,gettime(),": Error : files to analyse not specified"
  stop
endif

call get_command(commandline)
ALLOCATE (options(Narg))
if (narg .gt. 2) then
  do i=1,Narg-2
    call get_command_argument(i,options(i))
    if (options(i) .eq. "-n") then
      normalise = .true.
    endif
  enddo
endif

print *,gettime(),": command line: ",trim(commandline)

call get_command_argument(narg-1,spectrumfile)
call get_command_argument(narg,linelistfile)

! read in spectrum to fit and line list

call readspectrum(spectrumfile, realspec, spectrumlength, fittedspectrum)
call readlinelist(linelistfile, referencelinelist, nlines, linedata, fittedlines, realspec)
!call readfiles(spectrumfile,linelistfile,realspec,referencelinelist,spectrumlength, nlines, linedata, fittedspectrum, fittedlines)

! then subtract the continuum

print *,gettime(),": fitting continuum"
call fit_continuum(realspec,spectrumlength, continuum)

! now do the fitting
! first get guesses for the redshift and resolution

redshiftguess=1.0001
resolutionguess=4800.
tolerance=1.0

print *,gettime(),": fitting ",nlines," lines"
print *
print *,"Best fitting model parameters:       Resolution    Redshift    RMS min      RMS max"
call fit(realspec, referencelinelist, redshiftguess, resolutionguess, fittedspectrum, fittedlines, tolerance)

! then again with tighter tolerance

redshiftguess=fittedlines%redshift
resolutionguess=fittedlines%resolution
tolerance=0.05

print *,gettime(),": fitting ",nlines," lines"
print *
print *,"Best fitting model parameters:       Resolution    Redshift    RMS min      RMS max"
call fit(realspec, referencelinelist, redshiftguess, resolutionguess, fittedspectrum, fittedlines, tolerance)

! calculate the uncertainties

print *
print *,gettime(),": estimating uncertainties"
call get_uncertainties(fittedspectrum, realspec, fittedlines)

!write out line fluxes of best fitting spectrum

print *,gettime(),": writing output files ",trim(spectrumfile),"_lines.tex and ",trim(spectrumfile),"_fit"

!normalise Hb to 100 if requested

if (normalise) then
print *,gettime(),": normalising to Hb=100"
  do i=1,nlines
    if (fittedlines%wavelength(i) .eq. 4861.33) then
      normalisation = 100./gaussianflux(fittedlines%peak(i),(fittedlines%wavelength(i)/fittedlines%resolution))
      exit
    endif
  enddo
else
  normalisation = 1.D0
endif

open(100,file=trim(spectrumfile)//"_lines.tex")
open(101,file=trim(spectrumfile)//"_neat_input")
write(100,*) "Observed wavelength & Rest wavelength & Flux & Uncertainty & Ion & Multiplet & Lower term & Upper term & g_1 & g_2 \\"
write(101,*) "#Obs. wlen.  Rest wlen.   Flux   Uncertainty"
do i=1,nlines
  if (fittedlines%uncertainty(i) .gt. 3.0) then
    write (100,"(F7.2,' & ',F7.2,' & ',F12.3,' & ',F12.3,A85)") fittedlines%wavelength(i)*fittedlines%redshift,fittedlines%wavelength(i),normalisation*gaussianflux(fittedlines%peak(i),(fittedlines%wavelength(i)/fittedlines%resolution)), normalisation*gaussianflux(fittedlines%peak(i),(fittedlines%wavelength(i)/fittedlines%resolution))/fittedlines%uncertainty(i), linedata(i)
    write (101,"(F7.2, 3(F12.3))") fittedlines%wavelength(i),normalisation*gaussianflux(fittedlines%peak(i),(fittedlines%wavelength(i)/fittedlines%resolution)),normalisation*gaussianflux(fittedlines%peak(i),(fittedlines%wavelength(i)/fittedlines%resolution))/fittedlines%uncertainty(i), fittedlines%resolution
  end if
end do
close(101)
close(100)

! write out fit

open(100,file=trim(spectrumfile)//"_fit")

write (100,*) """wavelength""  ""fitted spectrum""  ""cont-subbed orig"" ""continuum""  ""residuals"""
do i=1,spectrumlength
  write(100,*) fittedspectrum(i)%wavelength,fittedspectrum(i)%flux, realspec(i)%flux, continuum(i)%flux, realspec(i)%flux - fittedspectrum(i)%flux
end do

close(100)
close(999) !temp XXXX

print *,gettime(),": all done"

end program alfa
