program alfa

use mod_readfiles
use mod_routines
use mod_types
use mod_quicksort
use mod_continuum
use mod_fit
use mod_uncertainties

implicit none
integer :: I, spectrumlength, nlines
character (len=512) :: spectrumfile,linelistfile

type(linelist) :: referencelinelist
type(linelist), dimension(:),allocatable :: population
type(spectrum), dimension(:,:), allocatable :: synthspec
type(spectrum), dimension(:), allocatable :: realspec
type(spectrum), dimension(:), allocatable :: continuum

real, dimension(:), allocatable :: rms

!temp XXXX

open (101,file="intermediate",status="replace")

! random seed

call init_random_seed()

! read in spectrum to fit and line list

call get_command_argument(1,spectrumfile)
call get_command_argument(2,linelistfile)

call readfiles(spectrumfile,linelistfile,realspec,referencelinelist,spectrumlength, nlines)

! then subtract the continuum

call fit_continuum(realspec,spectrumlength, continuum)

! now do the fitting

call fit(realspec, referencelinelist, population, synthspec, rms)

! calculate the uncertainties

call get_uncertainties(synthspec, realspec, population, rms)

!write out line fluxes of best fitting spectrum

open(100,file="outputlines")
write(100,*) """observed wavelength""  ""rest wavelength""  ""flux""  ""uncertainty"""
do i=1,nlines
  if (population(minloc(rms,1))%uncertainty(i) .gt. 1.0) then
    write (100,"(F7.2,2X,F7.2,2X,ES12.3,2X,ES12.3)") population(1)%wavelength(i)*population(1)%redshift,population(1)%wavelength(i),gaussianflux(population(minloc(rms,1))%peak(i),(population(minloc(rms,1))%wavelength(i)/population(minloc(rms,1))%resolution)), gaussianflux(population(minloc(rms,1))%peak(i),(population(minloc(rms,1))%wavelength(i)/population(minloc(rms,1))%resolution))/population(minloc(rms,1))%uncertainty(i)!, population(minloc(rms,1))%peak(i),population(minloc(rms,1))%width
  end if
end do
close(100)

! write out fit

open(100,file="outputfit")

write (100,*) """wavelength""  ""fitted spectrum""  ""cont-subbed orig"" ""continuum""  ""residuals"""
do i=1,spectrumlength
  write(100,*) synthspec(i,minloc(rms,1))%wavelength,synthspec(i,minloc(rms,1))%flux, realspec(i)%flux, continuum(i)%flux, realspec(i)%flux - synthspec(i,minloc(rms,1))%flux
end do

close(100)
close(101) !temp XXXX

end program alfa 
