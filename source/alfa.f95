program alfa

use mod_readfiles
use mod_routines
use mod_types
use mod_quicksort
use mod_continuum

implicit none
integer :: I, lineid, popnumber, gencount, spectrumlength, nlines
integer :: loc1, loc2
integer :: wlength
real :: random
character*512 :: spectrumfile,linelistfile

integer :: popsize, generations

real :: pressure ! the fraction of candidate spectra which breed in every generation

type(linelist) :: referencelinelist
type(linelist), dimension(:),allocatable :: population
type(linelist), dimension(:),allocatable ::  breed
type(spectrum), dimension(:,:), allocatable :: synthspec
type(spectrum), dimension(:), allocatable :: realspec
type(spectrum), dimension(:), allocatable :: continuum

real, dimension(:), allocatable :: rms

!temp XXXX
  open (101,file="intermediate",status="replace")

! initialise stuff for genetics

  popsize=100
  generations=1000

  pressure=0.05 !pressure * popsize needs to be an integer

  allocate(rms(popsize))

! random seed

  call init_random_seed()

!read files

! read in spectrum to fit and line list

  call get_command_argument(1,spectrumfile)
  call get_command_argument(2,linelistfile)

  call readfiles(spectrumfile,linelistfile,synthspec,realspec,referencelinelist,popsize,spectrumlength, nlines)

! then subtract the continuum

  call fit_continuum(realspec,spectrumlength, continuum)

!allocate some more arrays

  allocate(breed(int(popsize*pressure)))
  allocate(population(popsize))

  do i=1,int(popsize*pressure)
    allocate (breed(i)%peak(nlines))
    allocate (breed(i)%wavelength(nlines))
  end do
  do i=1,popsize
    allocate (population(i)%peak(nlines))
    allocate (population(i)%wavelength(nlines))
  end do

! now create population of synthetic spectra
! todo, make sure no lines are included which are outside the wavelength range
! of the observations

do popnumber=1,popsize
  population(popnumber)%wavelength = referencelinelist%wavelength
  population(popnumber)%peak=10.0
  population(popnumber)%width=7.5
  population(popnumber)%redshift=1.0
end do

do gencount=1,generations

!reset stuff to zero before doing calculations

synthspec%flux=0.D0
rms=0.D0

do popnumber=1,popsize

!calculate synthetic spectra - reset to 0 before synthesizing
!line fluxes are calculated within 7 sigma of mean

  do wlength=1,spectrumlength 
    do lineid=1,nlines
      if (abs(population(popnumber)%redshift*population(popnumber)%wavelength(lineid) - synthspec(wlength,popnumber)%wavelength) .lt. (7*population(popnumber)%width)) then
        synthspec(wlength,popnumber)%flux = synthspec(wlength,popnumber)%flux + &
        &gaussian(synthspec(wlength,popnumber)%wavelength,&
        &population(popnumber)%peak(lineid),population(popnumber)%redshift*population(popnumber)%wavelength(lineid), population(popnumber)%width)
      endif
    end do
  end do

  !now calculate RMS for the "models"

  rms(popnumber)=sum((synthspec(:,popnumber)%flux-realspec(:)%flux)**2)

end do

  !next, cream off the well performing models - put the population member with the lowest RMS into the breed array, replace the RMS with something very high so that it doesn't get copied twice, repeat until a fraction equal to the pressure factor have been selected

  do i=1,int(popsize*pressure) 
    breed(i) = population(minloc(rms,1))
    rms(minloc(rms,1))=1.e10
  end do

!then, "breed" pairs
!random approach will mean that some models have no offspring while others might
!have many.  Alternative approach could be to breed all adjacent pairs so that
!every model generates one offspring.

  if (gencount .ne. generations) then
    do i=1,popsize 
      call random_number(random)
      loc1=int(popsize*random*pressure)+1
      call random_number(random)
      loc2=int(popsize*random*pressure)+1 
      population(i)%peak=(breed(loc1)%peak + breed(loc2)%peak)/2.0
      population(i)%width=(breed(loc1)%width + breed(loc2)%width)/2.0
      population(i)%redshift=(breed(loc1)%redshift + breed(loc2)%redshift)/2.0
    end do
    !then, "mutate"
    do popnumber=1,popsize ! mutation of line width
      population(popnumber)%width = population(popnumber)%width * mutation()
if (population(popnumber)%width .lt. 0.5) then
  population(popnumber)%width = 0.5
endif
      population(popnumber)%redshift = population(popnumber)%redshift * ((99999.+mutation())/100000.)
      do lineid=1,nlines !mutation of line fluxes
        population(popnumber)%peak(lineid) = population(popnumber)%peak(lineid) * mutation()
      enddo
    enddo
  endif

  if (mod(gencount,generations/10) .eq.0 .or. gencount.eq.1) then
    print *,gettime()," : completed ",100*gencount/generations, "%", population(minloc(rms,1))%width,minval(rms,1), 3.e5*(population(minloc(rms,1))%redshift-1)
    do i=1,spectrumlength
      write (101,*) synthspec(i,minloc(rms,1))%wavelength,synthspec(i,minloc(rms,1))%flux
    enddo
    write (101,*)
  endif

!  if (mod(gencount,10) .eq. 0) then
!    do i=1,spectrumlength
!      print *,synthspec(i,minloc(rms,1))%wavelength,synthspec(i,minloc(rms,1))%flux,realspec(i)%flux
!    end do
!    print*
!  endif

!todo, work out why rms doesn't keep decreasing

end do
!write out line fluxes of best fitting spectrum

do i=1,nlines
  print *,population(1)%wavelength(i),gaussianflux(population(minloc(rms,1))%peak(i),population(minloc(rms,1))%width), population(minloc(rms,1))%peak(i),population(minloc(rms,1))%width
end do

open(100,file="outputfit")

write (100,*) """wavelength""  ""fitted spectrum""  ""cont-subbed orig"" ""continuum""  ""residuals"""
do i=1,spectrumlength
  write(100,*) synthspec(i,minloc(rms,1))%wavelength,synthspec(i,minloc(rms,1))%flux, realspec(i)%flux, continuum(i)%flux, realspec(i)%flux - synthspec(i,minloc(rms,1))%flux
end do

close(100)
close(101) !temp XXXX

end program alfa 
