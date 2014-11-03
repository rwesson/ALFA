module mod_fit
use mod_routines
use mod_types

contains
subroutine fit(realspec, referencelinelist, population, synthspec, rms)

implicit none

type(linelist) :: referencelinelist
type(linelist), dimension(:),allocatable :: population
type(linelist), dimension(:),allocatable ::  breed
type(spectrum), dimension(:,:), allocatable :: synthspec
type(spectrum), dimension(:) :: realspec
integer :: popsize, i, wlength, spectrumlength, generations, lineid, loc1, loc2, nlines, gencount, popnumber
real, dimension(:), allocatable :: rms
real :: random, pressure

!initialisation

  popsize=50
  generations=1000
  pressure=0.1 !pressure * popsize needs to be an integer

  nlines=size(referencelinelist%wavelength)
  spectrumlength=size(realspec%wavelength)

  call init_random_seed()

!allocate arrays

  allocate(synthspec(spectrumlength,popsize))
  allocate(rms(popsize))
  allocate(breed(int(popsize*pressure)))
  allocate(population(popsize))

  do i=1,int(popsize*pressure)
    allocate (breed(i)%peak(nlines))
    allocate (breed(i)%wavelength(nlines))
    allocate (breed(i)%uncertainty(nlines))
  end do
  do i=1,popsize
    allocate (population(i)%peak(nlines))
    allocate (population(i)%wavelength(nlines))
    allocate (population(i)%uncertainty(nlines))
  end do

! now create population of synthetic spectra
! todo, make sure no lines are included which are outside the wavelength range
! of the observations

do i=1,popsize
  synthspec(:,i)%wavelength=realspec%wavelength
end do

do popnumber=1,popsize
  population(popnumber)%wavelength = referencelinelist%wavelength
  population(popnumber)%peak=10.0
  population(popnumber)%resolution=3000.
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
        if (abs(population(popnumber)%redshift*population(popnumber)%wavelength(lineid) - synthspec(wlength,popnumber)%wavelength) .lt.  (5*population(popnumber)%wavelength(lineid)/population(popnumber)%resolution)) then
          synthspec(wlength,popnumber)%flux = synthspec(wlength,popnumber)%flux + &
          &gaussian(synthspec(wlength,popnumber)%wavelength,&
          &population(popnumber)%peak(lineid),population(popnumber)%redshift*population(popnumber)%wavelength(lineid), population(popnumber)%wavelength(lineid)/population(popnumber)%resolution)
        endif
      end do
    end do
  
    !now calculate RMS for the "models"
  
    rms(popnumber)=sum((synthspec(:,popnumber)%flux-realspec(:)%flux)**2)/spectrumlength
  
  end do
  
    !next, cream off the well performing models - put the population
    !member with the lowest RMS into the breed array, replace the RMS with
    !something very high so that it doesn't get copied twice, repeat until
    !a fraction equal to the pressure factor have been selected
  
    do i=1,int(popsize*pressure) 
      breed(i) = population(minloc(rms,1))
      rms(minloc(rms,1))=1.e10
    end do
  
  !then, "breed" pairs
  !random approach will mean that some models have no offspring while
  !others might
  !have many.  Alternative approach could be to breed all adjacent pairs
  !so that
  !every model generates one offspring.
  
    if (gencount .ne. generations) then
      do i=1,popsize 
        call random_number(random)
        loc1=int(popsize*random*pressure)+1
        call random_number(random)
        loc2=int(popsize*random*pressure)+1 
        population(i)%peak=(breed(loc1)%peak + breed(loc2)%peak)/2.0
        population(i)%resolution=(breed(loc1)%resolution + breed(loc2)%resolution)/2.0
        population(i)%redshift=(breed(loc1)%redshift + breed(loc2)%redshift)/2.0
      end do
      !then, "mutate"
      do popnumber=1,popsize ! mutation of spectral resolution
        population(popnumber)%resolution = population(popnumber)%resolution * mutation()
        if (population(popnumber)%resolution .gt. 20000.) then !this condition may not always be necessary
          population(popnumber)%resolution = 20000.
        endif
        population(popnumber)%redshift = population(popnumber)%redshift * ((999.+mutation())/1000.)
        do lineid=1,nlines !mutation of line fluxes
          population(popnumber)%peak(lineid) = population(popnumber)%peak(lineid) * mutation()
        enddo
      enddo
    endif
  
    if (mod(gencount,generations/10) .eq.0 .or. gencount.eq.1) then
      print *,gettime()," : completed ",100*gencount/generations, "%", population(minloc(rms,1))%resolution,minval(rms,1), 3.e5*(population(minloc(rms,1))%redshift-1)
      do i=1,spectrumlength
        write (101,*) synthspec(i,minloc(rms,1))%wavelength,synthspec(i,minloc(rms,1))%flux
      enddo
      write (101,*)
    endif
  
  !  if (mod(gencount,10) .eq. 0) then
  !    do i=1,spectrumlength
  !      print
  !      *,synthspec(i,minloc(rms,1))%wavelength,synthspec(i,minloc(rms,1))%flux,realspec(i)%flux
  !    end do
  !    print*
  !  endif
  

  
  end do
 
end subroutine fit
end module mod_fit 
