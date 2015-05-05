module mod_fit
use mod_routines
use mod_types

contains
subroutine fit(realspec, referencelinelist, population, synthspec, rms, redshiftguess, resolutionguess)

implicit none

type(linelist) :: referencelinelist
type(linelist), dimension(:),allocatable :: population
type(linelist), dimension(:),allocatable ::  breed
type(spectrum), dimension(:,:), allocatable :: synthspec
type(spectrum), dimension(:) :: realspec
integer :: popsize, i, specpoint, spectrumlength, lineid, loc1, loc2, nlines, gencount, popnumber
real, dimension(:), allocatable :: rms
real :: random, pressure, convergence, oldrms
real :: resolutionguess, redshiftguess
real :: tmpvar !XXX
!initialisation

  popsize=50
  pressure=0.1 !pressure * popsize needs to be an integer
  convergence=0.0 !new rms / old rms

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
  population(popnumber)%peak=referencelinelist%peak
  population(popnumber)%resolution=resolutionguess
  population(popnumber)%redshift=redshiftguess
end do

! iterate until rms changes by less than 1 percent

gencount=1

do while (gencount .lt. 1001)
!do while (convergence .lt. 0.99999)

  if (gencount.eq.1) then
   oldrms=1.e30
  else
    oldrms=minval(rms,1)
  endif

!reset stuff to zero before doing calculations

  synthspec%flux=0.D0
  rms=0.D0

  do popnumber=1,popsize
  
  !calculate synthetic spectra - reset to 0 before synthesizing
  !line fluxes are calculated within 5 sigma of mean
  
    do lineid=1,nlines
      where (abs(population(popnumber)%redshift*population(popnumber)%wavelength(lineid) - synthspec(:,popnumber)%wavelength) .lt. (5*population(popnumber)%wavelength(lineid)/population(popnumber)%resolution))
          synthspec(:,popnumber)%flux = synthspec(:,popnumber)%flux + &
          &population(popnumber)%peak(lineid)*exp((-(synthspec(:,popnumber)%wavelength-population(popnumber)%redshift*population(popnumber)%wavelength(lineid))**2)/(2*(population(popnumber)%wavelength(lineid)/population(popnumber)%resolution)**2))

      end where
    end do
  
    !now calculate RMS for the "models"
  
    rms(popnumber)=sum((synthspec(:,popnumber)%flux-realspec(:)%flux)**2)/spectrumlength
  
  end do
  
    !next, cream off the well performing models - put the population
    !member with the lowest RMS into the breed array, replace the RMS with
    !something very high so that it doesn't get copied twice, repeat until
    !a fraction equal to the pressure factor have been selected
tmpvar = maxval(rms,1) !XXX
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
      if (population(popnumber)%resolution .lt. 3000.) then !this condition may not always be necessary
        population(popnumber)%resolution = 7000.
      endif
      if (population(popnumber)%resolution .gt. 12000.) then !this condition may not always be necessary
        population(popnumber)%resolution = 9000.
      endif
      population(popnumber)%redshift = population(popnumber)%redshift * ((9999.+mutation())/10000.)
      if (abs(population(popnumber)%redshift) .gt. 1.0) population(popnumber)%redshift = 0.9996
      do lineid=1,nlines !mutation of line fluxes
        population(popnumber)%peak(lineid) = population(popnumber)%peak(lineid) * mutation()
      enddo
    enddo
  
  if (mod(gencount,100) .eq.0 .or. gencount.eq.1) then
      print "(X,A,A,i5,A,4(X,F12.3))",gettime()," : ",gencount, " generations  ", population(minloc(rms,1))%resolution, 3.e5*(population(minloc(rms,1))%redshift-1), minval(rms,1), tmpvar
! temporary
      do i=1,spectrumlength
        write (999,*) synthspec(i,minloc(rms,1))%wavelength,synthspec(i,minloc(rms,1))%flux,synthspec(i,maxloc(rms,1))%flux
      enddo
      write (999,*)
! end of temporariness
    endif
  
  !  if (mod(gencount,10) .eq. 0) then
  !    do i=1,spectrumlength
  !      print
  !      *,synthspec(i,minloc(rms,1))%wavelength,synthspec(i,minloc(rms,1))%flux,realspec(i)%flux
  !    end do
  !    print*
  !  endif
  
  gencount=gencount+1
  convergence=minval(rms,1)/oldrms
  if (convergence .gt. 1.0) then
    convergence = 1./convergence
  endif
!print *,convergence,gencount, oldrms
  end do
      print "(X,A,A,i5,A,4(X,F12.3))",gettime()," : ",gencount, " generations  ", population(minloc(rms,1))%resolution, 3.e5*(population(minloc(rms,1))%redshift-1), minval(rms,1), tmpvar
 
end subroutine fit
end module mod_fit 
