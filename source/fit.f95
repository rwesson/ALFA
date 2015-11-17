module mod_fit
use mod_routines
use mod_types

contains
subroutine fit(inputspectrum, referencelinelist, redshiftguess, resolutionguess, fittedlines, redshifttolerance, resolutiontolerance)

implicit none

type(linelist), dimension(:), allocatable :: referencelinelist, fittedlines
type(linelist), dimension(:,:), allocatable :: population
type(linelist), dimension(:,:), allocatable ::  breed
type(spectrum), dimension(:,:), allocatable :: synthspec
type(spectrum), dimension(:) :: inputspectrum
integer :: popsize, i, spectrumlength, lineid, loc1, loc2, nlines, gencount, generations, popnumber
real, dimension(:), allocatable :: rms
real :: random, pressure
real :: resolutionguess, redshiftguess, redshifttolerance, resolutiontolerance

!initialisation

  popsize=50
  pressure=0.1 !pressure * popsize needs to be an integer
  generations=500

  nlines=size(referencelinelist%wavelength)
  spectrumlength=size(inputspectrum%wavelength)

  call init_random_seed()

!allocate arrays

  allocate(synthspec(spectrumlength,popsize))
  allocate(rms(popsize))
  allocate(breed(int(popsize*pressure),nlines))
  allocate(population(popsize,nlines))

! now create population of synthetic spectra
! todo, make sure no lines are included which are outside the wavelength range
! of the observations

  do i=1,popsize
    synthspec(:,i)%wavelength=inputspectrum%wavelength
  enddo

  do popnumber=1,popsize
    population(popnumber,:)%wavelength = referencelinelist%wavelength
    population(popnumber,:)%peak=referencelinelist%peak
    population(popnumber,:)%resolution=resolutionguess
    population(popnumber,:)%redshift=redshiftguess
    population(popnumber,:)%linedata=referencelinelist%linedata
  enddo

! evolve

  do gencount=1,generations

!reset stuff to zero before doing calculations

    synthspec%flux=0.D0
    rms=0.D0

    do popnumber=1,popsize

!calculate synthetic spectra

    call makespectrum(population(popnumber,:),synthspec(:,popnumber))

!now calculate RMS for the "models"

    rms(popnumber)=sum((synthspec(:,popnumber)%flux-inputspectrum(:)%flux)**2)/spectrumlength

    enddo

    !if that was the last generation then exit

    if (gencount .eq. generations) then
      !copy fit results into arrays to return
      fittedlines = population(minloc(rms,1),:)
      !deallocate arrays
      deallocate(synthspec)
      deallocate(rms)
      deallocate(breed)
      deallocate(population)
      exit
    endif

    !next, cream off the well performing models - put the population
    !member with the lowest RMS into the breed array, replace the RMS with
    !something very high so that it doesn't get copied twice, repeat until
    !a fraction equal to the pressure factor have been selected

    do i=1,int(popsize*pressure)
      breed(i,:) = population(minloc(rms,1),:)
      rms(minloc(rms,1))=1.e10
    enddo

  !then, "breed" pairs
  !random approach will mean that some models have no offspring while others might have many.
  !Alternative approach could be to breed all adjacent pairs so that every model generates one offspring.

    do i=1,popsize
      call random_number(random)
      loc1=int(popsize*random*pressure)+1
      call random_number(random)
      loc2=int(popsize*random*pressure)+1
      population(i,:)%peak=(breed(loc1,:)%peak + breed(loc2,:)%peak)/2.0
      population(i,:)%resolution=(breed(loc1,:)%resolution + breed(loc2,:)%resolution)/2.0
      population(i,:)%redshift=(breed(loc1,:)%redshift + breed(loc2,:)%redshift)/2.0
    enddo
    !then, mutate
    do popnumber=1,popsize ! mutation of spectral resolution
      population(popnumber,:)%resolution = population(popnumber,:)%resolution * mutation()
      if (abs(population(popnumber,1)%resolution-resolutionguess) .gt. resolutiontolerance) then
        population(popnumber,:)%resolution = resolutionguess
      endif
      population(popnumber,:)%redshift = population(popnumber,:)%redshift * ((9999.+mutation())/10000.)
      if (abs(population(popnumber,1)%redshift-redshiftguess) .gt. redshifttolerance) then
        population(popnumber,:)%redshift = redshiftguess
      endif
      do lineid=1,nlines !mutation of line fluxes
        population(popnumber,lineid)%peak = population(popnumber,lineid)%peak * mutation()
      enddo
    enddo

  enddo

end subroutine fit
end module mod_fit
