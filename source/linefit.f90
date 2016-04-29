module mod_fit
use mod_routines
use mod_types

contains
subroutine fit(inputspectrum, redshiftguess, resolutionguess, fittedlines, redshifttolerance, resolutiontolerance, generations, popsize, pressure)

implicit none

type(linelist), dimension(:), allocatable :: fittedlines
type(linelist), dimension(:,:), allocatable :: population
type(linelist), dimension(:,:), allocatable ::  breed
type(spectrum), dimension(:,:), allocatable :: synthspec
type(spectrum), dimension(:) :: inputspectrum
integer :: popsize, i, spectrumlength, lineid, loc1, loc2, nlines, gencount, generations, popnumber
real, dimension(:), allocatable :: sumsquares
real :: random, pressure
real :: resolutionguess, redshiftguess, redshifttolerance, resolutiontolerance
real :: scalefactor

#ifdef CO
  print *,"subroutine: fit"
#endif

!initialisation

  nlines=size(fittedlines%wavelength)
  spectrumlength=size(inputspectrum%wavelength)

  call init_random_seed()

  scalefactor=1.d0
  if (maxval(inputspectrum%flux) .lt. 0.01) then
    !hacky fix, something goes wrong with very small numbers like fluxes in units of erg/cm2/s/A
    scalefactor = 1./maxval(inputspectrum%flux)
  else
  endif
  inputspectrum%flux = inputspectrum%flux * scalefactor

!allocate arrays

  allocate(synthspec(spectrumlength,popsize))
  allocate(sumsquares(popsize))
  allocate(breed(int(popsize*pressure),nlines))
  allocate(population(popsize,nlines))

! now create population of synthetic spectra
! todo, make sure no lines are included which are outside the wavelength range
! of the observations

  do i=1,popsize
    synthspec(:,i)%wavelength=inputspectrum%wavelength
  enddo

  do popnumber=1,popsize
    population(popnumber,:)%wavelength = fittedlines%wavelength
    population(popnumber,:)%peak=fittedlines%peak
    population(popnumber,:)%resolution=resolutionguess
    population(popnumber,:)%redshift=redshiftguess
    population(popnumber,:)%linedata=fittedlines%linedata
  enddo

! evolve

  do gencount=1,generations

!reset stuff to zero before doing calculations

    synthspec%flux=0.D0
    sumsquares=0.D0

    do popnumber=1,popsize

!calculate synthetic spectra

    call makespectrum(population(popnumber,:),synthspec(:,popnumber))

!now calculate sum of squares for the "models"

    sumsquares(popnumber)=sum((synthspec(:,popnumber)%flux-inputspectrum(:)%flux)**2)

    enddo

    !if that was the last generation then exit

    if (gencount .eq. generations) then
      !copy fit results into arrays to return
      fittedlines = population(minloc(sumsquares,1),:)
      !scale
      fittedlines%peak = fittedlines%peak / scalefactor
      !deallocate arrays
      deallocate(synthspec)
      deallocate(sumsquares)
      deallocate(breed)
      deallocate(population)
      exit
    endif

    !next, cream off the well performing models - put the population
    !member with the lowest sum of squares into the breed array, replace the sum of squares with
    !something very high so that it doesn't get copied twice, repeat until
    !a fraction equal to the pressure factor have been selected
    !best fitting member of generation is put into slot 1 of population array, to be retained unaltered in the next generation

    population(1,:)=population(minloc(sumsquares,1),:)
    sumsquares(minloc(sumsquares,1))=1.e30

    do i=1,int(popsize*pressure)
      breed(i,:) = population(minloc(sumsquares,1),:)
      sumsquares(minloc(sumsquares,1))=1.e20
    enddo

  !then, "breed" pairs
  !random approach will mean that some models have no offspring while others might have many.
  !Alternative approach could be to breed all adjacent pairs so that every model generates one offspring.

    do i=2,popsize
      call random_number(random)
      loc1=int(popsize*random*pressure)+1
      call random_number(random)
      loc2=int(popsize*random*pressure)+1
      population(i,:)%peak=(breed(loc1,:)%peak + breed(loc2,:)%peak)/2.0
      population(i,:)%resolution=(breed(loc1,:)%resolution + breed(loc2,:)%resolution)/2.0
      population(i,:)%redshift=(breed(loc1,:)%redshift + breed(loc2,:)%redshift)/2.0
    enddo
    !then, mutate
    do popnumber=2,popsize ! mutation of spectral resolution
      population(popnumber,:)%resolution = population(popnumber,:)%resolution * mutation()
      if (abs(population(popnumber,1)%resolution-resolutionguess) .gt. resolutiontolerance) then
        population(popnumber,:)%resolution = resolutionguess
      endif
!mutation returns a value between 0 and 2.
!useful values of redshift change are between -vtol/c and +vtol/c
!(mutation-1)*vtol/c is the useful variation
      population(popnumber,:)%redshift = population(popnumber,:)%redshift + ((mutation()-1.)*redshifttolerance)
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
