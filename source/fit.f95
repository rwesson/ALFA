module mod_fit
use mod_routines
use mod_types

contains
subroutine fit(inputspectrum, referencelinelist, redshiftguess, resolutionguess, fittedspectrum, fittedlines, redshifttolerance, resolutiontolerance)

implicit none

type(linelist), dimension(:), allocatable :: referencelinelist, fittedlines
type(linelist), dimension(:,:), allocatable :: population
type(linelist), dimension(:,:), allocatable ::  breed
type(spectrum), dimension(:,:), allocatable :: synthspec
type(spectrum), dimension(:) :: inputspectrum, fittedspectrum
integer :: popsize, i, spectrumlength, lineid, loc1, loc2, nlines, gencount, generations, popnumber
real, dimension(:), allocatable :: rms
real :: random, pressure, convergence, oldrms
real :: resolutionguess, redshiftguess, redshifttolerance, resolutiontolerance

!initialisation

  popsize=50
  pressure=0.1 !pressure * popsize needs to be an integer
  convergence=0.0 !new rms / old rms
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

! iterate until rms changes by less than 1 percent

  gencount=1

  do while (gencount .le. generations)
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
        where (abs(population(popnumber,lineid)%redshift*population(popnumber,lineid)%wavelength - synthspec(:,popnumber)%wavelength) .lt. (5*population(popnumber,lineid)%wavelength/population(popnumber,lineid)%resolution))
          synthspec(:,popnumber)%flux = synthspec(:,popnumber)%flux + &
          &population(popnumber,lineid)%peak*exp((-(synthspec(:,popnumber)%wavelength-population(popnumber,lineid)%redshift*population(popnumber,lineid)%wavelength)**2)/(2*(population(popnumber,lineid)%wavelength/population(popnumber,lineid)%resolution)**2))

        end where
      enddo

    !now calculate RMS for the "models"

      rms(popnumber)=sum((synthspec(:,popnumber)%flux-inputspectrum(:)%flux)**2)/spectrumlength

    enddo

    !if that was the last generation then exit before mutating

    if (gencount .eq. generations) exit

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

!    if (mod(gencount,100) .eq.0 .or. gencount.eq.1) then
!      print "(X,A,A,i5,A,4(X,F12.3))",gettime()," : ",gencount, " generations  ", population(minloc(rms,1),1)%resolution, 3.e5*(population(minloc(rms,1),1)%redshift-1), minval(rms,1), tmpvar
!    endif

    gencount=gencount+1
    convergence=minval(rms,1)/oldrms
    if (convergence .gt. 1.0) then
      convergence = 1./convergence
    endif

  enddo

!copy fit results into arrays to return

  fittedspectrum = synthspec(:,minloc(rms,1))
  fittedlines = population(minloc(rms,1),:)

!deallocate arrays

  deallocate(synthspec)
  deallocate(rms)
  deallocate(breed)
  deallocate(population)

end subroutine fit
end module mod_fit
