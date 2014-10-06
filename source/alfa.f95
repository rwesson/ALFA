program alfa

use mod_types
use mod_quicksort

implicit none
integer :: I, J, lineid, popnumber, gencount, IO, spectrumlength, nlines
integer :: loc1, loc2
integer :: wlength
real :: temp1, temp2, gaussian, gaussianflux, random, continuumtemp, mutation
character*512 :: filename
character*10 :: gettime

integer :: popsize, generations

real :: pressure ! the fraction of candidate spectra which breed in every generation

type(linelist) :: referencelinelist
type(linelist), dimension(:),allocatable :: population
type(linelist), dimension(:),allocatable ::  breed
type(spectrum), dimension(:,:), allocatable :: synthspec
type(spectrum), dimension(:), allocatable :: realspec
type(spectrum), dimension(20) :: spectrumchunk
type(spectrum), dimension(:), allocatable :: continuum

real :: null

real, dimension(:), allocatable :: rms

logical :: file_exists

!temp XXXX
open (101,file="intermediate",status="replace")

! initialise stuff for genetics

popsize=100
generations=1000

pressure=0.25 !pressure * popsize needs to be an integer

allocate(rms(popsize))

rms=0.D0
! random seed

call init_random_seed()

!read files

! read in spectrum to fit

  call get_command_argument(1,filename)

  if (trim(filename)=="") then
    print *,gettime()," : error: No input spectrum specified"
    stop
  endif

  inquire(file=filename, exist=file_exists) ! see if the input file is present

  if (.not. file_exists) then
    print *,gettime()," : error: input spectrum ",trim(filename)," does not exist"
    stop
  else
    print *, gettime()," : reading in spectrum ",trim(filename)
    I = 0
    OPEN(199, file=filename, iostat=IO, status='old')
      DO WHILE (IO >= 0)
      READ(199,*,end=112) null
      I = I + 1
    END DO
    112 spectrumlength=I
  endif

  !then allocate and read

  allocate (synthspec(spectrumlength,popsize))
  allocate (realspec(spectrumlength))
  allocate (continuum(spectrumlength))

  REWIND (199)
  DO I=1,spectrumlength
    READ(199,*) temp1, temp2
    synthspec(i,:)%wavelength = temp1
    realspec(i)%wavelength = temp1
    realspec(i)%flux = temp2
  END DO
  CLOSE(199)

  ! read in template spectrum

  call get_command_argument(2,filename)

  if (trim(filename)=="") then
    print *,gettime()," : error: No line catalogue specified"
    stop
  endif

  inquire(file=filename, exist=file_exists) ! see if the input file is present

  if (.not. file_exists) then
    print *,gettime()," : error: line catalogue ",trim(filename)," does not exist"
    stop
  else
    print *,gettime()," : fitting spectrum using line catalogue ",trim(filename)
    I = 0
    OPEN(199, file=filename, iostat=IO, status='old')
    DO WHILE (IO >= 0)
      READ(199,*,end=110) null
      if (null .ge. minval(realspec%wavelength) .and. null .le. maxval(realspec%wavelength)) then
      !only read in lines that lie within the observed wavelength range
        I = I + 1
      endif
    END DO
    110     nlines=I
  endif

  if (nlines .eq. 0) then
    print *,gettime()," : error : Line catalogue does not overlap with input spectrum"
    stop
  endif

!then allocate and read

  allocate(referencelinelist%peak(nlines))
  allocate(referencelinelist%wavelength(nlines))

  REWIND (199)
  I=1
  do while (i .le. nlines)
    READ(199,*) temp1
    if (temp1 .ge. minval(realspec%wavelength)) then
      referencelinelist%wavelength(i) = temp1
      i=i+1
    endif
    if (temp1 .ge.  maxval(realspec%wavelength)) then
      exit
    endif
  END DO
  CLOSE(199)

! First, subtract the continuum
! in 20-element chunks of the spectrum, calculate the mean of the lowest 5 points.
! Take this as the continuum

! todo: spline fit to points
! or at least linear interpolation

continuum%wavelength = realspec%wavelength

do i=11,spectrumlength-10,20
  spectrumchunk = realspec(i-10:i+9)
  do j=1,15
    spectrumchunk(maxloc(spectrumchunk%flux))%flux = 0.0
  enddo
  continuumtemp = sum(spectrumchunk%flux)/5.
  do j=-10,9
     continuum(i+j)%flux = continuumtemp
  enddo
enddo

realspec%flux = realspec%flux - continuum%flux

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
  population(popnumber)%width=1.5
  population(popnumber)%redshift=1.0
end do

do gencount=1,generations

!reset stuff to zero before doing calculations

rms = 0.D0
synthspec%flux=0.D0

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

  !now calculate "RMS" for the "models"

  do wlength=1,spectrumlength
    rms(popnumber)=rms(popnumber)+(((synthspec(wlength,popnumber)%flux-realspec(wlength)%flux)**2))
  end do

end do

  !next, cream off the well performing models - put the population member with the lowest RMS into the breed array, replace the RMS with something very high so that it doesn't get copied twice, repeat until a fraction equal to the pressure factor have been selected

  do i=1,int(popsize*pressure) 
    breed(i) = population(minloc(rms,1))
    rms(minloc(rms,1))=1.e10
  end do
!print *,population(minloc(rms,1))%peak(1)
!XXXX
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
      population(popnumber)%redshift = population(popnumber)%redshift * ((999999.+mutation())/1000000.)
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

real function gaussian(x,a,b,c)
!return the value of a gaussian function with parameters a, b, and c, at a value
!of x
  implicit none
  real :: x,a,b,c

  gaussian = a*exp((-(x-b)**2)/(2*c**2))
  return

end function gaussian

real function gaussianflux(a,c)
!return the integral of the gaussian, equal to a*c*(2*pi**0.5)
  implicit none
  real :: a,c,pi

  pi=3.14159265359
  gaussianflux = a*c*(2*pi)**0.5
  return

end function gaussianflux

character*10 function gettime()

character*10 :: time

  call DATE_AND_TIME(TIME=time)
  gettime = time(1:2)//":"//time(3:4)//":"//time(5:6)
  return

end function gettime

SUBROUTINE init_random_seed()
  INTEGER :: i, n, clock
  INTEGER, DIMENSION(:), ALLOCATABLE :: seed

  n=20
  i=n
  CALL RANDOM_SEED(size = n)
  ALLOCATE(seed(n))

  CALL SYSTEM_CLOCK(COUNT=clock)

  seed = clock + 37 * (/ (i - 1, i = 1, n) /)
  CALL RANDOM_SEED(PUT = seed)

  DEALLOCATE(seed)
END SUBROUTINE

real function mutation()
  implicit none
  real :: random

  mutation=1.0

  call random_number(random)
  if (random .le. 0.05) then
    mutation=random/0.05
  elseif (random .ge. 0.95) then
    mutation=2+((random-1)/0.05)
  endif

  return

end function mutation
