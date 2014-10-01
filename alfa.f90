program alfa

use mod_quicksort

implicit none
integer :: I, J, lineid, popnumber, gencount, IO, spectrumlength, nlines, counter
integer :: loc1, loc2, loc
integer :: wlength
real :: temp1, temp2, gaussian, gaussianflux, random, redshift, continuumtemp
character*512 :: filename
character*10 :: gettime

integer :: popsize, generations

real :: pressure ! the fraction of candidate spectra which breed in every generation
real :: mutationrate ! in each generation, this fraction of lines will be doubled, and this fraction of lines will be halved

type spectrum
  real :: wavelength
  real :: flux
end type

type templateline
  real :: wavelength
  real :: width
  real :: peak
end type

type(templateline), dimension(:),allocatable :: linelist
type(templateline), dimension(:,:),allocatable :: population
type(templateline), dimension(:,:),allocatable ::  breed
type(spectrum), dimension(:,:), allocatable :: synthspec
type(spectrum), dimension(:), allocatable :: realspec
type(spectrum), dimension(20) :: spectrumchunk
type(spectrum), dimension(:), allocatable :: continuum
type(spectrum), dimension(1) :: bestfit

real :: weightfactor

real :: null

real, dimension(:), allocatable :: rms,medianrms

logical :: file_exists

! initialise stuff for genetics

popsize=50
generations=500

pressure=0.5 !pressure * popsize needs to be an integer
mutationrate=0.1 !mutation rate * 3 needs to be less than one

redshift=0.0

allocate(rms(popsize))

rms=0.D0
! random seed

call init_random_seed()

!read files

! read in spectrum to fit

        call get_command_argument(1,filename)

        if (trim(filename)=="") then
          print *,gettime(),": error: No input spectrum specified"
          stop
        endif

        inquire(file=filename, exist=file_exists) ! see if the input file is present

        if (.not. file_exists) then
          print *,gettime(),": error: input spectrum ",trim(filename)," does not exist"
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
          print *,gettime(),": error: No line catalogue specified"
          stop
        endif

        inquire(file=filename, exist=file_exists) ! see if the input file is present

        if (.not. file_exists) then
          print *,gettime(),": error: line catalogue ",trim(filename)," does not exist"
          stop
        else
          print *,gettime()," : fitting spectrum using line catalogue ",trim(filename)
          I = 0
          OPEN(199, file=filename, iostat=IO, status='old')
    DO WHILE (IO >= 0)
      READ(199,*,end=110) null
      I = I + 1
    END DO
  110     nlines=I
endif

!then allocate and read
  allocate (linelist(nlines))

  REWIND (199)
  DO I=1,nlines
    READ(199,*) temp1, temp2
    linelist(i)%wavelength = temp1
    linelist(i)%peak = temp2
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

weightfactor=2.0*maxval(realspec%flux)

!allocate some more arrays

  allocate (breed(nlines,int(popsize*pressure)))
  allocate (population(nlines,popsize))

! now create population of synthetic spectra
! todo, make sure no lines are included which are outside the wavelength range
! of the observations

do lineid=1,nlines
  population(lineid,:)%wavelength = linelist(lineid)%wavelength
end do

population(:,:)%peak = 0.01
population(:,:)%width = 7.0

do gencount=1,generations

!reset stuff to zero before doing calculations

rms = 0.D0
synthspec%flux=0.D0

do popnumber=1,popsize

!calculate synthetic spectra - reset to 0 before synthesizing
!line fluxes are calculated within 7 sigma of mean

  do wlength=1,spectrumlength 
    do lineid=1,nlines
      if (abs(population(lineid,popnumber)%wavelength - synthspec(wlength,popnumber)%wavelength) .lt. (7*population(lineid,popnumber)%width)) then
        synthspec(wlength,popnumber)%flux = synthspec(wlength,popnumber)%flux + &
        &gaussian(synthspec(wlength,popnumber)%wavelength,&
        &population(lineid,popnumber)%peak,population(lineid,popnumber)%wavelength, population(lineid,popnumber)%width)
      endif
    end do
  end do

  !now calculate "RMS" for the "models"

  do wlength=1,spectrumlength
!least absolute difference
!rms(popnumber) = rms(popnumber) + abs(synthspec(wlength,popnumber)%flux-realspec(wlength)%flux)
!weighted rms
!       rms(popnumber)=rms(popnumber)+(((synthspec(wlength,popnumber)%flux-realspec(wlength)%flux)**2) &
!&* (1-(realspec(wlength)%flux/weightfactor)))
!sum of ratios
!       rms(popnumber) = rms(popnumber)+(synthspec(wlength,popnumber)%flux/realspec(wlength)%flux)
!simple rms
    rms(popnumber)=rms(popnumber)+(((synthspec(wlength,popnumber)%flux-realspec(wlength)%flux)**2))
  end do

end do
!print *, minval(rms,1), maxval(rms,1)

  !next, cream off the well performing models - put the population member with the lowest RMS into the breed array, replace the RMS with something very high so that it doesn't get copied twice, repeat until a fraction equal to the pressure factor have been selected

  do i=1,int(popsize*pressure) 
    breed(:,i) = population(:,minloc(rms,1))
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
      population(:,i)%peak=(breed(:,loc1)%peak + breed(:,loc2)%peak)/2.0
      population(:,i)%width=(breed(:,loc1)%width + breed(:,loc2)%width)/2.0
    end do
    !then, "mutate"
  
    do popnumber=1,popsize ! mutation of line width
      if (random .le. (0.5*mutationrate)) then
        population(:,popnumber)%width = population(:,popnumber)%width * 0.5
      elseif (random .gt. (0.5*mutationrate) .and. random .le. mutationrate) then
          population(:,popnumber)%width = population(:,popnumber)%width * 2.0
      elseif (random .gt. mutationrate .and. random .le. (2.0*mutationrate)) then
          population(:,popnumber)%width = population(:,popnumber)%width * 0.95
      elseif (random .gt. (2.0*mutationrate).and. random .le. (3.0*mutationrate)) then
          population(:,popnumber)%width = population(:,popnumber)%width *1.05
      endif
  
      do lineid=1,nlines !mutation of line fluxes
        call random_number(random)
  !continuous mutation:
  !if (random .le. 0.10) then 
  !population(lineid,popnumber)%peak = population(lineid,popnumber)%peak * (1+(10*(0.10-random)))
  !elseif (random .ge. 0.90) then
  !population(lineid,popnumber)%peak = population(lineid,popnumber)%peak * (1-(10*(random-0.90)))
  !endif
  !population(lineid,popnumber)%peak = population(lineid,popnumber)%peak * &
  !& 1+(10.*(0.5-random)**7.)
  !discrete mutation:
        if (random .le. (0.5*mutationrate)) then 
          population(lineid,popnumber)%peak = population(lineid,popnumber)%peak * 0.5
        elseif (random .gt. (0.5*mutationrate) .and. random .le. mutationrate) then
          population(lineid,popnumber)%peak = population(lineid,popnumber)%peak * 2.0
        elseif (random .gt. mutationrate .and. random .le. (2.0*mutationrate)) then
          population(lineid,popnumber)%peak = population(lineid,popnumber)%peak * 0.95
        elseif (random .gt. (2.0*mutationrate).and. random .le. (3.0*mutationrate)) then
          population(lineid,popnumber)%peak = population(lineid,popnumber)%peak * 1.05
        endif
      enddo
    enddo
  
  endif
  
  if (mod(gencount,generations/10) .eq.0) then
    print *,gettime()," : completed ",100*gencount/generations, "%"
    print *,gettime()," : line width = ",population(i,minloc(rms,1))%width
    print *,gettime()," : min rms = ",minval(rms,1)
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

!do i=1,nlines
!  print *,"#",population(i,1)%wavelength,gaussianflux(population(i,minloc(rms,1))%peak,linewidth), &
!    & gaussianflux(population(i,maxloc(rms,1))%peak,linewidth)
!end do
open(100,file="outputfit")

write (100,*) "#wavelength    fitted spectrum     cont-subbed orig   continuum    residuals"
do i=1,spectrumlength
  write(100,*) synthspec(i,minloc(rms,1))%wavelength,synthspec(i,minloc(rms,1))%flux, realspec(i)%flux, continuum(i)%flux, realspec(i)%flux - synthspec(i,minloc(rms,1))%flux
end do

close(100)

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
!return the integral of the gaussian, equal to a*c*(pi**0.5)
  implicit none
  real :: a,c,pi

  pi=3.14159265359
  gaussianflux = a*c*(pi**0.5)
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
