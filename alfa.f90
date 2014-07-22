program alfa

implicit none
integer :: I, lineid, popnumber, gencount, IO, spectrumlength, nlines, counter
integer :: loc1, loc2, loc
integer :: wlength
real :: temp1, temp2, gaussian, gaussianflux, random, linewidth, redshift
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
  real :: peak
end type

type(templateline), dimension(:),allocatable :: linelist
type(templateline), dimension(:,:),allocatable :: population
type(templateline), dimension(:,:),allocatable ::  breed
type(spectrum), dimension(:,:), allocatable :: synthspec
type(spectrum), dimension(:), allocatable :: realspec
type(spectrum), dimension(1) :: bestfit
real, dimension(:), allocatable :: continuum !linear continuum for each population member
real, dimension(:), allocatable :: breedcontinuum

real :: weightfactor

real :: null

real, dimension(:), allocatable :: rms,medianrms

logical :: file_exists

!read files

! read in spectrum to fit

call get_command_argument(1,filename)
filename=trim(filename)

if (filename=="") then
        print *,gettime(),": error: No input spectrum specified"
        stop
endif

inquire(file=filename, exist=file_exists) ! see if the input file is present

if (.not. file_exists) then
       print *,gettime(),": error: input spectrum ",filename," does not exist"
       stop
endif

print *, gettime()," : reading in spectrum ",filename

        I = 0
        OPEN(199, file=filename, iostat=IO, status='old')
                DO WHILE (IO >= 0)
                        READ(199,*,end=112) null
                        I = I + 1
                END DO
        112 spectrumlength=I

!then allocate and read

        allocate (synthspec(spectrumlength,popsize))
        allocate (realspec(spectrumlength))

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
filename=trim(filename)
if (filename=="") then
        print *,gettime(),": error: No line catalogue specified"
        stop
endif

inquire(file=filename, exist=file_exists) ! see if the input file is present

if (.not. file_exists) then
       print *,gettime(),": error: line catalogue ",filename," does not exist"
       stop
endif

print *,gettime()," : fitting spectrum using line catalogue ",filename

        I = 0
        OPEN(199, file=filename, iostat=IO, status='old')
          DO WHILE (IO >= 0)
            READ(199,*,end=110) null
            I = I + 1
          END DO

110     nlines=I

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
! in chunks of the spectrum, calculate the mean of the lowest 5% of the fluxes.
! Take this as the continuum

! do in chunks of 1/10 of the wavelength range



! initialise stuff for genetics

popsize=100
generations=10000

pressure=0.1 !pressure * popsize needs to be an integer
mutationrate=0.1 !mutation rate * 3 needs to be less than one

linewidth=0.7
redshift=0.0

allocate(continuum(popsize))
allocate(rms(popsize))

rms=0.D0
! random seed

call init_random_seed()

! read in spectrum to fit


call get_command_argument(1,filename)
filename=trim(filename)

if (filename=="") then
        print *,gettime(),": error: No input spectrum specified"
        stop
endif

inquire(file=filename, exist=file_exists) ! see if the input file is present

if (.not. file_exists) then
       print *,gettime(),": error: input spectrum ",filename," does not exist"
       stop
endif

print *, gettime()," : reading in spectrum ",filename

        I = 0
        OPEN(199, file=filename, iostat=IO, status='old')
                DO WHILE (IO >= 0)
                        READ(199,*,end=112) null
                        I = I + 1
                END DO
        112 spectrumlength=I

!then allocate and read

        allocate (synthspec(spectrumlength,popsize)) 
        allocate (realspec(spectrumlength)) 

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
filename=trim(filename)
if (filename=="") then
        print *,gettime(),": error: No line catalogue specified"
        stop
endif

inquire(file=filename, exist=file_exists) ! see if the input file is present

if (.not. file_exists) then
       print *,gettime(),": error: line catalogue ",filename," does not exist"
       stop
endif

print *,gettime()," : fitting spectrum using line catalogue ",filename

        I = 0
        OPEN(199, file=filename, iostat=IO, status='old')
          DO WHILE (IO >= 0)
            READ(199,*,end=110) null
            I = I + 1
          END DO

110     nlines=I

!then allocate and read
        allocate (linelist(nlines))

        REWIND (199)
        DO I=1,nlines
                READ(199,*) temp1, temp2
                linelist(i)%wavelength = temp1
                linelist(i)%peak = temp2
        END DO
        CLOSE(199)

!        110 PRINT "(X,A11,I4,A15,I4,A9)","read in ", I - 1," lines (out of
!        ",listlength," in file)"


weightfactor=2.0*maxval(realspec%flux)

!allocate some more arrays

        allocate (breed(nlines,int(popsize*pressure)))
        allocate (breedcontinuum(int(popsize*pressure)))
        allocate (population(nlines,popsize))

! now create population of synthetic spectra

do lineid=1,nlines
  population(lineid,:)%wavelength = linelist(lineid)%wavelength
  population(lineid,:)%peak = 0.01!0.1+random*0.5 !start small and definitely not really close to zero
end do
continuum(:)=0.05 !random small value for continuum

do gencount=1,generations

!reset stuff to zero before doing calculations

rms = 0.D0
synthspec%flux=0.D0

do popnumber=1,popsize

!calculate synthetic spectra - reset to 0 before synthesizing
!line fluxes are calculated within 7 sigma of mean

  do wlength=1,spectrumlength 
    do lineid=1,nlines
      if (abs(population(lineid,popnumber)%wavelength - synthspec(wlength,popnumber)%wavelength) .lt. (7*linewidth)) then
        synthspec(wlength,popnumber)%flux = synthspec(wlength,popnumber)%flux + &
        &gaussian(synthspec(wlength,popnumber)%wavelength,&
        &population(lineid,popnumber)%peak,population(lineid,popnumber)%wavelength, linewidth) 
      endif
    end do
  end do

        !add continuum to line fluxes

  synthspec(:,popnumber)%flux = synthspec(:,popnumber)%flux + continuum(popnumber)

        !now calculate "RMS" for the "models".  reset to zero first

  rms(popnumber)=0.D0

  do wlength=1,spectrumlength
!least absolute difference
!rms(popnumber) = rms(popnumber) + abs(synthspec(wlength,popnumber)%flux-realspec(wlength)%flux)
!weighted rms
!       rms(popnumber)=rms(popnumber)+(((synthspec(wlength,popnumber)%flux-realspec(wlength)%flux)**2) &
!&* (1-(realspec(wlength)%flux/weightfactor)))
!sum of ratios
!             rms(popnumber) = rms(popnumber)+(synthspec(wlength,popnumber)%flux/realspec(wlength)%flux)
!simple rms
    rms(popnumber)=rms(popnumber)+(((synthspec(wlength,popnumber)%flux-realspec(wlength)%flux)**2))
  end do

end do

        !next, cream off the well performing models - put the population member with the lowest RMS into the breed array, replace the RMS with something very high so that it doesn't get copied twice, repeat until a fraction equal to the pressure factor have been selected

        do i=1,int(popsize*pressure) 
          breed(:,i) = population(:,minloc(rms,1))
          breedcontinuum(i) = continuum(minloc(rms,1))
          rms(minloc(rms,1))=1.e10
        end do

!then, "breed" pairs
!breed line fluxes and continuum!
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
          continuum(i)=(breedcontinuum(loc1) + breedcontinuum(loc2))/2.0
        end do
        !then, "mutate"

       do popnumber=1,popsize
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


         call random_number(random) !mutation of continuum
!if (random .le. 0.10) then
!continuum(popnumber) = continuum(popnumber) * (1+(10*(0.10-random)))
!elseif (random .ge. 0.90) then
!continuum(popnumber) = continuum(popnumber) * (1-(10*(random-0.90)))
!endif

         if (random .le. (0.5*mutationrate)) then
           continuum(popnumber) = continuum(popnumber) * 0.5
         elseif (random .gt. (0.5*mutationrate) .and. random .le. mutationrate) then
           continuum(popnumber) = continuum(popnumber) * 2.0
         elseif (random .gt. mutationrate .and. random .le. (2.0*mutationrate)) then
           continuum(popnumber) = continuum(popnumber) * 0.95
         elseif (random .gt. (2.0*mutationrate).and. random .le. (3.0*mutationrate)) then
           continuum(popnumber) = continuum(popnumber) * 1.05
         endif
       enddo

     endif

end do
!write out line fluxes of best fitting spectrum

!do i=1,nlines
!  print *,"#",population(i,1)%wavelength,gaussianflux(population(i,minloc(rms,1))%peak,linewidth), &
!    & gaussianflux(population(i,maxloc(rms,1))%peak,linewidth)
!end do
do i=1,spectrumlength
!  print *,i
  print *,synthspec(i,minloc(rms,1))%wavelength,synthspec(i,minloc(rms,1))%flux
end do

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
