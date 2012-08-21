program alfa

implicit none
integer :: I, lineid, popnumber, gencount, IO, spectrumlength, nlines, counter
integer :: loc1, loc2, loc
integer :: wlength
double precision :: temp1, temp2, gaussian, random
character*512 :: filename

integer :: popsize, generations

double precision :: pressure ! the fraction of candidate spectra which breed in every generation
double precision :: mutationrate ! in each generation, this fraction of lines will be doubled, and this fraction of lines will be halved

type spectrum
  double precision :: wavelength
  double precision :: flux
end type

type templateline
  double precision :: wavelength
  double precision :: peak
end type

type(templateline), dimension(:),allocatable :: linelist
type(templateline), dimension(:,:),allocatable :: population
type(templateline), dimension(:,:),allocatable ::  breed
type(spectrum), dimension(:,:), allocatable :: synthspec
type(spectrum), dimension(:), allocatable :: realspec
double precision, dimension(:), allocatable :: continuum !linear continuum for each population member

CHARACTER*1 :: null

double precision, dimension(:), allocatable :: rms,medianrms

! initialise stuff

popsize=100
generations=100

pressure=0.1 !pressure * popsize needs to be an integer
mutationrate=0.1 !mutation rate * 3 needs to be less than one

allocate(continuum(popsize))
allocate(rms(popsize))

rms=0.D0
! random seed

call init_random_seed()

! read in template spectrum

! TODO replace with reading of filename from command line
filename=trim("peaks2.txt")

        I = 0
        OPEN(199, file=filename, iostat=IO, status='old') 
                DO WHILE (IO >= 0)
                        READ(199,*,end=111) null
                        I = I + 1
                END DO
        111 print *
        nlines=I

!then allocate and read
        allocate (linelist(nlines))

        REWIND (199)
        DO I=1,nlines
                READ(199,*,end=110) temp1, temp2
                linelist(i)%wavelength = temp1
                linelist(i)%peak = temp2 
        END DO
        CLOSE(199)

110 print *
!        110 PRINT "(X,A11,I4,A15,I4,A9)","read in ", I - 1," lines (out of ",listlength," in file)"

! read in spectrum to fit
! TODO replace with reading of filename from command line

filename="ngc6543_4553A_short.bdf.ascii"

        I = 0
        OPEN(199, file=filename, iostat=IO, status='old')
                DO WHILE (IO >= 0)
                        READ(199,*,end=112) null
                        I = I + 1
                END DO
        112 print *
        spectrumlength=I

!then allocate and read

        allocate (synthspec(spectrumlength,popsize))
        allocate (breed(nlines,popsize/2))
        allocate (realspec(spectrumlength))
        allocate (population(nlines,popsize))

        REWIND (199)
        DO I=1,spectrumlength
                READ(199,*,end=113) temp1, temp2
                synthspec(i,:)%wavelength = temp1 
                realspec(i)%wavelength = temp1
                realspec(i)%flux = temp2
        END DO
        CLOSE(199)

113 print *

! now create population of synthetic spectra

do gencount=1,generations

!reset stuff to zero before doing calculations

rms = 0.D0
synthspec%flux=0.D0

        do popnumber=1,popsize

        !randomize fluxes for first generation

            if (gencount .eq. 1) then
                do lineid=1,nlines
                  call RANDOM_NUMBER(random)
                  population(lineid,popnumber)%wavelength = linelist(lineid)%wavelength
!                  population(lineid,popnumber)%peak = (0.9+(0.2*random))*linelist(lineid)%peak !randomize to within +-10% of known flux
                  population(lineid,popnumber)%peak = random*0.1 !start small
                end do
                call RANDOM_NUMBER(random)
                continuum(popnumber)=random*0.05 !random small value for continuum
            endif

        !calculate synthetic spectra - reset to 0 before synthesizing

            do wlength=1,spectrumlength 
              do lineid=1,nlines 
                synthspec(wlength,popnumber)%flux = synthspec(wlength,popnumber)%flux + &
                &gaussian(synthspec(wlength,popnumber)%wavelength,&
                &population(lineid,popnumber)%peak,population(lineid,popnumber)%wavelength, dble(0.7)) 
              end do
            end do

        !add continuum to line fluxes

            synthspec(:,popnumber)%flux = synthspec(:,popnumber)%flux + continuum(popnumber)

        !now calculate RMS for the "models".  reset to zero first

            rms(popnumber)=0.D0

            do wlength=1,spectrumlength
              rms(popnumber)=rms(popnumber)+((synthspec(wlength,popnumber)%flux-realspec(wlength)%flux)**2)
            end do
            rms(popnumber)=(rms(popnumber)/spectrumlength)**0.5

        end do

        !next, cream off the well performing models - put the population member with the lowest RMS into the breed array, replace the RMS with something very high so that it doesn't get copied twice, repeat eg 500 times to get the best half of the models
PRINT *,minval(rms),maxval(rms)

if (gencount.eq.generations) then
  PRINT *
  do i=1,spectrumlength
    print *,synthspec(i,minloc(rms,1))%wavelength,synthspec(i,minloc(rms,1))%flux
  end do
end if

!print *,gencount,minval(rms)
        do i=1,int(popsize*pressure) 
          breed(:,i) = population(:,minloc(rms,1))
          rms(minloc(rms,1))=1.e10
        end do

        !then, "breed" pairs

        do i=1,popsize 
          call random_number(random)
          loc1=int(popsize*random*pressure)
          call random_number(random)
          loc2=int(popsize*random*pressure)
          population(:,i)%peak=(breed(:,loc1)%peak + breed(:,loc2)%peak)/2.0 
        end do
        !then, "mutate"

       do popnumber=1,popsize
         do lineid=1,nlines !mutation of line fluxes
           call random_number(random)
!continuous mutation:
!population(lineid,popnumber)%peak = population(lineid,popnumber)%peak * &
!& 1+(1000.*(0.5-random)**7.)
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

!if (gencount .gt. 500) then !only start mutating continuum once lines are stable
         call random_number(random) !mutation of continuum
         if (random .le. (0.5*mutationrate)) then
           continuum(popnumber) = continuum(popnumber) * 0.5
         elseif (random .gt. (0.5*mutationrate) .and. random .le. mutationrate) then
           continuum(popnumber) = continuum(popnumber) * 2.0
         elseif (random .gt. mutationrate .and. random .le. (2.0*mutationrate)) then
           continuum(popnumber) = continuum(popnumber) * 0.95
         elseif (random .gt. (2.0*mutationrate).and. random .le. (3.0*mutationrate)) then
           continuum(popnumber) = continuum(popnumber) * 1.05
         endif
!endif
       enddo

end do

end program alfa

double precision function gaussian(x,a,b,c)
!return the value of a gaussian function with parameters a, b, and c, at a value
!of x
  implicit none
  double precision :: x,a,b,c

  gaussian = a*exp((-(x-b)**2)/(2*c**2))
  return

end function 

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
