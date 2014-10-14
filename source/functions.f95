module mod_routines
contains

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
! write out the time that the function was called
  implicit none
  character*10 :: time

  call DATE_AND_TIME(TIME=time)
  gettime = time(1:2)//":"//time(3:4)//":"//time(5:6)
  return

end function gettime

SUBROUTINE init_random_seed()
! seed the random number generation
  implicit none
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
end module mod_routines
