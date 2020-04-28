!Copyright (C) 2013- Roger Wesson
!Free under the terms of the GNU General Public License v3

module mod_functions
contains

real function gaussianflux(a,c)
!return the integral of the gaussian, equal to a*c*(2*pi**0.5)
  implicit none
  real :: a,c,pi

  pi=3.14159265359
  gaussianflux = a*c*(2*pi)**0.5
  return

end function gaussianflux

character(len=11) function gettime()
implicit none
character(len=10) :: time
character(len=11), save :: oldtime

!debugging
#ifdef CO
        !print *,"function: gettime"
#endif

  call date_and_time(TIME=time)
  gettime = time(1:2)//":"//time(3:4)//":"//time(5:6)//" : "
  if (gettime .eq. oldtime) then
    gettime = "          "
  else
    oldtime = gettime
  endif
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
!generate a random number with a 90% chance of being 1., 10% chance of being between 0 and 2
  implicit none
  real :: random

  mutation=1.0

  call random_number(random)
  if (random .le. 0.05) then
    mutation=20.*random
  elseif (random .ge. 0.95) then
    mutation=2+(20.*(random-1))
  endif

  return

end function mutation

subroutine makespectrum(lines,spec)
!synthesizes a spectrum from fitted parameters
!calculate fluxes within 5 sigma of the mean
  use mod_types
  implicit none
  integer :: i
  type(spectrum), dimension(:) :: spec
  type(linelist), dimension(:), intent(in) :: lines
  real :: sigma

#ifdef CO
!not so useful, gets called hundreds of times
!  print *,"subroutine: makespectrum"
#endif

  do i=1,size(lines)
    sigma=lines(i)%wavelength/lines(i)%resolution
    where (abs(lines(i)%redshift*lines(i)%wavelength - spec%wavelength) .lt. 5.*sigma)
      spec%flux = spec%flux + &
      &lines(i)%peak*exp((-(spec%wavelength-lines(i)%redshift*lines(i)%wavelength)**2)/(2*sigma**2))
    end where
  enddo

end subroutine makespectrum

end module mod_functions
