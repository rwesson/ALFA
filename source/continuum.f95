module mod_continuum
use mod_types
use mod_quicksort

contains

subroutine fit_continuum(realspec, spectrumlength, continuum)

implicit none
type(spectrum), dimension(:), allocatable :: realspec
type(spectrum), dimension(:), allocatable :: continuum
real, dimension(51) :: spectrumchunk

integer :: i, j, spectrumlength

! in 20-element chunks of the spectrum, calculate the mean of the lowest 5 points.
! Take this as the continuum

! todo: spline fit to points
! or at least linear interpolation

  allocate(continuum(spectrumlength))

  continuum%wavelength = realspec%wavelength
  continuum%flux=0.D0

  do i=26,spectrumlength-25
    spectrumchunk = realspec(i-25:i+25)%flux
    call qsort(spectrumchunk)
    continuum(i)%flux = sum(spectrumchunk(1:10))/10
  enddo
!fill in the ends

  continuum(1:25)%flux=continuum(26)%flux
  continuum(spectrumlength-25:spectrumlength)%flux=continuum(spectrumlength-26)%flux

  realspec%flux = realspec%flux - continuum%flux

end subroutine fit_continuum

end module mod_continuum 
