module mod_continuum
use mod_types
contains

subroutine fit_continuum(realspec, spectrumlength, continuum)

implicit none
type(spectrum), dimension(:), allocatable :: realspec
type(spectrum), dimension(20) :: spectrumchunk
type(spectrum), dimension(:), allocatable :: continuum
integer :: i, j, spectrumlength

! in 20-element chunks of the spectrum, calculate the mean of the lowest 5 points.
! Take this as the continuum

! todo: spline fit to points
! or at least linear interpolation

  allocate(continuum(spectrumlength))

  continuum%wavelength = realspec%wavelength
  continuum%flux=0.D0

  do i=11,spectrumlength-10
    spectrumchunk = realspec(i-10:i+10)
    do j=1,15
      spectrumchunk(maxloc(spectrumchunk%flux))%flux = 0.0
    enddo
    continuum(i)%flux = sum(spectrumchunk%flux)/5
  enddo
!fill in the ends

  continuum(1:10)%flux=continuum(11)%flux
  continuum(spectrumlength-9:spectrumlength)%flux=continuum(spectrumlength-10)%flux

  realspec%flux = realspec%flux - continuum%flux

end subroutine fit_continuum

end module mod_continuum 
