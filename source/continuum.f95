module mod_continuum
use mod_types
use mod_quicksort

contains

subroutine fit_continuum(realspec, spectrumlength, continuum)

implicit none
type(spectrum), dimension(:), allocatable :: realspec
type(spectrum), dimension(:), allocatable :: continuum
real, dimension(101) :: spectrumchunk
integer :: i, spectrumlength

! take the 25th percentile value of 101-element chunks of the spectrum

  allocate(continuum(spectrumlength))

  continuum%wavelength = realspec%wavelength
  continuum%flux=0.D0

  do i=51,spectrumlength-50
    spectrumchunk = realspec(i-50:i+50)%flux
    call qsort(spectrumchunk)
    continuum(i)%flux = spectrumchunk(25)
  enddo

!fill in the ends

  continuum(1:50)%flux=continuum(51)%flux
  continuum(spectrumlength-50:spectrumlength)%flux=continuum(spectrumlength-51)%flux

  realspec%flux = realspec%flux - continuum%flux

end subroutine fit_continuum

end module mod_continuum
