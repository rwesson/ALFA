!Copyright (C) 2013- Roger Wesson
!Free under the terms of the GNU General Public License v3

module mod_continuum
use mod_types
use mod_quicksort

contains

subroutine fit_continuum(realspec, spectrumlength, continuum)

implicit none
type(spectrum), dimension(:), allocatable :: realspec
type(spectrum), dimension(:), allocatable :: continuum
real, dimension(:), allocatable :: spectrumchunk
integer :: i, spectrumlength
integer :: window,halfwindow

#ifdef CO
  print *,"subroutine: fit_continuum"
#endif

! take the 25th percentile value of chunks of the spectrum defined by window

  window=21 ! must be odd number, add check for this
  halfwindow=window/2
  !note that in integer maths only integer part is taken. ie 101/2 (=50.5) = 50
  allocate(continuum(spectrumlength))
  allocate(spectrumchunk(window))

  continuum%wavelength = realspec%wavelength
  continuum%flux=0.D0

  do i=halfwindow+1,spectrumlength-halfwindow
    spectrumchunk = realspec(i-halfwindow:i+halfwindow)%flux
    call qsort(spectrumchunk)
    continuum(i)%flux = spectrumchunk(halfwindow/2)
  enddo

!fill in the ends

  continuum(1:halfwindow)%flux=continuum(halfwindow+1)%flux
  continuum(spectrumlength-halfwindow:spectrumlength)%flux=continuum(spectrumlength-halfwindow-1)%flux

  realspec%flux = realspec%flux - continuum%flux

end subroutine fit_continuum

end module mod_continuum
