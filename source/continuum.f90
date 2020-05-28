!Copyright (C) 2013- Roger Wesson
!Free under the terms of the GNU General Public License v3

module mod_continuum
use mod_types
use mod_quicksort
use mod_globals

contains

subroutine fit_continuum(realspec,continuum)

implicit none
type(spectrum), dimension(:), allocatable :: realspec,continuum
real, dimension(:), allocatable :: spectrumchunk
integer :: i,halfwindow

#ifdef CO
  print *,"subroutine: fit_continuum"
#endif

  allocate(continuum(spectrumlength))
  continuum%wavelength = realspec%wavelength
  continuum%flux=0.D0

  if (subtractcontinuum) then
! take the 25th percentile value of chunks of the spectrum defined by window
    halfwindow=continuumwindow/2
! note that in integer maths only integer part is taken. ie 101/2 (=50.5) = 50
    allocate(spectrumchunk(continuumwindow))

    do i=halfwindow+1,spectrumlength-halfwindow
      spectrumchunk = realspec(i-halfwindow:i+halfwindow)%flux
      call qsort(spectrumchunk)
      continuum(i)%flux = spectrumchunk(halfwindow/2)
    enddo

!fill in the ends

    continuum(1:halfwindow)%flux=continuum(halfwindow+1)%flux
    continuum(spectrumlength-halfwindow:spectrumlength)%flux=continuum(spectrumlength-halfwindow-1)%flux

    realspec%flux = realspec%flux - continuum%flux

  endif

end subroutine fit_continuum

end module mod_continuum
