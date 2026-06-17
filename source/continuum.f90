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
integer :: i,halfwindow,cindex1,cindex2,cindex_median
real :: mean,std

#ifdef CO
  print *,"subroutine: fit_continuum"
#endif

! todo: better treatment of ends
! track running sum and sum of squares, adding new element and subtracting old each time

  allocate(continuum(spectrumlength))
  continuum%wavelength = realspec%wavelength
  continuum%flux=0.D0

  if (subtractcontinuum) then
! note that in integer maths only integer part is taken. ie 101/2 (=50.5) = 50
    halfwindow=continuumwindow/2
    allocate(spectrumchunk(continuumwindow))

     do i=halfwindow+1,spectrumlength-halfwindow
! get chunk of spectrum
       spectrumchunk = realspec(i-halfwindow:i+halfwindow)%flux
! get mean and standard deviation of chunk
      mean=sum(spectrumchunk)/size(spectrumchunk)
      std=(sum((spectrumchunk-mean)**2)/size(spectrumchunk))**0.5
! sort array
       call qsort(spectrumchunk)
! find locations where values are within 3std of the mean
      cindex1=minloc(abs(mean-(spectrumchunk+3*std)),1)
      cindex2=minloc(abs(mean-(spectrumchunk-3*std)),1)
! get the median. note than integer division only keeps integer part
      cindex_median=cindex1+(cindex2-cindex1)/2
      if (mod(cindex2 - cindex1,2).eq.1) then ! mean of two central elements
        continuum(i)%flux = (spectrumchunk(cindex_median)+spectrumchunk(cindex_median+1))*0.5
      else
        continuum(i)%flux = spectrumchunk(cindex_median)
      endif
     enddo

!fill in the ends

    continuum(1:halfwindow)%flux=continuum(halfwindow+1)%flux
    continuum(spectrumlength-halfwindow:spectrumlength)%flux=continuum(spectrumlength-halfwindow-1)%flux

    realspec%flux = realspec%flux - continuum%flux

  endif

end subroutine fit_continuum

end module mod_continuum
