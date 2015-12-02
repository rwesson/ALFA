program cont
use mod_quicksort
implicit none

character(len=50) :: filename
character(len=1) :: z
integer :: i, io, spectrumlength, begin, end, window
type spectrum
  real(kind=dp) :: wavelength
  real(kind=dp) :: flux
  real(kind=dp) :: difference
  real(kind=dp) :: continuum
end type
type(spectrum), dimension(:), allocatable :: input
type(spectrum), dimension(101) :: chunk
real(kind=dp), dimension(:), allocatable :: temporary
real(kind=dp) :: median

window=25

!take first differences of each point to points either side
!distribution of these should be mostly continuum with noise, with large values at lines
!select the lowest ones, smooth somehow

call get_command_argument(1,filename)
if (filename .eq. "") then
  print *,"file?"
  read(5,*) filename
endif

I = 0
OPEN(199, file=filename, iostat=IO, status='old')
DO WHILE (IO >= 0)
  READ(199,*,end=112) z
  I = I + 1
END DO
112 spectrumlength=I

!then allocate and read

allocate (input(spectrumlength))

REWIND (199)
DO I=1,spectrumlength
  READ(199,*) input(i)%wavelength,input(i)%flux
  input(i)%continuum = 0.d0
END DO
CLOSE(199)

!now get first differences

do i=2,spectrumlength-1
  input(i)%difference = abs(input(i)%flux - input(i-1)%flux) + abs(input(i+1)%flux - input(i)%flux)
end do

!in 100 unit windows, take continuum as average of points with less than median 1st difference

do i=1,spectrumlength

  if (i-window .lt. 1) then
    begin=1
  else
    begin=i-window
  endif

  if (i+window .gt. spectrumlength) then
    end=spectrumlength
  else
    end=i+window
  endif

  chunk = input(begin:end)

  allocate(temporary(end-begin+1))
  temporary=chunk%difference
  call qsort(temporary)
  median=temporary(int(real(size(temporary)/2)))

  deallocate(temporary)
  where (chunk%difference .gt. median)
    chunk%continuum=0.d0
  elsewhere
    chunk%continuum=chunk%flux
  endwhere

  input(i)%continuum=sum(chunk%continuum)/count(chunk%continuum.gt.0.d0)

enddo

open(123, file="continuum.dat")
do i=1,spectrumlength
  write (123,*) input(i)%wavelength,input(i)%flux,input(i)%continuum
enddo

end program cont
