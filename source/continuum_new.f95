program cont
use mod_quicksort
implicit none

character(len=50) :: filename
character(len=1) :: z
integer :: i, j, k, io, spectrumlength, begin, end
type spectrum
  real(kind=dp) :: wavelength
  real(kind=dp) :: flux
  real(kind=dp) :: difference
  real(kind=dp) :: continuum
  real(kind=dp) :: continuum_smoothed
end type
type(spectrum), dimension(:), allocatable :: input
type(spectrum), dimension(101) :: chunk
real(kind=dp), dimension(:), allocatable :: temporary
real(kind=dp) :: median, nzbefore, nzafter
real(kind=dp), dimension(101) :: sgcoeffs

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

!find median

allocate(temporary(spectrumlength))
temporary=input%difference
call qsort(temporary)
median=temporary(int(real(size(temporary)/4)))
deallocate(temporary)

!continuum fluxes are defined as any where the first difference is less than the median
!fill in the non-continuum fluxes with the previous value
!this should be replaced with proper interpolation

nzbefore=input(1)%flux
i=1

do while (i .lt. spectrumlength)
  if (input(i)%difference .gt. median) then
!find the next place where it's below the median
    do j=1,spectrumlength-i
      if (input(i+j)%difference .lt. median) then
        nzafter=input(i+j)%flux
        exit
      endif
    enddo
!now we have the points we can linearly interpolate
    do k=1,j
      input(i+k-1)%continuum=nzbefore+(nzafter-nzbefore)*dble(k)/dble(j)
!input (i+k-1)%continuum=0.d0
    end do
    i=i+k !skip ahead to next continuum point
  else
    input(i)%continuum=input(i)%flux
    nzbefore=input(i)%flux
    i=i+1
  endif
enddo

!next, smooth the data

!savitsky-golay in 101 unit window
!get coefficients. sum of coeffs=1

!do i=-50,50
!  sgcoeffs(i+51)=(3*101**2-7-(20*dble(i)**2))/(101*(101**2-4)) * (3./4.)
!enddo
!print *,sgcoeffs
!do i=51,spectrumlength-50
!  chunk=input(i-50:i+50)
!  input(i)%continuum_smoothed=sum(sgcoeffs*chunk%continuum)
!enddo

!least squares fit to chunk, use value at centre point as smoothed continuum

!average in 101-unit window excluding zero values
!ends take the average in a truncated window

do i=1,spectrumlength

  if (i-50 .lt. 1) then
    begin=1
  else
    begin=i-50
  endif

  if (i+50 .gt. spectrumlength) then
    end=spectrumlength
  else
    end=i+50
  endif

  chunk%continuum=0.d0
  chunk = input(begin:end)
  input(i)%continuum_smoothed=sum(chunk%continuum)/count(chunk%continuum.gt.0.d0)
enddo

open(123, file="continuum.dat")
do i=1,spectrumlength
  write (123,*) input(i)%wavelength,input(i)%flux,input(i)%continuum,input(i)%continuum_smoothed
enddo

end program cont
