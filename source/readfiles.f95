module mod_readfiles
use mod_types
use mod_routines

contains

subroutine readfiles(spectrumfile,linelistfile,realspec,referencelinelist,spectrumlength,nlines,linedata,fittedspectrum, fittedlines)
implicit none

character (len=512) :: spectrumfile, linelistfile
character (len=85) :: linedatainput
integer :: i
real :: input1, input2
integer :: io, nlines, spectrumlength
logical :: file_exists

type(linelist) :: referencelinelist, fittedlines
type(spectrum), dimension(:), allocatable :: realspec, fittedspectrum
character(len=85), dimension(:), allocatable :: linedata

! read in spectrum to fit

  if (trim(spectrumfile)=="") then
    print *,gettime(),": error: No input spectrum specified"
    stop
  endif

  inquire(file=spectrumfile, exist=file_exists) ! see if the input file is present

  if (.not. file_exists) then
    print *,gettime(),": error: input spectrum ",trim(spectrumfile)," does not exist"
    stop
  else
    print *,gettime(),": reading in spectrum ",trim(spectrumfile)
    I = 0
    OPEN(199, file=spectrumfile, iostat=IO, status='old')
      DO WHILE (IO >= 0)
      READ(199,*,end=112) input1
      I = I + 1
    END DO
    112 spectrumlength=I
  endif

  !then allocate and read

  allocate (realspec(spectrumlength))
  allocate (fittedspectrum(spectrumlength))

  REWIND (199)
  DO I=1,spectrumlength
    READ(199,*) input1, input2
    realspec(i)%wavelength = input1
    realspec(i)%flux = input2
  END DO
  CLOSE(199)

  if (trim(linelistfile)=="") then
    print *,gettime(),": error: No line catalogue specified"
    stop
  endif

  inquire(file=linelistfile, exist=file_exists) ! see if the input file is present

  if (.not. file_exists) then
    print *,gettime(),": error: line catalogue ",trim(linelistfile)," does not exist"
    stop
  else
    print *,gettime(),": reading in line catalogue ",trim(linelistfile)
    I = 0
    OPEN(199, file=linelistfile, iostat=IO, status='old')
    DO WHILE (IO >= 0)
      READ(199,*,end=110) input1
      if (input1 .ge. minval(realspec%wavelength) .and. input1 .le.  maxval(realspec%wavelength)) then
      !only read in lines that lie within the observed wavelength range
        I = I + 1
      endif
    END DO
    110     nlines=I
  endif

  if (nlines .eq. 0) then
    print *,gettime(),": error : Line catalogue does not overlap with input spectrum"
    stop
  endif

!then allocate and read

  allocate(referencelinelist%peak(nlines))
  allocate(referencelinelist%wavelength(nlines))
  allocate(fittedlines%peak(nlines))
  allocate(fittedlines%wavelength(nlines))
  allocate(linedata(nlines))

  REWIND (199)
  I=1
  do while (i .le. nlines)
    READ(199,'(F7.3,A)') input1, linedatainput
    if (input1 .ge. minval(realspec%wavelength)) then
      referencelinelist%wavelength(i) = input1
      ! seed the initial guess for line peaks with the nearest observed flux to
      ! the line centre
      referencelinelist%peak(i)=realspec(minloc((realspec%wavelength-input1)**2,1))%flux
      linedata(i) = linedatainput
      i=i+1
    endif
    if (input1 .ge.  maxval(realspec%wavelength)) then
      exit
    endif
  END DO
  CLOSE(199)

end subroutine readfiles
end module mod_readfiles
