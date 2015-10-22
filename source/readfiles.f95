module mod_readfiles
use mod_types
use mod_routines

contains

subroutine readspectrum(spectrumfile, realspec, spectrumlength, fittedspectrum)

  implicit none
  character (len=512) :: spectrumfile
  integer :: i
  real :: input1, input2
  integer :: io, spectrumlength
  logical :: file_exists

  type(spectrum), dimension(:), allocatable :: realspec, fittedspectrum

  !read in spectrum to fit

  if (trim(spectrumfile)=="") then
    print *,gettime(),": error: No input spectrum specified"
    stop
  endif

  inquire(file=spectrumfile, exist=file_exists) ! see if the input file is present

  if (.not. file_exists) then
    print *,gettime(),": error: input spectrum ",trim(spectrumfile)," does not exist"
    stop
  else
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
    realspec(i)%uncertainty = 0.d0
  END DO
  CLOSE(199)

  fittedspectrum%wavelength=realspec%wavelength
  fittedspectrum%flux=0.d0

end subroutine readspectrum

subroutine readlinelist(linelistfile,referencelinelist,nlines,fittedlines, wavelength1, wavelength2)
!this subroutine reads in the line catalogue
! - linelistfile is the name of the line catalogue file
! - referencelinelist is the array into which the values are read
! - nlines is the number of lines successfully read into the array
! - fittedlines is an array with the line wavelengths, created for subsequent use in the fit subroutine
! - wavelength1 and wavelength2 are the shortest and longest wavelengths to be considered - lines outside them will not be read in from the catalogue.  realspec is unaffected.  This is so that lines being fitted are completely within the spectrum chunk
  implicit none
  character (len=512) :: linelistfile
  character (len=85) :: linedatainput
  integer :: i
  real :: input1, wavelength1, wavelength2
  integer :: io, nlines
  logical :: file_exists

  type(linelist), dimension(:), allocatable :: referencelinelist, fittedlines

  ! deallocate if necessary
  if (allocated(referencelinelist)) deallocate(referencelinelist)
  if (allocated(fittedlines)) deallocate(fittedlines)

  if (trim(linelistfile)=="") then
    print *,gettime(),": error: No line catalogue specified"
    stop
  endif

  inquire(file=linelistfile, exist=file_exists) ! see if the input file is present

  if (.not. file_exists) then
    print *,gettime(),": error: line catalogue ",trim(linelistfile)," does not exist"
    stop
  else
    I = 0
    OPEN(199, file=linelistfile, iostat=IO, status='old')
    DO WHILE (IO >= 0)
      READ(199,*,end=110) input1
      if (input1 .ge. wavelength1 .and. input1 .le.  wavelength2) then
      !only read in lines that lie within the observed wavelength range
        I = I + 1
      endif
    END DO
    110     nlines=I
  endif

!then allocate and read

  allocate(referencelinelist(nlines))
  allocate(fittedlines(nlines))

  REWIND (199)
  I=1
  do while (i .le. nlines)
    READ(199,'(F7.3,A)') input1, linedatainput
    if (input1 .ge. wavelength1 .and. input1 .le. wavelength2) then
      referencelinelist(i)%wavelength = input1
      referencelinelist(i)%peak=1000.
!formerly abs(realspec(minloc((realspec%wavelength-input1)**2,1))%flux) but in case of weak lines near to negative flux values this prevented them being fitted. it makes a negligible difference to the running time
      referencelinelist(i)%linedata = linedatainput
      i=i+1
    endif
    if (input1 .ge.  wavelength2) then
      exit
    endif
  END DO
  CLOSE(199)

end subroutine readlinelist
end module mod_readfiles
