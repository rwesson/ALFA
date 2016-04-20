module mod_readfiles
use mod_types
use mod_routines

contains

subroutine getfiletype(spectrumfile, filetype, dimensions, axes)
!this subroutine determines what type of file the input file is
!input is the file name
!output is the file type, indicated by:
! 1 - 1d fits
! 2 - 2d fits
! 3 - 3dfits
! 4 - 1d ascii
!dimensions and lengths of axes are returned if it's a FITS file

  implicit none
  character (len=512) :: spectrumfile
  integer :: filetype, dimensions
  integer, dimension(:), allocatable :: axes
  logical :: file_exists

  !cfitsio variables

  integer :: status,unit,readwrite,blocksize,hdutype
  integer :: group,firstpix
  real :: nullval
  logical :: anynull
  integer :: alloc_err
  real :: wavelength, dispersion

  !read in spectrum to fit

  if (trim(spectrumfile)=="") then
    print *,gettime(),": error: No input spectrum specified"
    stop
  endif

  inquire(file=spectrumfile, exist=file_exists) ! see if the input file is present

  if (.not. file_exists) then
    print *,gettime(),": error: input spectrum ",trim(spectrumfile)," does not exist"
    stop
  endif

  if (index(spectrumfile,".fit").gt.0 .or. index(spectrumfile,".FIT").gt.0) then

    status=0
    !  Get an unused Logical Unit Number to use to open the FITS file.
    call ftgiou(unit,status)
    !  Open the FITS file

    readwrite=0
    call ftopen(unit,spectrumfile,readwrite,blocksize,status)

    ! get number of axes
    dimensions=0
    call ftgidm(unit,dimensions,status)
    do while (dimensions .eq. 0) ! if no axes found in first extension, advance andcheck again
      call ftmrhd(unit,1,hdutype,status)
      call ftgidm(unit,dimensions,status)
    end do
    if (dimensions .eq. 0) then ! still no axes found
      print *,gettime(),": error : no axes found in ",trim(spectrumfile)
      stop
    elseif (dimensions .gt. 3) then ! can't imagine what a 4D fits file would actually be, but alfa definitely can't handle it
      print *,gettime(),": error : more than 3 axes found in ",trim(spectrumfile)
      stop

    endif

    ! now get the dimensions of the axis

    allocate(axes(dimensions))
    call ftgisz(unit,dimensions,axes,status)

  else ! not FITS file, assume ascii

    filetype=4

  endif

end subroutine getfiletype

subroutine readspectrum(spectrumfile, realspec, spectrumlength, fittedspectrum)

  implicit none
  character (len=512) :: spectrumfile
  integer :: i
  real :: input1, input2
  integer :: io, spectrumlength
  logical :: file_exists

  type(spectrum), dimension(:), allocatable :: realspec, fittedspectrum

  !cfitsio variables

  integer :: status,unit,readwrite,blocksize,nfound,hdutype
  integer, dimension(:), allocatable :: naxes
  integer :: group,firstpix
  real :: nullval
  logical :: anynull
  integer :: alloc_err
  real :: wavelength, dispersion

  !read in spectrum to fit

  if (trim(spectrumfile)=="") then
    print *,gettime(),": error: No input spectrum specified"
    stop
  endif

  inquire(file=spectrumfile, exist=file_exists) ! see if the input file is present

  if (.not. file_exists) then
    print *,gettime(),": error: input spectrum ",trim(spectrumfile)," does not exist"
    stop
  endif

  if (index(spectrumfile,".fit").gt.0 .or. index(spectrumfile,".FIT").gt.0) then

    status=0
    !  Get an unused Logical Unit Number to use to open the FITS file.
    call ftgiou(unit,status)
    !  Open the FITS file

    readwrite=0
    call ftopen(unit,spectrumfile,readwrite,blocksize,status)

    ! check we have 1 axis
    nfound=0
    call ftgidm(unit,nfound,status)

    do while (nfound .eq. 0) ! if no axes found in first extension, advance and check again
      call ftmrhd(unit,1,hdutype,status)
      call ftgidm(unit,nfound,status)
    end do

    if (nfound .eq. 0) then ! still no axes found
      print *,gettime(),": error : no axes found in ",trim(spectrumfile)
      stop
    endif

    ! now get the dimensions of the axis

    allocate(naxes(nfound))
    call ftgisz(unit,nfound,naxes,status)

    group=1
    firstpix=1
    nullval=-999

    spectrumlength=naxes(1)

    if (nfound .eq. 1) then
      allocate(realspec(spectrumlength), stat=alloc_err)
      if (alloc_err .eq. 0) print *,gettime(), ": reading file into memory"
    elseif (nfound .eq. 2) then
      print *,gettime(),": 2D FITS file : please use alfarss"
      print *,gettime(),": dimensions seem to be ",naxes(1)," x ",naxes(2)
      stop
    elseif (nfound .eq. 3) then
      print *,gettime(),": 3D FITS file : please use alfacube"
      print *,gettime(),": dimensions seem to be ",naxes(1)," x ",naxes(2)," x ",naxes(3)
      stop
    endif

    ! find wavelength dispersion

    call ftgkye(unit,"CRVAL1",wavelength,"",status)
    call ftgkye(unit,"CD1_1",dispersion,"",status)

    ! calculate wavelength values

    do i=1,spectrumlength
      realspec(i)%wavelength = wavelength+(i-1)*dispersion
    enddo

    ! read spectrum into memory

    call ftgpve(unit,group,firstpix,spectrumlength,nullval,realspec%flux,anynull,status)

    if (status .eq. 0) then
      print "(X,A,A,I7,A)",gettime(), ": successfully read ",spectrumlength," pixels into memory"
    else
      print *,gettime(), ": couldn't read file into memory"
      stop
    endif

    ! close file

    call ftclos(unit, status)
    call ftfiou(unit, status)

  else ! read plain text

    I = 0
    OPEN(199, file=spectrumfile, iostat=IO, status='old')
      DO WHILE (IO >= 0)
      READ(199,*,end=112) input1
      I = I + 1
    END DO
    112 spectrumlength=I

    !then allocate and read

    allocate (realspec(spectrumlength))

    REWIND (199)
    DO I=1,spectrumlength
      READ(199,*) input1, input2
      realspec(i)%wavelength = input1
      realspec(i)%flux = input2
    END DO
    CLOSE(199)

  endif

  realspec%uncertainty = 0.d0

  allocate (fittedspectrum(spectrumlength))
  fittedspectrum%wavelength=realspec%wavelength
  fittedspectrum%flux=0.d0

end subroutine readspectrum

subroutine readlinelist(linelistfile,referencelinelist,nlines,wavelength1, wavelength2)
!this subroutine reads in the line catalogue
! - linelistfile is the name of the line catalogue file
! - referencelinelist is the array into which the values are read
! - nlines is the number of lines successfully read into the array
! - wavelength1 and wavelength2 define the range of wavelenths to be read in

  implicit none
  character (len=512) :: linelistfile
  character (len=85) :: linedatainput
  integer :: i
  real :: input1, wavelength1, wavelength2
  integer :: io, nlines
  logical :: file_exists

  type(linelist), dimension(:), allocatable :: referencelinelist

  ! deallocate if necessary
  if (allocated(referencelinelist)) deallocate(referencelinelist)

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

subroutine selectlines(referencelinelist,wavelength1,wavelength2,fittedlines,nlines)
!creates the array fittedlines, which contains the lines from referencelinelist which lie between wavelength1 and wavelength2

  implicit none
  real :: wavelength1, wavelength2
  integer :: startloc,nlines

  type(linelist), dimension(:), allocatable :: referencelinelist, fittedlines

!deallocate if necessary
  if (allocated(fittedlines)) deallocate(fittedlines)

!then copy the relevant lines
  nlines=count(referencelinelist%wavelength.gt.wavelength1 .and. referencelinelist%wavelength.le.wavelength2)
  if (nlines .gt. 0) then
    allocate(fittedlines(nlines))
    startloc=minloc(referencelinelist%wavelength-wavelength1,1,referencelinelist%wavelength-wavelength1.gt.0)
    fittedlines=referencelinelist(startloc:startloc+nlines-1)
  endif

end subroutine selectlines
end module mod_readfiles
