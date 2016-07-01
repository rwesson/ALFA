!Copyright (C) 2013- Roger Wesson
!Free under the terms of the GNU General Public License v3

module mod_readfiles
use mod_types
use mod_routines

contains

subroutine getfiletype(spectrumfile, filetype, dimensions, axes, wavelength, dispersion, referencepixel)
!this subroutine determines what type of file the input file is
!input is the file name
!output is the file type, indicated by:
! 1 - 1d fits
! 2 - 2d fits
! 3 - 3dfits
! 4 - 1d ascii
!dimensions and lengths of axes are returned if it's a FITS file

  implicit none
  character (len=*) :: spectrumfile
  integer :: filetype, dimensions, referencepixel
  integer, dimension(:), allocatable :: axes

  !cfitsio variables

  integer :: status,unit,readwrite,blocksize,hdutype
  real :: wavelength, dispersion

#ifdef CO
  print *,"subroutine: getfiletype"
#endif

  referencepixel=1

  !check if it's a fits file

  if (index(spectrumfile,".fit").gt.0 .or. index(spectrumfile,".FIT").gt.0) then !read header

    status=0
    !  Get an unused Logical Unit Number to use to open the FITS file.
    call ftgiou(unit,status)
    !  Open the FITS file

    readwrite=0
    call ftopen(unit,spectrumfile,readwrite,blocksize,status)

    ! get number of axes
    dimensions=0
    call ftgidm(unit,dimensions,status)
    do while (dimensions .eq. 0) ! if no axes found in first extension, advance and check again
      call ftmrhd(unit,1,hdutype,status)
      call ftgidm(unit,dimensions,status)
    enddo
    if (dimensions .eq. 0) then ! still no axes found
      print *,gettime(),"error : no axes found in ",trim(spectrumfile)
      call exit(1)
    elseif (dimensions .gt. 3) then ! can't imagine what a 4D fits file would actually be, but alfa definitely can't handle it
      print *,gettime(),"error : more than 3 axes found in ",trim(spectrumfile)
      call exit(1)
    endif

    ! now get the dimensions of the axis

    allocate(axes(dimensions))
    call ftgisz(unit,dimensions,axes,status)
    filetype=dimensions

    ! get wavelength and dispersion
    ! todo: make this more robust

    if (dimensions .lt. 3) then
      call ftgkye(unit,"CRVAL1",wavelength,"",status)
      call ftgkyj(unit,"CRPIX1",referencepixel,"",status)
      call ftgkye(unit,"CDELT1",dispersion,"",status)
      if (status.ne.0) then
        status=0
        call ftgkye(unit,"CD1_1",dispersion,"",status)
      endif
    else
      call ftgkye(unit,"CRVAL3",wavelength,"",status)
      call ftgkyj(unit,"CRPIX3",referencepixel,"",status)
      call ftgkye(unit,"CDELT3",dispersion,"",status)
      if (status.ne.0) then
        status=0
        call ftgkye(unit,"CD3_3",dispersion,"",status)
      endif
    endif

    ! close file

    call ftclos(unit, status)
    call ftfiou(unit, status)

  else ! not FITS file, assume ascii

    filetype=4

  endif

end subroutine getfiletype

subroutine readascii(spectrumfile, realspec, spectrumlength, fittedspectrum)
! read in a plain text file

  implicit none
  character (len=512) :: spectrumfile
  integer :: i
  real :: input1, input2
  integer :: io, spectrumlength

  type(spectrum), dimension(:), allocatable :: realspec, fittedspectrum

#ifdef CO
  print *,"subroutine: readascii"
#endif

  i = 0
  open(199, file=spectrumfile, iostat=IO, status='old')
    do while (IO >= 0)
    read(199,*,end=112) input1
    i = i + 1
  enddo
  112 spectrumlength=i

  !then allocate and read

  allocate (realspec(spectrumlength))

  rewind (199)
  do i=1,spectrumlength
    read(199,*) input1, input2
    realspec(i)%wavelength = input1
    realspec(i)%flux = input2
  enddo
  close(199)

  realspec%uncertainty = 0.d0

  allocate (fittedspectrum(spectrumlength))
  fittedspectrum%wavelength=realspec%wavelength
  fittedspectrum%flux=0.d0

end subroutine readascii

subroutine read1dfits(spectrumfile, realspec, spectrumlength, fittedspectrum, wavelength, dispersion, referencepixel)
! read in a 1D fits file

  implicit none
  character (len=512) :: spectrumfile
  integer :: i
  integer :: spectrumlength, referencepixel

  type(spectrum), dimension(:), allocatable :: realspec, fittedspectrum

  !cfitsio variables

  integer :: status,unit,readwrite,blocksize
  integer :: group
  real :: nullval
  logical :: anynull
  real :: wavelength, dispersion

#ifdef CO
  print *,"subroutine: read1dfits"
#endif

  allocate(realspec(spectrumlength))
  allocate(fittedspectrum(spectrumlength))

  status=0
  !  Get an unused Logical Unit Number to use to open the FITS file.
  call ftgiou(unit,status)
  !  Open the FITS file

  readwrite=0
  call ftopen(unit,spectrumfile,readwrite,blocksize,status)

  group=1
  nullval=-999

  ! calculate wavelength values

  do i=1,spectrumlength
    realspec(i)%wavelength = wavelength+(i-referencepixel)*dispersion
  enddo

  ! read spectrum into memory

  call ftgpve(unit,group,1,spectrumlength,nullval,realspec%flux,anynull,status)
!todo: report null values?
  if (status .eq. 0) then
    print "(X,A,A,I7,A)",gettime(),"read 1D fits file with ",spectrumlength," data points into memory."
  else
    print *,gettime(),"couldn't read file into memory"
    call exit(1)
  endif

  ! close file

  call ftclos(unit, status)
  call ftfiou(unit, status)

  realspec%uncertainty = 0.d0
  fittedspectrum%wavelength=realspec%wavelength
  fittedspectrum%flux=0.d0

end subroutine read1dfits

subroutine read2dfits(spectrumfile, rssdata, dimensions, axes)
!read a 2D FITS file.

  implicit none
  character(len=512) :: spectrumfile
  real, dimension(:,:), allocatable :: rssdata
  logical :: anynull
  integer :: alloc_err
  integer :: dimensions
  integer, dimension(:) :: axes
  !cfitsio variables

  integer :: status,unit,readwrite,blocksize,hdutype
  integer :: group
  real :: nullval

#ifdef CO
  print *,"subroutine: read2dfits"
#endif

  status=0
  !  Get an unused Logical Unit Number to use to open the FITS file.
  call ftgiou(unit,status)
  !  Open the FITS file
  readwrite=0
  call ftopen(unit,spectrumfile,readwrite,blocksize,status)

  group=1
  nullval=-999

  allocate(rssdata(axes(1),axes(2)), stat=alloc_err)
  if (alloc_err .eq. 0) print *,gettime(),"reading RSS file into memory"

! advance to first HDU containing axes

  dimensions=0
  call ftgidm(unit,dimensions,status)
  do while (dimensions .eq. 0) ! if no axes found in first extension, advance andcheck again
    call ftmrhd(unit,1,hdutype,status)
    call ftgidm(unit,dimensions,status)
  enddo

! read RSS file into memory

  status=0
  call ftg2de(unit,group,nullval,axes(1),axes(1),axes(2),rssdata,anynull,status)
!todo: report null values?
  if (status .eq. 0) then
    print "(X,A,A,I7,A)",gettime(),"read ",axes(2)," rows into memory."
  else
    print *,gettime(),"couldn't read RSS file into memory"
    print *,"error code ",status
    call exit(1)
  endif

end subroutine read2dfits

subroutine read3dfits(spectrumfile, cubedata, dimensions, axes)
!read a FITS cube.

  implicit none
  character(len=512) :: spectrumfile
  real, dimension(:,:,:), allocatable :: cubedata
  logical :: anynull
  integer :: alloc_err
  integer :: dimensions
  integer, dimension(:) :: axes
  !cfitsio variables

  integer :: status,unit,readwrite,blocksize,hdutype
  integer :: group
  real :: nullval

#ifdef CO
  print *,"subroutine: read3dfits"
#endif

  status=0
  !  Get an unused Logical Unit Number to use to open the FITS file.
  call ftgiou(unit,status)
  !  Open the FITS file
  readwrite=0
  call ftopen(unit,spectrumfile,readwrite,blocksize,status)

  group=1
  nullval=-999

  allocate(cubedata(axes(1),axes(2),axes(3)), stat=alloc_err)
  if (alloc_err .eq. 0) print *,gettime(),"reading data cube into memory"

! advance to first HDU containing axes

  dimensions=0
  call ftgidm(unit,dimensions,status)
  do while (dimensions .eq. 0) ! if no axes found in first extension, advance andcheck again
    call ftmrhd(unit,1,hdutype,status)
    call ftgidm(unit,dimensions,status)
  enddo

! read cube file into memory

  status=0
  call ftg3de(unit,group,nullval,axes(1),axes(2),axes(1),axes(2),axes(3),cubedata,anynull,status)
!todo: report null values?
  if (status .eq. 0) then
    print "(X,A,A,I7,A)",gettime(),"read ",axes(1)*axes(2)," pixels into memory."
  else
    print *,gettime(),"couldn't read cube into memory"
    call exit(1)
  endif

end subroutine read3dfits

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

#ifdef CO
  print *,"subroutine: readlinelist"
#endif

  ! deallocate if necessary
  if (allocated(referencelinelist)) deallocate(referencelinelist)

  if (trim(linelistfile)=="") then
    print *,gettime(),"error: No line catalogue specified"
    call exit(1)
  endif

  inquire(file=linelistfile, exist=file_exists) ! see if the input file is present

  if (.not. file_exists) then
    print *,gettime(),"error: line catalogue ",trim(linelistfile)," does not exist"
    call exit(1)
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

#ifdef CO
  print *,"subroutine: selectlines: ",wavelength1,wavelength2
#endif

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
