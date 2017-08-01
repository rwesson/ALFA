!Copyright (C) 2013- Roger Wesson
!Free under the terms of the GNU General Public License v3

module mod_readfiles
use mod_types
use mod_routines

contains

subroutine readdata(spectrumfile, spectrum_1d, spectrum_2d, spectrum_3d, wavelengths, wavelengthscaling, axes, rebinfactor)
!take the filename, check if it's FITS or plain text
!if FITS, then read the necessary keywords to set the wavelength scale, allocate the data array according to the number of dimensions found, and fill it
!if plain text, read two columns for wavelength and flux, return.

  implicit none
  character (len=*), intent(in) :: spectrumfile !input file name
  real, dimension(:), allocatable :: wavelengths !wavelength array

  real, dimension(:), allocatable :: spectrum_1d !array for 1d data
  real, dimension(:,:), allocatable :: spectrum_2d !array for 2d data
  real, dimension(:,:,:), allocatable :: spectrum_3d !array for 3d data

  character(len=1) :: checkrow ! for use in ignoring comment rows
  character(len=512) :: rowdata ! rows read into this variable then split into wavelength and flux

  real :: wavelength, dispersion, referencepixel
  real :: wavelengthscaling !factor to convert wavelengths into Angstroms
  logical :: loglambda !is the spectrum logarithmically sampled?
  integer :: dimensions !number of dimensions
  integer, dimension(:), allocatable :: axes !number of pixels in each dimension
  integer :: i, io !counter and io status for file reading
  integer, intent(in) :: rebinfactor

  !cfitsio variables

  integer :: status,unit,readwrite,blocksize,hdutype,group,numberofhdus,nrows
  character(len=80) :: key_cunit, key_ctype, key_crpix, key_crval, key_cdelt, key_cd
  character(len=80) :: cunit,ctype
  real :: nullval
  logical :: anynull,datafound

#ifdef CO
  print *,"subroutine: readdata"
#endif

  datafound=.false.

  cunit=""
  ctype=""

  !is the file a FITS file?
  !if it contains the string .fit or .FIT, assume that it is.

  if (index(spectrumfile,".fit").gt.2 .or. index(spectrumfile,".FIT").gt.2) then !filename contains at least 1 character followed by .fit or .FIT
    print *,gettime(),"this looks like a FITS file (.fit or .FIT appear in the the filename)"
    status=0
    !  Get an unused Logical Unit Number to use to open the FITS file.
    call ftgiou(unit,status)
    !  Open the FITS file
    readwrite=0
    call ftopen(unit,trim(spectrumfile),readwrite,blocksize,status)

    if (status .ne. 0) then
      print *,gettime(),"error: couldn't open FITS file ",trim(spectrumfile),". CFITSIO error code was ",status
      call exit(1)
    endif

    ! find number of HDUs

    call ftthdu(unit,numberofhdus,status)

    print *,gettime(),"FITS file has ",numberofhdus," extensions"

    ! get number of axes
    dimensions=0

    ! examine HDUs. if dimensions are found, treat as image. otherwise, look for table data

    do i=1,numberofhdus
      call ftghdt(unit,hdutype,status) ! get header type
      if (hdutype.eq.0) then
        print *,gettime(),"extension ",i," contains image data."
        call ftgidm(unit,dimensions,status)
        if (dimensions .ne. 0) then !extension has dimensions, now check if they look like actual data
          allocate(axes(dimensions))
          call ftgisz(unit,dimensions,axes,status)
          if (any(axes.gt.1)) then !at least one axis has more than one data point.
            print *,"reading in data"
            exit ! found dimensions in this HDU so leave the loop
          endif
          deallocate(axes)
        endif
        print *,gettime(),"no axes found, trying next extension"
      endif

      if (hdutype.eq.1.or.hdutype.eq.2) then
        if (hdutype.eq.1) print *,gettime(),"extension ",i," contains an ASCII table"
        if (hdutype.eq.2) print *,gettime(),"extension ",i," contains a binary table"
        call ftgnrw(unit,nrows,status) !
        print *,gettime(),"table contains ",nrows," rows"

! check for crappy ESO format which puts everything in one row
! look for header keyword nelem

        if (nrows.eq.1) then
          call ftgkyj(unit,"NELEM",nrows,"",status)
          print *,gettime(),"only one row so looked for NELEM keyword instead and found ",nrows," rows"
        endif

        allocate(spectrum_1d(nrows))
        allocate(wavelengths(nrows))
! assume first column is wavelength, second is flux. todo: might not always be so. either intelligently handle it within code or allow user to specify.
        call ftgcve(unit,1,1,1,nrows,0.,wavelengths,anynull,status)
        call ftgcve(unit,2,1,1,nrows,0.,spectrum_1d,anynull,status)
        datafound=.true.

    ! get units of wavelength
    ! current assumption is it will be A or nm

        if (wavelengthscaling .ne. 1.d0) then
          print *,gettime(),"  wavelength units: set by user. Angstroms per wavelength unit = ",wavelengthscaling
        else
          key_cunit="TUNIT1"
          call ftgkys(unit,key_cunit,cunit,"",status)
          if (trim(cunit) .eq. "nm" .or. trim(cunit) .eq. "NM") then
            wavelengthscaling=10.d0 ! convert to Angstroms if it's in nm
            print *,gettime(),"  wavelength units: nm.  Will convert to A."
          elseif (trim(cunit).eq."Angstrom" .or. trim(cunit).eq."Angstroms") then
            print *,gettime(),"  wavelength units: Angstroms"
            wavelengthscaling = 1.d0
          else
           print *,gettime(),"  wavelength units: not recognised - will assume A.  Set the --wavelength-scaling if this is not correct"
            wavelengthscaling = 1.d0
          endif
        endif

        goto 888 ! skip the subsequent file reading routines. todo: change this to a part of the if condition
        exit
      endif

      call ftmrhd(unit,1,hdutype,status) ! otherwise, advance to next HDU and search there
    enddo

    if (dimensions .eq. 0 .and. .not. datafound) then ! no axes found, no table data read
      print *,gettime(),"error : no axes found in ",trim(spectrumfile)
      print *,gettime(),"        (number of extensions searched: ",numberofhdus,")"
      call exit(1)
    elseif (dimensions .gt. 3) then ! can't imagine what a 4D fits file would actually be, but alfa definitely can't handle it
      print *,gettime(),"error : more than 3 axes found in ",trim(spectrumfile)
      call exit(1)
    endif

    print *,gettime(),"  number of dimensions: ",dimensions

    ! now get the dimensions of the axis

!    allocate(axes(dimensions)) already allocated
    call ftgisz(unit,dimensions,axes,status)

    ! set up array for wavelengths

    if (dimensions.eq.3) then
      allocate(wavelengths(axes(3)))
    else
      allocate(wavelengths(axes(1)))
    endif

    ! get wavelength, dispersion and reference pixel
    ! set which FITS keywords we are looking for depending on which axis represents wavelength
    ! this will be axis 3 for cubes, axis 1 otherwise

    status=0

    if (dimensions .lt. 3) then

      key_crval="CRVAL1"
      key_crpix="CRPIX1"
      key_ctype="CTYPE1"
      key_cunit="CUNIT1"
      key_cdelt="CDELT1"
      key_cd   ="CD1_1"

    else

      key_crval="CRVAL3"
      key_crpix="CRPIX3"
      key_ctype="CTYPE3"
      key_cunit="CUNIT3"
      key_cdelt="CDELT3"
      key_cd   ="CD3_3"

    endif

    call ftgkye(unit,key_crval,wavelength,"",status)
    if (status .ne. 0) then
      print *,gettime(),"error: couldn't find wavelength value at reference pixel - no keyword ",trim(key_crval),"."
      call exit(1)
    endif

    print *,gettime(),"  wavelength at reference pixel: ",wavelength

    call ftgkye(unit,key_crpix,referencepixel,"",status)
    if (status .ne. 0) then
      print *,gettime(),"warning: couldn't find reference pixel - no keyword ",trim(key_crpix),". Setting to 1.0"
      referencepixel=1.0
      status=0
    endif

    print *,gettime(),"  reference pixel: ",referencepixel

    call ftgkye(unit,key_cdelt,dispersion,"",status)
    if (status.ne.0) then
      status=0
      call ftgkye(unit,key_cd,dispersion,"",status)
        if (status .ne. 0) then
          print *,gettime(),"error: couldn't find wavelength dispersion - no keyword ",trim(key_cdelt)," or ",trim(key_cd),"."
          call exit(1)
        endif
    endif

    print *,gettime(),"  wavelength dispersion: ",dispersion

    ! check if the wavelength axis is log-sampled
    call ftgkey(unit,key_ctype,ctype,"",status)
    if (index(ctype,"-LOG").gt.0) then
      loglambda = .true.
      print *,gettime(),"  sampling: logarithmic"
    else
      loglambda = .false.
      print *,gettime(),"  sampling: linear"
    endif

    ! get units of wavelength
    ! current assumption is it will be A or nm

    if (wavelengthscaling .ne. 1.d0) then
      print *,gettime(),"  wavelength units: set by user. Angstroms per wavelength unit = ",wavelengthscaling
    else
      call ftgkey(unit,key_cunit,cunit,"",status)
      if (trim(cunit) .eq. "nm" .or. trim(cunit) .eq. "NM") then
        wavelengthscaling=10.d0 ! convert to Angstroms if it's in nm
        print *,gettime(),"  wavelength units: nm.  Will convert to A."
      elseif (trim(cunit).eq."Angstrom" .or. trim(cunit).eq."Angstroms") then
        print *,gettime(),"  wavelength units: Angstroms"
        wavelengthscaling = 1.d0
      else
        print *,gettime(),"  wavelength units: not recognised - will assume A.  Set the --wavelength-scaling if this is not correct"
        wavelengthscaling = 1.d0
      endif
    endif

    !now we have the information we need to read in the data

    group=1
    nullval=-1.e-27

    if (dimensions.eq.1) then

      allocate(spectrum_1d(axes(1)))

      status = 0
      call ftgpve(unit,group,1,axes(1),nullval,spectrum_1d,anynull,status)
      !todo: report null values?
      if (status .eq. 0) then
        print "(X,A,A,I7,A)",gettime(),"read 1D fits file with ",axes(1)," data points into memory."
      else
        print *,gettime(),"couldn't read file into memory"
        call exit(1)
      endif

    elseif (dimensions.eq.2) then

      allocate(spectrum_2d(axes(1),axes(2)))

      status=0
      call ftg2de(unit,group,nullval,axes(1),axes(1),axes(2),spectrum_2d,anynull,status)
    !todo: report null values?
      if (status .eq. 0) then
        print "(X,A,A,I7,A)",gettime(),"read ",axes(2)," rows into memory."
      else
        print *,gettime(),"couldn't read RSS file into memory"
        print *,"error code ",status
        call exit(1)
      endif

    elseif (dimensions.eq.3) then

      allocate(spectrum_3d(axes(1),axes(2),axes(3)))

      status=0
      call ftg3de(unit,group,nullval,axes(1),axes(2),axes(1),axes(2),axes(3),spectrum_3d,anynull,status)
    !todo: report null values?
      if (status .eq. 0) then
        print "(X,A,A,I7,A)",gettime(),"read ",axes(1)*axes(2)," pixels into memory."
      else
        print *,gettime(),"couldn't read cube into memory"
        call exit(1)
      endif

    else

      print *,gettime(),"More than 3 dimensions.  ALFA cannot comprehend that yet, sorry."
      call exit(1)

    endif

    ! calculate wavelength array

    if (loglambda) then !log-sampled case
      do i=1,size(wavelengths)
        wavelengths(i) = wavelength*exp((i-referencepixel)*dispersion/wavelength)
      enddo
    else !linear case
      do i=1,size(wavelengths) 
        wavelengths(i) = (wavelength+(i-referencepixel)*dispersion)
      enddo
    endif

  else ! end of FITS file loop. if we are here, assume file is 1D ascii with wavelength and flux columns.

    print *,gettime(),"filename does not contain .fit or .FIT, so assuming plain text format"

    allocate(axes(1))

    !get number of lines

    i = 0
    open(199, file=spectrumfile, iostat=IO, status='old')
      do while (IO >= 0)
      read(199,*,end=112) checkrow
      if (checkrow .ne. "#") i = i + 1
    enddo
    112 axes(1) = i
    print *,gettime(),"  number of data points: ",axes(1)

    !then allocate and read

    allocate(spectrum_1d(axes(1)))
    allocate(wavelengths(axes(1)))

    rewind (199)
    i=1
    do while (i<=axes(1))
      read(199,"(A)") rowdata
      if (index(rowdata,"#") .ne. 1) then
        read (rowdata,*) wavelengths(i), spectrum_1d(i)
        i=i+1
      endif
    enddo
    close(199)

  endif

888 continue ! fix this ugly code at some point

! wavelength scaling

  if (wavelengthscaling .eq. 0.d0) then
    print *,gettime(),"  wavelength units: assumed to be Angstroms"
    wavelengthscaling = 1.d0
  else
    print *,gettime(),"  wavelength units: set by user. Angstroms per wavelength unit = ",wavelengthscaling
  endif

  wavelengths = wavelengths * wavelengthscaling

  print *,gettime(),"wavelength range: ",wavelengths(1),wavelengths(size(wavelengths))
  print *

! rebinning

  if (rebinfactor .gt. 1) then
    if (allocated(spectrum_1d)) then
      print *,gettime(),"rebinned spectrum by factor of ",rebinfactor
      call rebinarray(wavelengths,rebinfactor)
      call rebinarray(spectrum_1d,rebinfactor)
    else
      print *,gettime(),"warning: can't rebin 2d or 3d spectra yet, left as is"
    endif
  endif

! todo: make this work for 2d and 3d data

end subroutine readdata

subroutine rebinarray(array,rebinfactor)

  implicit none
  real, dimension(:), allocatable :: array, array_temp
  integer :: rebinfactor, i, newsize, endbit, s1, s2

  newsize=size(array)/rebinfactor
  endbit=mod(size(array),rebinfactor)

  if (endbit .gt. 0) then
    newsize=newsize+1
  endif

  allocate(array_temp(newsize))

  do i = 1,newsize
    s1=(i-1)*rebinfactor + 1
    s2=i*rebinfactor
    array_temp(i)=sum(array(s1:s2)) / rebinfactor
  enddo

  if (endbit .ne. 0) then ! fill in last element
    array_temp(newsize) = sum(array(size(array)-endbit+1:size(array)))/endbit
  endif

  deallocate(array)
  allocate(array(newsize))
  array=array_temp

end subroutine rebinarray

subroutine readlinelist(linelistfile,referencelinelist,nlines,wavelength1,wavelength2,exclusions)
!this subroutine reads in the line catalogue
! - linelistfile is the name of the line catalogue file
! - referencelinelist is the array into which the values are read
! - nlines is the number of lines successfully read into the array
! - wavelength1 and wavelength2 define the range of wavelenths to be read in
! - exclusions is an array of wavelengths to be skipped if found in the input catalogue

  implicit none
  character (len=512) :: linelistfile
  character (len=85) :: linedatainput
  character (len=2) :: informatnumber
  character (len=20) :: informat
  integer :: i
  real :: input1, wavelength1, wavelength2
  integer :: io, nlines
  logical :: file_exists
  real, dimension(:) :: exclusions

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

  if (.not. file_exists) then ! try in default directory

    inquire(file=PREFIX//"/share/alfa/"//linelistfile, exist=file_exists)
    if (.not. file_exists) then
      print *,gettime(),"error: line catalogue not found: ",trim(linelistfile)," does not exist in current directory or in ",PREFIX,"/share/alfa"
      call exit(1)
    else
      linelistfile=PREFIX//"/share/alfa/"//trim(linelistfile)
    endif

  endif

  I = 0
  OPEN(199, file=linelistfile, iostat=IO, status='old')
  DO WHILE (IO >= 0)
    READ(199,*,end=110) input1
    if (input1 .ge. wavelength1 .and. input1 .le.  wavelength2 .and. .not. (any(exclusions.eq.input1))) then
    !only read in lines that lie within the observed wavelength range, and are not in the line exclusions array
      I = I + 1
    endif
  END DO
  110     nlines=I

!determine the format to use. ceiling of log gives number of sig figs before dp, add 3 for dp and 2 sig figs after

  write (informatnumber,"(I2)") ceiling(log10(max(wavelength1,wavelength2)))+3
  informat="(F"//trim(adjustl(informatnumber))//".2,A)"

!then allocate and read

  allocate(referencelinelist(nlines))

  REWIND (199)
  I=1
  do while (i .le. nlines)
    READ(199,informat) input1, linedatainput
    if (input1 .ge. wavelength1 .and. input1 .le. wavelength2 .and. .not. (any(exclusions.eq.input1))) then
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
