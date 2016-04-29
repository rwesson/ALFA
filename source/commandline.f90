module mod_commandline
use mod_routines

contains

subroutine readcommandline(commandline,normalise,normalisation,redshiftguess,resolutionguess,vtol1,vtol2,rtol1,rtol2,baddata,pressure,spectrumfile,outputdirectory,skylinelistfile,stronglinelistfile,deeplinelistfile,generations,popsize,subtractsky,resolution_estimated,file_exists,imagesection)

  implicit none

  logical :: normalise
  real :: normalisation,redshiftguess,resolutionguess,vtol1,vtol2,rtol1,rtol2,baddata,pressure,c
  character(len=2048) :: commandline
  character(len=2048), dimension(:), allocatable :: options
  character(len=512),intent(out) :: spectrumfile,outputdirectory,skylinelistfile,stronglinelistfile,deeplinelistfile
  character(len=32) :: imagesection
  integer,intent(out) :: generations,popsize
  integer :: Narg,nargused,i
  logical,intent(out) :: subtractsky,resolution_estimated,file_exists

  c=299792.458 !km/s

  spectrumfile=""

  narg = 0
  nargused = 0 !to count options specified
  narg = IARGC() !count input arguments

  if (narg .eq. 0) then
    print *,"Usage: alfa [options] [file]"
    print *,"  [file] is an ascii file with columns for wavelength and flux"
    print *,"  or a FITS file with 1, 2 or 3 dimensions, containing spectra."
    print *,"  see the man page or online documentation for details of the options"
    stop
  endif

  call get_command(commandline)
  allocate (options(Narg))
  options=""
  print *,gettime(),": command line: ",trim(commandline)

  do i=1,Narg
    call get_command_argument(i,options(i))
  enddo

  do i=1,narg

    if ((trim(options(i))=="-n" .or. trim(options(i))=="--normalise")) then
      if ((i+1) .le. Narg) then
        read (options(i+1),*) normalisation
        normalise=.true.
        options(i:i+1)=""
      else
        print *,gettime(),": error: no value specified for ",trim(options(i))
        stop
      endif
    endif

    if ((trim(options(i))=="-vg" .or. trim(options(i))=="--velocity-guess")) then
      if ((i+1) .le. Narg) then
        read (options(i+1),*) redshiftguess
        options(i:i+1)=""
      else
        print *,gettime(),": error: no value specified for ",trim(options(i))
        stop
      endif
    endif

    if ((trim(options(i))=="-rg" .or. trim(options(i))=="--resolution-guess")) then
      if ((i+1) .le. Narg) then
        read (options(i+1),*) resolutionguess
        resolution_estimated=.true.
        options(i:i+1)=""
        if (resolutionguess .lt. 0.) then
          print *,gettime(),": error: invalid value given for resolution guess"
          stop
        endif
      else
        print *,gettime(),": error: no value specified for ",trim(options(i))
        stop
      endif
    endif

    if ((trim(options(i))=="-vtol1" .or. trim(options(i))=="--velocity-tolerance-1")) then
      if ((i+1) .le. Narg) then
        read (options(i+1),*) vtol1
        vtol1 = vtol1/c
        options(i:i+1)=""
        if (vtol1 .lt. 0.) then
          print *,gettime(),": error: invalid value given for vtol1"
          stop
        endif
      else
        print *,gettime(),": error: no value specified for ",trim(options(i))
        stop
      endif
    endif
    if ((trim(options(i))=="-vtol2" .or. trim(options(i))=="--velocity-tolerance-2")) then
      if ((i+1) .le. Narg) then
        read (options(i+1),*) vtol2
        vtol2 = vtol2/c
        options(i:i+1)=""
        if (vtol2 .lt. 0.) then
          print *,gettime(),": error: invalid value given for vtol2"
          stop
        endif
      else
        print *,gettime(),": error: no value specified for ",trim(options(i))
        stop
      endif
    endif

    if ((trim(options(i))=="-rtol1" .or. trim(options(i))=="--resolution-tolerance-1")) then
      if ((i+1) .le. Narg) then
        read (options(i+1),*) rtol1
        options(i:i+1)=""
        if (rtol1 .lt. 0.) then
          print *,gettime(),": error: invalid value given for rtol1"
          stop
        endif
      else
        print *,gettime(),": error: no value specified for ",trim(options(i))
        stop
      endif
    endif

    if ((trim(options(i))=="-rtol2" .or. trim(options(i))=="--resolution-tolerance-2")) then
      if ((i+1) .le. Narg) then
        read (options(i+1),*) rtol2
        options(i:i+1)=""
        if (rtol1 .lt. 0.) then
          print *,gettime(),": error: invalid value given for rtol1"
          stop
        endif
      else
        print *,gettime(),": error: no value specified for ",trim(options(i))
        stop
      endif
    endif

    if (trim(options(i))=="-ss" .or. trim(options(i))=="--subtract-sky") then
      subtractsky=.true.
      options(i)=""
    endif

    if (trim(options(i))=="-b" .or. trim(options(i))=="--bad-data") then
      if ((i+1) .le. Narg) then
        read (options(i+1),*) baddata
        options(i:i+1)=""
      else
        print *,gettime(),": error: no value specified for ",trim(options(i))
        stop
      endif
    endif

    if ((trim(options(i))=="-o" .or. trim(options(i))=="--output-dir")) then
      if ((i+1) .le. Narg) then
        read (options(i+1),"(A)") outputdirectory
        outputdirectory=trim(outputdirectory)//"/"
        inquire(file=trim(outputdirectory), exist=file_exists) ! trailing slash ensures it's looking for a directory
      else
        print *,gettime(),": error: no value specified for ",trim(options(i))
        stop
      endif
      if (.not. file_exists) then
        print *,gettime(),": error: output directory does not exist"
        stop
      endif
      options(i:i+1)=""
    endif

    if (trim(options(i))=="-skycat") then
      if ((i+1) .le. Narg) then
        read (options(i+1),"(A)") skylinelistfile
        options(i:i+1)=""
      else
        print *,gettime(),": error: no value specified for ",trim(options(i))
        stop
      endif
    endif

    if (trim(options(i))=="-strongcat") then
      if ((i+1) .le. Narg) then
        read (options(i+1),"(A)") stronglinelistfile
        options(i:i+1)=""
      else
        print *,gettime(),": error: no value specified for ",trim(options(i))
        stop
      endif
    endif

    if (trim(options(i))=="-deepcat") then
      if ((i+1) .le. Narg) then
        read (options(i+1),"(A)") deeplinelistfile
        options(i:i+1)=""
      else
        print *,gettime(),": error: no value specified for ",trim(options(i))
        stop
      endif
    endif

    if ((trim(options(i))=="-g" .or. trim(options(i))=="--generations")) then
      if ((i+1) .le. Narg) then
        read (options(i+1),*) generations
        options(i:i+1)=""
        if (generations .lt. 1) then
          print *,gettime(),": error: invalid value given for generations"
          stop
        endif
      else
        print *,gettime(),": error: no value specified for ",trim(options(i))
        stop
      endif
    endif

    if ((trim(options(i))=="-ps" .or. trim(options(i))=="--populationsize")) then
      if ((i+1) .le. Narg) then
        read (options(i+1),*) popsize
        options(i:i+1)=""
        if (popsize .lt. 1) then
          print *,gettime(),": error: invalid value given for popsize"
          stop
        endif
      else
        print *,gettime(),": error: no value specified for ",trim(options(i))
        stop
      endif
    endif

    if ((trim(options(i))=="-pr" .or. trim(options(i))=="--pressure")) then
      if ((i+1) .le. Narg) then
        read (options(i+1),*) pressure
        options(i:i+1)=""
        if (pressure .lt. 0.d0 .or. pressure .gt. 1.d0) then
          print *,"error: pressure must be between 0 and 1"
          stop
        endif
      else
        print *,gettime(),": error: no value specified for ",trim(options(i))
        stop
      endif
    endif
  ! to implement:
  !   continuum window and percentile
  enddo

  nargused=narg-count(options.ne."")

  do i=1,narg
    if (len(trim(options(i))).gt.0) then
      spectrumfile=options(i)
      exit
    endif
  enddo

!check that an input file was specified and that no unrecognised options are present

  if (len(trim(spectrumfile)).eq.0) then
    print *,gettime(),": error: no input file specified"
    stop
  elseif (narg - nargused .gt. 1) then
    print *,gettime(),": error: some input options were not recognised:"
    do i=1,narg
      if (len(trim(options(i))).gt.0) then
        print *,trim(options(i))
      endif
    enddo
    stop
  endif

!deal with image sections

  if (index(spectrumfile,"[") .gt. 0) then !image section specified
    imagesection=spectrumfile(index(spectrumfile,"["):)
    spectrumfile=spectrumfile(1:index(spectrumfile,"[")-1)
  endif

!check if input file exists

  inquire(file=spectrumfile, exist=file_exists) ! see if the input file is present

  if (.not. file_exists) then
    print *,gettime(),": error: input spectrum ",trim(spectrumfile)," does not exist"
    stop
  endif

  deallocate(options)

!display the settings

  print *,gettime(),": ALFA is running with the following settings:"
  if (.not.normalise) then
    print *,"            normalisation:                    using measured value of Hb"
  else
  if (normalisation.eq.0.d0) then
    print *,"            normalisation:                    no normalisation"
  else
    print *,"            normalisation:                    to Hb=",normalisation
  endif
  endif
  print *,"            spectrum fitted if max value >    ",baddata
  print *,"            velocity guess:                   ",redshiftguess
  print *,"            resolution guess:                 ",resolutionguess
  print *,"            first pass velocity tolerance:    ",vtol1*c
  print *,"            second pass velocity tolerance:   ",vtol2*c
  print *,"            first pass resolution tolerance:  ",rtol1
  print *,"            second pass resolution tolerance: ",rtol2
  if (subtractsky) then
  print *,"            sky line fitting:                 enabled"
  print *,"            sky line catalogue:               ",trim(skylinelistfile)
  else
  print *,"            sky line fitting:                 disabled"
  endif
  print *,"            strong line catalogue:            ",trim(stronglinelistfile)
  print *,"            deep line catalogue:              ",trim(deeplinelistfile)
  print *,"            number of generations:            ",generations
  print *,"            population size:                  ",popsize
  print *,"            pressure factor:                  ",pressure
  print *,"            output directory:                 ",trim(outputdirectory)

end subroutine readcommandline

end module mod_commandline
