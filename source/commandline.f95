! command line reading routine
! included in both alfa.f95 and alfa_cube.f95

call get_command(commandline)
allocate (options(Narg))

do i=1,Narg
  call get_command_argument(i,options(i))
enddo

do i=1,narg
  if ((trim(options(i))=="-n" .or. trim(options(i))=="--normalise") .and. (i+1) .le. Narg) then
    read (options(i+1),*) normalisation
    normalise=.true.
    nargused = nargused + 2
  endif
  if ((trim(options(i))=="-vg" .or. trim(options(i))=="--velocity-guess") .and. (i+1) .le. Narg) then
    read (options(i+1),*) redshiftguess
    nargused = nargused + 2
  endif
  if ((trim(options(i))=="-rg" .or. trim(options(i))=="--resolution-guess") .and. (i+1) .le. Narg) then
    read (options(i+1),*) resolutionguess
    resolution_estimated=.true.
    nargused = nargused + 2
  endif
  if ((trim(options(i))=="-vtol1" .or. trim(options(i))=="--velocity-tolerance-1") .and. (i+1) .le. Narg) then
    read (options(i+1),*) vtol1
    vtol1 = vtol1/c
    nargused = nargused + 2
  endif
  if ((trim(options(i))=="-vtol2" .or. trim(options(i))=="--velocity-tolerance-2") .and. (i+1) .le. Narg) then
    read (options(i+1),*) vtol2
    vtol2 = vtol2/c
    nargused = nargused + 2
  endif
  if ((trim(options(i))=="-rtol1" .or. trim(options(i))=="--resolution-tolerance-1") .and. (i+1) .le. Narg) then
    read (options(i+1),*) rtol1
    nargused = nargused + 2
  endif
  if ((trim(options(i))=="-rtol2" .or. trim(options(i))=="--resolution-tolerance-2") .and. (i+1) .le. Narg) then
    read (options(i+1),*) rtol2
    nargused = nargused + 2
  endif
  if ((trim(options(i))=="-ss" .or. trim(options(i))=="--subtract-sky")) then
    subtractsky=.true.
    nargused = nargused + 1
  endif
  if ((trim(options(i))=="-o" .or. trim(options(i))=="--output-dir")) then
    read (options(i+1),*) outputdirectory
    outputdirectory=trim(outputdirectory)//"/"
    inquire(file=trim(outputdirectory), exist=file_exists) ! trailing slash ensures it's looking for a directory
    if (.not. file_exists) then
      print *,gettime(),": error: output directory does not exist"
      stop
    endif
    nargused = nargused + 2
  endif

! to implement:
!   continuum window and percentile
!   no. of generations, population size, pressure
!   linelists
enddo

if (narg - nargused .eq. 0) then
  print *,gettime(),": error: no input file specified"
  stop
elseif (narg - nargused .gt. 1) then
  print *,gettime(),": warning: some input options were not recognised"
else
  call get_command_argument(narg,spectrumfile)
  spectrumfile=trim(spectrumfile)
endif

deallocate(options) 
