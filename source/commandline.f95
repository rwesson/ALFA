module mod_commandline
use mod_routines

contains

subroutine readcommandline(commandline,normalise,normalisation,redshiftguess,resolutionguess,vtol1,vtol2,rtol1,rtol2,baddata,pressure,spectrumfile,outputdirectory,skylinelistfile,stronglinelistfile,deeplinelistfile,generations,popsize,subtractsky,resolution_estimated,file_exists)

  implicit none

  logical :: normalise
  real :: normalisation,redshiftguess,resolutionguess,vtol1,vtol2,rtol1,rtol2,baddata,pressure,c
  character(len=2048) :: commandline
  character(len=2048), dimension(:), allocatable :: options
  character(len=512),intent(out) :: spectrumfile,outputdirectory,skylinelistfile,stronglinelistfile,deeplinelistfile
  integer,intent(out) :: generations,popsize
  integer :: Narg,nargused,i
  logical,intent(out) :: subtractsky,resolution_estimated,file_exists

  c=299792.458 !km/s

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
  print *,gettime(),": command line: ",trim(commandline)

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
    if (trim(options(i))=="-ss" .or. trim(options(i))=="--subtract-sky") then
      subtractsky=.true.
      nargused = nargused + 1
    endif
    if (trim(options(i))=="-b" .or. trim(options(i))=="--bad-data") then
      read (options(i+1),*) baddata
      nargused = nargused + 2
    endif
    if ((trim(options(i))=="-o" .or. trim(options(i))=="--output-dir") .and. (i+1) .le. Narg) then
      read (options(i+1),"(A)") outputdirectory
      outputdirectory=trim(outputdirectory)//"/"
      inquire(file=trim(outputdirectory), exist=file_exists) ! trailing slash ensures it's looking for a directory
      if (.not. file_exists) then
        print *,gettime(),": error: output directory does not exist"
        stop
      endif
      nargused = nargused + 2
    endif
    if (trim(options(i))=="-skycat" .and. (i+1) .le. Narg) then
      read (options(i+1),"(A)") skylinelistfile
      nargused = nargused + 2
    endif
    if (trim(options(i))=="-strongcat" .and. (i+1) .le. Narg) then
      read (options(i+1),"(A)") stronglinelistfile
      nargused = nargused + 2
    endif
    if (trim(options(i))=="-deepcat" .and. (i+1) .le. Narg) then
      read (options(i+1),"(A)") deeplinelistfile
      nargused = nargused + 2
    endif
    if ((trim(options(i))=="-g" .or. trim(options(i))=="--generations") .and. (i+1) .le. Narg) then
      read (options(i+1),*) generations
      nargused = nargused + 2
    endif
    if ((trim(options(i))=="-ps" .or. trim(options(i))=="--populationsize") .and. (i+1) .le. Narg) then
      read (options(i+1),*) popsize
      nargused = nargused + 2
    endif
    if ((trim(options(i))=="-pr" .or. trim(options(i))=="--pressure") .and. (i+1) .le. Narg) then
      read (options(i+1),*) pressure
      nargused = nargused + 2
    endif
  ! to implement:
  !   continuum window and percentile
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
