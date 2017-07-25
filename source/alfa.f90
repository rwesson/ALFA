program alfa
!ALFA, the Automated Line Fitting Algorithm
!Copyright (C) 2013- Roger Wesson

!This program is free software: you can redistribute it and/or modify
!it under the terms of the GNU General Public License as published by
!the Free Software Foundation, either version 3 of the License, or
!(at your option) any later version.

!This program is distributed in the hope that it will be useful,
!but WITHOUT ANY WARRANTY; without even the implied warranty of
!MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!GNU General Public License for more details.

!You should have received a copy of the GNU General Public License
!along with this program.  If not, see <http://www.gnu.org/licenses/>.

use mod_readfiles
use mod_routines
use mod_types
use mod_continuum
use mod_fit
use mod_uncertainties
use mod_commandline

implicit none
integer :: I, spectrumlength, nlines, linearraypos, totallines, startpos, endpos
real :: startwlen, endwlen
character (len=512) :: spectrumfile,stronglinelistfile,deeplinelistfile,skylinelistfile,outputdirectory,outputbasename
character (len=32) :: imagesection

type(linelist), dimension(:), allocatable :: skylines_catalogue, stronglines_catalogue, deeplines_catalogue
type(linelist), dimension(:), allocatable :: fittedlines, fittedlines_section, skylines, skylines_section
type(spectrum), dimension(:), allocatable :: realspec, fittedspectrum, spectrumchunk, skyspectrum, continuum, stronglines

real :: baddata
integer :: cube_i, cube_j, rss_i
integer, dimension(:), allocatable :: axes
real, dimension(:), allocatable :: wavelengths !wavelength array
real, dimension(:), allocatable :: spectrum_1d
real, dimension(:,:), allocatable :: spectrum_2d
real, dimension(:,:,:), allocatable :: spectrum_3d
real :: minimumwavelength,maximumwavelength ! limits of spectrum, to be passed to catalogue reading subroutines
real :: wavelengthscaling ! default is 1.0 which is for wavelengths in A.  subroutine getfiletype checks if FITS header indicates units are nm, and sets to 10.0 if so.  todo: check for other units

CHARACTER(len=2048) :: commandline

real :: redshiftguess, redshiftguess_initial, redshiftguess_overall ! redshiftguess_initial is the user-specified value, used in the initial fit which determines redshiftguess_overall.  that is then used in the chunks to find redshift
real :: resolutionguess, resolutionguess_initial ! resolutionguess_initial is the user-specified value, used in the initial fit to determine resolutionguess.
real :: vtol1, vtol2, rtol1, rtol2
real :: blendpeak
real :: normalisation, hbetaflux
real :: c
integer :: linelocation, overlap
integer :: generations, popsize
real :: pressure

logical :: normalise=.false. !false means spectrum normalised to whatever H beta is detected, true means spectrum normalised to user specified value
logical :: resolution_estimated=.false. !true means user specified a value, false means estimate from sampling
logical :: subtractsky=.false. !attempt to fit night sky emission lines
logical :: upperlimits=.false. !if true, code reports 3 sigma limit for undetected lines
logical :: file_exists

logical :: collapse !if true, 2D or 3D data is summed into a single spectrum
logical :: messages

real, dimension(:), allocatable :: exclusions ! lines to be ignored when reading in catalogues
real :: detectionlimit ! sigma required to consider a line detected. default is 3.0

character(len=12) :: fluxformat !for writing out the line list
character(len=4),dimension(2) :: filenameformat !variable format to give suitable file names for 2D and 3D outputs

integer :: rebinfactor

! openmp variables

integer :: tid, omp_get_thread_num, omp_get_num_threads

! switches for writing out continuum fluxes

integer :: writeb1, writeb2, writep1, writep2

c=299792.458 !km/s
!default values in absence of user specificed guess
redshiftguess=0.0 !km/s
redshiftguess_initial=0.0
resolutionguess=0.0 !lambda/deltalambda, determined assuming nyquist sampling if not specified
resolutionguess_initial=0.0
rtol1=0.d0 !variation allowed in resolution on first pass.  determined later, either from user input, or to be equal to resolution guess.
rtol2=500. !second pass
vtol1=0.003 !variation allowed in velocity (expressed as redshift) on first pass. 0.003 = 900 km/s
vtol2=0.0002 !second pass. 0.0002 = 60 km/s
baddata=0.d0
wavelengthscaling=1.d0
detectionlimit=1.0
rebinfactor=1

stronglinelistfile=trim(PREFIX)//"/share/alfa/strong.cat"
deeplinelistfile=trim(PREFIX)//"/share/alfa/deep.cat"
skylinelistfile=trim(PREFIX)//"/share/alfa/sky.cat"

outputdirectory="./"
imagesection=""

collapse=.false.
messages=.false.

popsize=30
pressure=0.3 !pressure * popsize needs to be an integer
generations=500

! start

print *,"ALFA, the Automated Line Fitting Algorithm"
print *,"version ",VERSION," (jwst mod)"

print *
print *,gettime(),"starting code"

! random seed

call init_random_seed()

! read command line

call readcommandline(commandline,normalise,normalisation,redshiftguess_initial,resolutionguess_initial,vtol1,vtol2,rtol1,rtol2,baddata,pressure,spectrumfile,outputdirectory,skylinelistfile,stronglinelistfile,deeplinelistfile,generations,popsize,subtractsky,resolution_estimated,file_exists,imagesection,upperlimits,wavelengthscaling,collapse,exclusions,detectionlimit,rebinfactor)

! convert from velocity to redshift

redshiftguess_initial=1.+(redshiftguess_initial/c)

! read in spectrum to fit and line list

print *
print *,gettime(),"reading in file ",trim(spectrumfile),":"

!call subroutine to read in the data.  input is filename, output is 3D array containing data, length of dimensions dependent on whether file was 1D, 2D or 3D.
call readdata(trim(spectrumfile)//trim(imagesection), spectrum_1d, spectrum_2d, spectrum_3d, wavelengths, wavelengthscaling, axes, rebinfactor)

minimumwavelength = wavelengths(1)
maximumwavelength = wavelengths(size(wavelengths))
spectrumlength = size(wavelengths)

! collapse data to 1D if requested
! apply baddata to sum only pixels with good data in
i=0 ! to count the number of pixels summed

if (collapse) then
  if (allocated(spectrum_2d)) then

    print *,gettime(),"summing rows in 2D spectrum"
    allocate(spectrum_1d(spectrumlength))
    spectrum_1d = 0.d0

    do rss_i = 1,axes(2)
      if (maxval(spectrum_2d(:,rss_i)).gt. baddata) then
        spectrum_1d = spectrum_1d + spectrum_2d(:,rss_i)
        i=i+1
      endif
    enddo

    print *,gettime(),i," of ",axes(2)," rows co-added"

    deallocate(spectrum_2d)

  elseif (allocated(spectrum_3d)) then

    print *,gettime(),"summing pixels in 3D spectrum"
    allocate(spectrum_1d(spectrumlength))
    spectrum_1d = 0.d0

    do cube_i = 1,axes(1)
      do cube_j = 1,axes(2)
        if (maxval(spectrum_3d(cube_i,cube_j,:)) .gt. baddata) then
          spectrum_1d = spectrum_1d + spectrum_3d(cube_i,cube_j,:)
          i=i+1
        endif
      enddo
    enddo

    print *,gettime(),i," of ",axes(1)*axes(2)," pixels co-added"

    deallocate(spectrum_3d)

  endif
endif

!now we have all the flux values in an array and can fit the spectra

!read in the line catalogues

print *,gettime(),"reading in line catalogues"
call readlinelist(skylinelistfile, skylines_catalogue, nlines,minimumwavelength,maximumwavelength,exclusions)
call readlinelist(stronglinelistfile, stronglines_catalogue, nlines,minimumwavelength,maximumwavelength,exclusions)
call readlinelist(deeplinelistfile, deeplines_catalogue, nlines,minimumwavelength,maximumwavelength,exclusions)

if (allocated(spectrum_1d)) then !1d spectrum

  allocate(realspec(size(wavelengths)))
  allocate(fittedspectrum(spectrumlength))

  realspec%wavelength = wavelengths
  realspec%flux = spectrum_1d
  realspec%uncertainty = 0.d0
  redshiftguess=redshiftguess_initial
  resolutionguess=resolutionguess_initial

  fittedspectrum%wavelength=realspec%wavelength
  fittedspectrum%flux=0.d0

  if (maxval(realspec%flux) .lt. baddata) then
    print *,gettime(),"no good data in spectrum (all fluxes are less than ",baddata,")"
    call exit(1)
  endif
  messages=.true.

  tid=0
  write (outputbasename,"(A)") spectrumfile(index(spectrumfile,"/",back=.true.)+1:len(trim(spectrumfile)))
  include "spectralfit.f90"

elseif (allocated(spectrum_2d)) then !fit 2D data

  write (filenameformat,"(A,I1,A,I1)") "I",floor(log10(real(axes(2))))+1,".",floor(log10(real(axes(2))))+1

!$OMP PARALLEL private(outputbasename,realspec,fittedspectrum,spectrumlength,continuum,nlines,spectrumchunk,linearraypos,overlap,startpos,startwlen,endpos,endwlen,skylines,skylines_section,stronglines,fittedlines,fittedlines_section,blendpeak,hbetaflux,totallines,skyspectrum,redshiftguess_overall,rss_i,tid,redshiftguess,resolutionguess,file_exists,writeb1,writeb2,writep1,writep2) shared(skylines_catalogue,stronglines_catalogue,deeplines_catalogue,axes,spectrumfile,redshiftguess_initial,resolutionguess_initial,outputdirectory)
!$OMP MASTER
  if (omp_get_num_threads().gt.1) then
    print "(X,A9,X,A,I2,A)",gettime(),"using ",omp_get_num_threads()," processors"
  endif
!$OMP END MASTER

!$OMP DO schedule(dynamic)
  do rss_i=1,axes(2) ! wavelength is on axis 1, row position is on axis 2, so this loops over rows

    tid=OMP_GET_THREAD_NUM()

    write (outputbasename,"(A,"//filenameformat(1)//")") spectrumfile(index(spectrumfile,"/",back=.true.)+1:len(trim(spectrumfile)))//"_row_",rss_i

    allocate(realspec(axes(1)))
    spectrumlength=axes(1)
    realspec%wavelength = wavelengths
    realspec%flux=spectrum_2d(:,rss_i)
    redshiftguess=redshiftguess_initial
    resolutionguess=resolutionguess_initial

!check for valid data
!ultra crude at the moment

    inquire(file=trim(outputdirectory)//trim(outputbasename)//"_lines", exist=file_exists)

    if (maxval(realspec%flux) .lt. baddata) then
      print "(X,A,A,I2,A,"//filenameformat(1)//",A,ES10.2)",gettime(),"(thread ",tid,") : skipped row  ",rss_i,": all values below ",baddata
      deallocate(realspec)
      cycle
    elseif (file_exists) then
      print "(X,A,A,I2,A,"//filenameformat(1)//",A)",gettime(),"(thread ",tid,") : skipped row  ",rss_i,": already fitted"
      deallocate(realspec)
      cycle
    endif

    print *,gettime(),"fitting row ",rss_i

    allocate (fittedspectrum(spectrumlength))
    fittedspectrum%wavelength=realspec%wavelength
    fittedspectrum%flux=0.d0

!now do the fitting

    include "spectralfit.f90"

!deallocate arrays ready for the next pixel
    deallocate(realspec)
    deallocate(fittedspectrum)
    deallocate(continuum)
    if (allocated(skyspectrum)) deallocate(skyspectrum)

    print "(X,A,A,I2,A,"//filenameformat(1)//",A,F6.1,A,F6.0)",gettime(),"(thread ",tid,") : finished row ",rss_i,". approx velocity and resolution ",c*(redshiftguess_overall-1.d0),", ",fittedlines(1)%resolution

  enddo

!$OMP END DO
!$OMP END PARALLEL

  print *,gettime(),"all processing finished"

  deallocate(spectrum_2d)

elseif (allocated(spectrum_3d)) then !fit 3D data

  print *,gettime(),"processing cube"
!$OMP PARALLEL private(outputbasename,realspec,fittedspectrum,spectrumlength,continuum,nlines,spectrumchunk,linearraypos,overlap,startpos,startwlen,endpos,endwlen,skylines,skylines_section,stronglines,fittedlines,fittedlines_section,blendpeak,hbetaflux,totallines,skyspectrum,redshiftguess_overall,cube_i,cube_j,tid,redshiftguess,resolutionguess,file_exists,writeb1,writeb2,writep1,writep2) shared(skylines_catalogue,stronglines_catalogue,deeplines_catalogue,axes,spectrumfile,redshiftguess_initial,resolutionguess_initial,outputdirectory)

!$OMP MASTER
  if (omp_get_num_threads().gt.1) then
    print "(X,A9,X,A,I2,A)",gettime(),"using ",omp_get_num_threads()," processors"
  endif
!$OMP END MASTER

  write (filenameformat(1),"(A,I1,A,I1)") "I",floor(log10(real(axes(1))))+1,".",floor(log10(real(axes(1))))+1
  write (filenameformat(2),"(A,I1,A,I1)") "I",floor(log10(real(axes(2))))+1,".",floor(log10(real(axes(2))))+1

!$OMP DO schedule(dynamic)
  do cube_i=1,axes(1)
    do cube_j=1,axes(2)

      tid=OMP_GET_THREAD_NUM()

      write (outputbasename,"(A,"//filenameformat(1)//",A1,"//filenameformat(2)//")") spectrumfile(index(spectrumfile,"/",back=.true.)+1:len(trim(spectrumfile)))//"_pixel_",cube_i,"_",cube_j
      allocate(realspec(axes(3)))
      spectrumlength=axes(3)
      realspec%flux=spectrum_3d(cube_i,cube_j,:)
      realspec%wavelength=wavelengths
      redshiftguess=redshiftguess_initial
      resolutionguess=resolutionguess_initial

!check for valid data
!ultra crude at the moment

      inquire(file=trim(outputdirectory)//trim(outputbasename)//"_lines", exist=file_exists)

      if (maxval(realspec%flux) .lt. baddata) then
        print "(X,A,A,I2,A,"//filenameformat(1)//",A,"//filenameformat(2)//",A,ES10.2)",gettime(),"(thread ",tid,") : skipped pixel  ",cube_i,",",cube_j,": all values below ",baddata
        deallocate(realspec)
        cycle
      elseif (file_exists) then
        print "(X,A,A,I2,A,"//filenameformat(1)//",A,"//filenameformat(2)//",A)",gettime(),"(thread ",tid,") : skipped pixel  ",cube_i,",",cube_j,": already fitted"
        deallocate(realspec)
        cycle
      endif

      allocate (fittedspectrum(spectrumlength))
      fittedspectrum%wavelength=realspec%wavelength
      fittedspectrum%flux=0.d0

!now do the fitting

      include "spectralfit.f90"

!deallocate arrays ready for the next pixel

      deallocate(realspec)
      deallocate(fittedspectrum)
      deallocate(continuum)
      if (allocated(skyspectrum)) deallocate(skyspectrum)

      print "(X,A,A,I2,A,"//filenameformat(1)//",A,"//filenameformat(2)//",A,F6.1,A,F5.0)",gettime(),"(thread ",tid,") : finished pixel ",cube_i,",",cube_j,". approx velocity and resolution ",c*(redshiftguess_overall-1.d0),", ",fittedlines(1)%resolution

    enddo
  enddo

!$OMP END DO
!$OMP END PARALLEL

  print *,gettime(),"all processing finished"

  deallocate(spectrum_3d)

endif

!free memory

if (allocated(spectrum_2d)) deallocate(spectrum_2d)
if (allocated(spectrum_3d)) deallocate(spectrum_3d)
if (allocated(axes)) deallocate(axes)
if (allocated(wavelengths)) deallocate(wavelengths)
if (allocated(skylines_catalogue)) deallocate(skylines_catalogue)
if (allocated(stronglines_catalogue)) deallocate(stronglines_catalogue)
if (allocated(deeplines_catalogue)) deallocate(deeplines_catalogue)

print *,gettime(),"all done"
print *

end program alfa
