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

integer :: filetype, dimensions
real :: wavelength, dispersion, referencepixel, baddata
logical :: loglambda
integer :: cube_i, cube_j, cube_k, rss_i, rss_k
integer, dimension(:), allocatable :: axes
real, dimension(:,:), allocatable :: rssdata
real, dimension(:,:,:), allocatable :: cubedata
real :: minimumwavelength,maximumwavelength ! limits of spectrum, to be passed to catalogue reading subroutines

CHARACTER(len=2048) :: commandline

real :: redshiftguess, resolutionguess, redshiftguess_overall
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

logical :: messages

character(len=12) :: fluxformat !for writing out the line list

! openmp variables

integer :: tid, omp_get_thread_num, omp_get_num_threads

c=299792.458 !km/s
!default values in absence of user specificed guess
redshiftguess=0.0 !km/s
resolutionguess=0.0 !lambda/deltalambda, determined assuming nyquist sampling if not specified
rtol1=0.d0 !variation allowed in resolution on first pass.  determined later, either from user input, or to be equal to resolution guess.
rtol2=500. !second pass
vtol1=0.003 !variation allowed in velocity (expressed as redshift) on first pass. 0.003 = 900 km/s
vtol2=0.0002 !second pass. 0.0002 = 60 km/s
baddata=0.d0

stronglinelistfile=trim(PREFIX)//"/share/alfa/strong.cat"
deeplinelistfile=trim(PREFIX)//"/share/alfa/deep.cat"
skylinelistfile=trim(PREFIX)//"/share/alfa/sky.cat"

outputdirectory="./"
imagesection=""

messages=.false.

popsize=30
pressure=0.3 !pressure * popsize needs to be an integer
generations=500

! start

print *,"ALFA, the Automated Line Fitting Algorithm"
if (len(VERSION).gt.0) then
  print *,"version ",VERSION
else
  print *,"version 0.98"
endif

print *
print *,gettime(),"starting code"

! random seed

call init_random_seed()

! read command line

call readcommandline(commandline,normalise,normalisation,redshiftguess,resolutionguess,vtol1,vtol2,rtol1,rtol2,baddata,pressure,spectrumfile,outputdirectory,skylinelistfile,stronglinelistfile,deeplinelistfile,generations,popsize,subtractsky,resolution_estimated,file_exists,imagesection,upperlimits)

! convert from velocity to redshift

redshiftguess=1.+(redshiftguess/c)

! read in spectrum to fit and line list

print *,gettime(),"reading in file ",trim(spectrumfile)

!call subroutine to determine whether it's 1D, 2D or 3D fits, or ascii, or none of the above
call getfiletype(trim(spectrumfile)//imagesection,filetype,dimensions,axes,wavelength,dispersion,referencepixel,loglambda)

if (filetype.eq.1) then !1d fits file

  spectrumlength=axes(1)
  call read1dfits(spectrumfile, realspec, spectrumlength, fittedspectrum, wavelength, dispersion, referencepixel, loglambda)
  minimumwavelength=realspec(1)%wavelength
  maximumwavelength=realspec(spectrumlength)%wavelength
  if (maxval(realspec%flux) .lt. baddata) then
    print *,gettime(),"no good data in spectrum (all fluxes are less than ",baddata,")"
    call exit(1)
  endif
  messages=.true.

elseif (filetype .eq. 2) then !2d fits file

  call read2dfits(trim(spectrumfile)//imagesection, rssdata, dimensions, axes)

  if (loglambda) then
    minimumwavelength=wavelength*exp((1-referencepixel)*dispersion/wavelength)
    maximumwavelength=wavelength*exp((axes(1)-referencepixel)*dispersion/wavelength)
  else
    minimumwavelength=wavelength
    maximumwavelength=wavelength+(axes(1)-1)*dispersion
  endif

elseif (filetype .eq. 3) then !3d fits file

  call read3dfits(trim(spectrumfile)//imagesection, cubedata, dimensions, axes)

  if (loglambda) then
    minimumwavelength=wavelength*exp((1-referencepixel)*dispersion/wavelength)
    maximumwavelength=wavelength*exp((axes(3)-referencepixel)*dispersion/wavelength)
  else
    minimumwavelength=wavelength
    maximumwavelength=wavelength+(axes(3)-1)*dispersion
  endif

elseif (filetype .eq. 4) then !1d ascii file
  call readascii(spectrumfile, realspec, spectrumlength, fittedspectrum)
  minimumwavelength=realspec(1)%wavelength
  maximumwavelength=realspec(spectrumlength)%wavelength
  if (maxval(realspec%flux) .lt. baddata) then
    print *,gettime(),"no good data in spectrum (all fluxes are less than ",baddata,")"
    call exit(1)
  endif
  messages=.true.
else
  !not recognised, stop
  print *,"unrecognised file"
  call exit(1)
endif

print *,gettime(),"wavelength range ",minimumwavelength," to ",maximumwavelength,"(log: ",loglambda,")" ! PmW
if (loglambda) print *,gettime(),"warning: uncertainty estimation is not reliable for log-sampled spectra"

!read in catalogues

print *,gettime(),"reading in line catalogues"
call readlinelist(skylinelistfile, skylines_catalogue, nlines,minimumwavelength,maximumwavelength)
call readlinelist(stronglinelistfile, stronglines_catalogue, nlines,minimumwavelength,maximumwavelength)
call readlinelist(deeplinelistfile, deeplines_catalogue, nlines,minimumwavelength,maximumwavelength)

if (filetype .eq. 1 .or. filetype .eq. 4) then !fit 1D data
  tid=0
  write (outputbasename,"(A)") spectrumfile(index(spectrumfile,"/",back=.true.)+1:len(trim(spectrumfile)))
  include "spectralfit.f90"
elseif (filetype .eq. 2) then !fit 2D data

!$OMP PARALLEL private(outputbasename,realspec,fittedspectrum,spectrumlength,continuum,nlines,spectrumchunk,linearraypos,overlap,startpos,startwlen,endpos,endwlen,skylines,skylines_section,stronglines,fittedlines,fittedlines_section,blendpeak,hbetaflux,totallines,skyspectrum,redshiftguess_overall,rss_i,tid) firstprivate(redshiftguess,resolutionguess) shared(skylines_catalogue,stronglines_catalogue,deeplines_catalogue,axes,spectrumfile)
!$OMP MASTER
  if (omp_get_num_threads().gt.1) then
    print "(X,A9,X,A,I2,A)",gettime(),"using ",omp_get_num_threads()," processors"
  endif
!$OMP END MASTER

!$OMP DO schedule(dynamic)
  do rss_i=1,axes(2)

    tid=OMP_GET_THREAD_NUM()

    write (outputbasename,"(A,I5.5)") spectrumfile(index(spectrumfile,"/",back=.true.)+1:len(trim(spectrumfile)))//"_row_",rss_i
    allocate(realspec(axes(1)))
    spectrumlength=axes(1)
    realspec%flux=rssdata(:,rss_i)

    if (loglambda) then
      do rss_k=1,axes(1)
        realspec(rss_k)%wavelength=wavelength*exp((rss_k-referencepixel)*dispersion/wavelength)
      enddo
    else
      do rss_k=1,axes(1)
        realspec(rss_k)%wavelength=wavelength+(rss_k-referencepixel)*dispersion
      enddo
    endif

!check for valid data
!ultra crude at the moment

    inquire(file=trim(outputdirectory)//trim(outputbasename)//"_lines", exist=file_exists)

    if (maxval(realspec%flux) .lt. baddata .or. file_exists) then
      print "(X,A,A,I2,A,I5.5,A,I5.5)",gettime(),"(thread ",tid,") : skipped row  ",rss_i
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

    print "(X,A,A,I2,A,I5.5,A,I5.5)",gettime(),"(thread ",tid,") : finished row ",rss_i

  enddo

!$OMP END DO
!$OMP END PARALLEL

  print *,gettime(),"all processing finished"

  deallocate(rssdata)

elseif (filetype .eq. 3) then !fit 3D data
!process cube
  print *,gettime(),"processing cube"
!$OMP PARALLEL private(outputbasename,realspec,fittedspectrum,spectrumlength,continuum,nlines,spectrumchunk,linearraypos,overlap,startpos,startwlen,endpos,endwlen,skylines,skylines_section,stronglines,fittedlines,fittedlines_section,blendpeak,hbetaflux,totallines,skyspectrum,redshiftguess_overall,cube_i,cube_j,tid) firstprivate(redshiftguess,resolutionguess) shared(skylines_catalogue,stronglines_catalogue,deeplines_catalogue,axes,spectrumfile)

!$OMP MASTER
  if (omp_get_num_threads().gt.1) then
    print "(X,A9,X,A,I2,A)",gettime(),"using ",omp_get_num_threads()," processors"
  endif
!$OMP END MASTER

!$OMP DO schedule(dynamic)
  do cube_i=1,axes(1)
    do cube_j=1,axes(2)

      tid=OMP_GET_THREAD_NUM()

      write (outputbasename,"(A,I3.3,A1,I3.3)") spectrumfile(index(spectrumfile,"/",back=.true.)+1:len(trim(spectrumfile)))//"_pixel_",cube_i,"_",cube_j
      allocate(realspec(axes(3)))
      spectrumlength=axes(3)
      realspec%flux=cubedata(cube_i,cube_j,:)

      if (loglambda) then
        do cube_k=1,axes(3)
          realspec(cube_k)%wavelength=wavelength*exp((cube_k-referencepixel)*dispersion/wavelength)
        enddo
      else
        do cube_k=1,axes(3)
          realspec(cube_k)%wavelength=wavelength+(cube_k-referencepixel)*dispersion
        enddo
      endif

!check for valid data
!ultra crude at the moment

      inquire(file=trim(outputdirectory)//trim(outputbasename)//"_lines", exist=file_exists)

      if (maxval(realspec%flux) .lt. baddata .or. file_exists) then
        print "(X,A,A,I2,A,I3.3,A,I3.3)",gettime(),"(thread ",tid,") : skipped pixel  ",cube_i,",",cube_j
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

      print "(X,A,A,I2,A,I3.3,A,I3.3)",gettime(),"(thread ",tid,") : finished pixel ",cube_i,",",cube_j

    enddo
  enddo

!$OMP END DO
!$OMP END PARALLEL

  print *,gettime(),"all processing finished"

  deallocate(cubedata)

endif

!free memory

if (allocated(rssdata)) deallocate(rssdata)
if (allocated(cubedata)) deallocate(cubedata)

print *,gettime(),"all done"
print *

end program alfa
