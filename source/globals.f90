module mod_globals

use mod_types

implicit none
integer :: I, spectrumlength, nlines, linearraypos, totallines, startpos, endpos
real :: startwlen, endwlen
character (len=512) :: spectrumfile,stronglinelistfile,deeplinelistfile,skylinelistfile,outputdirectory,outputbasename
character (len=32) :: imagesection

type(linelist), dimension(:), allocatable :: skylines_catalogue, stronglines_catalogue, deeplines_catalogue
type(linelist), dimension(:), allocatable :: fittedlines, fittedlines_section, skylines, skylines_section
type(spectrum), dimension(:), allocatable :: realspec, fittedspectrum, spectrumchunk, skyspectrum, continuum, stronglines, maskedspectrum

real :: baddata
integer :: bdcount
integer :: cube_i, cube_j, rss_i
integer, dimension(:), allocatable :: axes !number of pixels in each dimension
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
logical :: subtractcontinuum=.true. !can be set to false if continuum is already subtracted
logical :: file_exists

logical :: collapse !if true, 2D or 3D data is summed into a single spectrum
logical :: messages

real, dimension(:), allocatable :: exclusions ! lines to be ignored when reading in catalogues
real :: detectionlimit ! sigma required to consider a line detected. default is 3.0

character(len=12) :: fluxformat !for writing out the line list
character(len=4),dimension(2) :: filenameformat !variable format to give suitable file names for 2D and 3D outputs
character(len=5) :: outputformat !text, csv, fits. json? xml?

integer :: rebinfactor,continuumwindow
integer :: tablewavelengthcolumn,tablefluxcolumn

! switches for writing out continuum fluxes

integer :: writeb1, writeb2, writep1, writep2

! openmp

integer :: threadnumber

end module mod_globals
