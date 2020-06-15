!Copyright (C) 2013- Roger Wesson
!Free under the terms of the GNU General Public License v3
module mod_globals

use mod_types

implicit none
integer :: I, spectrumlength
real :: startwlen, endwlen
character (len=512) :: spectrumfile,stronglinelistfile,deeplinelistfile,skylinelistfile,outputdirectory
character (len=32) :: imagesection

type(linelist), dimension(:), allocatable :: skylines_catalogue, stronglines_catalogue, deeplines_catalogue

real :: baddata
integer :: bdcount
integer :: cube_i, cube_j, rss_i

CHARACTER(len=2048) :: commandline

real :: vtol1, vtol2, rtol1, rtol2
real :: blendpeak
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
logical :: clobber !true to skip fitted spectra or pixels, false to overwrite

real, dimension(:), allocatable :: spectrum_1d
real, dimension(:,:), allocatable :: spectrum_2d
real, dimension(:,:,:), allocatable :: spectrum_3d

real, dimension(:), allocatable :: exclusions ! lines to be ignored when reading in catalogues
real :: detectionlimit ! sigma required to consider a line detected. default is 3.0
real :: minimumwavelength,maximumwavelength ! limits of spectrum, to be passed to catalogue reading subroutines
real :: wavelengthscaling ! default is 1.0 which is for wavelengths in A.  subroutine getfiletype checks if FITS header indicates units are nm, and sets to 10.0 if so.  todo: check for other units
integer, dimension(:), allocatable :: axes !number of pixels in each dimension
real, dimension(:), allocatable :: wavelengths !wavelength array

character(len=12) :: fluxformat !for writing out the line list
character(len=4),dimension(2) :: filenameformat !variable format to give suitable file names for 2D and 3D outputs
character(len=5) :: outputformat !text, csv, fits. json? xml?

integer :: rebinfactor,continuumwindow
integer :: tablewavelengthcolumn,tablefluxcolumn

! switches for writing out continuum fluxes

integer :: writeb1, writeb2, writep1, writep2

end module mod_globals
