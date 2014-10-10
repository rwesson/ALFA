module mod_uncertainties
use mod_types
use mod_quicksort

contains

subroutine get_uncertainties(synthspec, realspec, population, rms)
implicit none

real, dimension(:) :: rms
type(spectrum), dimension(:,:), allocatable :: synthspec
type(spectrum), dimension(:), allocatable :: realspec
type(linelist), dimension(:), allocatable :: population
real, dimension(:), allocatable :: residuals, medianresiduals
real, dimension(20) :: spectrumchunk
real :: wavelengthsampling
integer :: i, uncertaintywavelengthindex

allocate(residuals(size(realspec)))
allocate(medianresiduals(size(realspec)))

! median filter the residuals

residuals=realspec(:)%flux - synthspec(:,minloc(rms,1))%flux

do i=11,size(realspec)-10
  spectrumchunk=residuals(i-10:i+10)
  call qsort(spectrumchunk)
  medianresiduals(i)=spectrumchunk(11)
end do

medianresiduals(1:10)=medianresiduals(11)
medianresiduals(size(medianresiduals)-10:size(medianresiduals))=medianresiduals(size(medianresiduals)-11)

! calculate rms in moving 20 unit window

do i=11,size(realspec)-10
  spectrumchunk=medianresiduals(i-10:i+10)
  realspec(i)%uncertainty=(sum(spectrumchunk**2))**0.5
end do

realspec(1:10)%uncertainty=realspec(11)%uncertainty
realspec(size(realspec%uncertainty)-10:size(realspec%uncertainty))%uncertainty=realspec(size(realspec%uncertainty)-11)%uncertainty

! determine uncertainty for each line from ratio of peak flux to rms at wavelength of line

wavelengthsampling=realspec(2)%wavelength - realspec(1)%wavelength

do i=1,size(population(1)%uncertainty)
  uncertaintywavelengthindex=minloc(abs(realspec%wavelength-population(minloc(rms,1))%wavelength(i)),1)
  population(minloc(rms,1))%uncertainty(i)=0.67*(population(minloc(rms,1))%width/wavelengthsampling)**0.5&
  &*population(minloc(rms,1))%peak(i)&
  &/realspec(uncertaintywavelengthindex)%uncertainty
end do

end subroutine get_uncertainties

end module mod_uncertainties
