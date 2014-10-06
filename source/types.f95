module mod_types

type spectrum
  real :: wavelength
  real :: flux
end type

type linelist
  real :: redshift
  real :: width
  real, allocatable :: wavelength(:)
  real, allocatable :: peak(:)
end type

end module mod_types
