module mod_types

type spectrum
  real :: wavelength
  real :: flux
  real :: uncertainty
end type

type linelist
  real :: redshift
  real :: resolution
  real, allocatable :: wavelength(:)
  real, allocatable :: peak(:)
  real, allocatable :: uncertainty(:)
end type

end module mod_types
