module mod_types

type spectrum
  real :: wavelength
  real :: flux
  real :: uncertainty
end type

type linelist
  real :: wavelength
  real :: peak
  real :: uncertainty
  real :: redshift
  real :: resolution
  character(len=85) :: linedata
end type

end module mod_types
