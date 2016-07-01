!Copyright (C) 2013- Roger Wesson
!Free under the terms of the GNU General Public License v3

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
  integer :: blended
end type

end module mod_types
