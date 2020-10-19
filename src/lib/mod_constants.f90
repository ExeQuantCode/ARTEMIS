MODULE constants
  implicit none
  real, parameter, public :: k_b = 1.3806503e-23
  real, parameter, public :: hbar = 1.05457148e-34
  real, parameter, public :: h = 6.626068e-34
  real, parameter, public :: atomic_mass=1.67262158e-27
  real, parameter, public :: avogadros=6.022e23
  real, parameter, public :: bohrtoang=0.529177249
  integer, parameter, public :: real12 = Selected_real_kind(15,307)
  double precision, parameter, public :: pi = 4.D0*atan(1.D0)
  double precision, parameter, public :: INF = huge(0.D0)
  integer, public :: ierror = -1
end MODULE constants
