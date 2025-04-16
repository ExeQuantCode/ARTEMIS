module artemis__constants
  implicit none
  integer, parameter, public :: real32 = Selected_real_kind(6,37)
  real(real32), parameter, public :: k_b = 1.3806503e-23_real32
  real(real32), parameter, public :: hbar = 1.05457148e-34_real32
  real(real32), parameter, public :: h = 6.626068e-34_real32
  real(real32), parameter, public :: atomic_mass=1.67262158e-27_real32
  real(real32), parameter, public :: avogadros=6.022e23_real32
  real(real32), parameter, public :: bohrtoang=0.529177249_real32
  real(real32), parameter, public :: pi = 4._real32*atan(1._real32)
  real(real32), parameter, public :: INF = huge(0._real32)
  integer, public :: ierror = -1
  real(real32), parameter, public :: tolerance = 1.E-6_real32
end MODULE artemis__constants
