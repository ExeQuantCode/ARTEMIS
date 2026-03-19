module artemis__constants
  !! Module containing ARTEMIS physical and numerical constants.
  !!
  !! Re-exports fundamental constants from coreutils__const and defines
  !! ARTEMIS-specific constants such as error state and numerical tolerance.
  use coreutils__kind, only: real32
  use coreutils__const, only: pi, INF, k_b, hbar, h, &
       atomic_mass, avogadros, bohrtoang
  implicit none

  public :: real32, pi, INF, k_b, hbar, h, atomic_mass, avogadros, bohrtoang

  integer, public :: ierror = -1
  !! Error state flag. Default: -1 (no error).
  real(real32), parameter, public :: tolerance = 1.E-6_real32
  !! Numerical tolerance for floating-point comparisons.

end module artemis__constants
