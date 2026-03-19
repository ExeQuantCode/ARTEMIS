module artemis__constants
  use coreutils__kind, only: real32
  use coreutils__const, only: pi, INF, k_b, hbar, h, &
       atomic_mass, avogadros, bohrtoang
  implicit none

  public :: real32, pi, INF, k_b, hbar, h, atomic_mass, avogadros, bohrtoang

  integer, public :: ierror = -1
  real(real32), parameter, public :: tolerance = 1.E-6_real32
end MODULE artemis__constants
