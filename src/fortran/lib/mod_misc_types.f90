module artemis__misc_types
  !! Module containing custom derived types for ARTEMIS
  use artemis__constants, only: real32
  implicit none


  private

  public :: latmatch_type
  public :: tol_type


  type latmatch_type
     integer :: nfit
     logical :: lreduced
     character(1) :: abc(3)=(/'a','b','c'/)

     integer, dimension(2) :: axes
     integer, allocatable, dimension(:,:,:) :: tf1,tf2
     real(real32), allocatable, dimension(:,:) :: tol
     real(real32), dimension(3,3) :: lat1,lat2
  end type latmatch_type

  type tol_type
     integer :: maxsize,maxfit,nstore
     real(real32) :: maxlen=20._real32
     real(real32) :: maxarea=400._real32
     real(real32) :: vec,ang,area
     real(real32) :: ang_weight = 10._real32
     real(real32) :: area_weight = 100._real32
  end type tol_type

end module artemis__misc_types