! Module artemis defined in file ../fortran/artemis.f90

subroutine f90wrap_get_suppress_warnings(suppress_warnings)
  use artemis, only: artemis__suppress_warnings
  implicit none
  logical, intent(out) :: suppress_warnings

  suppress_warnings = artemis__suppress_warnings
end subroutine f90wrap_get_suppress_warnings

subroutine f90wrap_set_suppress_warnings(suppress_warnings)
  use artemis, only: artemis__suppress_warnings
  implicit none
  logical, intent(in) :: suppress_warnings

  artemis__suppress_warnings = suppress_warnings
end subroutine f90wrap_set_suppress_warnings

subroutine f90wrap_get_interface_translations(structure, t1, t2)
  use artemis, only: get_interface_translations
  use atomstruc, only: basis_type
  implicit none
  type basis_type_ptr_type
  type(basis_type), pointer :: p => NULL()
  end type basis_type_ptr_type
  integer, intent(in) :: structure(2)
  type(basis_type_ptr_type) :: structure_ptr
  real(4), intent(out) :: t1(3)
  real(4), intent(out) :: t2(3)

  structure_ptr = transfer(structure, structure_ptr)
  call get_interface_translations(structure_ptr%p, t1, t2)
end subroutine f90wrap_get_interface_translations

! End of module artemis defined in file ../fortran/artemis.f90
