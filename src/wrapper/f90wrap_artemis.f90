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

subroutine f90wrap_basis_type__crystals_equivalent__binding__basis_type(this, other, &
     tol, exact, allow_translation, has_plane, plane_normal, equivalent)
  use artemis__crystal_compare, only: crystals_equivalent
  use atomstruc, only: basis_type
  implicit none
  type basis_type_ptr_type
  type(basis_type), pointer :: p => NULL()
  end type basis_type_ptr_type
  integer, intent(in) :: this(2)
  integer, intent(in) :: other(2)
  type(basis_type_ptr_type) :: this_ptr
  type(basis_type_ptr_type) :: other_ptr
  real(4), intent(in) :: tol
  logical, intent(in) :: exact
  logical, intent(in) :: allow_translation
  logical, intent(in) :: has_plane
  real(4), intent(in) :: plane_normal(3)
  logical, intent(out) :: equivalent

  this_ptr = transfer(this, this_ptr)
  other_ptr = transfer(other, other_ptr)

  if (has_plane) then
     equivalent = crystals_equivalent(this_ptr%p, other_ptr%p, tol, exact=exact, &
          allow_translation=allow_translation, plane_normal=plane_normal)
  else
     equivalent = crystals_equivalent(this_ptr%p, other_ptr%p, tol, exact=exact, &
          allow_translation=allow_translation)
  end if
end subroutine f90wrap_basis_type__crystals_equivalent__binding__basis_type

! End of module artemis defined in file ../fortran/artemis.f90
