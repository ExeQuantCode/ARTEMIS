module artemis__structure_cache
  use artemis__geom_rw, only: basis_type
  implicit none

  private
  public :: store_last_generated_structures, retrieve_last_generated_structures

  type(basis_type), allocatable, dimension(:), save :: cached_structures

contains

  subroutine store_last_generated_structures(structures)
    implicit none
    type(basis_type), intent(in), allocatable :: structures(:)
    if (allocated(cached_structures)) deallocate(cached_structures)
    allocate(cached_structures(size(structures)))
    cached_structures = structures
  end subroutine store_last_generated_structures

  function retrieve_last_generated_structures() result(structures)
    implicit none
    type(basis_type), allocatable :: structures(:)
    if (.not.allocated(cached_structures)) then
        allocate(structures(0))
    else
        allocate(structures(size(cached_structures)))
        structures = cached_structures
    end if
  end function retrieve_last_generated_structures

end module artemis__structure_cache
