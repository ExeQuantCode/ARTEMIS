module artemis__structure_cache
  !! Module for caching generated structures.
  !!
  !! Provides save/retrieve operations for the most recently generated
  !! interface structures, enabling retrieval after generation.
  use artemis__geom_rw, only: basis_type
  implicit none

  private
  public :: store_last_generated_structures, retrieve_last_generated_structures

  type(basis_type), allocatable, dimension(:), save :: cached_structures
  !! Array of cached structures from the last generation run.

contains

!###############################################################################
  subroutine store_last_generated_structures(structures)
    !! Store generated structures in the module-level cache.
    implicit none

    ! Arguments
    type(basis_type), intent(in), allocatable :: structures(:)
    !! Array of structures to cache.

    if (allocated(cached_structures)) deallocate(cached_structures)
    allocate(cached_structures(size(structures)))
    cached_structures = structures
  end subroutine store_last_generated_structures
!###############################################################################


!###############################################################################
  function retrieve_last_generated_structures() result(structures)
    !! Retrieve cached structures from the last generation run.
    implicit none

    ! Local variables
    type(basis_type), allocatable :: structures(:)
    !! Returned array of cached structures.

    if (.not.allocated(cached_structures)) then
        allocate(structures(0))
    else
        allocate(structures(size(cached_structures)))
        structures = cached_structures
    end if
  end function retrieve_last_generated_structures
!###############################################################################

end module artemis__structure_cache
