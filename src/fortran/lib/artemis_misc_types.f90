module artemis__misc_types
  !! Module containing custom derived types for ARTEMIS
  use artemis__constants, only: real32, pi
  use coreutils__string, only: to_lower
  use atomstruc, only: basis_type, geom_write
  use artemis__geom_utils, only: MATNORM
  implicit none


  private

  public :: struc_data_type
  public :: latmatch_type
  public :: tol_type
  public :: abstract_artemis_generator_type


  type struc_data_type
     !! Type for storing structure data associated with a generated interface.
     integer :: match_idx = 0
     !! Index of the lattice match.
     integer :: match_and_term_idx = 0
     !! Combined match and termination index.
     integer :: shift_idx = 0
     !! Index of the interface shift.
     integer :: swap_idx  = 0
     !! Index of the atomic swap configuration.
     logical :: from_pricel_lw = .false.
     !! Whether the lower slab originates from a primitive cell.
     logical :: from_pricel_up = .false.
     !! Whether the upper slab originates from a primitive cell.
     integer, dimension(2) :: term_lw_idx = 0
     !! Lower slab termination indices.
     integer, dimension(2) :: term_up_idx = 0
     !! Upper slab termination indices.
     real(real32), dimension(4) :: term_lw_bounds = 0._real32
     !! Lower slab termination bounds.
     real(real32), dimension(4) :: term_up_bounds = 0._real32
     !! Upper slab termination bounds.
     integer, dimension(2) :: term_lw_natom = 0
     !! Number of atoms in the lower slab termination.
     integer, dimension(2) :: term_up_natom = 0
     !! Number of atoms in the upper slab termination.
     integer, dimension(3,3) :: transform_lw = 0
     !! Transformation matrix for the lower slab.
     integer, dimension(3,3) :: transform_up = 0
     !! Transformation matrix for the upper slab.
     real(real32) :: approx_thickness_lw = 0._real32
     !! Approximate thickness of the lower slab.
     real(real32) :: approx_thickness_up = 0._real32
     !! Approximate thickness of the upper slab.
     real(real32), dimension(3) :: mismatch
     !! Lattice mismatch vector.
     real(real32), dimension(3) :: shift = 0._real32
     !! Interface shift vector.
     real(real32) :: swap_density = 0._real32
     !! Atomic swap density.
     real(real32), dimension(2) :: approx_eff_swap_conc = 0._real32
     !! Approximate effective swap concentration.

  end type struc_data_type

  interface struc_data_type
     module function init_struc_data_type( &
          match_idx, &
          match_and_term_idx, &
          from_pricel_lw, from_pricel_up, &
          term_lw_idx, term_up_idx, &
          term_lw_bounds, term_up_bounds, &
          term_lw_natom, term_up_natom, &
          transform_lw, transform_up, &
          approx_thickness_lw, approx_thickness_up, &
          mismatch, &
          shift_idx, shift, &
          swap_idx, swap_density, approx_eff_swap_conc &
     ) result(output)
       integer, intent(in) :: match_idx
       integer, intent(in) :: match_and_term_idx
       logical, intent(in) :: from_pricel_lw, from_pricel_up
       integer, dimension(2), intent(in) :: term_lw_idx, term_up_idx
       real(real32), dimension(4), intent(in) :: term_lw_bounds, term_up_bounds
       integer, dimension(2), intent(in) :: term_lw_natom, term_up_natom
       integer, dimension(3,3), intent(in) :: transform_lw, transform_up
       real(real32), intent(in) :: approx_thickness_lw, approx_thickness_up
       real(real32), dimension(3), intent(in) :: mismatch
       integer, intent(in), optional :: shift_idx
       real(real32), dimension(3), intent(in), optional :: shift
       integer, intent(in), optional :: swap_idx
       real(real32), intent(in), optional :: swap_density
       real(real32), dimension(2), intent(in), optional :: approx_eff_swap_conc
       type(struc_data_type) :: output
     end function init_struc_data_type
  end interface struc_data_type

  type latmatch_type
     !! Type for storing lattice matching data.
     integer :: nfit
     !! Number of fitted matches found.
     integer :: max_num_matches = 5
     !! Maximum number of lattice matches to store.
     logical :: reduce = .false.
     !! Whether to reduce matches.
     logical :: reduced = .false.
     !! Whether matches have been reduced.
     character(1) :: abc(3)= [ 'a', 'b', 'c' ]
     !! Lattice vector labels.

     integer, dimension(2) :: axes
     !! Axes to constrain for matching.
     integer, allocatable, dimension(:,:,:) :: tf1,tf2
     !! Transformation matrices for lower and upper lattices.
     real(real32), allocatable, dimension(:,:) :: tol
     !! Tolerance values for each match.
     real(real32), dimension(3,3) :: lat1,lat2
     !! Normalised lower and upper lattice matrices.
   contains
     procedure, pass(this) :: init => latmatch_init
     procedure, pass(this) :: constrain_axes
  end type latmatch_type

  type tol_type
     !! Type for storing lattice matching tolerances.
     integer :: maxfit = 100
     !! Maximum number of fits to evaluate.
     integer :: maxsize = 10
     !! Maximum supercell size multiplier.
     real(real32) :: maxlen  = 20._real32
     !! Maximum lattice vector length.
     real(real32) :: maxarea = 400._real32
     !! Maximum supercell area.
     real(real32) :: vec  = 5._real32 / 100._real32
     !! Vector length tolerance (fractional).
     real(real32) :: ang  = 1._real32 * pi / 180._real32
     !! Angle tolerance (radians).
     real(real32) :: area = 10._real32 / 100._real32
     !! Area tolerance (fractional).
     real(real32) :: ang_weight  = 10._real32
     !! Weight factor for angle mismatch.
     real(real32) :: area_weight = 100._real32
     !! Weight factor for area mismatch.
  end type tol_type

  type :: abstract_artemis_generator_type
     !! Abstract base type for ARTEMIS structure generators.
     integer :: num_structures = 0
     !! Number of generated structures.
     integer :: max_num_structures = 100
     !! Maximum number of structures to generate.
     
     integer :: axis = 3
     !! Axis along which to align the slab/interface normal vector.

     real(real32) :: vacuum_gap = 14._real32
     !! Vacuum thickness in Å.

     real(real32) :: warning_min_bond = 1.5_real32
     !! Minimum bond length to trigger a warning in Å.

     type(basis_type), dimension(:), allocatable :: structures
     !! Array of generated structures.
   contains
     procedure, pass(this) :: write_structures
     procedure, pass(this) :: get_structures
     procedure, pass(this) :: set_structures
  end type abstract_artemis_generator_type


contains
  
!###############################################################################
  module function init_struc_data_type( &
       match_idx, &
       match_and_term_idx, &
       from_pricel_lw, from_pricel_up, &
       term_lw_idx, term_up_idx, &
       term_lw_bounds, term_up_bounds, &
       term_lw_natom, term_up_natom, &
       transform_lw, transform_up, &
       approx_thickness_lw, approx_thickness_up, &
       mismatch, &
       shift_idx, shift, &
       swap_idx, swap_density, approx_eff_swap_conc &
  ) result(output)
    !! Initialise a struc_data_type instance.
    implicit none

    ! Arguments
    integer, intent(in) :: match_idx
    !! Index of the lattice match.
    integer, intent(in) :: match_and_term_idx
    !! Combined match and termination index.
    logical, intent(in) :: from_pricel_lw, from_pricel_up
    !! Whether the lower/upper slab originates from a primitive cell.
    integer, dimension(2), intent(in) :: term_lw_idx, term_up_idx
    !! Lower/upper slab termination indices.
    real(real32), dimension(4), intent(in) :: term_lw_bounds, term_up_bounds
    !! Lower/upper slab termination bounds.
    integer, dimension(2), intent(in) :: term_lw_natom, term_up_natom
    !! Number of atoms in the lower/upper slab termination.
    integer, dimension(3,3), intent(in) :: transform_lw, transform_up
    !! Transformation matrices for the lower/upper slab.
    real(real32), intent(in) :: approx_thickness_lw, approx_thickness_up
    !! Approximate thickness of the lower/upper slab.
    real(real32), dimension(3), intent(in) :: mismatch
    !! Lattice mismatch vector.
    integer, intent(in), optional :: shift_idx
    !! Index of the interface shift.
    real(real32), dimension(3), intent(in), optional :: shift
    !! Interface shift vector.
    integer, intent(in), optional :: swap_idx
    !! Index of the atomic swap configuration.
    real(real32), intent(in), optional :: swap_density
    !! Atomic swap density.
    real(real32), dimension(2), intent(in), optional :: approx_eff_swap_conc
    !! Approximate effective swap concentration.
    type(struc_data_type) :: output
    !! Initialised structure data instance.

    output%match_idx = match_idx
    output%match_and_term_idx = match_and_term_idx
    output%from_pricel_lw = from_pricel_lw
    output%from_pricel_up = from_pricel_up
    output%term_lw_idx = term_lw_idx
    output%term_up_idx = term_up_idx
    output%term_lw_bounds = term_lw_bounds
    output%term_up_bounds = term_up_bounds
    output%term_lw_natom = term_lw_natom
    output%term_up_natom = term_up_natom
    output%transform_lw = transform_lw
    output%transform_up = transform_up
    output%approx_thickness_lw = approx_thickness_lw
    output%approx_thickness_up = approx_thickness_up
    output%mismatch = mismatch

    if(present(shift)) output%shift = shift
    if(present(shift_idx)) output%shift_idx = shift_idx

    if(present(swap_idx)) output%swap_idx = swap_idx
    if(present(swap_density)) output%swap_density = swap_density
    if(present(approx_eff_swap_conc)) output%approx_eff_swap_conc = approx_eff_swap_conc

  end function init_struc_data_type
!###############################################################################


!###############################################################################
  subroutine latmatch_init( &
       this, tol, lattice_lw, lattice_up, max_num_matches, reduce_matches &
  )
    !! Initialise a latmatch_type instance.
    implicit none

    ! Arguments
    class(latmatch_type), intent(inout) :: this
    !! Instance of latmatch_type to initialise.
    type(tol_type), intent(in) :: tol
    !! Lattice matching tolerances.
    integer, intent(in) :: max_num_matches
    !! Maximum number of lattice matches.
    real(real32), dimension(3,3), intent(in) :: lattice_lw,lattice_up
    !! Lower and upper lattice matrices.
    logical, intent(in) :: reduce_matches
    !! Whether to reduce matches.

    this%max_num_matches = max_num_matches
    allocate(this%tf1(this%max_num_matches,3,3))
    allocate(this%tf2(this%max_num_matches,3,3))
    allocate(this%tol(this%max_num_matches,3))

    this%tol(:,:) = huge(0._real32)
    this%lat1 = MATNORM(lattice_lw)
    this%lat2 = MATNORM(lattice_up)

    this%reduce = reduce_matches

  end subroutine latmatch_init
!###############################################################################


!###############################################################################
  subroutine constrain_axes(this, miller_lw, miller_up, verbose)
    !! Constrain lattice matching axes based on Miller indices.
    implicit none

    ! Arguments
    class(latmatch_type), intent(inout) :: this
    !! Instance of latmatch_type.
    integer, dimension(3), intent(in) :: miller_lw, miller_up
    !! Miller indices for the lower and upper slabs.
    integer, intent(in) :: verbose
    !! Verbosity level.


    if(all(miller_lw.eq.0))then
       this%axes(1) = 3
       if(verbose.gt.0) write(*,*) &
            "Finding matches for all possible lower planes."
    else
       this%axes(1) = 2
       if(verbose.gt.0) write(*,*) "Finding matches for the lower ab plane."
    end if

    if(all(miller_up.eq.0))then
       this%axes(2) = 3
       if(verbose.gt.0) write(*,*) &
            "Finding matches for all possible upper planes."
    else
       this%axes(2) = 2
       if(verbose.gt.0) write(*,*) "Finding matches for the upper ab plane."
    end if

  end subroutine constrain_axes
!###############################################################################


!###############################################################################
  subroutine write_structures( &
       this, directory, prefix &
  )
    !! Write the generated terminations to file
    implicit none
   
    ! Arguments
    class(abstract_artemis_generator_type), intent(in) :: this
    !! Instance of artemis generator type
    character(len=*), intent(in) :: directory
    !! Directory to write the files to
    character(len=*), intent(in), optional :: prefix
    !! Prefix for the output files
   
    ! Local variables
    integer :: i
    !! Loop variable
    integer :: unit
    !! File unit number
    character(len=256) :: filename, filename_template
    !! File name for the output files
    character(len=:), allocatable :: prefix_
    !! Prefix for the output files



    if(trim(directory).ne."") then
       call system('mkdir -p '//trim(adjustl(directory)))
    end if

    filename_template = "POSCAR"
    if(present(prefix)) then
       prefix_ = trim(to_lower(prefix))
       filename_template = trim(filename_template) // "_" // trim(prefix_)
    end if
    if(allocated(this%structures))then
       do i = 1, size(this%structures)
          write(filename,'(A,I0)') trim(filename_template), i
          if(trim(directory).ne."") then
             filename = trim(directory) // "/" // trim(filename)
          end if
          open(newunit=unit,file=filename)
          call geom_write(unit, this%structures(i))
          close(unit)
       end do
    else
       write(0,'(1X,"No structures to write.")')
    end if
   
  end subroutine write_structures
!###############################################################################


!###############################################################################
  function get_structures(this) result(structures)
    !! Get the generated structures.
    implicit none
    ! Arguments
    class(abstract_artemis_generator_type), intent(in) :: this
    !! Instance of the artemis generator.
    type(basis_type), dimension(:), allocatable :: structures
    !! Generated structures.

    structures = this%structures
  end function get_structures
!###############################################################################


!###############################################################################
  subroutine set_structures(this, structures)
    !! Set the generated structures.
    implicit none
    ! Arguments
    class(abstract_artemis_generator_type), intent(inout) :: this
    !! Instance of the artemis generator.
    type(basis_type), dimension(:), allocatable :: structures
    !! Generated structures.

    this%structures = structures
    this%num_structures = size(structures)
  end subroutine set_structures
!###############################################################################

end module artemis__misc_types
