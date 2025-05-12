!!!#############################################################################
!!! INTERFACES CARD SUBROUTINES
!!! Code written by Ned Thaddeus Taylor and Isiah Edward Mikel Rudkin
!!! Code part of the ARTEMIS group (Hepplestone research group).
!!! Think Hepplestone, think HRG.
!!!#############################################################################
module artemis__generator
  use artemis__constants,     only: real32, pi
  use artemis__misc,          only: to_lower, to_upper
  use artemis__misc_types,    only: abstract_artemis_generator_type, &
       latmatch_type, tol_type, struc_data_type
  use artemis__geom_rw,       only: basis_type
  use lat_compare,            only: lattice_matching, cyc_lat1
  use artemis__io_utils,      only: err_abort, print_warning, stop_program
  use artemis__io_utils_extd, only: err_abort_print_struc
  use misc_linalg,            only: uvec,modu,get_area,inverse,cross
  use artemis__interface_identifier,   only: intf_info_type,&
       get_interface,get_layered_axis,gen_DON
  use artemis__geom_utils,    only: planecutter, primitive_lat, ortho_axis,&
       shift_region, set_vacuum, transformer, shifter, reducer, &
       get_min_bulk_bond, get_min_bond, get_shortest_bond, bond_type, &
       share_strain, MATNORM, basis_stack, compare_stoichiometry, &
       get_primitive_cell
  use artemis__sym,           only: confine_type, gldfnd
  use artemis__terminations,  only: get_termination_info, term_arr_type, &
       set_layer_tol, build_slab_supercell, cut_slab_to_height
  use swapping,               only: rand_swapper
  use shifting !!! CHANGE TO SHIFTER?
  implicit none


  private

  public :: artemis_generator_type


  type, extends(abstract_artemis_generator_type) :: artemis_generator_type
    !! Interface generator type
    type(basis_type) :: structure_lw, structure_up
    !! Lower and upper bulk structures
    real(real32), dimension(:), allocatable :: elastic_constants_lw, elastic_constants_up
    !! Elastic constants for the lower and upper bulk structures
    logical :: use_pricel_lw = .true., use_pricel_up = .true.
    !! Use primitive cell for lower and upper bulk structures
    logical :: require_stoichiometry_lw = .false., &
         require_stoichiometry_up = .false.
    !! Boolean whether to require stoichiometry for the lower and upper bulk structures
    
    integer, dimension(3) :: miller_lw = [ 0, 0, 0 ], miller_up = [ 0, 0, 0 ]
    !! Miller indices for the lower and upper bulk structures
    logical :: is_layered_lw = .false., is_layered_up = .false.
    !! Boolean whether the lower and upper bulk structures are layered
    logical :: ludef_is_layered_lw = .false., ludef_is_layered_up = .false.
    !! Boolean whether the user defined whether to use layered structures

    integer :: shift_method = 4
    !! Shift method
    integer :: num_shifts = 5
    !! Number of shifts per lattice match
    real(real32), dimension(:,:), allocatable :: shifts
    !! Shift values
    real(real32) :: interface_depth = 1.5_real32
    !! Interface depth
    real(real32) :: separation_scale = 1._real32
    !! Separation scale
    integer :: depth_method = 0
    !! Method for determining the depth to which consider atoms from interface

    type(struc_data_type), dimension(:), allocatable :: structure_data
    !! Structure data

    integer :: swap_method = 0
    !! Swap method
    integer :: num_swaps = 0
    !! Number of swaps per shifted interface
    real(real32) :: swap_density = 5.E-2_real32
    !! Swap density
    real(real32) :: swap_depth = 3._real32
    !! Swap depth
    real(real32) :: swap_sigma = -1._real32
    !! Swap sigma
    logical :: require_mirror_swaps = .true.
    !! Require mirror swaps

    integer :: match_method = 0
    !! Match method
    integer :: max_num_matches = 5
    !! Maximum number of matches
    integer :: max_num_terms = 5
    !! Maximum number of terminations
    integer :: max_num_planes = 10
    !! Maximum number of planes

    logical :: compensate_normal = .true.
    !! Compensate mismatch strain by adjusting the axes parallel to the interface normal vector
    !! Compensate = false = strained
    !! Compensate = true = relaxed (compensate for interfacial strain by extending/compressing)
    
    real(real32) :: bondlength_cutoff = 6._real32
    !! Maximum bond length cutoff for the bulk structures
    real(real32), dimension(2) :: layer_separation_cutoff = 1._real32
    !! Minimum separation between layers

    type(tol_type) :: tolerance
    !! Tolerance structure
    real(real32) :: tol_sym = 1.E-6_real32
    !! Tolerance for symmetry operations

   contains
    procedure, pass(this) :: get_all_structures_data
    !! Get the structure data for all structures
    procedure, pass(this) :: get_structure_data
    !! Get the structure data for a specific structure
    procedure, pass(this) :: get_all_structures_mismatch
    !! Get the mismatch data for all structures
    procedure, pass(this) :: get_structure_mismatch
    !! Get the mismatch data for a specific structure
    procedure, pass(this) :: get_all_structures_transform
    !! Get the structure data for a specific structure
    procedure, pass(this) :: get_structure_transform
    !! Get the structure data for a specific structure
    procedure, pass(this) :: get_all_structures_shift
    !! Get the shifts for all structures
    procedure, pass(this) :: get_structure_shift
    !! Get the shifts for a specific structure

    procedure, pass(this) :: write_match_and_term_data
    !! Write the match and termination data to a file
    procedure, pass(this) :: write_shift_data
    !! Write the shift data to a file

    procedure, pass(this) :: set_tolerance
    !! Set tolerance for identifying good lattice matches
    procedure, pass(this) :: set_shift_method
    !! Set the shift method and associated data
    procedure, pass(this) :: set_swap_method
    !! Set the swap method and associated data
    procedure, pass(this) :: set_match_method
    !! Set the lattice match method and associated data

    procedure, pass(this) :: set_materials
    !! Set the input materials for the interface generator
    procedure, pass(this) :: set_surface_properties
    !! Set the surface properties for the interface generator
    procedure, pass(this) :: reset_is_layered_lw
    !! Reset the is_layered flags for the lower bulk structure
    procedure, pass(this) :: reset_is_layered_up
    !! Reset the is_layered flags for the upper bulk structure

    procedure, pass(this) :: get_terminations
    !! Return the terminations for structure
    procedure, pass(this) :: get_interface_location
    !! Get the interface location for the given structure

    procedure, pass(this) :: generate => generate_interfaces
    !! Generate interfaces from two bulk structures
    procedure, pass(this) :: restart => generate_interfaces_from_existing
    !! Generate interfaces from existing bulk structures
    procedure, pass(this) :: generate_perturbations => generate_shifts_and_swaps
    !! Generate perturbations for the given basis
  end type artemis_generator_type

contains

!###############################################################################
  function get_all_structures_data(this) result(output)
    !! Get the structure data for all structures
    implicit none

    ! Arguments
    class(artemis_generator_type), intent(in) :: this
    !! Instance of artemis generator type

    type(struc_data_type), dimension(this%num_structures) :: output
    !! Structure data

    ! Local variables
    integer :: i

    do i = 1, this%num_structures
       output(i) = this%structure_data(i)
    end do

  end function get_all_structures_data
!###############################################################################


!###############################################################################
  function get_structure_data(this, idx) result(output)
    !! Get the structure data for a specific structure
    implicit none

    ! Arguments
    class(artemis_generator_type), intent(in) :: this
    !! Instance of artemis generator type
    integer, intent(in) :: idx
    !! Index of the structure

    type(struc_data_type) :: output
    !! Structure data

    output = this%structure_data(idx)

  end function get_structure_data
!###############################################################################


!###############################################################################
  function get_all_structures_mismatch(this) result(output)
    !! Get the mismatch data for all structures
    implicit none

    ! Arguments
    class(artemis_generator_type), intent(in) :: this
    !! Instance of artemis generator type

    real(real32), dimension(3,this%num_structures) :: output
    !! Mismatch data

    ! Local variables
    integer :: i

    do i = 1, this%num_structures
       output(:,i) = this%structure_data(i)%mismatch
    end do

  end function get_all_structures_mismatch
!###############################################################################


!###############################################################################
  function get_structure_mismatch(this, idx) result(output)
    !! Get the mismatch data for a specific structure
    implicit none

    ! Arguments
    class(artemis_generator_type), intent(in) :: this
    !! Instance of artemis generator type
    integer, intent(in) :: idx
    !! Index of the structure
      
    real(real32), dimension(3) :: output
    !! Mismatch data

    output = this%structure_data(idx)%mismatch

  end function get_structure_mismatch
!###############################################################################


!###############################################################################
  function get_all_structures_transform(this) result(output)
    !! Get the structure data for a specific structure
    implicit none

    ! Arguments
    class(artemis_generator_type), intent(in) :: this
    !! Instance of artemis generator type

    integer, dimension(3,3,2,this%num_structures) :: output
    !! Transformation data

    ! Local variables
    integer :: i
    ! Loop over all structures

    do i = 1, this%num_structures
       output(:,:,1,i) = this%structure_data(i)%transform_lw
       output(:,:,2,i) = this%structure_data(i)%transform_up
    end do

  end function get_all_structures_transform
!###############################################################################


!###############################################################################
  function get_structure_transform(this, idx) result(output)
    !! Get the structure data for a specific structure
    implicit none

    ! Arguments
    class(artemis_generator_type), intent(in) :: this
    !! Instance of artemis generator type
    integer, intent(in) :: idx
    !! Index of the structure

    integer, dimension(3,3,2) :: output
    !! Transformation data

    output(:,:,1) = this%structure_data(idx)%transform_lw
    output(:,:,2) = this%structure_data(idx)%transform_up

  end function get_structure_transform
!###############################################################################


!###############################################################################
  function get_all_structures_shift(this) result(output)
    !! Get the shifts for all structures
    implicit none

    ! Arguments
    class(artemis_generator_type), intent(in) :: this
    !! Instance of artemis generator type

    real(real32), dimension(3,this%num_structures) :: output
    !! Shift data

    ! Local variables
    integer :: i

    do i = 1, this%num_structures
       output(:,i) = this%structure_data(i)%shift
    end do

  end function get_all_structures_shift
!###############################################################################


!###############################################################################
   function get_structure_shift(this, idx) result(output)
      !! Get the shifts for a specific structure
      implicit none
   
      ! Arguments
      class(artemis_generator_type), intent(in) :: this
      !! Instance of artemis generator type
      integer, intent(in) :: idx
      !! Index of the structure
   
      real(real32), dimension(3) :: output
      !! Shift data
   
      output = this%structure_data(idx)%shift
   
   end function get_structure_shift
!###############################################################################


!###############################################################################
  subroutine set_tolerance( &
       this, &
       tolerance, &
       vector_mismatch, angle_mismatch, area_mismatch, &
       max_length, max_area, max_fit, max_extension, &
       angle_weight, area_weight &
  )
    !! Set tolerance for the best match
    implicit none

    ! Arguments
    class(artemis_generator_type), intent(inout) :: this
    !! Instance of artemis generator type
    type(tol_type), intent(in), optional :: tolerance
    !! Tolerance structure
    real(real32), intent(in), optional :: vector_mismatch
    !! Tolerance for the vector mismatch
    real(real32), intent(in), optional :: angle_mismatch
    !! Tolerance for the angle mismatch
    real(real32), intent(in), optional :: area_mismatch
    !! Tolerance for the area mismatch
    real(real32), intent(in), optional :: max_length
    !! Maximum allowed length of a lattice vector
    real(real32), intent(in), optional :: max_area
    !! Maximum allowed area parallel to the surface
    integer, intent(in), optional :: max_fit
    !! Maximum allowed number of matches for each individial ... ???? area mapped out on a plane
    integer, intent(in), optional :: max_extension
    !! Maximum allowed integer extension of each lattice vector
    real(real32), intent(in), optional :: angle_weight
    !! Importance weighting of angle mismatch
    real(real32), intent(in), optional :: area_weight
    !! Importance weighting of area mismatch

    if(present(tolerance))then
       this%tolerance = tolerance
    else
       if(present(vector_mismatch)) this%tolerance%vec = vector_mismatch
       if(present(angle_mismatch)) this%tolerance%ang = angle_mismatch
       if(present(area_mismatch)) this%tolerance%area = area_mismatch
       if(present(max_length)) this%tolerance%maxlen = max_length
       if(present(max_area)) this%tolerance%maxarea = max_area
       if(present(max_fit)) this%tolerance%maxfit = max_fit
       if(present(max_extension)) this%tolerance%maxsize = max_extension
       if(present(angle_weight)) this%tolerance%ang_weight = angle_weight
       if(present(area_weight)) this%tolerance%area_weight = area_weight
    end if

    !!! TOLERANCE EXPECTED IN FRACTIONS OF Å, radians, and Å^2

  end subroutine set_tolerance
!###############################################################################


!###############################################################################
  subroutine set_shift_method( &
       this, &
       method, num_shifts, shifts, &
       interface_depth, separation_scale, depth_method, &
       bondlength_cutoff &
  )
    !! Set the shift method
    implicit none

    ! Arguments
    class(artemis_generator_type), intent(inout) :: this
    !! Instance of artemis generator type
    integer, intent(in), optional :: method
    !! Shift method
    integer, intent(in), optional :: num_shifts
    !! Number of shifts
    real(real32), dimension(..), intent(in), optional :: shifts
    !! Shift values
    real(real32), intent(in), optional :: interface_depth
    !! Interface depth
    real(real32), intent(in), optional :: separation_scale
    !! Separation scale
    integer, intent(in), optional :: depth_method
    !! Method for determining the depth to which consider atoms from interface
    real(real32), intent(in), optional :: bondlength_cutoff
    !! Bond length cutoff for the bulk structures

    ! Local variables
    character(len=256) :: err_msg

    if(present(method)) this%shift_method = method
    if(present(num_shifts)) this%num_shifts = num_shifts
    if(present(interface_depth)) this%interface_depth = interface_depth
    if(present(separation_scale)) this%separation_scale = separation_scale
    if(present(depth_method)) this%depth_method = depth_method
    if(present(bondlength_cutoff)) this%bondlength_cutoff = bondlength_cutoff
    if(present(shifts)) then
       if(allocated(this%shifts)) deallocate(this%shifts)
       select rank(shifts)
       rank(0)
          allocate(this%shifts(1,3))
          this%shifts(1,this%axis) = shifts
       rank(1)
          allocate(this%shifts(1,3))
          select case(size(shifts,dim=1))
          case(1)
             this%shifts(1,this%axis) = shifts(1)
          case(3)
             this%shifts(1,:) = shifts
          case default
             ! check if length of shifts is divisible by 3
             if(mod(size(shifts,dim=1),3).eq.0) then
               allocate(this%shifts(size(shifts,dim=1)/3,3))
               this%shifts = reshape(shifts, [ size(shifts,dim=1)/3,3 ])
             else
                write(err_msg,'(A,I0,A)') &
                     "The shifts vector has ", size(shifts, dim=1), &
                     " components. It should have 1 or 3."
                call stop_program(trim(err_msg))
                return
             end if
          end select
       rank(2)
          select case(size(shifts,dim=2))
          case(1)
             allocate(this%shifts(size(shifts,1),3))
             this%shifts(:,3) = shifts(:,1)
          case(3)
             allocate(this%shifts(size(shifts,1),3))
             this%shifts = shifts
          case default
             write(err_msg,'(A,I0,A)') &
                  "The shifts argument was improperly defined."
             call stop_program(trim(err_msg))
             return
          end select
       rank default
          write(err_msg,'(A,I0,A)') &
               "The shifts vector has ", size(shifts, dim=1), &
               " components. It should have 1, 2, or 3."
          call stop_program(trim(err_msg))
          return
       end select
    else
       if(allocated(this%shifts)) deallocate(this%shifts)
       allocate(this%shifts(1,3), source = -1._real32)
    end if

  end subroutine set_shift_method
!###############################################################################


!###############################################################################
  subroutine set_swap_method( &
       this, method, num_swaps, swap_density, swap_depth, swap_sigma, &
       require_mirror_swaps &
  )
    !! Set the swap method
    implicit none

    ! Arguments
    class(artemis_generator_type), intent(inout) :: this
    !! Instance of artemis generator type
    integer, intent(in), optional :: method
    !! Swap method
    integer, intent(in), optional :: num_swaps
    !! Number of swaps
    real(real32), intent(in), optional :: swap_density
    !! Swap density
    real(real32), intent(in), optional :: swap_depth
    !! Swap depth
    real(real32), intent(in), optional :: swap_sigma
    !! Swap sigma
    logical, intent(in), optional :: require_mirror_swaps
    !! Require mirror swaps

    if(present(method)) this%swap_method = method
    if(present(num_swaps)) this%num_swaps = num_swaps
    if(present(swap_density)) this%swap_density = swap_density
    if(present(swap_depth)) this%swap_depth = swap_depth
    if(present(swap_sigma)) this%swap_sigma = swap_sigma
    if(present(require_mirror_swaps)) &
         this%require_mirror_swaps = require_mirror_swaps

  end subroutine set_swap_method
!###############################################################################


!###############################################################################
  subroutine set_match_method( &
       this, method, max_num_matches, max_num_terms, max_num_planes, &
       compensate_normal &
  )
    !! Set the lattice match method
    implicit none
   
    ! Arguments
    class(artemis_generator_type), intent(inout) :: this
    !! Instance of artemis generator type
    integer, intent(in), optional :: method
    !! Match method
    integer, intent(in), optional :: max_num_matches
    !! Maximum number of matches
    integer, intent(in), optional :: max_num_terms
    !! Maximum number of terminations
    integer, intent(in), optional :: max_num_planes
    !! Maximum number of planes
    logical, intent(in), optional :: compensate_normal
    !! Compensate mismatch strain by adjusting the axes parallel to the interface normal vector
   
    if(present(method)) this%match_method = method
    if(present(max_num_matches)) this%max_num_matches = max_num_matches
    if(present(max_num_terms)) this%max_num_terms = max_num_terms
    if(present(max_num_planes)) this%max_num_planes = max_num_planes
    if(present(compensate_normal)) this%compensate_normal = compensate_normal
   
  end subroutine set_match_method
!###############################################################################


!###############################################################################
  subroutine set_materials( &
       this, structure_lw, structure_up, &
       elastic_constants_lw, elastic_constants_up, &
       use_pricel_lw, use_pricel_up &
  )
    !! Set the materials for the interface generator
    implicit none

    ! Arguments
    class(artemis_generator_type), intent(inout) :: this
    !! Instance of artemis generator type
    type(basis_type), intent(in), optional :: structure_lw
    !! Lower bulk structure
    type(basis_type), intent(in), optional :: structure_up
    !! Upper bulk structure
    real(real32), dimension(:), intent(in), optional :: elastic_constants_lw
    !! Elastic constants for the lower bulk structure
    real(real32), dimension(:), intent(in), optional :: elastic_constants_up
    !! Elastic constants for the upper bulk structure
    logical, intent(in), optional :: use_pricel_lw
    !! Use primitive cell for lower bulk structure
    logical, intent(in), optional :: use_pricel_up


    if(present(structure_lw))then
       if(structure_lw%natom.gt.0) call this%structure_lw%copy(structure_lw, length=4)
    end if
    if(present(structure_up))then
       if(structure_up%natom.gt.0) call this%structure_up%copy(structure_up, length=4)
    end if

    !---------------------------------------------------------------------------
    ! Handle the elastic constants
    !---------------------------------------------------------------------------
    if(present(elastic_constants_lw))then
       if(allocated(this%elastic_constants_lw)) deallocate(this%elastic_constants_lw)
       allocate(this%elastic_constants_lw(size(elastic_constants_lw)))
       this%elastic_constants_lw = elastic_constants_lw
    end if
    if(present(elastic_constants_up))then
       if(allocated(this%elastic_constants_up)) deallocate(this%elastic_constants_up)
       allocate(this%elastic_constants_up(size(elastic_constants_up)))
       this%elastic_constants_up = elastic_constants_up
    end if

    if(present(use_pricel_lw)) this%use_pricel_lw = use_pricel_lw
    if(present(use_pricel_up)) this%use_pricel_up = use_pricel_up


  end subroutine set_materials
!###############################################################################


!###############################################################################
  subroutine set_surface_properties( &
       this, &
       miller_lw, miller_up, &
       is_layered_lw, is_layered_up, &
       require_stoichiometry_lw, require_stoichiometry_up, &
       layer_separation_cutoff_lw, layer_separation_cutoff_up, &
       layer_separation_cutoff, &
       vacuum_gap &
  )
    !! Set the surface properties for the interface generator
    implicit none
   
    ! Arguments
    class(artemis_generator_type), intent(inout) :: this
    !! Instance of artemis generator type
    integer, dimension(3), intent(in), optional :: miller_lw
    !! Miller indices for the lower bulk structure
    integer, dimension(3), intent(in), optional :: miller_up
    !! Miller indices for the upper bulk structure
   
    logical, intent(in), optional :: is_layered_lw
    !! Boolean whether the lower bulk structure is layered
    logical, intent(in), optional :: is_layered_up
    !! Boolean whether the upper bulk structure is layered

    logical, intent(in), optional :: require_stoichiometry_lw
    !! Boolean whether to require stoichiometry for the lower bulk structure
    logical, intent(in), optional :: require_stoichiometry_up
    !! Boolean whether to require stoichiometry for the upper bulk structure

    real(real32), intent(in), optional :: layer_separation_cutoff_lw
    !! Layer separation cutoff for the lower bulk structure
    real(real32), intent(in), optional :: layer_separation_cutoff_up
    !! Layer separation cutoff for the upper bulk structure
    real(real32), dimension(..), intent(in), optional :: layer_separation_cutoff
    !! Layer separation cutoff

    real(real32), intent(in), optional :: vacuum_gap
    !! Vacuum gap for termination generator

    ! Local variables
    character(len=256) :: err_msg
    !! Error message


    if(present(miller_lw)) this%miller_lw = miller_lw
    if(present(miller_up)) this%miller_up = miller_up

    if(present(is_layered_lw))then
       this%is_layered_lw = is_layered_lw
       this%ludef_is_layered_lw = .true.
    end if
    if(present(is_layered_up))then
       this%is_layered_up = is_layered_up
       this%ludef_is_layered_up = .true.
    end if

    if(present(require_stoichiometry_lw)) &
         this%require_stoichiometry_lw = require_stoichiometry_lw
    if(present(require_stoichiometry_up)) &
         this%require_stoichiometry_up = require_stoichiometry_up

    if(present(vacuum_gap)) this%vacuum_gap = vacuum_gap

    if(present(layer_separation_cutoff_lw)) &
         this%layer_separation_cutoff(1) = layer_separation_cutoff_lw
    if(present(layer_separation_cutoff_up)) &
         this%layer_separation_cutoff(2) = layer_separation_cutoff_up

    if( ( present(layer_separation_cutoff_lw) .or. &
         present(layer_separation_cutoff_up) ) .and. &
         present(layer_separation_cutoff) ) then
       write(err_msg,'(A)') &
            "The layer separation cutoff is defined in two ways. Please use only one."
       call stop_program(trim(err_msg))
       return
    elseif(present(layer_separation_cutoff))then
       select rank(layer_separation_cutoff)
       rank(0)
          this%layer_separation_cutoff(:) = layer_separation_cutoff
       rank(1)
          select case(size(layer_separation_cutoff,dim=1))
          case(1)
             this%layer_separation_cutoff = layer_separation_cutoff(1)
          case(2)
             this%layer_separation_cutoff = layer_separation_cutoff
          case default
             write(err_msg,'(A,I0,A)') &
                  "The layer separation cutoff vector has ", &
                  size(layer_separation_cutoff,dim=1), &
                  " components. It should have 1 or 2."
             call stop_program(trim(err_msg))
             return
          end select
       rank default
          write(err_msg,'(A,I0,A)') &
               "The layer separation cutoff only accepts rank 0 or 1."
          call stop_program(trim(err_msg))
          return
       end select
    end if
    if(any(this%layer_separation_cutoff.lt.1.E-2_real32))then
       write(err_msg,'(A,I0,A)') &
            "A layer separation this small is not realistic: ", &
            this%layer_separation_cutoff
       call stop_program(trim(err_msg))
       return
    end if
   
   end subroutine set_surface_properties
!###############################################################################


!###############################################################################
  subroutine reset_is_layered_lw(this)
   !! Reset the is_layered flags
   implicit none

   ! Arguments
   class(artemis_generator_type), intent(inout) :: this
   !! Instance of artemis generator type

   this%is_layered_lw = .false.
   this%ludef_is_layered_lw = .false.

  end subroutine reset_is_layered_lw
!###############################################################################


!###############################################################################
  subroutine reset_is_layered_up(this)
   !! Reset the is_layered flags
   implicit none

   ! Arguments
   class(artemis_generator_type), intent(inout) :: this
   !! Instance of artemis generator type

   this%is_layered_up = .false.
   this%ludef_is_layered_up = .false.

  end subroutine reset_is_layered_up
!###############################################################################


!###############################################################################
  function get_terminations( &
       this, identifier, miller, surface, num_layers, thickness, &
       orthogonalise, normalise, break_on_fail, &
       print_termination_info, verbose, exit_code &
  ) result(output)
    !! Generate and prints terminations parallel to the supplied miller plane
    implicit none

    ! Arguments
    class(artemis_generator_type), intent(inout) :: this
    !! Instance of artemis generator type
    integer, intent(in) :: identifier
    !! Identifier for the material (1=lower, 2=upper)
    integer, dimension(3), intent(in), optional :: miller
    !! Miller plane
    integer, dimension(:), intent(in), optional :: surface
    !! Surface termination indices
    integer, intent(in), optional :: num_layers
    !! Number of layers in the slab
    real(real32), intent(in), optional :: thickness
    !! Thickness of the slab (in Å)
    logical, intent(in), optional :: orthogonalise
    !! Boolean whether to orthogonalise the lattice
    logical, intent(in), optional :: normalise
    !! Boolean whether to normalise the lattice and basis
    logical, intent(in), optional :: break_on_fail
    !! Boolean whether to break on failure
    logical, intent(in), optional :: print_termination_info
    !! Boolean whether to print termination information
    integer, intent(in), optional :: verbose
    !! Boolean whether to print verbose output
    integer, intent(out), optional :: exit_code
    !! Exit code for the program

    type(basis_type), dimension(:), allocatable :: output
    !! Output structures

    ! Local variables
    integer :: itmp1, iterm, term_start, term_end, term_step, i
    !! Termination loop variables
    integer :: num_cells, ntrans
    !! Number of cells in the slab
    integer :: num_structures
    !! Number of structures to be generated
    integer, dimension(2) :: surface_
    !! Surface termination indices
    integer, dimension(3) :: miller_
    !! Miller plane
    integer :: num_layers_
    !! Number of layers in the slab
    real(real32) :: height, thickness_
    !! Height of the slab
    logical :: lcycle
    !! Boolean whether to cycle through the slab
    type(basis_type) :: structure, structure_compare
    !! Temporary basis structures
    type(confine_type) :: confine
    !! Confine structure along the specified axis
    type(term_arr_type) :: term
    !! List of terminations
    real(real32), dimension(3,3) :: tfmat
    !! Transformation matrix
    logical :: orthogonalise_
    !! Boolean whether to orthogonalise the lattice
    logical :: normalise_
    !! Boolean whether to normalise the lattice
    logical :: break_on_fail_
    !! Boolean whether to break on failure
    logical :: print_termination_info_
    !! Boolean whether to print termination information


    real(real32) :: layer_sep
    character(len=2) :: prefix
    character(len=256) :: warn_msg, err_msg
    integer :: exit_code_
    !! Exit code for the program
    integer :: verbose_
    !! Verbosity level

    integer, allocatable, dimension(:,:,:) :: bas_map,t1bas_map
    real(real32), allocatable, dimension(:,:) :: trans


    !---------------------------------------------------------------------------
    ! Initialise variables
    !---------------------------------------------------------------------------
    exit_code_ = 0
    verbose_ = 0
    print_termination_info_ = .true.
    if(present(verbose)) verbose_ = verbose
    if(present(print_termination_info)) &
         print_termination_info_ = print_termination_info


    !---------------------------------------------------------------------------
    ! Handle identifier
    !---------------------------------------------------------------------------
    select case(identifier)
    case(1)
       call structure%copy(this%structure_lw, length=4)
       call structure_compare%copy(this%structure_lw, length=4)
       if(this%use_pricel_lw)then
          if(verbose_.gt.0) write(*,'(1X,"Using primitive cell for material")')
          call get_primitive_cell(structure, tol_sym=this%tol_sym)
       end if
       miller_ = this%miller_lw
       prefix = "lw"
       layer_sep = this%layer_separation_cutoff(1)
    case(2)
       call structure%copy(this%structure_up, length=4)
       call structure_compare%copy(this%structure_up, length=4)
       if(this%use_pricel_up)then
          if(verbose_.gt.0) write(*,'(1X,"Using primitive cell for material")')
          call get_primitive_cell(structure, tol_sym=this%tol_sym)
       end if
       miller_ = this%miller_up
       prefix = "up"
       layer_sep = this%layer_separation_cutoff(2)
    case default
       write(err_msg,'(A,I0,A)') &
            "The identifier for the material is not valid: ", identifier
       call stop_program(trim(err_msg))
       return
    end select
    ! check if the structures have anything (i.e. atoms) in them
    if(structure%natom.eq.0)then
       write(err_msg,'(A,I0,A)') &
            "The structure has ", structure%natom, &
            " atoms. It should have at least 1."
       call stop_program(trim(err_msg))
       return
    end if


    ! set thickness if provided by user
    thickness_ = -1._real32
    num_layers_ = 0
    if(present(num_layers)) num_layers_ = num_layers
    if(present(thickness)) thickness_ = thickness
    if(num_layers_.eq.0.and.abs(thickness_+1._real32).lt.1.E-6_real32)then
       thickness_ = 10._real32
    elseif(num_layers_.le.0.and.thickness_.le.0._real32)then
       write(err_msg,'(A,I0,A)') &
            "The number of layers for the material is ", &
            num_layers_, " and the thickness is ", thickness_, &
            " One of these must be greater than 0."
       call stop_program(trim(err_msg))
       exit_code_ = 1
       return
    end if


    !---------------------------------------------------------------------------
    ! Handle the miller plane
    !---------------------------------------------------------------------------
    if(present(miller)) miller_ = miller
    if(all(miller_.eq.0))then
       write(err_msg,'(A,I0,A)') &
            "The miller plane is not valid: ", identifier
       call stop_program(trim(err_msg))
       exit_code_ = 1
       return
    end if


    orthogonalise_ = .true.
    if(present(orthogonalise)) orthogonalise_ = orthogonalise
    break_on_fail_ = .false.
    if(present(break_on_fail)) break_on_fail_ = break_on_fail
    normalise_ = .true.
    if(present(normalise)) normalise_ = normalise
    surface_ = 0
    if(present(surface))then
       select case(size(surface,dim=1))
       case(1)
          surface_(:) = surface(1)
       case(2)
          surface_ = surface
       case default
          write(err_msg,'(A,I0,A)') &
               "The surface termination indices have ", size(surface,dim=1), &
               " components. It should have 1 or 2."
          exit_code_ = 1
          call stop_program( &
               trim(err_msg), &
               exit_code=exit_code_, &
               block_stop = present(exit_code) &
          )
          return
       end select
    end if

    !! copy lattice and basis for manipulating
    allocate(bas_map(structure%nspec,maxval(structure%spec(:)%num,dim=1),2))
    bas_map = -1


    if(verbose_.gt.0) write(*,'(1X,"Using supplied plane...")')
    tfmat = planecutter(structure%lat,real(miller_,real32))
    call transformer(structure,tfmat,bas_map)


    !---------------------------------------------------------------------------
    ! Finds smallest thickness of the slab and increases to ...
    ! ... user-defined thickness
    !---------------------------------------------------------------------------
    confine%l = .false.
    confine%axis = this%axis
    confine%laxis = .false.
    confine%laxis(this%axis) = .true.
    if(allocated(trans)) deallocate(trans)
    allocate(trans(minval(structure%spec(:)%num+2),3))
    call gldfnd(confine, structure, structure, trans, ntrans, this%tol_sym)
    tfmat(:,:) = 0._real32
    tfmat(1,1) = 1._real32
    tfmat(2,2) = 1._real32
    if(ntrans.eq.0)then
       tfmat(3,3) = 1._real32
    else
       itmp1=minloc(abs(trans(:ntrans,this%axis)),dim=1,&
            mask=abs(trans(:ntrans,this%axis)).gt.1.D-3/modu(structure%lat(this%axis,:)))
       tfmat(3,:) = trans(itmp1,:)
    end if
    if(all(abs(tfmat(3,:)).lt.1.E-5_real32)) tfmat(3,3) = 1._real32
    call transformer(structure,tfmat,bas_map)
    if(.not.compare_stoichiometry(structure,structure_compare))then
       write(err_msg,'(A,I0,A)') &
            "The transformed structure stoichiometry does not match the &
            &original structure."
       exit_code_ = 1
       call stop_program( &
            trim(err_msg), &
            exit_code=exit_code_, &
            block_stop = present(exit_code) &
       )
       return
    end if


    ! get the terminations
    term = get_termination_info( &
         structure, this%axis, &
         verbose = merge(1,verbose_,print_termination_info_), &
         tol_sym = this%tol_sym, &
         layer_sep = layer_sep, &
         exit_code = exit_code_ &
    )
    if(exit_code_.ne.0)then
       write(err_msg,'(A,I0,A)') &
            "The termination generator failed with exit code ", exit_code_
       if(break_on_fail_)then
          call stop_program(trim(err_msg))
          return
       end if
    end if
    if(term%nterm .eq. 0)then
       write(warn_msg, '(A,I0,1X,I0,1X,I0,A)') &
            "No terminations found for Miller plane (",miller_,")"
       call print_warning(trim(warn_msg))
       return
    end if

    ! determine tolerance for layer separations (termination tolerance)
    ! ... this is different from layer_sep
    call set_layer_tol(term)

    ! determine required extension and perform that
    call build_slab_supercell(structure, bas_map, term, surface_,&
         height, num_layers_, thickness_, num_cells,&
         term_start, term_end, term_step &
    )
    

    !---------------------------------------------------------------------------
    ! loop over terminations and write them
    !---------------------------------------------------------------------------
    num_structures = ( term_end - term_start ) / term_step + 1
    allocate(output(num_structures))
    do iterm = term_start, term_end, term_step
       i = ( iterm - term_start ) / term_step + 1 
       call output(i)%copy(structure, length=4)
       if(allocated(t1bas_map)) deallocate(t1bas_map)
       allocate(t1bas_map,source=bas_map)
       call cut_slab_to_height(output(i),bas_map,term,[iterm,surface_(2)],&
            thickness_, num_cells, num_layers_, height,&
            prefix, lcycle, orthogonalise_, this%vacuum_gap &
       )
       ! Normalise lattice
       !------------------------------------------------------------------------
       if(normalise_)then
          call reducer(output(i), verbose = verbose_)
          output(i)%lat = MATNORM(output(i)%lat)
       end if
    end do

   end function get_terminations
!###############################################################################


!###############################################################################
   function get_interface_location( &
       this, structure, axis, return_fractional, verbose, exit_code &
   ) result(output)
    !! Get the interface location for the given structure
    implicit none

    ! Arguments
    class(artemis_generator_type), intent(inout) :: this
    !! Instance of artemis generator type
    type(basis_type), intent(in) :: structure
    !! Atomic structure data
    integer, intent(in), optional :: axis
    !! Axis for the interface
    logical, intent(in), optional :: return_fractional
    !! Return the interface location in fractional coordinates
    integer, intent(in), optional :: verbose
    !! Verbosity level
    integer, intent(out), optional :: exit_code
    !! Exit code for the program

    type(intf_info_type) :: output
    !! Output interface location

    ! Local variables
    integer :: axis_
    !! Axis for the interface
    logical :: return_fractional_
    !! Return fractional coordinates
    integer :: exit_code_
    !! Exit code for the program

    axis_ = 0
    exit_code_ = 0
    return_fractional_ = .false.
    if(present(axis)) axis_ = axis
    if(present(return_fractional)) return_fractional_ = return_fractional

    output = get_interface(structure, axis_)

    if(return_fractional_)then
       output%loc = output%loc / modu(structure%lat(output%axis,:))
    end if

    if(present(exit_code)) exit_code = exit_code_

   end function get_interface_location
!###############################################################################


!###############################################################################
  subroutine generate_interfaces_from_existing( &
       this, structure, interface_location, &
       print_shift_info, seed, verbose, exit_code &
  )
    !! Generate swaps and shifts for an existing interface
    implicit none

    ! Arguments
    class(artemis_generator_type), intent(inout) :: this
    !! Instance of artemis generator type
    type(basis_type), intent(in) :: structure
    !! Atomic structure data
    real(real32), dimension(2), intent(in), optional :: interface_location
    !! Interface location
    logical, intent(in), optional :: print_shift_info
    !! Print shift information
    integer, intent(in), optional :: seed
    !! Random seed for generating random numbers
    integer, intent(in), optional :: verbose
    !! Verbosity level
    integer, intent(out), optional :: exit_code
    !! Exit code for the program

    ! Local variables
    integer :: is, ia, js, ja
    !! Loop variables
    real(real32) :: rtmp1,min_bond,min_bond1,min_bond2
    !! Minimum bond length
    type(intf_info_type) :: intf
    !! Interface information
    type(struc_data_type) :: struc_data
    !! Structure data
    real(real32), dimension(3) :: vtmp1
    !! Temporary vector
    logical :: print_shift_info_
    !! Print shift information
    integer :: num_seed
    !! Number of seeds for the random number generator.
    integer, dimension(:), allocatable :: seed_arr
    !! Array of seeds for the random number generator.

    type(bulk_DON_type), dimension(2) :: bulk_DON
    !! Distribution functions for the lower and upper bulk structures

    integer :: verbose_
    !! Verbosity level
    integer :: exit_code_
    !! Exit code for the program


    !---------------------------------------------------------------------------
    ! Initialise variables
    !---------------------------------------------------------------------------
    exit_code_ = 0
    verbose_ = 0
    if(present(verbose)) verbose_ = verbose


    !---------------------------------------------------------------------------
    ! Set the random seed
    !---------------------------------------------------------------------------
    if(present(seed))then
       call random_seed(size=num_seed)
       allocate(seed_arr(num_seed))
       seed_arr = seed
       call random_seed(put=seed_arr)
    else
       call random_seed(size=num_seed)
       allocate(seed_arr(num_seed))
       call random_seed(get=seed_arr)
    end if

    print_shift_info_ = .false.
    if(present(print_shift_info)) print_shift_info_ = print_shift_info

    if(.not.allocated(this%structures)) allocate(this%structures(0))


    min_bond1=huge(0._real32)
    min_bond2=huge(0._real32)
    if(present(interface_location))then
      intf%axis = this%axis
      intf%loc = interface_location
    else
       intf=get_interface(structure,this%axis)
       intf%loc=intf%loc/modu(structure%lat(intf%axis,:))
       if(verbose_.gt.0) write(*,*) "interface axis:",intf%axis
       if(verbose_.gt.0) write(*,*) "interface loc:",intf%loc
    end if
    specloop1: do is=1,structure%nspec
       atomloop1: do ia=1,structure%spec(is)%num

          specloop2: do js=1,structure%nspec
             atomloop2: do ja=1,structure%spec(js)%num
                if(is.eq.js.and.ia.eq.ja) cycle atomloop2
                if( &
                     ( structure%spec(is)%atom(ia,intf%axis).gt.intf%loc(1).and.&
                     structure%spec(is)%atom(ia,intf%axis).lt.intf%loc(2) ).and.&
                     ( structure%spec(js)%atom(ja,intf%axis).gt.intf%loc(1).and.&
                     structure%spec(js)%atom(ja,intf%axis).lt.intf%loc(2) ) )then
                   vtmp1 = (structure%spec(is)%atom(ia,:3)-structure%spec(js)%atom(ja,:3))
                   vtmp1 = matmul(vtmp1,structure%lat)
                   rtmp1 = modu(vtmp1)
                   if(rtmp1.lt.min_bond1) min_bond1 = rtmp1
                elseif( &
                     ( structure%spec(is)%atom(ia,intf%axis).lt.intf%loc(1).or.&
                     structure%spec(is)%atom(ia,intf%axis).gt.intf%loc(2) ).and.&
                     ( structure%spec(js)%atom(ja,intf%axis).lt.intf%loc(1).or.&
                     structure%spec(js)%atom(ja,intf%axis).gt.intf%loc(2) ) )then
                   vtmp1 = (structure%spec(is)%atom(ia,:3)-structure%spec(js)%atom(ja,:3))
                   vtmp1 = matmul(vtmp1,structure%lat)
                   rtmp1 = modu(vtmp1)
                   if(rtmp1.lt.min_bond2) min_bond2 = rtmp1
                end if

             end do atomloop2
          end do specloop2
    
       end do atomloop1
    end do specloop1

    min_bond = ( min_bond1 + min_bond2 ) / 2._real32
    if(verbose_.gt.0) write(*,'(1X,"Avg min bulk bond: ",F0.3," Å")') min_bond
    if(verbose_.gt.0) write(*,'(1X,"Trans-interfacial scaling factor:",F0.3)') this%separation_scale
    this%axis = intf%axis
    call this%generate_perturbations( &
         structure, intf%loc, &
         min_bond, bulk_DON, &
         struc_data, &
         print_shift_info_, seed_arr, verbose_, exit_code_ &
    )

    if(present(exit_code)) exit_code = exit_code_

  end subroutine generate_interfaces_from_existing
!###############################################################################


!###############################################################################
  subroutine generate_interfaces( &
       this, &
       surface_lw, surface_up, &
       thickness_lw, thickness_up, &
       num_layers_lw, num_layers_up, &
       reduce_matches, &
       print_lattice_match_info, print_termination_info, print_shift_info, &
       break_on_fail, &
       icheck_term_pair, interface_idx, &
       generate_structures, &
       seed, verbose, exit_code &
  )
    !! Generate interfaces from two bulk structures
    implicit none

    ! Arguments
    class(artemis_generator_type), intent(inout) :: this
    !! Instance of artemis generator type
    integer, intent(in), dimension(:), optional :: surface_lw
    !! Surface indices for the lower bulk structure
    integer, intent(in), dimension(:), optional :: surface_up
    !! Surface indices for the upper bulk structure
    real(real32), intent(in), optional :: thickness_lw
    !! Thickness of the lower slab
    real(real32), intent(in), optional :: thickness_up
    !! Thickness of the upper slab
    integer, intent(in), optional :: num_layers_lw
    !! Number of layers in the lower slab
    integer, intent(in), optional :: num_layers_up
    !! Number of layers in the upper slab
    logical, intent(in), optional :: reduce_matches
    !! Reduce lattice matches to their smallest cell (UNSTABLE)

    logical, intent(in), optional :: break_on_fail
    !! Break on failure
    logical, intent(in), optional :: print_lattice_match_info
    !! Print lattice match information
    logical, intent(in), optional :: print_termination_info
    !! Print termination information
    logical, intent(in), optional :: print_shift_info
    !! Print shift information
    integer, intent(in), optional :: icheck_term_pair
    !! Index of the lattice match to check
    integer, intent(in), optional :: interface_idx
    !! Index of the interface to output
    logical, intent(in), optional :: generate_structures
    !! Boolean whether to generate structures or just print information
    integer, intent(in), optional :: seed
    !! Random seed for generating random numbers
    integer, intent(in), optional :: verbose
    !! Verbosity level
    integer, intent(out), optional :: exit_code
    !! Exit code for the function

    ! Local variables
    real(real32) :: avg_min_bond
    !! Average minimum bond length

    type(basis_type) :: structure_lw, structure_up, supercell_lw, supercell_up
    !! Copy of the basis structures
    type(basis_type) :: slab_lw, slab_up
    !! Slab structures
    type(basis_type) :: intf_basis
    !! Interface structure
    character(len=256) :: err_msg
    !! Error message

    integer :: j, is, ia
    !! Loop indices
    integer :: unit
    !! Unit number for file I/O
    integer :: ifit, intf_start, intf_end
    !! Interface loop indices
    integer :: iterm_lw, term_lw_start_idx, term_lw_end_idx, term_lw_step
    !! Lower bulk termination loop indices
    integer :: iterm_up, term_up_start_idx, term_up_end_idx, term_up_step
    !! Upper bulk termination loop indices

    ! slab thickness variables
    integer :: num_cells_lw, num_cells_up
    !! Number of cells in the slab
    real(real32) :: height_lw, height_up
    !! Height of the slab
    real(real32) :: thickness_lw_, thickness_up_
    !! Thickness of the slab
    integer :: num_layers_lw_, num_layers_up_
    !! Number of layers in the slab

    integer, dimension(3) :: miller_lw, miller_up
    !! Miller indices for the lower and upper bulk structures
    integer, dimension(2) :: surface_lw_, surface_up_
    !! Surface indices for the lower and upper bulk structures
    logical :: ludef_surface_lw, ludef_surface_up
    !! Boolean whether surfaces are defined 
    logical :: lcycle
    !! Boolean whether to skip the cycle

    logical :: break_on_fail_
    !! Boolean whether to break on failure
    logical :: print_lattice_match_info_, print_termination_info_, &
         print_shift_info_
    !! Boolean whether to print lattice match, termination, and shift info
    integer :: num_seed
    !! Number of seeds for the random number generator.
    integer, dimension(:), allocatable :: seed_arr
    !! Array of seeds for the random number generator.
    integer :: icheck_term_pair_
    !! Index of the lattice match to check
    integer :: interface_idx_
    !! Index of the interface to output
    logical :: generate_structures_
    !! Boolean whether to generate structures or just print information


    type(struc_data_type) :: struc_data
    !! Structure data (i.e. mismatch, terminations, etc)
    character(len=256) :: filename
    !! Filename for error output data
    real(real32) :: rtmp1, bondlength
    !! Temporary variables

    integer :: ntrans, iunique, itmp1, num_structures_old
    integer :: layered_axis_lw, layered_axis_up
    type(confine_type) :: confine
    type(latmatch_type) :: SAV
    type(term_arr_type) :: lw_term, up_term
    integer, dimension(3) :: ivtmp1
    real(real32), dimension(2) :: intf_loc
    real(real32), dimension(3) :: init_offset
    logical :: reduce_matches_
    !! Boolean whether to reduce lattice matches to their smallest cell
    real(real32), dimension(3,3) :: tfmat
    !! Transformation matrix
    type(bulk_DON_type), dimension(2) :: bulk_DON
    !! Distribution functions for the lower and upper bulk structures
    integer, allocatable, dimension(:,:,:) :: lw_map, t1lw_map, t2lw_map
    integer, allocatable, dimension(:,:,:) :: up_map, t1up_map, t2up_map
    real(real32), allocatable, dimension(:,:) :: trans

    integer :: exit_code_
    !! Exit code for the function
    integer :: verbose_
    !! Verbosity level


    !---------------------------------------------------------------------------
    ! Initialise variables
    !---------------------------------------------------------------------------
    exit_code_ = 0
    verbose_ = 0
    reduce_matches_ = .false.
    if(present(verbose)) verbose_ = verbose
    if(present(reduce_matches)) reduce_matches_ = reduce_matches

    icheck_term_pair_ = -1; interface_idx_ = -1
    if(present(icheck_term_pair)) icheck_term_pair_ = icheck_term_pair
    if(present(interface_idx)) interface_idx_ = interface_idx

    break_on_fail_ = .false.
    if(present(break_on_fail)) break_on_fail_ = break_on_fail

    generate_structures_ = .true.
    if(present(generate_structures)) generate_structures_ = generate_structures

    print_lattice_match_info_ = .false.
    print_termination_info_ = .false.
    print_shift_info_ = .false.
    if(present(print_lattice_match_info)) &
         print_lattice_match_info_ = print_lattice_match_info
    if(present(print_termination_info)) &
         print_termination_info_ = print_termination_info
    if(present(print_shift_info)) print_shift_info_ = print_shift_info

    init_offset = [0._real32,0._real32,2._real32]
    if(.not.allocated(this%shifts)) call this%set_shift_method()


    !---------------------------------------------------------------------------
    ! Set the random seed
    !---------------------------------------------------------------------------
    if(present(seed))then
       call random_seed(size=num_seed)
       allocate(seed_arr(num_seed))
       seed_arr = seed
       call random_seed(put=seed_arr)
    else
       call random_seed(size=num_seed)
       allocate(seed_arr(num_seed))
       call random_seed(get=seed_arr)
    end if


    !---------------------------------------------------------------------------
    ! Check if the structures are valid
    !---------------------------------------------------------------------------
    ! check if the structures have anything (i.e. atoms) in them
    if(this%structure_lw%natom.eq.0)then
       write(err_msg,'(A,I0,A)') &
            "The lower structure has ", this%structure_lw%natom, &
            " atoms. It should have at least 1."
       exit_code_ = 1
       call stop_program( &
            trim(err_msg), &
            exit_code=exit_code_, &
            block_stop = present(exit_code) &
       )
       return
    end if
    if(this%structure_up%natom.eq.0)then
       write(err_msg,'(A,I0,A)') &
            "The upper structure has ", this%structure_lw%natom, &
            " atoms. It should have at least 1."
       exit_code_ = 1
       call stop_program( &
            trim(err_msg), &
            exit_code=exit_code_, &
            block_stop = present(exit_code) &
       )
       return
    end if
    call structure_lw%copy(this%structure_lw, length=4)
    call structure_up%copy(this%structure_up, length=4)
    if(.not.allocated(this%structures)) allocate(this%structures(0))
    if(.not.allocated(this%structure_data)) allocate(this%structure_data(0))


    !---------------------------------------------------------------------------
    ! Retrieve the primitive cells if necessary
    !---------------------------------------------------------------------------
    if(this%use_pricel_lw)then
       if(verbose_.gt.0) write(*,'(1X,"Using primitive cell for lower material")')
       call get_primitive_cell(structure_lw, tol_sym=this%tol_sym)
    else
       if(verbose_.gt.0) write(*,'(1X,"Using supplied cell for lower material")')
       call primitive_lat(structure_lw)
    end if
    if(this%use_pricel_up)then
       if(verbose_.gt.0) write(*,'(1X,"Using primitive cell for upper material")')
       call get_primitive_cell(structure_up, tol_sym=this%tol_sym)
    else
       if(verbose_.gt.0) write(*,'(1X,"Using supplied cell for upper material")')
       call primitive_lat(structure_up)
    end if


    !---------------------------------------------------------------------------
    ! Handle surface properties
    !---------------------------------------------------------------------------
    miller_lw = this%miller_lw
    miller_up = this%miller_up
    surface_lw_ = 0
    surface_up_ = 0
    if(present(surface_lw))then
       select case(size(surface_lw, dim=1))
       case(1)
          surface_lw_ = surface_lw(1)
       case(2)
          surface_lw_ = surface_lw
       case default
          write(err_msg,'(A,I0,A)') &
               "The surface vector for the lower material has ", &
               size(surface_lw, dim=1), " components. It should have 1 or 2."
          call stop_program(trim(err_msg))
          return
       end select
    end if
    if(present(surface_up))then
       select case(size(surface_up, dim=1))
       case(1)
          surface_up_ = surface_up(1)
       case(2)
          surface_up_ = surface_up
       case default
          write(err_msg,'(A,I0,A)') &
               "The surface vector for the upper material has ", &
               size(surface_up, dim=1), " components. It should have 1 or 2."
          call stop_program(trim(err_msg))
          return
       end select
    end if

    ludef_surface_lw = .false.
    ludef_surface_up = .false.
    if(all(surface_lw_.gt.0)) ludef_surface_lw = .true.
    if(all(surface_up_.gt.0)) ludef_surface_up = .true.

    thickness_lw_ = -1._real32
    thickness_up_ = -1._real32
    num_layers_lw_ = 0
    num_layers_up_ = 0
    if(present(num_layers_lw)) num_layers_lw_ = num_layers_lw
    if(present(num_layers_up)) num_layers_up_ = num_layers_up
    if(present(thickness_lw)) thickness_lw_ = thickness_lw
    if(present(thickness_up)) thickness_up_ = thickness_up
    if(num_layers_lw_.eq.0.and.abs(thickness_lw_+1._real32).lt.1.E-6_real32)then
       thickness_lw_ = 10._real32
    elseif(num_layers_lw_.le.0.and.thickness_lw_.le.0._real32)then
       write(err_msg,'(A,I0,A)') &
            "The number of layers for the lower material is ", &
            num_layers_lw_, " and the thickness is ", thickness_lw_, &
            " One of these must be greater than 0."
       call stop_program(trim(err_msg))
       return
    end if
    if(num_layers_up_.eq.0.and.abs(thickness_up_+1._real32).lt.1.E-6_real32)then
       thickness_up_ = 10._real32
    elseif(num_layers_up_.le.0.and.thickness_up_.le.0._real32)then
       write(err_msg,'(A,I0,A)') &
            "The number of layers for the upper material is ", &
            num_layers_up_, " and the thickness is ", thickness_up_, &
            " One of these must be greater than 0."
       call stop_program(trim(err_msg))
       return
    end if


    !---------------------------------------------------------------------------
    ! Get the average bond length
    !---------------------------------------------------------------------------
    avg_min_bond = &
         ( &
              get_min_bulk_bond(structure_lw) + &
              get_min_bulk_bond(structure_up) &
         ) / 2._real32
    if(verbose_.gt.0) write(*,'(1X,"Avg min bulk bond: ",F0.3," Å")') avg_min_bond
    if(verbose_.gt.0) write(*,'(1X,"Trans-interfacial scaling factor: ",F0.3)') this%separation_scale
    if(this%shift_method.eq.-1) this%num_shifts = 1
    

    !---------------------------------------------------------------------------
    ! Gets bulk distribution functions (i.e. densities of neighbours)
    ! ... if shift_method = 4
    !---------------------------------------------------------------------------
    allocate(lw_map(structure_lw%nspec,maxval(structure_lw%spec(:)%num,dim=1),2))
    allocate(up_map(structure_up%nspec,maxval(structure_up%spec(:)%num,dim=1),2))    
    if(this%shift_method.eq.4.or.this%shift_method.eq.0)then
       lw_map=0
       bulk_DON(1)%spec=gen_DON(structure_lw%lat,structure_lw,&
            dist_max=this%bondlength_cutoff,&
            scale_dist=.false.,&
            norm=.true. &
       )
       do is = 1, structure_lw%nspec
          if(all(abs(bulk_DON(1)%spec(is)%atom(:,:)).lt.1._real32))then
             bondlength = huge(0._real32)
             do ia = 1, structure_lw%spec(is)%num
                rtmp1 = modu(get_min_bond(structure_lw, is, ia))
                if(rtmp1.lt.bondlength) bondlength = rtmp1
                if(rtmp1.gt.this%bondlength_cutoff)then
                   write(filename,'("lw_DON_",I0,"_",I0,".dat")') is,ia
                   open(newunit=unit, file=filename)
                   do j=1,1000
                      write(unit,*) &
                           (j-1)*this%bondlength_cutoff/1000,&
                           bulk_DON(1)%spec(is)%atom(ia,j)
                   end do
                   close(unit)
                  end if
             end do
             if(bondlength.gt.this%bondlength_cutoff)then
                write(err_msg,'(A,I0,A,F0.3,A,F0.3)') &
                     "Minimum bondlength for species ", &
                     is, " in lower structure is ", bondlength, achar(10) // &
                     "To account for this, increase bondlength cutoff from ", &
                     this%bondlength_cutoff
                call stop_program(trim(err_msg))
             end if
             exit_code_ = 1
             return
          end if
       end do
       up_map=0
       bulk_DON(2)%spec=gen_DON(structure_up%lat,structure_up,&
            dist_max=this%bondlength_cutoff,&
            scale_dist=.false.,&
            norm=.true.)
       do is = 1, structure_up%nspec
          if(all(abs(bulk_DON(2)%spec(is)%atom(:,:)).lt.1._real32))then
             bondlength = huge(0._real32)
             do ia = 1, structure_up%spec(is)%num
                rtmp1 = modu(get_min_bond(structure_up, is, ia))
                if(rtmp1.lt.bondlength) bondlength = rtmp1
                if(rtmp1.gt.this%bondlength_cutoff)then
                   write(filename,'("up_DON_",I0,"_",I0,".dat")') is,ia
                   open(newunit=unit, file=filename)
                   do j=1,1000
                      write(unit,*) &
                           (j-1)*this%bondlength_cutoff/1000,&
                           bulk_DON(2)%spec(is)%atom(ia,j)
                   end do
                   close(unit)
                  end if
             end do
             if(bondlength.gt.this%bondlength_cutoff)then
                write(err_msg,'(A,I0,A,F0.3,A,F0.3)') &
                     "Minimum bondlength for species ", &
                     is, " in upper structure is ", bondlength, achar(10) // &
                     "To account for this, increase bondlength cutoff from ", &
                     this%bondlength_cutoff
                call stop_program(trim(err_msg))
             end if
             exit_code_ = 1
             return
          end if
       end do
    else
       lw_map=-1
       up_map=-1       
    end if


    !---------------------------------------------------------------------------
    ! Check whether system appears layered
    !---------------------------------------------------------------------------
    layered_axis_lw = get_layered_axis( structure_lw%lat, structure_lw )
    if(.not.this%is_layered_lw.and.layered_axis_lw.gt.0)then
       ivtmp1 = 0
       ivtmp1(layered_axis_lw)=1
       if(this%ludef_is_layered_lw)then
          write(err_msg,'("Lower crystal appears layered along axis ",I0,"\n&
               &Partial layer terminations will be generated\n&
               &We suggest using LW_MILLER =",3(1X,I1))') layered_axis_lw,ivtmp1
          call print_warning(trim(err_msg))
       else
          write(err_msg,'("Lower crystal has been identified as layered\nalong",3(1X,I1),"\n&
               &Confining crystal to this plane and\nstoichiometric terminations.\n&
               &If you don''t want this, set\nLW_LAYERED = .FALSE.")') &
               ivtmp1
          call print_warning(trim(err_msg))
          miller_lw=ivtmp1
          this%is_layered_lw=.true.
       end if
    elseif(this%is_layered_lw.and.layered_axis_lw.gt.0.and.all(miller_lw.eq.0))then
       miller_lw(layered_axis_lw)=1
    end if

    layered_axis_up = get_layered_axis( structure_up%lat, structure_up )
    if(.not.this%is_layered_up.and.layered_axis_up.gt.0)then
       ivtmp1=0
       ivtmp1(layered_axis_up)=1
       if(this%ludef_is_layered_up)then
          write(err_msg,'("Upper crystal appears layered along axis ",I0,"\n&
               &Partial layer terminations will be generated\n&
               &We suggest using UP_MILLER =",3(1X,I1))') layered_axis_up,ivtmp1
          call print_warning(trim(err_msg))
       else
          write(err_msg,'("Upper crystal has been identified as layered\nalong",3(1X,I1),"\n&
               &Confining crystal to this plane and\nstoichiometric terminations.\n&
               &If you don''t want this, set\nUP_LAYERED = .FALSE.")') &
               ivtmp1
          call print_warning(trim(err_msg))
          miller_up=ivtmp1
          this%is_layered_up=.true.
       end if
    elseif(this%is_layered_up.and.layered_axis_up.gt.0.and.all(miller_up.eq.0))then
       miller_up(layered_axis_up)=1
    end if


    !---------------------------------------------------------------------------
    ! Finds and stores the best matches between the materials
    !---------------------------------------------------------------------------
    num_structures_old = -1
    if(this%match_method.ne.0.and.(any(miller_lw.ne.0).or.any(miller_up.ne.0)))then
       call stop_program('Cannot use LW_MILLER or UP_MILLER with IMATCH>0')
       exit_code_ = 1
       return
    elseif(this%match_method.ne.0)then
       write(err_msg,'("&
            &IMATCH /= 0 methods are experimental and may\n&
            &not work as expected.\n&
            &They are not intended to be thorough searches.\n&
            &This method is not recommended unless you\n&
            &are clear on its intended use and\n&
            &limitations.&
       &")')
       call print_warning(trim(err_msg))
       tfmat = planecutter(structure_lw%lat,real(miller_lw,real32))
       call transformer(structure_lw,tfmat,lw_map)
    end if
    call SAV%init( &
         this%tolerance, structure_lw%lat, structure_up%lat, &
         this%max_num_matches, reduce_matches_ &
    )
    select case(this%match_method)
    case(0)
       call lattice_matching(&
            SAV, this%tolerance, &
            structure_lw, structure_up, &
            miller_lw = miller_lw, miller_up = miller_up, &
            max_num_planes = this%max_num_planes, &
            verbose = merge(1,verbose_,print_lattice_match_info_), &
            tol_sym = this%tol_sym &
       )
    case default
       call SAV%constrain_axes(miller_lw, miller_up, verbose = verbose_)
       call cyc_lat1(SAV, this%tolerance, this%match_method, verbose = verbose_)
    end select
    if(min(this%max_num_matches,SAV%nfit).eq.0)then
       write(err_msg,'("No matches found between the two structures")')
       call print_warning(trim(err_msg))
       return
    else
       if(verbose_.gt.0) write(*,'(1X,"Number of matches found: ",I0)')&
            min(this%max_num_matches,SAV%nfit)
    end if
    if(verbose_.gt.0) write(*,'(1X,"Maximum number of generated interfaces will be: ",I0)')&
         this%max_num_terms * this%num_shifts * this%max_num_matches
    if(.not.generate_structures_)then
       if(verbose_.gt.0) write(*,'(1X,"Told not to generate structures, just find matches.")')
       return
    end if

       
!!!-----------------------------------------------------------------------------
!!! Saves current directory and moves to new directory
!!!-----------------------------------------------------------------------------
    if(interface_idx_.gt.0)then
       intf_start=interface_idx_
       intf_end=interface_idx_
       if(verbose_.gt.0) write(*,'(1X,"Generating only interfaces for match ",I0)') interface_idx_
    else
       intf_start=1
       intf_end=min(this%max_num_matches,SAV%nfit)
    end if
    iunique=0
!!!-----------------------------------------------------------------------------
!!! Applies the best match transformations
!!!-----------------------------------------------------------------------------
    intf_loop: do ifit = intf_start, intf_end
       if(verbose_.gt.0) write(*,'("Fit number: ",I0)') ifit
       call supercell_lw%copy(structure_lw)
       call supercell_up%copy(structure_up)
       if(allocated(t1lw_map)) deallocate(t1lw_map)
       if(allocated(t1up_map)) deallocate(t1up_map)
       allocate(t1lw_map,source=lw_map)
       allocate(t1up_map,source=up_map)
       

       !!-----------------------------------------------------------------------
       !! Applies the best match transformations
       !!-----------------------------------------------------------------------
       call transformer(supercell_lw,real(SAV%tf1(ifit,:,:),real32),t1lw_map)
       call transformer(supercell_up,real(SAV%tf2(ifit,:,:),real32),t1up_map)


       !!-----------------------------------------------------------------------
       !! Determines the cell change for the upper lattice to get the new DON
       !!-----------------------------------------------------------------------
       if(this%shift_method.eq.4)then
          t1up_map=0 !TEMPORARY TO USE SUPERCELL DONS.
          !DONsupercell_up%lat = matmul(mtmp1,inverse(real(SAV%tf2(ifit,:,:),real32)))
          deallocate(bulk_DON(2)%spec)
          bulk_DON(2)%spec=gen_DON(supercell_up%lat,supercell_up,&
               dist_max=this%bondlength_cutoff,&
               scale_dist=.false.,&
               norm=.true.)
       end if


       !!-----------------------------------------------------------------------
       !! Finds smallest thickness of the lower slab and increases to ...
       !!user-defined thickness
       !! SHOULD MAKE IT LATER MAKE DIFFERENT SETS OF THICKNESSES
       !!-----------------------------------------------------------------------
       confine%l=.false.
       confine%axis=this%axis
       confine%laxis=.false.
       confine%laxis(this%axis)=.true.
       if(allocated(trans)) deallocate(trans)
       allocate(trans(minval(supercell_lw%spec(:)%num+2),3))
       call gldfnd(confine, supercell_lw, supercell_lw, trans, ntrans, this%tol_sym)
       tfmat(:,:)=0._real32
       tfmat(1,1)=1._real32
       tfmat(2,2)=1._real32
       if(ntrans.eq.0)then
          tfmat(3,3)=1._real32
       else
          itmp1=minloc(abs(trans(:ntrans,this%axis)),dim=1,&
               mask=abs(trans(:ntrans,this%axis)).gt.1.D-3/modu(supercell_lw%lat(this%axis,:)))
          tfmat(3,:)=trans(itmp1,:)
       end if
       if(all(abs(tfmat(3,:)).lt.1.E-5_real32)) tfmat(3,3) = 1._real32
       call transformer(supercell_lw,tfmat,t1lw_map)
       ! check the stoichiometry ratios are still maintained
       if(.not.compare_stoichiometry(structure_lw,supercell_lw))then
          write(0,'(1X,"ERROR: Internal error in generate_interfaces")')
          write(0,'(2X,"&
               &The gldfnd subroutine could not reproduce a valid primitive &
               &cell for the lower material on match ",I0)') ifit
          if(verbose_.gt.1)then
             call err_abort_print_struc(supercell_lw, "broken_primitive.vasp", &
              "Code exiting due to IPRINT = 1")
          end if
          write(0,'(2X,"Skipping this lattice match")')
          cycle intf_loop
       end if

       
       !!-----------------------------------------------------------------------
       !! Finds all terminations parallel to the surface plane
       !!-----------------------------------------------------------------------
       if(allocated(lw_term%arr)) deallocate(lw_term%arr)
       lw_term = get_termination_info( &
            supercell_lw, this%axis, &
            verbose = merge(1,verbose_,print_termination_info_), &
            tol_sym = this%tol_sym, &
            layer_sep = this%layer_separation_cutoff(1), &
            exit_code = exit_code_ &
       )
       if(exit_code_.ne.0)then
          write(err_msg,'(A,I0,A)') &
               "The termination generator failed with exit code ", exit_code_
          if(break_on_fail_)then
             call stop_program(trim(err_msg))
             return
          end if
       end if
       if(lw_term%nterm .eq. 0)then
          write(0,'("WARNING: &
               &No terminations found for lower material Miller plane &
               &(",3(1X,I0)," )")' &
          ) SAV%tf1(ifit,3,1:3)
          cycle intf_loop
       end if
       if(any(surface_lw_.gt.lw_term%nterm))then
          write(err_msg, '("surface_lw_ACE VALUES INVALID!\nOne or more value &
               &exceeds the maximum number of terminations in the &
               &structure.\n&
               &  Supplied values: ",I0,1X,I0,"\n&
               &  Maximum allowed: ",I0)') surface_lw_, lw_term%nterm
          call err_abort(trim(err_msg),fmtd=.true.)
          return
       end if


       !!-----------------------------------------------------------------------
       !! Sort out ladder rungs (checks whether the material is centrosymmetric)
       !!-----------------------------------------------------------------------
       !call setup_ladder(supercell_lw%lat,supercell_lw,this%axis,lw_term)
       if(sum(lw_term%arr(:)%natom)*lw_term%nstep.ne.supercell_lw%natom)then
          write(err_msg, '("Number of atoms in lower layers not correct: "&
               &I0,2X,I0)') sum(lw_term%arr(:)%natom)*lw_term%nstep,supercell_lw%natom
          call err_abort(trim(err_msg),fmtd=.true.)
          return
       end if
       call set_layer_tol(lw_term)


       !!-----------------------------------------------------------------------
       !! Defines height of lower slab from user-defined values
       !!-----------------------------------------------------------------------
       call build_slab_supercell(supercell_lw,t1lw_map,lw_term,surface_lw_,&
            height_lw,num_layers_lw_, thickness_lw_,num_cells_lw,&
            term_lw_start_idx,term_lw_end_idx,term_lw_step &
       )
       if(term_lw_end_idx.gt.this%max_num_terms) term_lw_end_idx = this%max_num_terms


       !!-----------------------------------------------------------------------
       !! Finds smallest thickness of the upper slab and increases to ...
       !! ... user-defined thickness
       !! SHOULD MAKE IT LATER MAKE DIFFERENT SETS OF THICKNESSES
       !!-----------------------------------------------------------------------
       deallocate(trans)
       allocate(trans(minval(supercell_up%spec(:)%num+2),3))
       call gldfnd(confine, supercell_up, supercell_up, trans,ntrans, this%tol_sym)
       tfmat(:,:)=0._real32
       tfmat(1,1)=1._real32
       tfmat(2,2)=1._real32
       if(ntrans.eq.0)then
          tfmat(3,3)=1._real32
       else
          itmp1=minloc(abs(trans(:ntrans,this%axis)),dim=1,&
               mask=abs(trans(:ntrans,this%axis)).gt.1.D-3/modu(supercell_lw%lat(this%axis,:)))
          tfmat(3,:)=trans(itmp1,:)
       end if
       if(all(abs(tfmat(3,:)).lt.1.E-5_real32)) tfmat(3,3) = 1._real32
       call transformer(supercell_up,tfmat,t1up_map)
       ! check the stoichiometry ratios are still maintained
       if(.not.compare_stoichiometry(structure_up,supercell_up))then
          write(0,'(1X,"ERROR: Internal error in generate_interfaces")')
          write(0,'(2X,"&
               &The gldfnd subroutine could not reproduce a valid primitive &
               &cell for the upper material on match ",I0)') ifit
          if(verbose_.gt.1)then
             call err_abort_print_struc(supercell_up, "broken_primitive.vasp", &
              "Code exiting due to IPRINT = 1")
          end if
          write(0,'(2X,"Skipping this lattice match")')
          cycle intf_loop
       end if

       
       !!-----------------------------------------------------------------------
       !! Finds all supercell_up%lat unique terminations parallel to the surface plane
       !!-----------------------------------------------------------------------
       if(allocated(up_term%arr)) deallocate(up_term%arr)
       up_term = get_termination_info( &
            supercell_up, this%axis, &
            verbose = merge(1,verbose_,print_termination_info_), &
            tol_sym = this%tol_sym, &
            layer_sep = this%layer_separation_cutoff(2), &
            exit_code = exit_code_ &
       )
       if(exit_code_.ne.0)then
          write(err_msg,'(A,I0,A)') &
               "The termination generator failed with exit code ", exit_code_
          if(break_on_fail_)then
             call stop_program(trim(err_msg))
             return
          end if
       end if
       if(up_term%nterm .eq. 0)then
          write(0,'("WARNING: &
               &No terminations found for upper material Miller plane &
               &(",3(1X,I0)," )")' &
          ) SAV%tf2(ifit,3,1:3)
          cycle intf_loop
       end if
       if(any(surface_up_.gt.up_term%nterm))then
          write(err_msg, '("surface_up_ACE VALUES INVALID!\nOne or more value &
               &exceeds the maximum number of terminations in the &
               &structure.\n&
               &  Supplied values: ",I0,1X,I0,"\n&
               &  Maximum allowed: ",I0)') surface_up_, up_term%nterm
          call err_abort(trim(err_msg),fmtd=.true.)
          return
       end if


       !!-----------------------------------------------------------------------
       !! Sort out ladder rungs (checks whether the material is centrosymmetric)
       !!-----------------------------------------------------------------------
       !call setup_ladder(supercell_up%lat,supercell_up,this%axis,up_term)
       if(sum(up_term%arr(:)%natom)*up_term%nstep.ne.supercell_up%natom)then
          write(err_msg, '("Number of atoms in upper layers not correct: "&
               &I0,2X,I0)') sum(up_term%arr(:)%natom)*up_term%nstep,supercell_up%natom
          call stop_program(trim(err_msg))
          return
       end if
       call set_layer_tol(up_term)


       !!-----------------------------------------------------------------------
       !! Defines height of upper slab from user-defined values
       !!-----------------------------------------------------------------------
       call build_slab_supercell(supercell_up,t1up_map,up_term,surface_up_,&
            height_up,num_layers_up_, thickness_up_, num_cells_up,&
            term_up_start_idx,term_up_end_idx,term_up_step &
       )
       if(term_up_end_idx.gt.this%max_num_terms) term_up_end_idx = this%max_num_terms


       !!-----------------------------------------------------------------------
       !! Print termination plane locations
       !!-----------------------------------------------------------------------
       if(verbose_.gt.0) write(*,'(1X,"Number of unique terminations: ",I0,2X,I0)') &
            lw_term%nterm,up_term%nterm

       !!-----------------------------------------------------------------------
       !! Cycle over terminations of both materials and generates interfaces ...
       !! ... composed of all of the possible combinations of the two
       !!-----------------------------------------------------------------------
       lw_term_loop: do iterm_lw = term_lw_start_idx, term_lw_end_idx, term_lw_step
          call slab_lw%copy(supercell_lw)
          if(allocated(t2lw_map)) deallocate(t2lw_map)
          allocate(t2lw_map,source=t1lw_map)
          !!--------------------------------------------------------------------
          !! Shifts lower material to specified termination
          !!--------------------------------------------------------------------
          call cut_slab_to_height(slab_lw,t2lw_map,lw_term,[iterm_lw,surface_lw_(2)],&
               thickness_lw_, num_cells_lw, num_layers_lw_, height_lw,&
               "lw",lcycle, &
               vacuum = this%vacuum_gap &
          )
          if(lcycle) cycle lw_term_loop

          
          !!--------------------------------------------------------------------
          !! Cycles over terminations of upper material
          !!--------------------------------------------------------------------
          up_term_loop: do iterm_up = term_up_start_idx, term_up_end_idx, term_up_step
             call slab_up%copy(supercell_up)
             if(allocated(t2up_map)) deallocate(t2up_map)
             allocate(t2up_map,source=t1up_map)
             call cut_slab_to_height(slab_up,t2up_map,up_term,[iterm_up,surface_up_(2)],&
                  thickness_up_, num_cells_up, num_layers_up_, height_up,&
                  "up",lcycle, &
                  vacuum = this%vacuum_gap &
             )
             if(lcycle) cycle up_term_loop

             
             !!-----------------------------------------------------------------
             !! Checks stoichiometry
             !!-----------------------------------------------------------------
             if(slab_lw%nspec.ne.structure_lw%nspec.or.any(&
                  (structure_lw%spec(1)%num*slab_lw%spec(:)%num)&
                  /slab_lw%spec(1)%num.ne.structure_lw%spec(:)%num))then
                write(*,'("WARNING: This lower surface termination is not &
                     &stoichiometric")')
                if(this%is_layered_lw)then
                   write(*,'(2X,"As lower structure is layered, stoichiometric &
                        &surfaces are required.")')
                   write(*,'(2X,"Skipping this termination...")')
                   cycle lw_term_loop
                elseif(this%require_stoichiometry_lw)then
                   write(*,'(2X,"Skipping this termination...")')
                   cycle lw_term_loop
                end if
             end if
             if(slab_up%nspec.ne.structure_up%nspec.or.any(&
                  (structure_up%spec(1)%num*slab_up%spec(:)%num)&
                  /slab_up%spec(1)%num.ne.structure_up%spec(:)%num))then
                write(*,'("WARNING: This upper surface termination is not &
                     &stoichiometric")')
                if(this%is_layered_up)then
                   write(*,'(2X,"As upper structure is layered, stoichiometric &
                        &surfaces are required.")')
                   write(*,'(2X,"Skipping this termination...")')
                   cycle up_term_loop
                elseif(this%require_stoichiometry_up)then
                   write(*,'(2X,"Skipping this termination...")')
                   cycle up_term_loop
                end if
             end if


             !------------------------------------------------------------------
             ! Use the bulk moduli to determine the strain sharing
             !------------------------------------------------------------------
             if(allocated(this%elastic_constants_lw).and. &
                allocated(this%elastic_constants_up))then
                select case(size(this%elastic_constants_lw))
                case(1)
                   if( abs(this%elastic_constants_lw(1)).gt.0.E0 .and. &
                         abs(this%elastic_constants_up(1)).gt.0.E0 &
                   )then
                      call share_strain(slab_lw%lat,slab_up%lat,&
                            this%elastic_constants_lw(1), &
                            this%elastic_constants_up(1), &
                            lcompensate = this%compensate_normal &
                      )
                   end if
                case default
                   write(err_msg,'("Elastic constants not yet set up to handle &
                        &the full tensor.")')
                   call stop_program(trim(err_msg))
                   exit_code_ = 1
                   return
                end select
             elseif(allocated(this%elastic_constants_lw).neqv. &
                   allocated(this%elastic_constants_up))then
                write(err_msg,'(A)') &
                     "Elastic constants not set up for both materials."
                call stop_program(trim(err_msg))
                exit_code_ = 1
                return
             end if

               

             !------------------------------------------------------------------
             ! Merge the two bases and lattices and define the interface loc
             !------------------------------------------------------------------
             intf_basis = basis_stack(&
                  basis1 = slab_lw, basis2 = slab_up, &
                  axis = this%axis, offset = init_offset(:), &
                  map1 = t2lw_map, map2 = t2up_map &
             )
             intf_loc(1) = ( modu(slab_lw%lat(this%axis,:)) + 0.5_real32*init_offset(this%axis) - &
                  this%vacuum_gap)/modu(intf_basis%lat(this%axis,:))
             intf_loc(2) = ( modu(slab_lw%lat(this%axis,:)) + modu(slab_up%lat(this%axis,:)) + &
                  1.5_real32*init_offset(this%axis) - 2._real32*this%vacuum_gap )/modu(intf_basis%lat(this%axis,:))
             if(verbose_.ge.1)then
                write(0,*) "interface:",intf_loc
                if(verbose_.eq.1.and.iunique.eq.icheck_term_pair_-1)then
                  !  call chdir(intf_dir)
                   call err_abort_print_struc(slab_lw,"lw_term.vasp",&
                        "",.false.)
                   call err_abort_print_struc(slab_up,"up_term.vasp",&
                        "As IPRINT = 1 and ICHECK has been set, &
                        &code is now exiting...")
                elseif(verbose_.eq.2.and.iunique.eq.icheck_term_pair_-1)then
                  !  call chdir(intf_dir)
                   call err_abort_print_struc(intf_basis,"test_intf.vasp",&
                        "As IPRINT = 2 and ICHECK has been set, &
                        &code is now exiting...")
                end if
             end if



             !------------------------------------------------------------------
             ! Saves current directory and moves to new directory
             !------------------------------------------------------------------
             if(this%num_structures.gt.num_structures_old) iunique = iunique + 1
             num_structures_old = this%num_structures

             
             !------------------------------------------------------------------
             ! Write information of current match to file in save directory
             !------------------------------------------------------------------
            !  call output_intf_data(SAV, ifit, lw_term, iterm_lw, up_term, iterm_up,&
            !       this%use_pricel_lw, this%use_pricel_up)
             struc_data = struc_data_type( &
                  match_idx = ifit, &
                  match_and_term_idx = iunique, &
                  from_pricel_lw = this%use_pricel_lw, &
                  from_pricel_up = this%use_pricel_up, &
                  term_lw_idx = [iterm_lw,max(surface_lw_(2),iterm_lw)], &
                  term_up_idx = [iterm_up,max(surface_up_(2),iterm_up)], &
                  term_lw_bounds = [ lw_term%arr(iterm_lw)%hmin, &
                                     lw_term%arr(iterm_lw)%hmax, &
                                     lw_term%arr(max(surface_lw_(2),iterm_lw))%hmin, &
                                     lw_term%arr(max(surface_lw_(2),iterm_lw))%hmax &
                  ], &
                  term_up_bounds = [ up_term%arr(iterm_up)%hmin, &
                                     up_term%arr(iterm_up)%hmax, &
                                     up_term%arr(max(surface_up_(2),iterm_up))%hmin, &
                                     up_term%arr(max(surface_up_(2),iterm_up))%hmax &
                  ], &
                  term_lw_natom = [ lw_term%arr(iterm_lw)%natom, &
                       lw_term%arr(max(surface_lw_(2),iterm_lw))%natom &
                  ], &
                  term_up_natom = [ up_term%arr(iterm_up)%natom, &
                       up_term%arr(max(surface_up_(2),iterm_up))%natom &
                  ], &
                  approx_thickness_lw = max(thickness_lw_,height_lw), &
                  approx_thickness_up = max(thickness_up_,height_up), &
                  transform_lw = SAV%tf1(ifit,:,:), &
                  transform_up = SAV%tf2(ifit,:,:), &
                  mismatch = SAV%tol(ifit,:3) &
             )


             !------------------------------------------------------------------
             ! Generate shifts and swaps and prints the subsequent structures
             !------------------------------------------------------------------
             call this%generate_perturbations( &
                  intf_basis, intf_loc, avg_min_bond, &
                  bulk_DON, &
                  struc_data, &
                  print_shift_info_, &
                  seed_arr, &
                  verbose_, &
                  exit_code_, &
                  t2lw_map &
             )

             if(this%num_structures.ge.this%max_num_structures) exit intf_loop

             if(ludef_surface_up) exit up_term_loop
          end do up_term_loop
          if(ludef_surface_lw) exit lw_term_loop
       end do lw_term_loop

    end do intf_loop

    if(present(exit_code)) exit_code = exit_code_

  end subroutine generate_interfaces
!###############################################################################


!!!#############################################################################
!!! Takes input interface structure and generates a set of shifts and swaps.
!!!#############################################################################
!!! ISWAP METHOD NOT YET SET UP
  subroutine generate_shifts_and_swaps( &
       this, basis, intf_loc, bond, bulk_DON, struc_data, print_shift_info, &
       seed_arr, verbose, exit_code, map &
  )
    implicit none
    class(artemis_generator_type), intent(inout) :: this
    type(basis_type), intent(in) :: basis
    real(real32), dimension(2), intent(in) :: intf_loc
    real(real32), intent(in) :: bond
    type(bulk_DON_type), dimension(2), intent(in) :: bulk_DON
    !! Distribution functions for the lower and upper bulk structures
    type(struc_data_type), intent(in) :: struc_data
    logical, intent(in) :: print_shift_info
    integer, dimension(:), intent(in) :: seed_arr
    integer, intent(in) :: verbose
    integer, intent(inout) :: exit_code
    integer, dimension(:,:,:), optional, intent(in) :: map

    integer :: iaxis,k,l
    integer :: ngen_swaps,nswaps_per_cell
    real(real32) :: rtmp1
    type(basis_type) :: tbas
    type(bond_type) :: min_bond
    type(struc_data_type) :: struc_data_shift
    type(struc_data_type), dimension(:), allocatable :: struc_data_swaps
    character(len=256) :: err_msg
    integer, dimension(3) :: abc
    real(real32), dimension(3) :: toffset
    type(basis_type), allocatable, dimension(:) :: basis_arr
    real(real32), allocatable, dimension(:,:) :: output_shifts



!!!-----------------------------------------------------------------------------
!!! Sets up shift axis
!!!-----------------------------------------------------------------------------
    abc = [ 1, 2, 3 ]
    abc = cshift(abc,this%axis)


!!!-----------------------------------------------------------------------------
!!! Generates sets of shifts based on shift version
!!!-----------------------------------------------------------------------------
    if(this%shift_method.eq.0.or.this%shift_method.eq.1) allocate(output_shifts(this%num_shifts,3))
    select case(this%shift_method)
    case(1)
       output_shifts(1,:3)=0._real32
       do k=2,this%num_shifts
          do iaxis = 1, 2
             call random_number(output_shifts(k,iaxis))
          end do
       end do
    case(2)
       output_shifts = get_fit_shifts(&
            lat=basis%lat,bas=basis,&
            bond=bond,&
            axis=this%axis,&
            intf_loc=intf_loc,&
            depth=this%interface_depth,&
            nstore=this%num_shifts)
    case(3)
       output_shifts = get_descriptive_shifts(&
            lat=basis%lat,bas=basis,&
            bond=bond,&
            axis=this%axis,&
            intf_loc=intf_loc,&
            depth=this%interface_depth, &
            c_scale=this%separation_scale,&
            nstore=this%num_shifts,lprint=print_shift_info)
    case(4)
       if(present(map))then
          output_shifts = get_shifts_DON(&
               bas=basis,&
               axis=this%axis,&
               intf_loc=intf_loc,&
               nstore=this%num_shifts, &
               c_scale=this%separation_scale, &
               offset=this%shifts(1,:3),&
               verbose=merge(1,verbose,print_shift_info), &
               bulk_DON=bulk_DON,bulk_map=map,&
               max_bondlength=this%bondlength_cutoff,&
               tol_sym=this%tol_sym)
       else
          output_shifts = get_shifts_DON(&
               bas=basis,&
               axis=this%axis,&
               intf_loc=intf_loc,&
               nstore=this%num_shifts, &
               c_scale=this%separation_scale, &
               offset=this%shifts(1,:3),&
               verbose=merge(1,verbose,print_shift_info), &
               max_bondlength=this%bondlength_cutoff,&
               tol_sym=this%tol_sym)
       end if
       if(size(output_shifts(:,1)).eq.0)then
          write(0,'(2X,"No shifts were identified with ISHIFT = 4 for this lattice match")')
          write(0,'(2X,"We suggest increasing MBOND_MAXLEN to find shifts")')
          write(0,'("Skipping interface...")')
          return
       end if
    case default
       if(.not.allocated(output_shifts)) allocate(output_shifts(1,3))
       output_shifts(:,:) = this%shifts
       do iaxis = 1, 2
          output_shifts(1,iaxis) = output_shifts(1,iaxis)!/modu(lat(iaxis,:))
       end do
    end select
    if(this%shift_method.gt.0)then
       output_shifts(:,this%axis) = output_shifts(:,this%axis) * modu(basis%lat(this%axis,:))
    end if


!!!-----------------------------------------------------------------------------
!!! Prints number of shifts to terminal
!!!-----------------------------------------------------------------------------
    if(verbose.gt.0) write(*,'(3X,"Number of unique shifts structures: ",I0)') size(output_shifts,1)


!!!-----------------------------------------------------------------------------
!!! Determines number of swaps across the interface
!!!-----------------------------------------------------------------------------
    nswaps_per_cell = nint(this%swap_density*get_area([basis%lat(abc(1),:)],[basis%lat(abc(2),:)]))
    if(this%swap_method.ne.0)then
       if(verbose.gt.0) write(*,&
            '(" Generating ",I0," swaps per structure ")') nswaps_per_cell
    end if


!!!-----------------------------------------------------------------------------
!!! Prints each unique shift structure
!!!-----------------------------------------------------------------------------
    shift_loop: do k = 1, size(output_shifts,1), 1
       call tbas%copy(basis)
       toffset=output_shifts(k,:3)
       do iaxis=1,2
          call shift_region(tbas,this%axis,&
               intf_loc(1),intf_loc(2),&
               shift_axis=iaxis,shift=toffset(iaxis),renorm=.true.)
       end do
       rtmp1=modu(tbas%lat(this%axis,:))
       call set_vacuum(&
            basis=tbas,&
            axis=this%axis,loc=maxval(intf_loc(:)),&
            vac=toffset(this%axis))
       rtmp1=minval(intf_loc(:))*rtmp1/modu(tbas%lat(this%axis,:))
       call set_vacuum(&
            basis=tbas,&
            axis=this%axis,loc=rtmp1,&
            vac=toffset(this%axis))
       min_bond = get_shortest_bond(tbas)
       if(min_bond%length.le.1.5_real32)then
          write(err_msg,'("Smallest bond in the interface structure is\nless than 1.5 Å.")')
          call print_warning(trim(err_msg))
          write(*,'(2X,"bond length: ",F9.6)') min_bond%length
          write(*,'(2X,"atom 1:",I4,2X,I4)') min_bond%atoms(1,:)
          write(*,'(2X,"atom 2:",I4,2X,I4)') min_bond%atoms(2,:)
       end if


       !!-----------------------------------------------------------------------
       !! Merges lower and upper materials
       !! Writes interfaces to output directories
       !!-----------------------------------------------------------------------
       struc_data_shift = struc_data
       struc_data_shift%shift_idx = k
       struc_data_shift%shift = toffset
       this%structures = [ this%structures, tbas ]
       this%num_structures = size(this%structures, dim = 1)
       this%structure_data = [ this%structure_data, struc_data_shift ]
       if(this%num_structures.ge.this%max_num_structures) return


       !!-----------------------------------------------------------------------
       !! Performs swaps within the shifted structures if requested
       !!-----------------------------------------------------------------------
       if_swap: if(this%swap_method.ne.0)then
          basis_arr = rand_swapper(tbas%lat,tbas,this%axis,this%swap_depth,&
               nswaps_per_cell,this%num_swaps,intf_loc,this%swap_method,&
               seed_arr, &
               tol_sym = this%tol_sym, &
               verbose = verbose, &
               sigma=this%swap_sigma, &
               require_mirror=this%require_mirror_swaps &
          )
          ngen_swaps = this%num_swaps
          LOOPswaps: do l=1,this%num_swaps
             if (basis_arr(l)%nspec.eq.0) then
                ngen_swaps = l - 1
                exit LOOPswaps
             end if
          end do LOOPswaps
          if(ngen_swaps.eq.0)then
             exit if_swap
          end if
          if(allocated(struc_data_swaps)) deallocate(struc_data_swaps)
          allocate(struc_data_swaps(ngen_swaps))
          do l=1,ngen_swaps
             struc_data_swaps(l) = struc_data_shift
             struc_data_swaps(l)%swap_idx = l
             struc_data_swaps(l)%swap_density = this%swap_density
             ! struc_data_swaps(l)%approx_eff_swap_conc = 
          end do
          this%structures = [ this%structures, basis_arr(1:ngen_swaps) ]
          this%structure_data = [ this%structure_data, struc_data_swaps ]
          deallocate(basis_arr)
       end if if_swap


    end do shift_loop

  end subroutine generate_shifts_and_swaps
!!!#############################################################################


!###############################################################################
  subroutine write_match_and_term_data(this, idx, directory, filename)
    !! This subroutine writes the match and termination data to a file
    implicit none

    ! Arguments
    class(artemis_generator_type), intent(in) :: this
    !! Instance of artemis generator type
    integer, intent(in) :: idx
    !! List of indices for the structures to be written
    character(len=*), intent(in) :: directory
    !! Directory where the files will be written
    character(len=*), intent(in) :: filename
    !! Name of the file to be written

    ! Local variables
    integer :: unit


    open(newunit=unit, file=trim(adjustl(directory))//"/"//trim(adjustl(filename)))
    associate( struc_data => this%structure_data(idx) )
       write(unit,'("Lower material primitive cell used: ",L1)') struc_data%from_pricel_lw
       write(unit,'("Upper material primitive cell used: ",L1)') struc_data%from_pricel_up
       write(unit,*)
       write(unit,'("Match and termination identifier: ",I0)') struc_data%match_and_term_idx
       write(unit,'("Lattice match: ",I0)') struc_data%match_idx
       write(unit,'((1X,3(3X,A1),3X,3(3X,A1)),3(/,2X,3(I3," "),3X,3(I3," ")))') &
            "a", "b", "c", "a", "b", "c", &
            struc_data%transform_lw(1,1:3), struc_data%transform_up(1,1:3), &
            struc_data%transform_lw(2,1:3), struc_data%transform_up(2,1:3), &
            struc_data%transform_lw(3,1:3), struc_data%transform_up(3,1:3)
       write(unit,'(" vector mismatch (%) = ",F0.9)') struc_data%mismatch(1)
       write(unit,'(" angle mismatch (°)  = ",F0.9)') struc_data%mismatch(2) * 180._real32 / pi
       write(unit,'(" area mismatch (%)   = ",F0.9)') struc_data%mismatch(3)
       write(unit,*)
       write(unit,'(" Lower crystal Miller plane: ",3(I3," "))') struc_data%transform_lw(3,1:3)
       write(unit,'(" Lower termination")')
       write(unit,'(1X,"Term.",3X,"Min layer loc",3X,"Max layer loc",3X,"no. atoms")')
       write(unit,'(1X,I3,8X,F7.5,9X,F7.5,8X,I3)') &
               struc_data%term_lw_idx(1), &
               struc_data%term_lw_bounds(1:2), &
               struc_data%term_lw_natom(1)
       write(unit,'(1X,I3,8X,F7.5,9X,F7.5,8X,I3)') &
               struc_data%term_lw_idx(2), &
               struc_data%term_lw_bounds(3:4), &
               struc_data%term_lw_natom(2)
       write(unit,*)
       write(unit,'(" Upper crystal Miller plane: ",3(I3," "))') struc_data%transform_up(3,1:3)
       write(unit,'(" Upper termination")')
       write(unit,'(1X,"Term.",3X,"Min layer loc",3X,"Max layer loc",3X,"no. atoms")')
       write(unit,'(1X,I3,8X,F7.5,9X,F7.5,8X,I3)') &
               struc_data%term_up_idx(1), &
               struc_data%term_up_bounds(1:2), &
               struc_data%term_up_natom(1)
       write(unit,'(1X,I3,8X,F7.5,9X,F7.5,8X,I3)') &
               struc_data%term_up_idx(2), &
               struc_data%term_up_bounds(3:4), &
               struc_data%term_up_natom(2)
       write(unit,*)
    end associate

    close(unit)

  end subroutine write_match_and_term_data
!###############################################################################


!###############################################################################
  subroutine write_shift_data(this, idx_list, directory, filename)
    !! This subroutine writes the shift data to a file
    implicit none

    ! Arguments
    class(artemis_generator_type), intent(in) :: this
    !! Instance of artemis generator type
    integer, dimension(:), intent(in) :: idx_list
    !! List of indices for the structures to be written
    character(len=*), intent(in) :: directory
    !! Directory where the files will be written
    character(len=*), intent(in) :: filename
    !! Name of the file to be written

    ! Local variables
    integer :: i
    integer :: unit


    open(newunit=unit, file=trim(adjustl(directory))//"/"//trim(adjustl(filename)))
    write(unit, &
         '("# shift_num    shift (a,b,c) units=(direct,direct,Å)")')
    do i = 1, size(idx_list), 1
       write(unit,'(2X,I0.2,15X,"(",2(" ",F9.6,", ")," ",F9.6," )")') &
            i, this%structure_data(idx_list(i))%shift
    end do
    close(unit)

  end subroutine write_shift_data
!###############################################################################

end module artemis__generator
