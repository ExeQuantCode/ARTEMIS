!!!#############################################################################
!!! INTERFACES CARD SUBROUTINES
!!! Code written by Ned Thaddeus Taylor and Isiah Edward Mikel Rudkin
!!! Code part of the ARTEMIS group (Hepplestone research group).
!!! Think Hepplestone, think HRG.
!!!#############################################################################
module artemis__interface_generator
  use artemis__constants,     only: real32, ierror, pi
  use artemis__misc,          only: to_lower,to_upper
  use artemis__misc_types,    only: abstract_artemis_generator_type, latmatch_type, tol_type
  use artemis__geom_rw,       only: basis_type,geom_write
  use lat_compare,            only: get_best_match
  use artemis__io_utils,      only: err_abort, print_warning
  use artemis__io_utils_extd, only: err_abort_print_struc
  use misc_linalg,            only: uvec,modu,get_area,inverse,cross
  use interface_identifier,   only: intf_info_type,&
       get_interface,get_layered_axis,gen_DON
  use edit_geom,              only: planecutter,primitive_lat,ortho_axis,&
       shift_region,set_vacuum,transformer,shifter,reducer,&
       get_min_bulk_bond,get_shortest_bond,bond_type,&
       share_strain, MATNORM, basis_stack, compare_stoichiometry
  use artemis__sym,           only: confine_type,gldfnd,&
       get_primitive_cell
  use artemis__terminations,  only: get_termination_info, term_arr_type, set_slab_height, set_layer_tol, build_slab
  use swapping,               only: rand_swapper
  use shifting !!! CHANGE TO SHIFTER?
  implicit none


  private

  public :: artemis_interface_generator_type


  type, extends(abstract_artemis_generator_type) :: artemis_interface_generator_type
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
    real(real32), dimension(:,:), allocatable :: shift_data
    !! Data of shifts for each interface, where index 1 is the interface number in structures

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
    integer :: max_num_matches = 5
    integer :: max_num_terms = 5
    integer :: max_num_planes = 10

    logical :: fix_normal = .true. !! compensate_strains_parallel = .true.
    !! Fix the lattice constants parallel to the interface normal vector
    !! Fix = true = strained
    !! Fix = false = relaxed (compensate for interfacial strain by extending/compressing)
    
    real(real32) :: bondlength_cutoff = 6._real32
    real(real32), dimension(2) :: layer_separation_cutoff = 1._real32

    type(tol_type) :: tolerance

   !  type(basis_type), dimension(:), allocatable :: term_structures_lw
   !  type(basis_type), dimension(:), allocatable :: term_structures_up
   contains
    procedure, pass(this) :: set_tolerance
    procedure, pass(this) :: set_shift_method
    procedure, pass(this) :: generate => generate_interfaces
    procedure, pass(this) :: restart => generate_intefaces_from_existing
    procedure, pass(this) :: generate_perturbations => generate_shifts_and_swaps
  end type artemis_interface_generator_type

contains

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
    class(artemis_interface_generator_type), intent(inout) :: this
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
       ! if(present(nstore)) this%tolerance%nstore = nstore
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
       interface_depth, separation_scale, depth_method &
  )
    !! Set the shift method
    implicit none

    ! Arguments
    class(artemis_interface_generator_type), intent(inout) :: this
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

    ! Local variables
    character(len=256) :: err_msg

    if(present(method)) this%shift_method = method
    if(present(num_shifts)) this%num_shifts = num_shifts
    if(present(interface_depth)) this%interface_depth = interface_depth
    if(present(separation_scale)) this%separation_scale = separation_scale
    if(present(depth_method)) this%depth_method = depth_method
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
                     "ERROR: The shifts vector has ", size(shifts, dim=1), &
                     " components. It should have 1 or 3."
                call err_abort(trim(err_msg),fmtd=.true.)
             end if
          end select
       rank(2)
          if(size(shifts,dim=2).eq.3) then
             allocate(this%shifts(size(shifts,1),3))
             this%shifts = shifts
          else
             write(err_msg,'(A,I0,A)') &
                  "ERROR: The shifts vector has ", size(shifts, dim=1), &
                  " components. It should have 3."
             call err_abort(trim(err_msg),fmtd=.true.)
          end if
       rank default
          write(err_msg,'(A,I0,A)') &
               "ERROR: The shifts vector has ", size(shifts, dim=1), &
               " components. It should have 1, 2, or 3."
          call err_abort(trim(err_msg),fmtd=.true.)
       end select
    else
       if(allocated(this%shifts)) deallocate(this%shifts)
       allocate(this%shifts(1,3), source = -1._real32)
    end if

  end subroutine set_shift_method
!###############################################################################


!###############################################################################
  subroutine generate_intefaces_from_existing(this, basis, interface_location, &
       print_shift_info, seed &
  )
    !! Generate interfaces for the given basis
    implicit none

    ! Arguments
    class(artemis_interface_generator_type), intent(inout) :: this
    !! Instance of artemis generator type
    type(basis_type), intent(in) :: basis
    !! Atomic structure data
    real(real32), dimension(2), intent(in), optional :: interface_location
    !! Interface location
    logical, intent(in), optional :: print_shift_info
    !! Print shift information
    integer, intent(in), optional :: seed
    !! Random seed for generating random numbers

    ! Local variables
    integer :: is,ia,js,ja
    !! Loop variables
    real(real32) :: dtmp1,min_bond,min_bond1,min_bond2
    !! Minimum bond length
    type(intf_info_type) :: intf
    !! Interface information
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
       intf=get_interface(basis%lat,basis,this%axis)
       intf%loc=intf%loc/modu(basis%lat(intf%axis,:))
       write(6,*) "interface axis:",intf%axis
       write(6,*) "interface loc:",intf%loc
       !! write interface location to a file for user to refer back to
       open(unit=10,file="interface_location.dat")
       write(10,'(1X,"AXIS = ",I0)') intf%axis
       write(10,'(1X,"INTF_LOC = ",2(2X,F9.6))') intf%loc
       close(10)
    end if
    specloop1: do is=1,basis%nspec
       atomloop1: do ia=1,basis%spec(is)%num

          specloop2: do js=1,basis%nspec
             atomloop2: do ja=1,basis%spec(js)%num
                if(is.eq.js.and.ia.eq.ja) cycle atomloop2
                if( &
                     ( basis%spec(is)%atom(ia,intf%axis).gt.intf%loc(1).and.&
                     basis%spec(is)%atom(ia,intf%axis).lt.intf%loc(2) ).and.&
                     ( basis%spec(js)%atom(ja,intf%axis).gt.intf%loc(1).and.&
                     basis%spec(js)%atom(ja,intf%axis).lt.intf%loc(2) ) )then
                   vtmp1 = (basis%spec(is)%atom(ia,:3)-basis%spec(js)%atom(ja,:3))
                   vtmp1 = matmul(vtmp1,basis%lat)
                   dtmp1 = modu(vtmp1)
                   if(dtmp1.lt.min_bond1) min_bond1 = dtmp1
                elseif( &
                     ( basis%spec(is)%atom(ia,intf%axis).lt.intf%loc(1).or.&
                     basis%spec(is)%atom(ia,intf%axis).gt.intf%loc(2) ).and.&
                     ( basis%spec(js)%atom(ja,intf%axis).lt.intf%loc(1).or.&
                     basis%spec(js)%atom(ja,intf%axis).gt.intf%loc(2) ) )then
                   vtmp1 = (basis%spec(is)%atom(ia,:3)-basis%spec(js)%atom(ja,:3))
                   vtmp1 = matmul(vtmp1,basis%lat)
                   dtmp1 = modu(vtmp1)
                   if(dtmp1.lt.min_bond2) min_bond2 = dtmp1
                end if

             end do atomloop2
          end do specloop2
    
       end do atomloop1
    end do specloop1

    min_bond = ( min_bond1 + min_bond2 ) / 2._real32
    write(6,'(1X,"Avg min bulk bond: ",F0.3," Å")') min_bond
    write(6,'(1X,"Trans-interfacial scaling factor:",F0.3)') this%separation_scale
    this%axis = intf%axis
    call this%generate_perturbations(basis, intf%loc, min_bond, bulk_DON, print_shift_info_, seed_arr)


  end subroutine generate_intefaces_from_existing
!###############################################################################


!###############################################################################
  subroutine generate_interfaces( &
       this, basis_lw, basis_up, &
       miller_lw, miller_up, &
       surface_lw, surface_up, &
       thickness_lw, thickness_up, &
       num_layers_lw, num_layers_up, &
       use_pricel_lw, use_pricel_up, &
       is_layered_lw, is_layered_up, &
       elastic_constants_lw, elastic_constants_up, &
       print_lattice_match_info, print_termination_info, print_shift_info, &
       break_on_fail, &
       icheck_match, interface_idx, &
       generate_structures, &
       seed &
  )
    !! Generate interfaces from two bulk structures
    implicit none

    ! Arguments
    class(artemis_interface_generator_type), intent(inout) :: this
    !! Instance of artemis generator type
    type(basis_type), intent(in) :: basis_lw
    !! Lower bulk structure
    type(basis_type), intent(in) :: basis_up
    !! Upper bulk structure
    integer, intent(in), optional :: miller_lw(3)
    !! Miller indices for the lower bulk structure
    integer, intent(in), optional :: miller_up(3)
    !! Miller indices for the upper bulk structure
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

    logical, intent(in), optional :: use_pricel_lw
    !! Use primitive cell for lower bulk structure
    logical, intent(in), optional :: use_pricel_up
    !! Use primitive cell for upper bulk structure
    logical, intent(in), optional :: is_layered_lw
    !! Boolean whether the lower bulk structure is layered
    logical, intent(in), optional :: is_layered_up
    !! Boolean whether the upper bulk structure is layered

    real(real32), dimension(:), intent(in), optional :: elastic_constants_lw
    !! Elastic constants for the lower bulk structure
    real(real32), dimension(:), intent(in), optional :: elastic_constants_up
    !! Elastic constants for the upper bulk structure

    logical, intent(in), optional :: break_on_fail
    !! Break on failure
    logical, intent(in), optional :: print_lattice_match_info
    !! Print lattice match information
    logical, intent(in), optional :: print_termination_info
    !! Print termination information
    logical, intent(in), optional :: print_shift_info
    !! Print shift information
    integer, intent(in), optional :: icheck_match
    !! Index of the lattice match to check
    integer, intent(in), optional :: interface_idx
    !! Index of the interface to output
    logical, intent(in), optional :: generate_structures
    !! Boolean whether to generate structures or just print information
    integer, intent(in), optional :: seed
    !! Random seed for generating random numbers

    ! Local variables
    real(real32) :: avg_min_bond
    !! Average minimum bond length

    type(basis_type) :: basis_lw_, basis_up_
    !! Temporary basis structures
    type(basis_type) :: supercell_lw, supercell_up
    !! Copy of the basis structures
    type(basis_type) :: slab_lw, slab_up
    !! Slab structures
    type(basis_type) :: intf_basis
    !! Interface structure
    character(len=256) :: err_msg
    !! Error message

    integer :: j
    !! Loop index
    integer :: ifit, intf_start, intf_end
    !! Interface loop indices
    integer :: iterm_lw, term_lw_start_idx, term_lw_end_idx, term_lw_step
    !! Lower bulk termination loop indices
    integer :: iterm_up, term_up_start_idx, term_up_end_idx, term_up_step
    !! Upper bulk termination loop indices

    ! slab thickness variables
    integer :: ncells_lw, ncells_up
    !! Number of cells in the slab
    real(real32) :: height_lw, height_up
    !! Height of the slab
    real(real32) :: thickness_lw_, thickness_up_
    !! Thickness of the slab
    integer :: num_layers_lw_, num_layers_up_
    !! Number of layers in the slab
    logical :: use_pricel_lw_, use_pricel_up_
    !! Use primitive cell for lower and upper bulk structures
    logical :: is_layered_lw_, is_layered_up_
    !! Boolean whether the bulk structures are layered
    logical :: ludef_is_layered_lw, ludef_is_layered_up
    !! Boolean whether the user defined whether to use layered structures

    integer, dimension(3) :: miller_lw_, miller_up_
    !! Miller indices for the lower and upper bulk structures
    integer, dimension(2) :: surface_lw_, surface_up_
    !! Surface indices for the lower and upper bulk structures
    logical :: ludef_surface_lw, ludef_surface_up
    !! Boolean whether surfaces are defined 
    logical :: lcycle
    !! Boolean whether to skip the cycle

    logical :: break_on_fail_
    !! Boolean whether to break on failure
    logical :: print_lattice_match_info_, print_termination_info_, print_shift_info_
    !! Boolean whether to print lattice match, termination, and shift information
    integer :: num_seed
    !! Number of seeds for the random number generator.
    integer, dimension(:), allocatable :: seed_arr
    !! Array of seeds for the random number generator.
    integer :: icheck_match_
    !! Index of the lattice match to check
    integer :: interface_idx_
    !! Index of the interface to output
    logical :: generate_structures_
    !! Boolean whether to generate structures or just print information

    real(real32), dimension(:), allocatable :: elastic_constants_lw_, elastic_constants_up_
    !! Elastic constants for the lower and upper bulk structures

    type(bulk_DON_type), dimension(2) :: bulk_DON
    !! Distribution functions for the lower and upper bulk structures

    integer :: ntrans,iunique,itmp1,old_intf
    integer :: layered_axis_lw,layered_axis_up
    character(3) :: abc
    character(1024) :: pwd,intf_dir,dirpath,msg
    type(confine_type) :: confine
    type(latmatch_type) :: SAV
    type(term_arr_type) :: lw_term,up_term
    integer, dimension(3) :: ivtmp1
    real(real32), dimension(2) :: intf_loc
    real(real32), dimension(3) :: init_offset = [0._real32,0._real32,2._real32]
    !real(real32), dimension(3,3) :: mtmp1,DONsupercell_up%lat
    real(real32), dimension(3,3) :: tfmat
    integer, allocatable, dimension(:,:,:) :: lw_map,t1lw_map,t2lw_map
    integer, allocatable, dimension(:,:,:) :: up_map,t1up_map,t2up_map
    real(real32), allocatable, dimension(:,:) :: trans


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
    icheck_match_ = -1; interface_idx_ = -1
    if(present(icheck_match)) icheck_match_ = icheck_match
    if(present(interface_idx)) interface_idx_ = interface_idx
    break_on_fail_ = .true.
    if(present(break_on_fail)) break_on_fail_ = break_on_fail
    generate_structures_ = .true.
    if(present(generate_structures)) generate_structures_ = generate_structures

    if(.not.allocated(this%shifts)) call this%set_shift_method()


    !---------------------------------------------------------------------------
    ! Handle the elastic constants
    !---------------------------------------------------------------------------
    if(present(elastic_constants_lw))then
       if(allocated(elastic_constants_lw_)) deallocate(elastic_constants_lw_)
       allocate(elastic_constants_lw_(size(elastic_constants_lw)))
       elastic_constants_lw_ = elastic_constants_lw
    else
       if(allocated(elastic_constants_lw_)) deallocate(elastic_constants_lw_)
       allocate(elastic_constants_lw_(1))
       elastic_constants_lw_ = 0._real32
    end if
    if(present(elastic_constants_up))then
       if(allocated(elastic_constants_up_)) deallocate(elastic_constants_up_)
       allocate(elastic_constants_up_(size(elastic_constants_up)))
       elastic_constants_up_ = elastic_constants_up
    else
       if(allocated(elastic_constants_up_)) deallocate(elastic_constants_up_)
       allocate(elastic_constants_up_(1))
       elastic_constants_up_ = 0._real32
    end if


!!!-----------------------------------------------------------------------------
!!! determines the primitive and niggli reduced cell for each bulk
!!!-----------------------------------------------------------------------------
    call basis_lw_%copy(basis_lw)
    call basis_up_%copy(basis_up)
    write(6,*)
    use_pricel_lw_ = .false.
    use_pricel_up_ = .false.
    if(present(use_pricel_lw)) use_pricel_lw_ = use_pricel_lw
    if(present(use_pricel_up)) use_pricel_up_ = use_pricel_up
    if(use_pricel_lw_)then
       write(6,'(1X,"Using primitive cell for lower material")')
       call get_primitive_cell(basis_lw_)
    else
       write(6,'(1X,"Using supplied cell for lower material")')
       call reducer(basis_lw_)
       basis_lw_%lat=primitive_lat(basis_lw_%lat)
    end if
    if(use_pricel_up_)then
       write(6,'(1X,"Using primitive cell for upper material")')
       call get_primitive_cell(basis_up_)
    else
       write(6,'(1X,"Using supplied cell for upper material")')
       call reducer(basis_up_)
       basis_up_%lat=primitive_lat(basis_up_%lat)
    end if
    write(6,*)


    surface_lw_ = 0
    surface_up_ = 0
    if(present(surface_lw))then
       select case(size(surface_lw, dim=1))
       case(1)
          surface_lw_ = surface_lw(1)
       case(2)
          surface_lw_ = surface_lw
       case default
          write(msg,'(A,I0,A)') &
               "ERROR: The surface vector for the lower material has ", &
               size(surface_lw, dim=1), " components. It should have 1 or 2."
          call err_abort(trim(msg),fmtd=.true.)
       end select
    end if
    if(present(surface_up))then
       select case(size(surface_up, dim=1))
       case(1)
          surface_up_ = surface_up(1)
       case(2)
          surface_up_ = surface_up
       case default
          write(msg,'(A,I0,A)') &
               "ERROR: The surface vector for the upper material has ", &
               size(surface_up, dim=1), " components. It should have 1 or 2."
          call err_abort(trim(msg),fmtd=.true.)
       end select
    end if

    ludef_surface_lw = .false.
    ludef_surface_up = .false.
    if(all(surface_lw_.gt.0)) ludef_surface_lw = .true.
    if(all(surface_up_.gt.0)) ludef_surface_up = .true.

    miller_lw_ = 0
    miller_up_ = 0
    if(present(miller_lw)) miller_lw_ = miller_lw
    if(present(miller_up)) miller_up_ = miller_up

    print_lattice_match_info_ = .false.
    print_termination_info_ = .false.
    print_shift_info_ = .false.
    if(present(print_lattice_match_info)) print_lattice_match_info_ = print_lattice_match_info
    if(present(print_termination_info)) print_termination_info_ = print_termination_info
    if(present(print_shift_info)) print_shift_info_ = print_shift_info
    
    if(.not.allocated(this%structures)) allocate(this%structures(0))

    thickness_lw_ = 10._real32
    thickness_up_ = 10._real32
    num_layers_lw_ = 0
    num_layers_up_ = 0
    if(present(num_layers_lw)) num_layers_lw_ = num_layers_lw
    if(present(num_layers_up)) num_layers_up_ = num_layers_up
    if(present(thickness_lw)) thickness_lw_ = thickness_lw
    if(present(thickness_up)) thickness_up_ = thickness_up
    if(num_layers_lw_.le.0.and.thickness_lw_.le.0._real32)then
       write(msg,'(A,I0,A)') &
            "ERROR: The number of layers for the lower material is ", &
            num_layers_lw_, " and the thickness is ", thickness_lw_, &
            " One of these must be greater than 0."
       call err_abort(trim(msg),fmtd=.true.)
    end if
    if(num_layers_up_.le.0.and.thickness_up_.le.0._real32)then
       write(msg,'(A,I0,A)') &
            "ERROR: The number of layers for the upper material is ", &
            num_layers_up_, " and the thickness is ", thickness_up_, &
            " One of these must be greater than 0."
       call err_abort(trim(msg),fmtd=.true.)
    end if


    
!!!-----------------------------------------------------------------------------
!!! investigates individual bulks and their bondlengths
!!!-----------------------------------------------------------------------------
    avg_min_bond = &
         ( get_min_bulk_bond(basis_lw_) + get_min_bulk_bond(basis_up_) )/2._real32
    write(6,'(1X,"Avg min bulk bond: ",F0.3," Å")') avg_min_bond
    write(6,'(1X,"Trans-interfacial scaling factor: ",F0.3)') this%separation_scale
    if(this%shift_method.eq.-1) this%num_shifts=1
    

!!!-----------------------------------------------------------------------------
!!! gets bulk DONs, if shift_method = 4
!!!-----------------------------------------------------------------------------
    allocate(lw_map(basis_lw_%nspec,maxval(basis_lw_%spec(:)%num,dim=1),2))
    allocate(up_map(basis_up_%nspec,maxval(basis_up_%spec(:)%num,dim=1),2))    
    if(this%shift_method.eq.4.or.this%shift_method.eq.0)then
       lw_map=0
       bulk_DON(1)%spec=gen_DON(basis_lw_%lat,basis_lw_,&
            dist_max=this%bondlength_cutoff,&
            scale_dist=.false.,&
            norm=.true.)
       if(all(abs(bulk_DON(1)%spec(1)%atom(:,:)).lt.1._real32))then
          open(unit=13,file="lw_DON.dat")
          do j=1,1000
             write(13,*) &
                  (j-1)*this%bondlength_cutoff/1000,&
                  bulk_DON(1)%spec(1)%atom(1,j)
          end do
          close(13)
          write(err_msg,'(A,F0.3,A)') &
               "The lower bulk DON identified no atoms within the bulk cutoff" //&
               &" distance (MAX_BONDLENGTH = ", this%bondlength_cutoff, " Å)." // achar(10) //&
               &" To proceed with the current shift method," //&
               &" increase MAX_BONDLENGTH."
          call err_abort(trim(err_msg),fmtd=.true.)
       end if
       !call exit()
       up_map=0
       bulk_DON(2)%spec=gen_DON(basis_up_%lat,basis_up_,&
            dist_max=this%bondlength_cutoff,&
            scale_dist=.false.,&
            norm=.true.)
       if(all(abs(bulk_DON(2)%spec(1)%atom(:,:)).lt.1._real32))then
          open(unit=13,file="up_DON.dat")
          do j=1,1000
             write(13,*) &
                  (j-1)*this%bondlength_cutoff/1000,&
                  bulk_DON(2)%spec(1)%atom(1,j)
          end do
          close(13)
          write(err_msg,'(A,F0.3,A)') &
               "The upper bulk DON identified no atoms within the bulk cutoff" //&
               &" distance (MAX_BONDLENGTH = ", this%bondlength_cutoff, " Å)." // achar(10) //&
               &" To proceed with the current shift method," //&
               &" increase MAX_BONDLENGTH."
          call err_abort(trim(err_msg),fmtd=.true.)
       end if
    else
       lw_map=-1
       up_map=-1       
    end if


!!!-----------------------------------------------------------------------------
!!! checks whether system appears layered
!!!-----------------------------------------------------------------------------
    if(present(is_layered_lw))then
       is_layered_lw_ = is_layered_lw
       ludef_is_layered_lw = .true.
    else
       is_layered_lw_ = .false.
       ludef_is_layered_lw = .false.
    end if
    if(present(is_layered_up))then
       is_layered_up_ = is_layered_up
       ludef_is_layered_up = .true.
    else
       is_layered_up_ = .false.
       ludef_is_layered_up = .false.
    end if



    layered_axis_lw=get_layered_axis(basis_lw_%lat,basis_lw_)
    if(.not.is_layered_lw_.and.layered_axis_lw.gt.0)then
       ivtmp1=0
       ivtmp1(layered_axis_lw)=1
       if(ludef_is_layered_lw)then
          write(msg,'("Lower crystal appears layered along axis ",I0,"\n&
               &Partial layer terminations will be generated\n&
               &We suggest using LW_MILLER =",3(1X,I1))') layered_axis_lw,ivtmp1
          call print_warning(trim(msg))
       else
          write(msg,'("Lower crystal has been identified as layered\nalong",3(1X,I1),"\n&
               &Confining crystal to this plane and\nstoichiometric terminations.\n&
               &If you don''t want this, set\nLW_LAYERED = .FALSE.")') &
               ivtmp1
          call print_warning(trim(msg))
          miller_lw_=ivtmp1
          is_layered_lw_=.true.
       end if
    elseif(is_layered_lw_.and.layered_axis_lw.gt.0.and.all(miller_lw_.eq.0))then
       miller_lw_(layered_axis_lw)=1
    end if

    layered_axis_up=get_layered_axis(basis_up_%lat,basis_up_)
    if(.not.is_layered_up_.and.layered_axis_up.gt.0)then
       ivtmp1=0
       ivtmp1(layered_axis_up)=1
       if(ludef_is_layered_up)then
          write(msg,'("Upper crystal appears layered along axis ",I0,"\n&
               &Partial layer terminations will be generated\n&
               &We suggest using UP_MILLER =",3(1X,I1))') layered_axis_up,ivtmp1
          call print_warning(trim(msg))
       else
          write(msg,'("Upper crystal has been identified as layered\nalong",3(1X,I1),"\n&
               &Confining crystal to this plane and\nstoichiometric terminations.\n&
               &If you don''t want this, set\nUP_LAYERED = .FALSE.")') &
               ivtmp1
          call print_warning(trim(msg))
          miller_up_=ivtmp1
          is_layered_up_=.true.
       end if
    elseif(is_layered_up_.and.layered_axis_up.gt.0.and.all(miller_up_.eq.0))then
       miller_up_(layered_axis_up)=1
    end if


!!!-----------------------------------------------------------------------------
!!! Finds and stores the best matches between the materials
!!!-----------------------------------------------------------------------------
   !  call getcwd(pwd)
    old_intf = -1
    abc="abc"
    if(any(miller_lw_.ne.0))then
       if(this%match_method.ne.0)then
          abc="ab"
          tfmat=planecutter(basis_lw_%lat,real(miller_lw_,real32))
          call transformer(basis_lw_,tfmat,lw_map)
          SAV=get_best_match(&
               this%tolerance,&
               basis_lw_,basis_up_,&
               trim(abc),"abc",print_lattice_match_info_,ierror,imatch=this%match_method)
       elseif(any(miller_up_.ne.0))then
          SAV=get_best_match(&
               this%tolerance,&
               basis_lw_,basis_up_,&
               trim(abc),"abc",print_lattice_match_info_,ierror,imatch=this%match_method,&
               plane1=miller_lw_,plane2=miller_up_,nmiller=this%max_num_planes)
       else
          SAV=get_best_match(&
               this%tolerance,&
               basis_lw_,basis_up_,&
               trim(abc),"abc",print_lattice_match_info_,ierror,imatch=this%match_method,&
               plane1=miller_lw_,nmiller=this%max_num_planes)
       end if
    elseif(any(miller_up_.ne.0))then
       SAV=get_best_match(&
            this%tolerance,&
            basis_lw_,basis_up_,&
            trim(abc),"abc",print_lattice_match_info_,ierror,imatch=this%match_method,&
            plane2=miller_up_,nmiller=this%max_num_planes)
    else
       SAV=get_best_match(&
            this%tolerance,&
            basis_lw_,basis_up_,&
            trim(abc),"abc",print_lattice_match_info_,ierror,imatch=this%match_method,&
            nmiller=this%max_num_planes)
    end if
    if(min(this%tolerance%nstore,SAV%nfit).eq.0)then
       write(0,'("No matches found.")')
       write(0,'("Exiting...")')
       call exit()
    else
       write(0,'(1X,"Number of matches found: ",I0)')&
            min(this%tolerance%nstore,SAV%nfit)
    end if
    write(6,'(1X,"Maximum number of generated interfaces will be: ",I0)')&
         this%max_num_terms*this%num_shifts*this%tolerance%nstore
    if(.not.generate_structures_)then
       write(0,'(1X,"Told not to generate interfaces, just find matches.")')
       write(0,'("Exiting...")')
       call exit()
    end if

       
!!!-----------------------------------------------------------------------------
!!! Saves current directory and moves to new directory
!!!-----------------------------------------------------------------------------
   !  call system('mkdir -p '//trim(adjustl(dirname)))
   !  call chdir(dirname)
   !  call getcwd(intf_dir)

    if(interface_idx_.gt.0)then
       intf_start=interface_idx_
       intf_end=interface_idx_
       write(6,'(1X,"Generating only interfaces for match ",I0)') interface_idx_
    else
       intf_start=1
       intf_end=min(this%tolerance%nstore,SAV%nfit)
    end if
    iunique=0
!!!-----------------------------------------------------------------------------
!!! Applies the best match transformations
!!!-----------------------------------------------------------------------------
    intf_loop: do ifit = intf_start, intf_end
       write(6,'("Fit number: ",I0)') ifit
       call supercell_lw%copy(basis_lw_)
       call supercell_up%copy(basis_up_)
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
          !! Issue with using this method when large deformations result in large
          !! angle changes. REMOVING IT FOR NOW AND RETURNING TO CALCULATING DONS
          !! FOR THE SUPERCELL.
          t1up_map=0 !TEMPORARY TO USE SUPERCELL DONS.
          !do i=1,2
          !   mtmp1(i,:) = &
          !        ( modu(lw_lat(i,:)) )*uvec(supercell_up%lat(i,:))
          !end do
          !mtmp1(3,:) = supercell_up%lat(3,:)
          !DONsupercell_up%lat = matmul(mtmp1,inverse(real(SAV%tf2(ifit,:,:),real32)))
          !if(ierror.eq.1)then
          !   write(0,*) "#####################################"
          !   write(0,*) "ifit", ifit
          !   write(0,*) "undeformed lattice"
          !   write(0,'(3(2X,F6.2))') (mtmp1(i,:),i=1,3)
          !   write(0,*)
          !   write(0,*) "deformed lattice"
          !   write(0,'(3(2X,F8.4))') (DONsupercell_up%lat(i,:),i=1,3)
          !   write(0,*)
          !end if
          deallocate(bulk_DON(2)%spec)
          bulk_DON(2)%spec=gen_DON(supercell_up%lat,supercell_up,&
               dist_max=this%bondlength_cutoff,&
               scale_dist=.false.,&
               norm=.true.)
          !call err_abort_print_struc(basis_up_,"bulk_up_term.vasp",&
          !     "",.false.)
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
       call gldfnd(confine,supercell_lw,supercell_lw,trans,ntrans)
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
       if(.not.compare_stoichiometry(basis_lw_,supercell_lw))then
          write(0,'(1X,"ERROR: Internal error in generate_interfaces")')
          write(0,'(2X,"The gldfnd subroutine could not reproduce a valid primitive cell for the lower material on match ",I0)') ifit
          if(ierror.eq.1)then
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
            lprint = print_termination_info_, layer_sep = this%layer_separation_cutoff(1), &
            break_on_fail = break_on_fail_ &
       )
       if(lw_term%nterm .eq. 0)then
          write(0,'("WARNING: &
               &No terminations found for lower material Miller plane &
               &(",3(1X,I0)," )")' &
          ) SAV%tf1(ifit,3,1:3)
          cycle intf_loop
       end if
       if(any(surface_lw_.gt.lw_term%nterm))then
          write(msg, '("surface_lw_ACE VALUES INVALID!\nOne or more value &
               &exceeds the maximum number of terminations in the &
               structure.\n&
               &  Supplied values: ",I0,1X,I0,"\n&
               &  Maximum allowed: ",I0)') surface_lw_, lw_term%nterm
          call err_abort(trim(msg),fmtd=.true.)
       end if


       !!-----------------------------------------------------------------------
       !! Sort out ladder rungs (checks whether the material is centrosymmetric)
       !!-----------------------------------------------------------------------
       !call setup_ladder(supercell_lw%lat,supercell_lw,this%axis,lw_term)
       if(sum(lw_term%arr(:)%natom)*lw_term%nstep.ne.supercell_lw%natom)then
          write(msg, '("ERROR: Number of atoms in lower layers not correct: "&
               &I0,2X,I0)') sum(lw_term%arr(:)%natom)*lw_term%nstep,supercell_lw%natom
          call err_abort(trim(msg),fmtd=.true.)
       end if
       call set_layer_tol(lw_term)


       !!-----------------------------------------------------------------------
       !! Defines height of lower slab from user-defined values
       !!-----------------------------------------------------------------------
       call set_slab_height(supercell_lw,t1lw_map,lw_term,surface_lw_,&
            height_lw,num_layers_lw_, thickness_lw_,ncells_lw,&
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
       call gldfnd(confine,supercell_up,supercell_up,trans,ntrans)
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
       if(.not.compare_stoichiometry(basis_up_,supercell_up))then
          write(0,'(1X,"ERROR: Internal error in generate_interfaces")')
          write(0,'(2X,"The gldfnd subroutine could not reproduce a valid primitive cell for the upper material on match ",I0)') ifit
          if(ierror.eq.1)then
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
            lprint = print_termination_info_, layer_sep = this%layer_separation_cutoff(2), &
            break_on_fail = break_on_fail_ &
       )
       if(up_term%nterm .eq. 0)then
          write(0,'("WARNING: &
               &No terminations found for upper material Miller plane &
               &(",3(1X,I0)," )")' &
          ) SAV%tf2(ifit,3,1:3)
          cycle intf_loop
       end if
       if(any(surface_up_.gt.up_term%nterm))then
          write(msg, '("surface_up_ACE VALUES INVALID!\nOne or more value &
               &exceeds the maximum number of terminations in the &
               structure.\n&
               &  Supplied values: ",I0,1X,I0,"\n&
               &  Maximum allowed: ",I0)') surface_up_, up_term%nterm
          call err_abort(trim(msg),fmtd=.true.)
       end if


       !!-----------------------------------------------------------------------
       !! Sort out ladder rungs (checks whether the material is centrosymmetric)
       !!-----------------------------------------------------------------------
       !call setup_ladder(supercell_up%lat,supercell_up,this%axis,up_term)
       if(sum(up_term%arr(:)%natom)*up_term%nstep.ne.supercell_up%natom)then
          write(msg, '("ERROR: Number of atoms in upper layers not correct: "&
               &I0,2X,I0)') sum(up_term%arr(:)%natom)*up_term%nstep,supercell_up%natom
          call err_abort(trim(msg),fmtd=.true.)
       end if
       call set_layer_tol(up_term)


       !!-----------------------------------------------------------------------
       !! Defines height of upper slab from user-defined values
       !!-----------------------------------------------------------------------
       call set_slab_height(supercell_up,t1up_map,up_term,surface_up_,&
            height_up,num_layers_up_, thickness_up_, ncells_up,&
            term_up_start_idx,term_up_end_idx,term_up_step &
       )
       if(term_up_end_idx.gt.this%max_num_terms) term_up_end_idx = this%max_num_terms


       !!-----------------------------------------------------------------------
       !! Print termination plane locations
       !!-----------------------------------------------------------------------
       write(6,'(1X,"Number of unique terminations: ",I0,2X,I0)') &
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
          call build_slab(slab_lw,t2lw_map,lw_term,[iterm_lw,surface_lw_(2)],&
               thickness_lw_, ncells_lw, num_layers_lw_, height_lw,&
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
             call build_slab(slab_up,t2up_map,up_term,[iterm_up,surface_up_(2)],&
                  thickness_up_, ncells_up, num_layers_up_, height_up,&
                  "up",lcycle, &
                  vacuum = this%vacuum_gap &
             )
             if(lcycle) cycle up_term_loop

             
             !!-----------------------------------------------------------------
             !! Checks stoichiometry
             !!-----------------------------------------------------------------
             if(slab_lw%nspec.ne.basis_lw_%nspec.or.any(&
                  (basis_lw_%spec(1)%num*slab_lw%spec(:)%num)&
                  /slab_lw%spec(1)%num.ne.basis_lw_%spec(:)%num))then
                write(6,'("WARNING: This lower surface termination is not &
                     &stoichiometric")')
                if(is_layered_lw_)then
                   write(6,'(2X,"As lower structure is layered, stoichiometric &
                        &surfaces are required.")')
                   write(6,'(2X,"Skipping this termination...")')
                   cycle lw_term_loop
                end if
             end if
             if(slab_up%nspec.ne.basis_up_%nspec.or.any(&
                  (basis_up_%spec(1)%num*slab_up%spec(:)%num)&
                  /slab_up%spec(1)%num.ne.basis_up_%spec(:)%num))then
                write(6,'("WARNING: This upper surface termination is not &
                     &stoichiometric")')
                if(is_layered_up_)then
                   write(6,'(2X,"As upper structure is layered, stoichiometric &
                        &surfaces are required.")')
                   write(6,'(2X,"Skipping this termination...")')
                   cycle up_term_loop
                end if
             end if


             !!-----------------------------------------------------------------
             !! Use the bulk moduli to determine the strain sharing
             !!-----------------------------------------------------------------
             if( all(abs(elastic_constants_lw_).gt.0.E0) .and. &
                  all(abs(elastic_constants_up_).gt.0.E0) &
             )then
                call share_strain(slab_lw%lat,slab_up%lat,&
                     elastic_constants_lw_(1), &
                     elastic_constants_up_(1), &
                     lcompensate = .not.this%fix_normal &
                )
             end if
             

             !!-----------------------------------------------------------------
             !! Merge the two bases and lattices and define the interface loc
             !!-----------------------------------------------------------------
             intf_basis = basis_stack(&
                  basis1 = slab_lw, basis2 = slab_up, &
                  axis = this%axis, offset = init_offset(:), &
                  map1 = t2lw_map, map2 = t2up_map &
             )
             intf_loc(1) = ( modu(slab_lw%lat(this%axis,:)) + 0.5_real32*init_offset(this%axis) - &
                  this%vacuum_gap)/modu(intf_basis%lat(this%axis,:))
             intf_loc(2) = ( modu(slab_lw%lat(this%axis,:)) + modu(slab_up%lat(this%axis,:)) + &
                  1.5_real32*init_offset(this%axis) - 2._real32*this%vacuum_gap )/modu(intf_basis%lat(this%axis,:))
             if(ierror.ge.1)then
                write(0,*) "interface:",intf_loc
                if(ierror.eq.1.and.iunique.eq.icheck_match_-1)then
                  !  call chdir(intf_dir)
                   call err_abort_print_struc(slab_lw,"lw_term.vasp",&
                        "",.false.)
                   call err_abort_print_struc(slab_up,"up_term.vasp",&
                        "As IPRINT = 1 and ICHECK has been set, &
                        &code is now exiting...")
                elseif(ierror.eq.2.and.iunique.eq.icheck_match_-1)then
                  !  call chdir(intf_dir)
                   call err_abort_print_struc(intf_basis,"test_intf.vasp",&
                        "As IPRINT = 2 and ICHECK has been set, &
                        &code is now exiting...")
                end if
             end if


             !!-----------------------------------------------------------------
             !! Saves current directory and moves to new directory
             !!-----------------------------------------------------------------
             if(this%num_structures.gt.old_intf)then
                iunique=iunique+1
               !  if(this%shift_method.gt.0.and.this%num_shifts.gt.1) &
               !       write(6,'(1X,"Generating shifts for unique interface ",&
               !       &I0,":")') iunique
               !  write(dirpath,'(A,I0.2)') trim(adjustl(subdir_prefix)),iunique
               !  call system('mkdir -p '//trim(adjustl(dirpath)))
             else
               !  write(dirpath,'(A,I0.2)') trim(adjustl(subdir_prefix)),iunique
             end if
            !  call chdir(dirpath)
             old_intf = this%num_structures

             
             !!-----------------------------------------------------------------
             !! Writes information of current match to file in save directory
             !!-----------------------------------------------------------------
             call  output_intf_data(SAV, ifit, lw_term, iterm_lw, up_term, iterm_up,&
                  use_pricel_lw_,use_pricel_up_)


             !!-----------------------------------------------------------------
             !! Generates shifts and swaps and prints the subsequent structures
             !!-----------------------------------------------------------------
             call this%generate_perturbations( &
                  intf_basis, intf_loc, avg_min_bond, &
                  bulk_DON, &
                  print_shift_info_, &
                  seed_arr, &
                  t2lw_map &
             )

             if(this%num_structures.ge.this%max_num_structures) exit intf_loop
             !call chdir(dirname)
            !  call chdir(intf_dir)

             if(ludef_surface_up) exit up_term_loop
          end do up_term_loop
          if(ludef_surface_lw) exit lw_term_loop
       end do lw_term_loop
       !!-----------------------------------------------------------------------
       !! Returns to working directory
       !!-----------------------------------------------------------------------
      !  call chdir(intf_dir)

    end do intf_loop

   !  call chdir(pwd)


    return
  end subroutine generate_interfaces
!###############################################################################


!!!#############################################################################
!!! Takes input interface structure and generates a set of shifts and swaps.
!!! Prints these new structures to POSCARs.
!!!#############################################################################
!!! ISWAP METHOD NOT YET SET UP
  subroutine generate_shifts_and_swaps( &
       this, basis, intf_loc, bond, bulk_DON, print_shift_info, seed_arr, map &
  )
    implicit none
    class(artemis_interface_generator_type), intent(inout) :: this
    type(basis_type), intent(in) :: basis
    real(real32), dimension(2), intent(in) :: intf_loc
    real(real32), intent(in) :: bond
    type(bulk_DON_type), dimension(2), intent(in) :: bulk_DON
    !! Distribution functions for the lower and upper bulk structures
    logical, intent(in) :: print_shift_info
    integer, dimension(:), intent(in) :: seed_arr
    integer, dimension(:,:,:), optional, intent(in) :: map

    integer :: shift_unit
    integer :: ounit,iaxis,k,l
    integer :: ngen_swaps,nswaps_per_cell
    real(real32) :: dtmp1
    type(basis_type) :: tbas
    type(bond_type) :: min_bond
    character(1024) :: filename,dirpath,pwd1,pwd2,msg
    integer, dimension(3) :: abc
    real(real32), dimension(3) :: toffset
    type(basis_type), allocatable, dimension(:) :: bas_arr
    real(real32), allocatable, dimension(:,:) :: output_shifts



!!!-----------------------------------------------------------------------------
!!! Sets up shift axis
!!!-----------------------------------------------------------------------------
    abc = [ 1, 2, 3 ]
    abc = cshift(abc,this%axis)


!!!-----------------------------------------------------------------------------
!!! Sets up and moves to appropriate directories
!!!-----------------------------------------------------------------------------
   !  call getcwd(pwd1)
   !  if(this%shift_method.gt.0.or.this%num_shifts.gt.1)then
   !     call system('mkdir -p '//trim(adjustl(shiftdir)))
   !     call chdir(shiftdir)
   !  end if
   !  call getcwd(pwd2)
   !  open(newunit=shift_unit,file="shift_vals.txt")
   !  write(shift_unit,&
   !       '("# interface_num    shift (a,b,c) units=(direct,direct,Å)")')


!!!-----------------------------------------------------------------------------
!!! Generates sets of shifts based on shift version
!!!-----------------------------------------------------------------------------
    if(this%shift_method.eq.0.or.this%shift_method.eq.1) allocate(output_shifts(this%num_shifts,3))
    select case(this%shift_method)
    case(1)
       output_shifts(1,:3)=0._real32
       do k=2,this%num_shifts
          do iaxis=1,2
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
               lprint=print_shift_info, &
               bulk_DON=bulk_DON,bulk_map=map,&
               max_bondlength=this%bondlength_cutoff)
       else
          output_shifts = get_shifts_DON(&
               bas=basis,&
               axis=this%axis,&
               intf_loc=intf_loc,&
               nstore=this%num_shifts, &
               c_scale=this%separation_scale, &
               offset=this%shifts(1,:3),&
               lprint=print_shift_info,&
               max_bondlength=this%bondlength_cutoff)
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
       output_shifts(:,this%axis) = output_shifts(:,this%axis)*modu(basis%lat(this%axis,:))
    end if


!!!-----------------------------------------------------------------------------
!!! Prints number of shifts to terminal
!!!-----------------------------------------------------------------------------
    write(6,'(3X,"Number of unique shifts structures: ",I0)') this%num_shifts


!!!-----------------------------------------------------------------------------
!!! Determines number of swaps across the interface
!!!-----------------------------------------------------------------------------
    nswaps_per_cell=nint(this%swap_density*get_area([basis%lat(abc(1),:)],[basis%lat(abc(2),:)]))
    if(this%swap_method.ne.0)then
       write(6,&
            '(" Generating ",I0," swaps per structure ")') nswaps_per_cell
    end if


!!!-----------------------------------------------------------------------------
!!! Prints each unique shift structure
!!!-----------------------------------------------------------------------------
    shift_loop: do k=1,this%num_shifts
       call tbas%copy(basis)
       toffset=output_shifts(k,:3)
       do iaxis=1,2
          call shift_region(tbas,this%axis,&
               intf_loc(1),intf_loc(2),&
               shift_axis=iaxis,shift=toffset(iaxis),renorm=.true.)
       end do
       dtmp1=modu(tbas%lat(this%axis,:))
       call set_vacuum(&
            basis=tbas,&
            axis=this%axis,loc=maxval(intf_loc(:)),&
            vac=toffset(this%axis))
       dtmp1=minval(intf_loc(:))*dtmp1/modu(tbas%lat(this%axis,:))
       call set_vacuum(&
            basis=tbas,&
            axis=this%axis,loc=dtmp1,&
            vac=toffset(this%axis))
       min_bond = get_shortest_bond(tbas)
       if(min_bond%length.le.1.5_real32)then
          write(msg,'("Smallest bond in the interface structure is\nless than 1.5 Å.")')
          call print_warning(trim(msg))
          write(6,'(2X,"bond length: ",F9.6)') min_bond%length
          write(6,'(2X,"atom 1:",I4,2X,I4)') min_bond%atoms(1,:)
          write(6,'(2X,"atom 2:",I4,2X,I4)') min_bond%atoms(2,:)
       end if


       !!-----------------------------------------------------------------------
       !! prints shift vector to shift_vals.txt
       !!-----------------------------------------------------------------------
      !  write(shift_unit,'(2X,I0.2,15X,"(",2(" ",F9.6,", ")," ",F9.6," )")') &
      !       k,toffset(:)


       !!-----------------------------------------------------------------------
       !! Merges lower and upper materials
       !! Writes interfaces to output directories
       !!-----------------------------------------------------------------------
      !  ounit=100+intf
      !  if(this%shift_method.gt.0.or.this%num_shifts.gt.1)then
      !     write(dirpath,'(A,I0.2)') trim(adjustl(subdir_prefix)),k
      !     call system('mkdir -p '//trim(adjustl(dirpath)))
      !     write(filename,'(A,"/",A)') trim(adjustl(dirpath)),trim(out_filename)
      !  else
      !     filename = trim(out_filename)
      !  end if
      !  write(6,'(2X,"Writing interface ",I0,"...")') intf
      !  open(unit=ounit,file=trim(adjustl(filename)))
      !  call geom_write(ounit,tbas)
      !  close(ounit)
       this%structures = [ this%structures, tbas ]
       this%num_structures = size(this%structures, dim = 1)
       if(this%num_structures.ge.this%max_num_structures) return


       !!-----------------------------------------------------------------------
       !! Performs swaps within the shifted structures if requested
       !!-----------------------------------------------------------------------
       if_swap: if(this%swap_method.ne.0)then
          bas_arr = rand_swapper(tbas%lat,tbas,this%axis,this%swap_depth,&
               nswaps_per_cell,this%num_swaps,intf_loc,this%swap_method,seed_arr,sigma=this%swap_sigma,&
               require_mirror=this%require_mirror_swaps)
          ngen_swaps = this%num_swaps
          LOOPswaps: do l=1,this%num_swaps
             if (bas_arr(l)%nspec.eq.0) then
                ngen_swaps = l - 1
                exit LOOPswaps
             end if
          end do LOOPswaps
          if(ngen_swaps.eq.0)then
             exit if_swap
          end if
         !  call chdir(dirpath)
         !  call system('mkdir -p '//trim(adjustl(swapdir)))
         !  call chdir(swapdir)
         !  write(6,'(3X,"Number of unique swap structures: ",I0)') ngen_swaps
          this%structures = [ this%structures, bas_arr(1:ngen_swaps) ]
         !  do l=1,ngen_swaps
         !     write(dirpath,'(A,I0.2)') trim(adjustl(subdir_prefix)),l
         !     call system('mkdir -p '//trim(adjustl(dirpath)))
         !     write(filename,'(A,"/",A)') &
         !          trim(adjustl(dirpath)),trim(out_filename)
         !     ounit=100+l
         !     write(6,'(3X,"Writing swap ",I0,"...")') l
         !     open(unit=ounit,file=trim(adjustl(filename)))
         !     call geom_write(ounit,bas_arr(l))
         !     close(ounit)
         !  end do
          deallocate(bas_arr)
         !  call chdir(pwd2)
       end if if_swap


    end do shift_loop
   !  call chdir(pwd1)
   !  close(unit=shift_unit)


  end subroutine generate_shifts_and_swaps
!!!#############################################################################


!!!#############################################################################
!!! write structure data in each structure directory
!!!#############################################################################
  subroutine output_intf_data(SAV, ifit, lw_term, term_lw_idx, up_term, term_up_idx, lw_pricel,up_pricel)
    implicit none
    integer :: unit

    integer, intent(in) :: ifit, term_lw_idx, term_up_idx
    logical, intent(in) :: lw_pricel,up_pricel
    type(term_arr_type), intent(in) :: lw_term, up_term
    type(latmatch_type), intent(in) :: SAV


    
    unit=99
    open(unit=unit, file="struc_dat.txt")
    write(unit,'("Lower material primitive cell used: ",L1)') lw_pricel
    write(unit,'("Upper material primitive cell used: ",L1)') lw_pricel
    write(unit,*)
    write(unit,'("Lattice match:")')
    write(unit,'((1X,3(3X,A1),3X,3(3X,A1)),3(/,2X,3(I3," "),3X,3(I3," ")))') &
         SAV%abc,SAV%abc,&
         SAV%tf1(ifit,1,1:3),SAV%tf2(ifit,1,1:3),&
         SAV%tf1(ifit,2,1:3),SAV%tf2(ifit,2,1:3),&
         SAV%tf1(ifit,3,1:3),SAV%tf2(ifit,3,1:3)
    write(unit,'(" vector mismatch (%) = ",F0.9)') SAV%tol(ifit,1)
    write(unit,'(" angle mismatch (°)  = ",F0.9)') SAV%tol(ifit,2)*180/pi
    write(unit,'(" area mismatch (%)   = ",F0.9)') SAV%tol(ifit,3)
    write(unit,*)
    write(unit,'(" Lower crystal Miller plane: ",3(I3," "))') SAV%tf1(ifit,3,1:3)
    write(unit,'(" Lower termination")')
    write(unit,'(1X,"Term.",3X,"Min layer loc",3X,"Max layer loc",3X,"no. atoms")')
    write(unit,'(1X,I3,8X,F7.5,9X,F7.5,8X,I3)') &
            term_lw_idx,lw_term%arr(term_lw_idx)%hmin,lw_term%arr(term_lw_idx)%hmax,lw_term%arr(term_lw_idx)%natom
    write(unit,*)
    write(unit,'(" Upper crystal Miller plane: ",3(I3," "))') SAV%tf2(ifit,3,1:3)
    write(unit,'(" Upper termination")')
    write(unit,'(1X,"Term.",3X,"Min layer loc",3X,"Max layer loc",3X,"no. atoms")')
    write(unit,'(1X,I3,8X,F7.5,9X,F7.5,8X,I3)') &
            term_up_idx,up_term%arr(term_up_idx)%hmin,up_term%arr(term_up_idx)%hmax,up_term%arr(term_up_idx)%natom
    write(unit,*)
    close(unit)
    
    return
  end subroutine output_intf_data
!!!#############################################################################

end module artemis__interface_generator
