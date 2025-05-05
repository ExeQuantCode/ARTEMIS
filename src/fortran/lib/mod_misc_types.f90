module artemis__misc_types
  !! Module containing custom derived types for ARTEMIS
  use artemis__constants, only: real32, pi
  use artemis__misc, only: to_lower
  use artemis__geom_rw, only: basis_type, geom_write
  use artemis__geom_utils, only: MATNORM
  implicit none


  private

  public :: struc_data_type
  public :: latmatch_type
  public :: tol_type
  public :: abstract_artemis_generator_type


  type struc_data_type
     integer :: match_idx = 0
     integer :: match_and_term_idx = 0
     integer :: shift_idx = 0
     integer :: swap_idx  = 0
     logical :: from_pricel_lw = .false.
     logical :: from_pricel_up = .false.
     integer, dimension(2) :: term_lw_idx = 0
     integer, dimension(2) :: term_up_idx = 0
     real(real32), dimension(4) :: term_lw_bounds = 0._real32
     real(real32), dimension(4) :: term_up_bounds = 0._real32
     integer, dimension(2) :: term_lw_natom = 0
     integer, dimension(2) :: term_up_natom = 0
     integer, dimension(3,3) :: transform_lw = 0
     integer, dimension(3,3) :: transform_up = 0
     real(real32) :: approx_thickness_lw = 0._real32
     real(real32) :: approx_thickness_up = 0._real32
     real(real32), dimension(3) :: mismatch
     real(real32), dimension(3) :: shift = 0._real32
     ! real(real32), dimension(:,:) :: swaps !!! UNSURE HOW TO DO THIS
     real(real32) :: swap_density = 0._real32
     real(real32), dimension(2) :: approx_eff_swap_conc = 0._real32

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
     integer :: nfit
     integer :: max_num_matches = 5
     logical :: reduce = .false.
     logical :: reduced = .false.
     character(1) :: abc(3)= [ 'a', 'b', 'c' ]

     integer, dimension(2) :: axes
     integer, allocatable, dimension(:,:,:) :: tf1,tf2
     real(real32), allocatable, dimension(:,:) :: tol
     real(real32), dimension(3,3) :: lat1,lat2
   contains
     procedure, pass(this) :: init => latmatch_init
     procedure, pass(this) :: constrain_axes
  end type latmatch_type

  type tol_type
     integer :: maxfit = 100
     integer :: maxsize = 10
     real(real32) :: maxlen  = 20._real32
     real(real32) :: maxarea = 400._real32
     real(real32) :: vec  = 5._real32 / 100._real32
     real(real32) :: ang  = 1._real32 * pi / 180._real32
     real(real32) :: area = 10._real32 / 100._real32
     real(real32) :: ang_weight  = 10._real32
     real(real32) :: area_weight = 100._real32
  end type tol_type

  type :: abstract_artemis_generator_type
     integer :: num_structures = 0
     integer :: max_num_structures = 100
     
     integer :: axis = 3
     !! Axis along which to align the slab/interface normal vector

     real(real32) :: vacuum_gap = 14._real32
     !! Vacuum thickness in Å

     type(basis_type), dimension(:), allocatable :: structures
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
    implicit none
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
    implicit none
    class(latmatch_type), intent(inout) :: this
    type(tol_type), intent(in) :: tol
    integer, intent(in) :: max_num_matches
    real(real32), dimension(3,3), intent(in) :: lattice_lw,lattice_up
    logical, intent(in) :: reduce_matches

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
    implicit none
    class(latmatch_type), intent(inout) :: this
    integer, dimension(3), intent(in) :: miller_lw, miller_up
    integer, intent(in) :: verbose


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
