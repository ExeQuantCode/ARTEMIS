module artemis__misc_types
  !! Module containing custom derived types for ARTEMIS
  use artemis__constants, only: real32, pi
  use artemis__misc, only: to_lower
  use artemis__geom_rw, only: basis_type, geom_write
  use artemis__geom_utils, only: MATNORM
  implicit none


  private

  public :: latmatch_type
  public :: tol_type
  public :: abstract_artemis_generator_type


  type latmatch_type
     integer :: nfit
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
     integer :: nstore = 5
     integer :: maxfit = 100
     integer :: maxsize = 10
     real(real32) :: maxlen=20._real32
     real(real32) :: maxarea=400._real32
     real(real32) :: vec = 5._real32 / 100._real32
     real(real32) :: ang = 1._real32 * pi / 180._real32
     real(real32) :: area = 10._real32 / 100._real32
     real(real32) :: ang_weight = 10._real32
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
  subroutine latmatch_init( &
       this, tol, lattice_lw, lattice_up, reduce_matches &
  )
    implicit none
    class(latmatch_type), intent(inout) :: this
    type(tol_type), intent(in) :: tol
    real(real32), dimension(3,3), intent(in) :: lattice_lw,lattice_up
    logical, intent(in) :: reduce_matches

    allocate(this%tf1(tol%nstore,3,3))
    allocate(this%tf2(tol%nstore,3,3))
    allocate(this%tol(tol%nstore,3))

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
    !! Instance of the raffle generator.
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
    !! Instance of the raffle generator.
    type(basis_type), dimension(:), allocatable :: structures
    !! Generated structures.

    this%structures = structures
    this%num_structures = size(structures)
  end subroutine set_structures
!###############################################################################

end module artemis__misc_types
