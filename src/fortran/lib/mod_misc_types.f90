module artemis__misc_types
  !! Module containing custom derived types for ARTEMIS
  use artemis__constants, only: real32, pi
  use artemis__misc, only: to_lower
  use artemis__geom_rw, only: basis_type, geom_write
  implicit none


  private

  public :: latmatch_type
  public :: tol_type
  public :: abstract_artemis_generator_type


  type latmatch_type
     integer :: nfit
     logical :: lreduced
     character(1) :: abc(3)= [ 'a', 'b', 'c' ]

     integer, dimension(2) :: axes
     integer, allocatable, dimension(:,:,:) :: tf1,tf2
     real(real32), allocatable, dimension(:,:) :: tol
     real(real32), dimension(3,3) :: lat1,lat2
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

     real(real32) :: tol_cart
     real(real32), dimension(3) :: tol_crys

     type(basis_type), dimension(:), allocatable :: structures
   contains
     procedure, pass(this) :: write_structures
     procedure, pass(this) :: get_structures
     procedure, pass(this) :: set_structures
  end type abstract_artemis_generator_type


contains
  
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
