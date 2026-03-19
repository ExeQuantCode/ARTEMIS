module artemis__io_utils_extd
  use artemis__misc, only: to_upper
  use artemis__io_utils, only: err_abort

  private

  public :: err_abort_print_struc
  public :: setup_input_fmt, setup_output_fmt



contains

!###############################################################################
  subroutine err_abort_print_struc(basis,filename,msg,lstop)
    !! Print structure to file and stops
    use artemis__geom_rw, only: basis_type, geom_write
    implicit none

    ! Arguments
    type(basis_type), intent(in) :: basis
    !! Structure to print
    character(len=*), intent(in) :: filename
    !! File name to print to
    character(len=*), intent(in) :: msg
    !! Message to print
    logical, intent(in), optional :: lstop
    !! Boolean whether to stop or not

    ! Local variables
    integer :: unit
    !! File unit

    
    open(newunit=unit,file=filename)
    call geom_write(unit, basis)
    close(unit)
    if(msg.ne.'') write(0,'(A)') trim(msg)
    if(present(lstop))then
       if(lstop) stop
    else
       stop
    end if

  end subroutine err_abort_print_struc
!###############################################################################
  

!###############################################################################
  subroutine setup_input_fmt(fmt)
    !! Set the structure file input format for the program
    use artemis__geom_rw, only : igeom_input
    implicit none

    ! Arguments
    character(len=*), intent(in) :: fmt
    !! Format of the input file

    ! Local variables
    character(len=:), allocatable :: form
    !! Formatted string for the input file


    allocate(character(len=len(trim(adjustl(fmt)))) ::  form)
    form = trim(adjustl(to_upper(fmt)))
    
    select case(form)
    case("VASP")
       write(*,*) "Input files will be VASP formatted"
       igeom_input=1
    case("CASTEP")
       write(*,*) "Input files will be CASTEP formatted"
       igeom_input=2
       !call err_abort('ERROR: ARTEMIS not yet set up for CASTEP')
    case("QE","QUANTUMESPRESSO")
       write(*,*) "Input files will be QuantumEspresso formatted"
       igeom_input=3
       !call err_abort('ERROR: ARTEMIS not yet set up for Quantum Espresso')
    case("CRYSTAL")
       write(*,*) "Input files will be CRYSTAL formatted"
       igeom_input=4
       call err_abort('ERROR: ARTEMIS not yet set up for CRYSTAL')
    end select

  end subroutine setup_input_fmt
!###############################################################################
  

!###############################################################################
  subroutine setup_output_fmt(fmt,out_filename)
    !! Set the structure file input format for the program
    use artemis__geom_rw, only : igeom_output
    implicit none

    ! Arguments
    character(len=*), intent(in) :: fmt
    !! Format of the output file
    character(len=*), intent(inout) :: out_filename
    !! File name to print to

    ! Local variables
    character(len=:), allocatable :: form
    !! Formatted string for the output file
    

    allocate(character(len=len(trim(adjustl(fmt)))) ::  form)
    form = trim(adjustl(to_upper(fmt)))
    
    select case(form)
    case("VASP")
       write(*,*) "Output files will be VASP formatted"
       if(out_filename.eq.'') out_filename="POSCAR"
       igeom_output=1
    case("CASTEP")
       write(*,*) "Output files will be CASTEP formatted"
       if(out_filename.eq.'') out_filename="struc.cell"
       igeom_output=2
       !call err_abort('ERROR: ARTEMIS not yet set up for CASTEP')
    case("QE","QUANTUMESPRESSO")
       write(*,*) "Output files will be QuantumEspresso formatted"
       if(out_filename.eq.'') out_filename="struc.geom"
       igeom_output=3
       !call err_abort('ERROR: ARTEMIS not yet set up for Quantum Espresso')
    case("CRYSTAL")
       write(*,*) "Output files will be CRYSTAL formatted"
       if(out_filename.eq.'') out_filename="INPUT_geom"
       igeom_output=4
       call err_abort('ERROR: ARTEMIS not yet set up for CRYSTAL')
    end select
    
  end subroutine setup_output_fmt
!###############################################################################

end module artemis__io_utils_extd