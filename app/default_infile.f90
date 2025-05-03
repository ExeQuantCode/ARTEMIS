module infile_print
  !!  This module contains a subroutine to print a default input file
  implicit none


  private

  public :: print_default_file



contains

!###############################################################################
  subroutine print_default_file(file)
    !! Print a default input file for the program
    implicit none

    ! Arguments
    character(*), intent(in), optional :: file
    !! The name of the file to print to (default = stdout)

    ! Local variables
    integer :: unit, status, i
    !! unit number, status, and loop counter
    logical :: exist
    !! logical variable to check if file exists
    character(len=16) :: buffer
    !! buffer for user input


    ! Check if file is present
    ! If not, use stdout
    unit = 6
    if(present(file))then
       ! check if file exists
       inquire(file=file,exist=exist)
       i = 0
       file_overwrite_check: do while(exist)
         i = i + 1
         if(i.gt. 10) then
            write(0,*) "Too many attempts to overwrite file. Exiting."
            return
         end if
         ! file exists, ask if overwrite
         if(i.eq.1) write(*,'("File ",A," already exists. ")',advance='no') trim(adjustl(file))
         write(*,'("Overwrite? (y/n) ")',advance='no')
         read(*,'(A)',iostat=status) buffer
         if(status .ne. 0) return
         buffer = trim(adjustl(buffer))
         select case(buffer(1:1))
         case('y','Y')
            ! overwrite
            write(*,'(" Overwriting file ",A)') trim(adjustl(file))
            exit file_overwrite_check
         case('n','N')
            ! do not overwrite, exit
            write(*,'(" Exiting without overwriting file ",A)') trim(adjustl(file))
            return
         case default
            ! invalid input, ask again
            write(0,'(" Invalid input. Please enter ''y'' or ''n''.")')
         end select
       end do file_overwrite_check
       open(newunit=unit,file=file,action='write')
    end if

    ! Print the default input file
    write(unit,'("SETTINGS")')
    write(unit,'(2X,"TASK        = 1")')
    write(unit,'(2X,"RESTART     = 0")')
    write(unit,'(2X,"STRUC1_FILE = POSCAR1  ! lower structure/interface structure")')
    write(unit,'(2X,"STRUC2_FILE = POSCAR2  ! upper structure (not used if RESTART > 0)")')
    write(unit,'(2X,"MASTER_DIR  = DINTERFACES")')
    write(unit,'(2X,"SUBDIR_PREFIX = D")')
    write(unit,'(2X,"IPRINT = 0")')
    write(unit,'(2X,"CLOCK =               ! taken from the time clock by default")')
    write(unit,'("END SETTINGS")')
    write(unit,*)
    write(unit,*)
    write(unit,'("CELL_EDITS")')
    write(unit,'(2X,"LSURF_GEN   = T")')
    write(unit,'(2X,"MILLER_PLANE  = 1 2 1")')
    write(unit,'(2X,"SLAB_THICKNESS = 6")')
    write(unit,'("END CELL_EDITS")')
    write(unit,*)
    write(unit,*)
    write(unit,'("INTERFACES")')
    write(unit,'(2X,"LGEN_INTERFACES = T   ! generate interfaces")')
    write(unit,'(2X,"IMATCH =  0           ! interface matching method")')
    write(unit,'(2X,"NINTF = 100           ! max number of interfaces")')
    write(unit,'(2X,"NMATCH = 5            ! max number of lattice matches")')
    write(unit,'(2X,"TOL_VEC = 5.D0        ! max vector tolerance (in percent %)")')
    write(unit,'(2X,"TOL_ANG = 1.D0        ! max angle tolerance (in degrees (°))")')
    write(unit,'(2X,"TOL_AREA = 10.D0      ! max area tolerance (in percent %)")')
    write(unit,'(2X,"TOL_MAXFIND = 100     ! max number of good fits to find per plane")')
    write(unit,'(2X,"TOL_MAXSIZE = 10      ! max increase of any lattice vector")')
    write(unit,'(2X,"LW_USE_PRICEL = T     ! extract and use the primitive cell of lower")')
    write(unit,'(2X,"UP_USE_PRICEL = T     ! extract and use the primitive cell of upper")')
    write(unit,*)
    write(unit,'(2X,"NMILLER = 10          ! number of Miller planes to consider")')
    write(unit,'(2X,"LW_MILLER =           ! written as a miller plane, e.g. 0 0 1")')
    write(unit,'(2X,"UP_MILLER =           ! written as a miller plane, e.g. 0 0 1")')
    write(unit,*)
    write(unit,'(2X,"LW_MIN_THICKNESS = 10 ! thickness of lower material (in Angstrom)")')
    write(unit,'(2X,"UP_MIN_THICKNESS = 10 ! thickness of upper material (in Angstrom)")')
    write(unit,'(2X,"NTERM = 5             ! max number of terminations per material per match")')
    write(unit,'(2X,"LW_SURFACE =          ! surface to force for interface generation")')
    write(unit,'(2X,"UP_SURFACE =          ! surface to force for interface generation")')
    write(unit,*)
    write(unit,'(2X,"SHIFTDIR =  DSHIFT    ! shift directory name")')
    write(unit,'(2X,"ISHIFT = 4            ! shifting method")')
    write(unit,'(2X,"NSHIFT = 5            ! number of shifts to apply")')
    write(unit,'(2X,"C_SCALE = 1.D0        ! interface-separation scaling factor")')
    write(unit,*)
    write(unit,'(2X,"SWAPDIR =  DSWAP      ! swap directory name")')
    write(unit,'(2X,"ISWAP = 0             ! swapping method")')
    write(unit,'(2X,"NSWAP = 5             ! number of swap structures generated per interface")')
    write(unit,'(2X,"SWAP_DENSITY = 5.D-2  ! intermixing area density")')
    write(unit,*)
    write(unit,'(2X,"LSURF_GEN      = F      ! generate surfaces of a plane")')
    write(unit,'(2X,"LPRINT_TERMS   = F      ! prints all found terminations")')
    write(unit,'(2X,"LPRINT_MATCHES = F    ! prints all found lattice matches")')
    write(unit,'("END INTERFACES")')
    write(unit,*)
    !write(unit,*)
    !write(unit,'("DEFECTS")')
    !write(unit,'("! NOT CURRENTLY IMPLEMENTED")')
    !write(unit,'("END DEFECTS")')


    if(present(file)) close(unit)

  end subroutine print_default_file
!###############################################################################

end module infile_print
