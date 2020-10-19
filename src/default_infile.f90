!!!#############################################################################
!!! module to write example input file
!!!#############################################################################
module infile_print
  implicit none


!!!updated 2020/02/26


contains
!!!#############################################################################
!!! print example.in
!!!#############################################################################
  subroutine print_default_file(file)
    implicit none
    integer :: UNIT
    character(*), optional :: file

    UNIT=0
    if(present(file))then
       UNIT=20
       open(unit=UNIT,file=file)
    end if

    write(UNIT,'("SETTINGS")')
    write(UNIT,'(2X,"TASK        = 1")')
    write(UNIT,'(2X,"RESTART     = 0")')
    write(UNIT,'(2X,"STRUC1_FILE = POSCAR1  ! lower structure/interface structure")')
    write(UNIT,'(2X,"STRUC2_FILE = POSCAR2  ! upper structure (not used if IRESTART > 0)")')
    write(UNIT,'(2X,"MASTER_DIR  = DINTERFACES")')
    write(UNIT,'(2X,"SUBDIR_PREFIX = D")')
    write(UNIT,'(2X,"IPRINT = 0")')
    write(UNIT,'(2X,"CLOCK =               ! taken from the time clock by default")')
    write(UNIT,'("END SETTINGS")')
    write(UNIT,*)
    write(UNIT,*)
    write(UNIT,'("CELL_EDITS")')
    write(UNIT,'(2X,"LSURF_GEN   = T")')
    write(UNIT,'(2X,"MILLER_PLANE  = 1 2 1")')
    write(UNIT,'(2X,"SLAB_THICKNESS = 6")')
    write(UNIT,'("END CELL_EDITS")')
    write(UNIT,*)
    write(UNIT,*)
    write(UNIT,'("INTERFACES")')
    write(UNIT,'(2X,"LGEN_INTERFACES = T   ! generate interfaces")')
    write(UNIT,'(2X,"IMATCH =  0           ! interface matching method")')
    write(UNIT,'(2X,"NINTF = 100           ! max number of interfaces")')
    write(UNIT,'(2X,"NMATCH = 5            ! max number of lattice matches")')
    write(UNIT,'(2X,"TOL_VEC = 5.D0        ! max vector tolerance (in percent %)")')
    write(UNIT,'(2X,"TOL_ANG = 1.D0        ! max angle tolerance (in degrees (°))")')
    write(UNIT,'(2X,"TOL_AREA = 10.D0      ! max area tolerance (in percent %)")')
    write(UNIT,'(2X,"TOL_MAXFIND = 100     ! max number of good fits to find per plane")')
    write(UNIT,'(2X,"TOL_MAXSIZE = 10      ! max increase of any lattice vector")')
    write(UNIT,'(2X,"NMILLER = 10          ! number of Miller planes to consider")')
    write(UNIT,'(2X,"LW_MILLER =           ! written as a miller plane, e.g. 0 0 1")')
    write(UNIT,'(2X,"UP_MILLER =           ! written as a miller plane, e.g. 0 0 1")')
    write(UNIT,*)
    write(UNIT,'(2X,"LW_SLAB_THICKNESS = 3 ! thickness of lower material")')
    write(UNIT,'(2X,"UP_SLAB_THICKNESS = 3 ! thickness of upper material")')
    write(UNIT,'(2X,"NTERM = 5             ! max number of terminations per material per match")')
    write(UNIT,'(2X,"LW_SURFACE =          ! surface to force for interface generation")')
    write(UNIT,'(2X,"UP_SURFACE =          ! surface to force for interface generation")')
    write(UNIT,*)
    write(UNIT,'(2X,"SHIFTDIR =  DSHIFT    ! shift directory name")')
    write(UNIT,'(2X,"ISHIFT = 4            ! shifting method")')
    write(UNIT,'(2X,"NSHIFT = 5            ! number of shifts to apply")')
    write(UNIT,'(2X,"C_SCALE = 1.D0        ! interface-separation scaling factor")')
    write(UNIT,*)
    write(UNIT,'(2X,"SWAPDIR =  DSWAP      ! swap directory name")')
    write(UNIT,'(2X,"ISWAP = 0             ! swapping method")')
    write(UNIT,'(2X,"NSWAP = 5             ! number of swap structures generated per interface")')
    write(UNIT,'(2X,"SWAP_DENSITY = 5.D-2  ! intermixing area density")')
    write(UNIT,*)
    write(UNIT,'(2X,"LSURF_GEN      = F      ! generate surfaces of a plane")')
    write(UNIT,'(2X,"LPRINT_TERMS   = F      ! prints all found terminations")')
    write(UNIT,'(2X,"LPRINT_MATCHES = F    ! prints all found lattice matches")')
    write(UNIT,'("END INTERFACES")')
    write(UNIT,*)
    !write(UNIT,*)
    !write(UNIT,'("DEFECTS")')
    !write(UNIT,'("! NOT CURRENTLY IMPLEMENTED")')
    !write(UNIT,'("END DEFECTS")')


    if(UNIT.ne.0) close(UNIT)


  end subroutine print_default_file
!!!#############################################################################


end module infile_print
