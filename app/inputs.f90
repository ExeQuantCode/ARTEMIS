!!!#############################################################################
!!! Module to define all global variables
!!! Code written by:
!!!    Ned Thaddeus Taylor
!!!    Francis Huw Davies
!!!    Isiah Edward Mikel Rudkin
!!! Code part of the ARTEMIS group
!!!#############################################################################
module inputs
  use artemis__constants, only: real32, pi
  use artemis__misc, only: flagmaker,file_check
  use artemis__geom_rw, only: basis_type,geom_read
  use artemis__io_utils, only: &
       artemis__version__, &
       print_warning, print_header, &
       err_abort
  use artemis__io_utils_extd, only: setup_input_fmt, setup_output_fmt
  use aspect, only: aspect_type, edit_structure
  use lat_compare, only: tol_type
  use infile_tools
  use infile_print
  implicit none


  integer :: max_num_matches, max_num_terms, max_num_planes
  !! Maximum number of matches, terminations and Miller planes for matching
  logical :: compensate_normal
  !! Boolean whether to compensate for mismatch strain by adjusting the
  !! interface normal axis
  integer :: match_method, shift_method, swap_method, depth_method
  !! Integer to determine which method to use for matching, shifting and swapping
  integer :: num_shifts
  !! Number of shifts to be generated per termination pair
  real(real32) :: interface_depth, separation_scale, bondlength_cutoff
  !! Interface depth, separation scale, and maximum bondlength considered for
  !! the shifting method
  real(real32), allocatable, dimension(:,:) :: shifts
  !! Array of shifts to be applied to the upper structure in the interface

  integer :: num_swaps
  !! Number of swaps to be generated per shift
  real(real32) :: swap_density, swap_sigma, swap_depth
  !! Swap density, swap sigma and swap depth for the swapping method
  logical :: require_mirror_swaps
  !! Boolean whether to require swaps to be mirrors on each interface

  logical :: reduce_matches
  !! Reduce lattice matches to their smallest cell (UNSTABLE)

  logical :: break_on_fail
  integer :: icheck_term_pair, interface_idx
  integer :: clock, verbose

  real(real32) :: vacuum_gap
  !! Vacuum gap (FOR SURFACE GENERATION ONLY)
  logical :: lortho
  !! Boolean whether to orthogonalise the lattice (FOR SURFACE GENERATION ONLY)
  integer :: max_num_structures
  !!! Maximum number of structures to be generated

  integer :: axis
  !! Integer to determine which axis to use for the interface

  type(tol_type) :: tolerance
  !! Tolerance settings for lattice matchings

  logical :: lw_use_pricel, up_use_pricel
  !! Boolean whether to use the primitive cell of the lower and upper
  logical :: lw_layered, up_layered
  !! Boolean whether the lower and upper structures are layered
  logical :: lw_require_stoich, up_require_stoich
  !! Boolean whether to require terminations of the lower and upper structures
  !! to be stoichiometrically equivalent to their provided structure

  integer :: nout,task,task_defect
  integer :: irestart
  integer :: lw_num_layers,up_num_layers
  real(real32) :: lw_thickness, up_thickness
  real(real32) :: lw_bulk_modulus, up_bulk_modulus
  real(real32) :: layer_sep,lw_layer_sep,up_layer_sep,tol_sym
  character(len=20) :: input_fmt,output_fmt
  character(200) :: struc1_file,struc2_file,out_filename
  character(100) :: dirname,shiftdir,swapdir,subdir_prefix
  logical :: lsurf_gen,lprint_matches,lprint_terms,lgen_interfaces,lprint_shifts
  logical :: lnorm_lat
  logical :: ludef_lw_layered,ludef_up_layered,ludef_axis
  logical :: lpresent_struc2
  type(basis_type) :: struc1_bas,struc2_bas
  type(aspect_type) :: edits
  integer, dimension(2) :: lw_surf,up_surf
  integer, dimension(3) :: lw_mplane,up_mplane
  integer, allocatable, dimension(:) :: seed
  real(real32), dimension(2) :: udef_intf_loc
  real(real32), dimension(3,3) :: struc1_lat,struc2_lat



contains
!!!#############################################################################
  subroutine set_global_vars()
    use mod_help, only: settings_help,cell_edits_help,interface_help
    implicit none
    integer :: GEOMunit,Reason
    integer :: i,j,n
    character(1024) :: buffer,flag,input_file
    logical :: skip,empty,lout_name
    integer, dimension(8) :: date_time_vals


!!!-----------------------------------------------------------------------------
!!! initialises variables
!!!-----------------------------------------------------------------------------
    Reason=0
    irestart=0
    input_file=""
    nout=10
    empty=.false.
    skip=.false.
    output_fmt="VASP"
    struc1_file="POSCAR"
    struc2_file=""
    out_filename=""
    dirname="DINTERFACES"
    shiftdir="DSHIFT"
    swapdir="DSWAP"
    subdir_prefix="D"
    n=1
    clock = 0
    verbose = 0
    allocate(seed(n))
    match_method = 0
    shift_method = 4
    depth_method = 0   !!! SWAP DEFAULT DEPTH METHOD !!!
    interface_depth = 1.5_real32
    layer_sep=1._real32
    lw_layer_sep=0._real32
    up_layer_sep=0._real32
    lortho = .true.
    lsurf_gen=.false.
    up_mplane=(/0,0,0/)
    lw_mplane=(/0,0,0/)
    axis=3
    lw_num_layers=0
    up_num_layers=0
    lw_thickness=-1._real32
    up_thickness=-1._real32
    vacuum_gap=14._real32
    lw_surf=0
    up_surf=0
    separation_scale = 1._real32
    bondlength_cutoff = 4._real32
    max_num_planes = 10
    num_shifts = 5
    max_num_terms = 5
    max_num_structures=100
    max_num_matches=5
    tolerance%maxlen=20._real32
    tolerance%maxarea=400._real32
    tolerance%maxfit=100
    tolerance%maxsize=10
    tolerance%vec=5._real32
    tolerance%ang=1._real32
    tolerance%area=10._real32
    lprint_terms=.false.
    lprint_shifts=.false.
    lprint_matches=.false.
    lgen_interfaces=.true.
    reduce_matches=.false.
    swap_method = 0
    num_swaps = 5
    swap_density = 5.E-2_real32
    swap_sigma = -1.0
    swap_depth = 3.0
    require_mirror_swaps = .true.
    icheck_term_pair=-1
    lw_layered=.false.
    up_layered=.false.
    ludef_lw_layered=.false.
    ludef_up_layered=.false.
    ludef_axis=.false.
    lnorm_lat=.true.
    lw_surf=0
    up_surf=0
    interface_idx=-1
    tol_sym = 1.E-6_real32
    udef_intf_loc = [ -1._real32, -1._real32 ]
    lw_use_pricel=.true.
    up_use_pricel=.true.

    lw_bulk_modulus=0.E0
    up_bulk_modulus=0.E0
    compensate_normal=.true.
    break_on_fail = .true.


!!!-----------------------------------------------------------------------------
!!! Reads flags and assigns to variables
!!!-----------------------------------------------------------------------------
    flagloop: do i=1,iargc()
       empty=.false.
       if (skip) then
          skip=.false.
          cycle flagloop
       end if
       call getarg(i,buffer)
       buffer=trim(buffer)
!!!------------------------------------------------------------------------
!!! FILE AND DIRECTORY FLAGS
!!!------------------------------------------------------------------------
       if(index(buffer,'-f').eq.1)then
          flag="-f"
          call flagmaker(buffer,flag,i,skip,empty)
          if(.not.empty)then
             read(buffer,'(A)') input_file
          else
             write(*,'("ERROR: No input filename supplied, but the flag ''-f'' was used")')
             infilename_do: do j=1,3
                write(*,'("Please supply an input filename:")')
                read(5,'(A)') input_file
                if(trim(input_file).ne.'')then
                   write(*,'("Input filename supplied")')
                   exit infilename_do
                else
                   write(*,'(1X,"Not a valid filename")')
                end if
                if(j.eq.3)then
                   call err_abort('ERROR: No valid input filename supplied\nExiting...',.true.)
                end if
             end do infilename_do
          end if
       elseif(index(buffer,'-i').eq.1)then
          flag="-i"
          call flagmaker(buffer,flag,i,skip,empty)
          if(.not.empty) read(buffer,'(A)') struc1_file
       elseif(index(buffer,'-I').eq.1)then
          flag="-I"
          call flagmaker(buffer,flag,i,skip,empty)
          if(.not.empty) read(buffer,'(A)') struc2_file
       elseif(index(buffer,'-o').eq.1)then
          flag="-o"
          call flagmaker(buffer,flag,i,skip,empty)
          if(.not.empty) read(buffer,'(A)') subdir_prefix
       elseif(index(buffer,'-D').eq.1)then
          flag="-D"
          call flagmaker(buffer,flag,i,skip,empty)
          if(.not.empty) read(buffer,'(A)') dirname
!!!------------------------------------------------------------------------
!!! JOB CONTROL FLAGS
!!!------------------------------------------------------------------------
       elseif(index(buffer,'-s').eq.1)then
          flag="-s"
          call flagmaker(buffer,flag,i,skip,empty)
          if(.not.empty) read(buffer,*) clock
       elseif(index(buffer,'-p').eq.1)then
          flag="-p"
          call flagmaker(buffer,flag,i,skip,empty)
          if(.not.empty) read(buffer,*) nout
!!!------------------------------------------------------------------------
!!! VERBOSE PRINTS
!!!------------------------------------------------------------------------
       elseif(index(buffer,'-d').eq.1)then
          flag="-d"
          call flagmaker(buffer,flag,i,skip,empty)
          if(empty)then
             call print_default_file()
          else
             call print_default_file(buffer)
          end if
          call exit()
       elseif(index(buffer,'-v').eq.1)then
          flag="-v"
          call flagmaker(buffer,flag,i,skip,empty)
          if(.not.empty) read(buffer,*) verbose
       elseif(index(buffer,'--version').eq.1)then
          flag="--version"
          write(*,'(1X,"ARTEMIS version: ",A)') trim(artemis__version__)
          stop
       elseif(index(buffer,'-h').eq.1.or.index(buffer,'--help').eq.1)then
          flag="--help"
          if(index(buffer,'-h').eq.1) flag="-h"
          call flagmaker(buffer,flag,i,skip,empty)
          if(empty)then
             write(*,'("Flags:")')
             write(*,'("-----------------FILE-NAME-FLAGS-----------------")')
             write(*,'(2X,"-f<STR>         : Input file name (Default = (empty)). (ALTERNATIVE TO FLAGS)")')
             write(*,'(2X,"-i<STR>         : Structure file 1 (Default = POSCAR)")')
             write(*,'(2X,"-I<STR>         : Structure file 2 (Default = (empty)")')
             write(*,'(2X,"-D<STR>         : Output directory name for generated structures (Default = DInterfaces)")')
             write(*,'(2X,"-o<STR>         : Subdirectory prefix (Default = D)")')
             write(*,'("--------------------JOB-FLAGS--------------------")')
             write(*,'(2X,"--restart       : Restart job from where left off (NOT YET IMPLEMENTED)")')
             write(*,'(2X,"--gen-surfaces  : Generates the surfaces and labels them (NOT YET IMPLEMENTED)")')
             write(*,'("------------------VERBOSE-FLAGS------------------")')
             write(*,'(2X,"--version       : Prints the version number")')
             write(*,'(2X,"-v<INT>         : Verbose printing type")')
             write(*,'(2X,"-d[STR]         : Print example input file (to file STR if present)")')
             write(*,'(2X,"-h|--help [tag] : Prints the help for flags and tags (describes [tag] if supplied)")')
             write(*,'(2X,"     ""     all : Prints a list of all input file tags")')
             write(*,'(2X,"--search <str>  : Searches the help for tags including the string <str>")')
          else
             write(*,*) 
             call settings_help(6,trim(adjustl(buffer)))
             call cell_edits_help(6,trim(adjustl(buffer)))
             call interface_help(6,trim(adjustl(buffer)))
             write(*,'("======================================")')
          end if
          stop
       elseif(index(buffer,'--search').eq.1)then
          flag="--search"
          write(*,*) 
          call flagmaker(buffer,flag,i,skip,empty)
          call settings_help(6,trim(adjustl(buffer)),search=.true.)
          call cell_edits_help(6,trim(adjustl(buffer)),search=.true.)
          call interface_help(6,trim(adjustl(buffer)),search=.true.)  
          write(*,'("======================================")')        
          stop
       end if
    end do flagloop


!!!-----------------------------------------------------------------------------
!!! print header
!!!-----------------------------------------------------------------------------
    call print_header(6)


!!!-----------------------------------------------------------------------------
!!! print execution date and time
!!!-----------------------------------------------------------------------------
    call date_and_time(values=date_time_vals)
    write(*,'(" executed on ",&
         &I4,".",I2.2,".",I2.2," at ",&
         &I0,":",I0,":",I0)')&
         date_time_vals(1:3),date_time_vals(5:7)


!!!-----------------------------------------------------------------------------
!!! check if input file was specified and read if true
!!!-----------------------------------------------------------------------------
    if(trim(input_file).ne."")then
       call read_input_file(input_file)
    end if
    if(trim(input_fmt).eq."") input_fmt=output_fmt
    call setup_input_fmt(input_fmt)
    if(task.eq.0.and.out_filename.ne.'')then
       lout_name = .true.
    else
       lout_name = .false.
    end if
    call setup_output_fmt(output_fmt,out_filename)
    if(lout_name)then
       if(out_filename.eq.struc1_file) out_filename=trim(out_filename)//"_out"
    end if



!!!-----------------------------------------------------------------------------
!!! readjust interface tolerances
!!!-----------------------------------------------------------------------------
    tolerance%vec=tolerance%vec/100._real32
    tolerance%ang=tolerance%ang*pi/180._real32
    tolerance%area=tolerance%area/100._real32


!!!-----------------------------------------------------------------------------
!!! set up random seed
!!!-----------------------------------------------------------------------------
    call random_seed(size=n)
    if(clock.eq.0) CALL SYSTEM_CLOCK(COUNT=clock)
    seed = clock + 37 * (/ (i - 1, i = 1, n) /)
    call random_seed(put=seed)

    write(*,'(1X,A,I0)') "clock seed: ",clock


!!!-----------------------------------------------------------------------------
!!! uses module to read vasp POSCAR 1 structure file
!!!-----------------------------------------------------------------------------
    GEOMunit=10
    call file_check(GEOMunit,struc1_file)
    call geom_read(GEOMunit,struc1_bas,4)
    close(GEOMunit)
    lpresent_struc2 = .false.
    !!--------------------------------------------------------------------------
    !! determines whether struc2_file is necessary
    !!--------------------------------------------------------------------------
    if( (irestart.eq.1.and.task.eq.1).or.&
         (lsurf_gen.and.task.eq.1.and.struc2_file.eq.'').or.&
         (task.eq.0.and.struc2_file.eq.'') )then
       write(*,'("2nd structure file not supplied")')
       write(*,'(2X,"As is not necessary for this run, skipping...")')
    elseif(struc2_file.eq.'')then
       call err_abort('ERROR: 2nd structure file not supplied\n&
            &  Supply a filename to the tag STRUC2_FILE in the SETTINGS card\n&
            &Exiting...',.true.)
    else
       !!-----------------------------------------------------------------------
       !! uses module to read vasp POSCAR 2 structure file
       !!-----------------------------------------------------------------------
       lpresent_struc2 = .true.
       GEOMunit=11
       call file_check(GEOMunit,struc2_file)
       call geom_read(GEOMunit,struc2_bas,4)
       close(GEOMunit)
    end if


!!!-----------------------------------------------------------------------------
!!! changes interface depth depending on IDEPTH method
!!!-----------------------------------------------------------------------------
    if(depth_method.eq.0) interface_depth=0._real32



!!!-----------------------------------------------------------------------------
!!! make the output directory
!!!-----------------------------------------------------------------------------
    if(task.ne.0.and..not.lsurf_gen)then
       call system('mkdir -p ' //adjustl(trim(dirname)))
       call write_settings(adjustl(trim(dirname)))
    end if

    write(*,'(A)') repeat("#",50)

    if(lw_thickness.gt.0._real32.and.lw_num_layers.gt.0)then
       write(0,'(1X,A)') "WARNING: SLAB THICKNESS AND NUMBER OF LAYERS BOTH DEFINED"
       write(0,'(1X,A)') "         SLAB THICKNESS OVERRIDES NUMBER OF LAYERS"
       lw_num_layers=0
    elseif(lw_thickness.le.0._real32.and.lw_num_layers.le.0)then
       lw_thickness = 10._real32
    end if
    if(up_thickness.gt.0._real32.and.up_num_layers.gt.0)then
       write(0,'(1X,A)') "WARNING: SLAB THICKNESS AND NUMBER OF LAYERS BOTH DEFINED"
       write(0,'(1X,A)') "         SLAB THICKNESS OVERRIDES NUMBER OF LAYERS"
       up_num_layers=0
    elseif(up_thickness.le.0._real32.and.up_num_layers.le.0)then
       up_thickness = 10._real32
    end if

    return
  end subroutine set_global_vars
!!!#############################################################################


!!!#############################################################################
!!! read input file to get variables
!!!#############################################################################
  subroutine read_input_file(file_name)
    implicit none
    integer :: Reason,unit,count
    character(*) :: file_name
    logical :: lskip
    integer, dimension(27)::readvar
    character(1024) :: buffer,cardname
    

    unit=20
    call file_check(unit,file_name)

    count=0
    readvar=0
    inread: do
       count=count+1
       read(unit,'(A)',iostat=Reason) buffer
       if(Reason.ne.0) exit
       call rm_comments(buffer,count)
       if(trim(buffer).eq.'') cycle inread
       !!---------------------------------------
       !! assignI for integers
       !! assignD for doubles
       !! assignS for strings
       !! assignL for logicals
       !!---------------------------------------
       lskip=.false.
       cardname=trim(adjustl(buffer))
       select case(cardname)
       case("SETTINGS")
          call read_card_settings(unit,count)
       case("CELL_EDITS")  !ASPECT
          if(task.ne.0) lskip=.true.
          call read_card_cell_edits(unit,count,lskip)
       case("INTERFACES")  !ARTEMIS
          if(task.ne.1) lskip=.true.
          call read_card_interfaces(unit,count,lskip)
       case("DEFECTS")     !ARTIE
          if(task.ne.2) lskip=.true.
          call read_card_defects(unit,count,lskip)
          !!DOPER AND SUB_AND_ROTATION

       case default
          write(0,'("NOTE: unable to associate statement on line ",I0," with a card")') count
       end select
    end do inread
    close(unit)


    return
  end subroutine read_input_file
!!!#############################################################################


!!!#############################################################################
!!! reads the settings card
!!!#############################################################################
  subroutine read_card_settings(unit,count,skip)
    implicit none
    integer :: Reason
    character(1024) :: buffer,tagname,jobname
    integer, intent(in) :: unit
    integer, intent(inout) :: count
    integer, dimension(11) :: readvar
    logical, optional, intent(in) :: skip


    task=1
    readvar=0
    settings_read: do 
       count=count+1
       read(unit,'(A)',iostat=Reason) buffer
       if(Reason.ne.0) exit
       call rm_comments(buffer,count)
       if(trim(buffer).eq.'') cycle settings_read
       if(index(trim(buffer),"END").ne.0.and.&
            index(trim(buffer),"SETTINGS").ne.0) exit settings_read
       if(present(skip))then
          if(skip) cycle
       end if
       tagname=trim(adjustl(buffer))
       if(scan(buffer,"=").ne.0) tagname=trim(tagname(:scan(tagname,"=")-1))
       select case(trim(tagname))
       case("TASK")
          if(verify(buffer,' 0123456789').ne.0) then
             call assign(buffer,     task, readvar(1))
          else
             call assign(buffer,  jobname, readvar(1))
             if(index(jobname,"cell_edits").ne.0)then
                task=0
             elseif(index(jobname,"interface").ne.0)then
                task=1
             elseif(index(jobname,"defects").ne.0)then
                task=2
             end if
          end if
       case("RESTART")
          call assign(buffer,irestart,     readvar(2))
       case("STRUC1_FILE")
          call assign(buffer,struc1_file,  readvar(3))
       case("STRUC2_FILE")
          call assign(buffer,struc2_file,  readvar(4))
       case("MASTER_DIR")
          call assign(buffer,dirname,      readvar(5))
       case("SUBDIR_PREFIX")
          call assign(buffer,subdir_prefix,readvar(6))
       case("IPRINT")
          call assign(buffer,verbose,      readvar(7))
       case("CLOCK")
          call assign(buffer,clock,        readvar(8))
       case("INPUT_FMT")
          call assign(buffer,input_fmt,    readvar(9))
       case("OUTPUT_FMT")
          call assign(buffer,output_fmt,   readvar(10))
       case("TOL_SYM")
          call assign(buffer,tol_sym,      readvar(11))
       case default
          write(*,'("NOTE: unable to assign variable on line ",I0)') count
       end select
    end do settings_read

    if(any(readvar.gt.1)) then
       write(0,*)
       write(0,'(A43)') '###############################'
       write(0,'(A43)') '##########   ERROR   ##########'
       write(0,'(A43)') '###############################'
       write(0,*)
       write(0,'(A)') ' ---       Error in subroutine "checkINPUT"       ---'
       write(0,'(A)') ' --- ERROR: same KEYWORD appears more than once   ---'
       call exit
    end if

    
    return
  end subroutine read_card_settings
!!!#############################################################################


!!!#############################################################################
!!! reads the cell_edits card
!!!#############################################################################
  subroutine read_card_cell_edits(unit,count,skip)
    implicit none
    integer :: Reason
    character(1024) :: buffer,tagname,store
    integer, intent(in) :: unit
    integer, intent(inout) :: count
    integer, dimension(15) :: readvar
    logical, optional, intent(in) :: skip
    character(len=6), dimension(4) :: &
         tag_list = ["axis  ","loc   ","val   ","bounds"]


    edits%nedits=0

    readvar=0
    cell_edits_read: do 
       count=count+1
       read(unit,'(A)',iostat=Reason) buffer
       if(Reason.ne.0) exit
       call rm_comments(buffer,count)
       if(trim(buffer).eq.'') cycle cell_edits_read
       if(index(trim(buffer),"END").ne.0.and.&
            index(trim(buffer),"CELL_EDITS").ne.0) exit cell_edits_read
       if(present(skip))then
          if(skip) cycle
       end if
       tagname=trim(adjustl(buffer))
       if(scan(buffer,"=").ne.0) tagname=trim(tagname(:scan(tagname,"=")-1))
       if(scan(trim(adjustl(tagname))," ").ne.0) read(tagname,*) tagname
       select case(trim(tagname))
       case("OUTPUT_FILE")
          call assign(buffer,out_filename,   readvar(1))
       case("LSURF_GEN")
          call assign(buffer,lsurf_gen,      readvar(2))
       case("MILLER_PLANE")
          call assign_vec(buffer,lw_mplane,  readvar(3))
       case("NUM_LAYERS", "SLAB_THICKNESS")
          if(index(buffer,"SLAB_THICKNESS").ne.0)then
             write(0,'(1X,A)') "WARNING: SLAB_THICKNESS is deprecated, use NUM_LAYERS instead"
          end if
          call assign(buffer,lw_num_layers,   readvar(4))
       case("SHIFT")
          edits%nedits=edits%nedits+1
          store=buffer(index(buffer,"SHIFT")+len("SHIFT"):)
          if(trim(store).eq.'')then
             call cat(unit=unit,end_string="END",end_string2="SHIFT",&
                  line=count,string=store,rm_cmt=.true.)
          end if
          edits%axis(edits%nedits)= nint( assign_list(store,tag_list,1) )
          edits%val(edits%nedits)= assign_list(store,tag_list,3)
          if(index(store,"bounds").eq.0)then
             readvar(5) = readvar(5) + 1
             edits%list(edits%nedits)=1
          else
             readvar(6) = readvar(6) + 1
             edits%list(edits%nedits)=4
             edits%bounds(edits%nedits,:)=assign_listvec(store,tag_list,4)
          end if
       case("VACUUM")
          edits%nedits=edits%nedits+1
          edits%list(edits%nedits)=2
          store=buffer(index(buffer,"VACUUM")+len("VACUUM"):)
          if(trim(store).eq.'')then
             readvar(7) = readvar(7) + 1
             call cat(unit=unit,end_string="END",end_string2="VACUUM",&
                  line=count,string=store,rm_cmt=.true.)
             edits%axis(edits%nedits)= nint( assign_list(store,tag_list,1) )
             edits%bounds(edits%nedits,1)=assign_list(store,tag_list,2)
             edits%val(edits%nedits)=assign_list(store,tag_list,3)
          else
             call assign(buffer, vacuum_gap,        readvar(7))
          end if
       case("TFMAT")
          readvar(8) = readvar(8) + 1
          edits%nedits=edits%nedits+1
          edits%list(edits%nedits)=3
          store=''
          call cat(unit=unit,end_string="END",end_string2="TFMAT",&
               line=count,string=store,rm_cmt=.true.)
          read(store,*) edits%tfmat(1,:),edits%tfmat(2,:),edits%tfmat(3,:)
       case("LAYER_SEP")
          call assign(buffer,layer_sep,        readvar(9))
       case("LORTHO")
          call assign(buffer,lortho,           readvar(10))
       case("SURFACE")
          call assign(buffer,store,            readvar(11))
          select case(icount(store))
          case(1)
             read(store,*) lw_surf(1)
             lw_surf(2) = lw_surf(1)
          case(2)
             read(store,*) lw_surf
          end select
       case("LNORM_LAT")
          call assign(buffer,lnorm_lat,         readvar(12))
       case("MIN_THICKNESS")
          call assign(buffer,lw_thickness,      readvar(13))
       case("USE_PRICEL")
          call assign(buffer,lw_use_pricel,     readvar(14))
       case("REQUIRE_STOICH")
          call assign(buffer,lw_require_stoich, readvar(15))
       case default
          write(*,'("NOTE: unable to assign variable on line ",I0)') count
       end select
    end do cell_edits_read

    if(any(readvar.gt.1)) then
       write(0,*)
       write(0,'(A43)') '###############################'
       write(0,'(A43)') '##########   ERROR   ##########'
       write(0,'(A43)') '###############################'
       write(0,*)
       write(0,'(A)') ' ---       Error in subroutine "checkINPUT"       ---'
       write(0,'(A)') ' --- ERROR: same KEYWORD appears more than once   ---'
       call exit
    end if


    return
  end subroutine read_card_cell_edits
!!!#############################################################################


!!!#############################################################################
!!! reads the interfaces card
!!!#############################################################################
  subroutine read_card_interfaces(unit,count,skip)
    implicit none
    integer :: Reason,j,iudef_nshift
    character(1024) :: store
    character(1024) :: buffer,tagname
    logical :: ludef_shifts, ludef_lw_layer_sep, ludef_up_layer_sep
    integer, intent(in) :: unit
    integer, intent(inout) :: count
    integer, dimension(59) :: readvar
    logical, optional, intent(in) :: skip


    ludef_shifts=.false.
    ludef_lw_layer_sep=.false.
    ludef_up_layer_sep=.false.
    readvar=0
    interfaces_read: do
       count=count+1
       read(unit,'(A)',iostat=Reason) buffer
       if(Reason.ne.0) exit
       call rm_comments(buffer,count)
       if(trim(buffer).eq.'') cycle interfaces_read
       if(index(trim(buffer),"END").ne.0.and.&
            index(trim(buffer),"INTERFACES").ne.0) exit interfaces_read
       if(present(skip))then
          if(skip) cycle
       end if
       tagname=trim(adjustl(buffer))
       if(scan(buffer,"=").ne.0) tagname=trim(tagname(:scan(tagname,"=")-1))
       select case(trim(tagname))
       case("LSURF_GEN")
          call assign(buffer,lsurf_gen,        readvar(1)) 
       case("AXIS")
          ludef_axis=.true.
          call assign(buffer,axis,             readvar(2))
       case("LW_NUM_LAYERS", "LW_SLAB_THICKNESS")
          if(index(buffer,"LW_SLAB_THICKNESS").ne.0)then
             write(0,'(1X,A)') "WARNING: LW_SLAB_THICKNESS is deprecated, use LW_NUM_LAYERS instead"
          end if
          call assign(buffer,lw_num_layers,     readvar(3))
       case("UP_NUM_LAYERS", "UP_SLAB_THICKNESS")
          if(index(buffer,"UP_SLAB_THICKNESS").ne.0)then
             write(0,'(1X,A)') "WARNING: UP_SLAB_THICKNESS is deprecated, use UP_NUM_LAYERS instead"
          end if
          call assign(buffer,up_num_layers,     readvar(4))
       case("LW_MILLER")
          call assign_vec(buffer,lw_mplane,    readvar(5))
       case("UP_MILLER")
          call assign_vec(buffer,up_mplane,    readvar(6))
       case("LW_SURFACE")
          call assign(buffer,store,            readvar(7))
          select case(icount(store))
          case(1)
             read(store,*) lw_surf(1)
             lw_surf(2) = lw_surf(1)
          case(2)
             read(store,*) lw_surf
          end select
       case("UP_SURFACE")
          call assign(buffer,store,            readvar(8))
          select case(icount(store))
          case(1)
             read(store,*) up_surf(1)
             up_surf(2) = up_surf(1)
          case(2)
             read(store,*) up_surf
          end select
       case("SHIFT")
          ludef_shifts=.true.
          iudef_nshift=0
          store=''
          store=buffer(index(buffer,"SHIFT")+len("SHIFT"):)
          if(trim(store).eq.'')then
             readvar(9) = readvar(9) + 1
             call cat(unit=unit,end_string="END",end_string2="SHIFT",&
                  line=iudef_nshift,string=store,rm_cmt=.true.)
             count=count+iudef_nshift
             iudef_nshift=iudef_nshift-1 !removes counting of ENDSHIFT line
             allocate(shifts(iudef_nshift,3))
             read(store,*) (shifts(j,:3),j=1,iudef_nshift)
          else
             call assign(buffer,store,         readvar(9))
             allocate(shifts(1,3))
             select case(icount(store))
             case(1)
                shifts(1,:)=0._real32
                read(store,*) shifts(1,3)
                iudef_nshift = 1
             case(3)
                read(store,*) shifts(1,:)
                if(all(shifts.ge.0._real32)) iudef_nshift=1
             case default
                call err_abort('ERROR: Invalid number of arguments provided to SHIFT&
                     &\nValid number of arguments is 1 or 3.&')
             end select
          end if
       case("NSHIFT")
          call assign(buffer,num_shifts,         readvar(10))
       case("NTERM")
          call assign(buffer,max_num_terms,      readvar(11))
       case("NMATCH")
          call assign(buffer,max_num_matches,    readvar(12))
       case("TOL_VEC")
          call assign(buffer,tolerance%vec,      readvar(13))
       case("TOL_ANG")
          call assign(buffer,tolerance%ang,      readvar(14))
       case("TOL_AREA")
          call assign(buffer,tolerance%area,     readvar(15))
       case("TOL_MAXFIND")
          call assign(buffer,tolerance%maxfit,   readvar(16))
       case("TOL_MAXSIZE")
          call assign(buffer,tolerance%maxsize,  readvar(17))
       case("LPRINT_MATCHES")
          call assign(buffer,lprint_matches,     readvar(18))
       case("LPRINT_TERMS")
          call assign(buffer,lprint_terms,       readvar(19))
       case("LGEN_INTERFACES")
          call assign(buffer,lgen_interfaces,    readvar(20))
       case("IMATCH")
          call assign(buffer,match_method,       readvar(21))
       case("ISHIFT")
          call assign(buffer,shift_method,       readvar(22))
       case("LREDUCE")
          call assign(buffer,reduce_matches,     readvar(23))
       case("LPRINT_SHIFTS")
          call assign(buffer,lprint_shifts,      readvar(24))
       case("C_SCALE")
          call assign(buffer,separation_scale,   readvar(25))
       case("INTF_DEPTH")
          call assign(buffer,interface_depth,    readvar(26))
          depth_method=0
       case("IDEPTH")
          call assign(buffer,depth_method,       readvar(27))
       case("NINTF")
          call assign(buffer,max_num_structures, readvar(28))
       case("ISWAP")
          call assign(buffer,swap_method,        readvar(29))
       case("NSWAP")
          call assign(buffer,num_swaps,          readvar(30))
       case("SWAP_DENSITY")
          call assign(buffer,swap_density,       readvar(31))
       case("SHIFTDIR")
          call assign(buffer,shiftdir,           readvar(32))
       case("SWAPDIR")
          call assign(buffer,swapdir,            readvar(33))
       case("ICHECK")
          call assign(buffer,icheck_term_pair,   readvar(34))
       case("NMILLER")
          call assign(buffer,max_num_planes,     readvar(35))
       case("MAXLEN")
          call assign(buffer,tolerance%maxlen,   readvar(36))
       case("MAXAREA")
          call assign(buffer,tolerance%maxarea,  readvar(37))
       case("LW_LAYERED")
          call assign(buffer,lw_layered,         readvar(38))
          ludef_lw_layered=.true.
       case("UP_LAYERED")
          call assign(buffer,up_layered,         readvar(39))
          ludef_up_layered=.true.
       case("IINTF")
          call assign(buffer,interface_idx,      readvar(40))
       case("LAYER_SEP")
          call assign(buffer,layer_sep,          readvar(41))
       case("LW_LAYER_SEP")
          call assign(buffer,lw_layer_sep,       readvar(42))
          ludef_lw_layer_sep=.true.
       case("UP_LAYER_SEP")
          call assign(buffer,up_layer_sep,       readvar(43))
          ludef_up_layer_sep=.true.
       case("MBOND_MAXLEN")
          call assign(buffer,bondlength_cutoff,  readvar(44))
       case("SWAP_SIGMA")
          call assign(buffer,swap_sigma,         readvar(45))
       case("SWAP_DEPTH")
          call assign(buffer,swap_depth,         readvar(46))
       case("INTF_LOC")
          call assign_vec(buffer,udef_intf_loc,  readvar(47))
       case("LMIRROR")
          call assign(buffer,require_mirror_swaps, readvar(48))
       case("LORTHO")
          call assign(buffer,lortho,             readvar(49))
       case("LW_USE_PRICEL")
          call assign(buffer,lw_use_pricel,      readvar(50))
       case("UP_USE_PRICEL")
          call assign(buffer,up_use_pricel,      readvar(51))
       case("LW_BULK_MODULUS")
          call assign(buffer,lw_bulk_modulus,    readvar(52))
       case("UP_BULK_MODULUS")
          call assign(buffer,up_bulk_modulus,    readvar(53))
       case("LC_FIX")
          call assign(buffer,compensate_normal,  readvar(54))
       case("LBREAK_ON_NO_TERM")
          call assign(buffer,break_on_fail,      readvar(55))
       case("LW_MIN_THICKNESS")
          call assign(buffer,lw_thickness,       readvar(56))
       case("UP_MIN_THICKNESS")
          call assign(buffer,up_thickness,       readvar(57))
       case("LW_REQUIRE_STOICH")
          call assign(buffer,lw_require_stoich,  readvar(58))
       case("UP_REQUIRE_STOICH")
          call assign(buffer,up_require_stoich,  readvar(59))
       case default
          write(0,'("NOTE: unable to assign variable on line ",I0)') count
       end select
    end do interfaces_read


    if(readvar(25).eq.0)then
       select case(shift_method)
       case(0,4)
          separation_scale = 1._real32
       end select
    end if


    if(ludef_shifts)then
       if(readvar(22).eq.1.and.shift_method.ne.0.and.all(shifts.ge.0._real32))then
          write(0,*) "ISHIFT = ",shift_method
          write(0,*) "SHIFT = ",shifts
          call err_abort('ERROR: Contradictory tags used (ISHIFT and SHIFT) &
               &\nNo free shifting directions available&
               &\nExiting...',.true.)
       elseif(readvar(22).eq.1.and.shift_method.ne.0.and.size(shifts(:,1),dim=1).gt.1)then
          call err_abort('ERROR: Contradictory tags used (ISHIFT and SHIFT) &
               &\nExiting...',.true.)
       elseif(all(shifts.ge.0._real32))then
          shift_method=0
          num_shifts=iudef_nshift
       end if
    else
       allocate(shifts(1,3))
       shifts(1,:)=(/-1._real32,-1._real32,-1._real32/)
    end if

    ! set lw_ and up_layer_sep if not defined
    if(.not.ludef_lw_layer_sep) lw_layer_sep = layer_sep
    if(.not.ludef_up_layer_sep) up_layer_sep = layer_sep


    if(any(readvar.gt.1)) then
       write(0,*)
       write(0,'(A43)') '###############################'
       write(0,'(A43)') '##########   ERROR   ##########'
       write(0,'(A43)') '###############################'
       write(0,*)
       write(0,'(A)') ' ---       Error in subroutine "checkINPUT"       ---'
       write(0,'(A)') ' --- ERROR: same KEYWORD appears more than once   ---'
       call exit
    end if


    return
  end subroutine read_card_interfaces
!!!#############################################################################


!!!#############################################################################
!!! reads the defects card
!!!#############################################################################
  subroutine read_card_defects(unit,count,skip)
    implicit none
    integer :: Reason
    character(1024) :: buffer,tagname
    integer, intent(in) :: unit
    integer, intent(inout) :: count
    integer, dimension(22) :: readvar
    logical, optional, intent(in) :: skip


    readvar=0
    defects_read: do
       count=count+1
       read(unit,'(A)',iostat=Reason) buffer
       if(Reason.ne.0) exit
       call rm_comments(buffer,count)
       if(trim(buffer).eq.'') cycle defects_read
       if(index(trim(buffer),"END").ne.0.and.&
            index(trim(buffer),"DEFECTS").ne.0) exit defects_read
       if(present(skip))then
          if(skip) cycle
       end if
       tagname=trim(adjustl(buffer))
       if(scan(buffer,"=").ne.0) tagname=trim(tagname(:scan(tagname,"=")-1))
      select case(trim(tagname))
      case("DEFECT_TASK")
         call assign(buffer,task_defect,readvar(1))
         !! defect task 1 = doper
         !! defect task 2 = molec rotater
      case default
         write(*,'("NOTE: unable to assign variable on line ",I0)') count
      end select
    end do defects_read

    if(any(readvar.gt.1)) then
       write(0,*)
       write(0,'(A43)') '###############################'
       write(0,'(A43)') '##########   ERROR   ##########'
       write(0,'(A43)') '###############################'
       write(0,*)
       write(0,'(A)') ' ---       Error in subroutine "checkINPUT"       ---'
       write(0,'(A)') ' --- ERROR: same KEYWORD appears more than once   ---'
       call exit
    end if


    return
  end subroutine read_card_defects
!!!#############################################################################


!!!############################################################################
!!! write settings file to the output folder
!!!############################################################################
  subroutine write_settings(dirname)
    implicit none
    integer :: UNIT
    character(*) :: dirname
    character(100) :: filename


    filename="settings.txt"
    UNIT=50
    open(unit,file=trim(dirname)//"/"//trim(filename))
!!! WRITE OUT THE SETTINGS CARD IN FULL !!!
    write(UNIT,'("SETTINGS")')
    write(UNIT,'(2X,"TASK    = ",I0)') task
    write(UNIT,'(2X,"RESTART = ",I0)') irestart
    write(UNIT,'(2X,"CLOCK   = ",I0)') clock
    write(UNIT,'("END SETTINGS")')
    write(UNIT,*)
    write(UNIT,*)

!!! WRITE OUT THE CELL_EDITS CARD IN FULL !!!
    if(task.eq.0)then
       write(UNIT,'("CELL_EDITS")')
       write(UNIT,'("END CELL_EDITS")')
       write(UNIT,*)
       write(UNIT,*)
!!! WRITE OUT THE INTERFACES CARD IN FULL !!!
    elseif(task.eq.1)then
       write(UNIT,'("INTERFACES")')
       write(UNIT,'(2X,"LGEN_INTERFACES = ",L)') lgen_interfaces
       write(UNIT,'(2X,"NINTF = ",I0)') max_num_structures
       write(UNIT,'(2X,"IMATCH = ",I0)') match_method
       write(UNIT,'(2X,"NMATCH = ",I0)') max_num_matches
       write(UNIT,'(2X,"TOL_VEC = ",F0.7)') tolerance%vec*100
       write(UNIT,'(2X,"TOL_ANG = ",F0.7)') tolerance%ang*360/(2*pi)
       write(UNIT,'(2X,"TOL_AREA = ",F0.7)') tolerance%area*100
       write(UNIT,*)
       write(UNIT,'(2X,"NMILLER  = ",3(I0,1X))') max_num_planes
       write(UNIT,'(2X,"LW_MILLER_PLANE  = ",3(I0,1X))') lw_mplane
       write(UNIT,'(2X,"UP_MILLER_PLANE  = ",3(I0,1X))') up_mplane
       write(UNIT,'(2X,"LW_SLAB_THICKNESS = ",I0)') lw_num_layers
       write(UNIT,'(2X,"UP_SLAB_THICKNESS = ",I0)') up_num_layers
       if(ludef_lw_layered) write(UNIT,'(2X,"LW_LAYERED = ",L)') lw_layered
       if(ludef_up_layered) write(UNIT,'(2X,"UP_LAYERED = ",L)') lw_layered
       write(UNIT,'(2X,"NTERM = ",I0)') max_num_terms
       write(UNIT,*)
       write(UNIT,'(2X,"ISHIFT = ",I0)') shift_method
       write(UNIT,'(2X,"NSHIFT = ",I0)') num_shifts
       write(UNIT,'(2X,"C_SCALE = ",F0.7)') separation_scale
       write(UNIT,*)
       write(UNIT,'(2X,"ISWAP = ",I0)') swap_method
       write(UNIT,'(2X,"NSWAP = ",I0)') num_swaps
       write(UNIT,'(2X,"SWAP_DENSITY = ",F0.5)') swap_density
       write(UNIT,*)
       write(UNIT,'(2X,"LSURF_GEN   = ",L1)') lsurf_gen
       write(UNIT,'("END INTERFACES")')
       write(UNIT,*)
       write(UNIT,*)
!!! WRITE OUT THE DEFECTS CARD IN FULL !!!
    elseif(task.eq.2)then
       write(UNIT,'("DEFECTS")')
       write(UNIT,'("END DEFECTS")')
       write(UNIT,*)
       write(UNIT,*)
    end if

    close(UNIT)

    return
  end subroutine write_settings
!!!############################################################################

end module inputs
