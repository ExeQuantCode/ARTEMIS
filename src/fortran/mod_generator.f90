!!!#############################################################################
!!! INTERFACES CARD SUBROUTINES
!!! Code written by Ned Thaddeus Taylor and Isiah Edward Mikel Rudkin
!!! Code part of the ARTEMIS group (Hepplestone research group).
!!! Think Hepplestone, think HRG.
!!!#############################################################################
module artemis__generator
  use artemis__constants, only: real32, ierror, pi
  use artemis__misc, only: to_lower,to_upper
  use artemis__geom_rw, only: basis_type,geom_write
  use lat_compare, only: get_best_match,latmatch_type,tol_type
  use artemis__io_utils, only: err_abort
  use artemis__io_utils_extd, only: err_abort_print_struc
  use misc_linalg,          only: uvec,modu,get_area,inverse,cross
  use inputs
  use interface_identifier, only: intf_info_type,&
       get_interface,get_layered_axis,gen_DON
  use edit_geom,            only: planecutter,primitive_lat,ortho_axis,&
       shift_region,set_vacuum,transformer,shifter,reducer,&
       get_min_bulk_bond,get_min_bond,get_shortest_bond,bond_type,&
       share_strain, MATNORM, basis_stack, compare_stoichiometry
  use artemis__sym,              only: confine_type,gldfnd,&
       get_primitive_cell
  use artemis__terminations, only: get_termination_info, term_arr_type, set_slab_height, set_layer_tol, build_slab
  use swapping,              only: rand_swapper
  use shifting !!! CHANGE TO SHIFTER?
  implicit none
  integer, private :: intf=0
  real(real32), private, parameter :: tmp_vac = 14._real32


  type(bulk_DON_type), dimension(2) :: bulk_DON

  type :: abstract_artemis_generator_type
     integer :: max_num_structures = 100

     real(real32) :: tol_cart
     real(real32), dimension(3) :: tol_crys

     type(basis_type), dimension(:), allocatable :: structures
   contains
     procedure, pass(this) :: write_structures
  end type abstract_artemis_generator_type



  type, extends(abstract_artemis_generator_type) :: artemis_termination_generator_type
    
    real(real32) :: layer_separation_cutoff = 1._real32

   contains
     procedure, pass(this) :: generate => generate_terminations
  end type artemis_termination_generator_type


  type, extends(abstract_artemis_generator_type) :: artemis_interface_generator_type
    integer :: match_method = 0
    integer :: max_num_matches = 5
    integer :: max_num_term = 5
    integer :: num_miller_planes = 10
    
    integer :: num_shifts = 5
    integer :: shift_method = 4
    real(real32) :: bondlength_cutoff = 6._real32
    real(real32), dimension(2) :: layer_separation_cutoff = 1._real32

    type(tol_type) :: tolerance

   !  type(basis_type), dimension(:), allocatable :: term_structures_lw
   !  type(basis_type), dimension(:), allocatable :: term_structures_up
   contains
    procedure, pass(this) :: set_tolerance
    procedure, pass(this) :: generate => generate_interfaces
    procedure, pass(this) :: restart => generate_intefaces_from_existing
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

    if(present(tolerance)) this%tolerance = tolerance

    if(present(vector_mismatch)) then
       this%tolerance%vec = vector_mismatch
    else
       this%tolerance%vec = 5._real32
    end if

    if(present(angle_mismatch)) then
       this%tolerance%ang = angle_mismatch
    else
       this%tolerance%ang = 5._real32
    end if

    if(present(area_mismatch)) then
       this%tolerance%area = area_mismatch
    else
       this%tolerance%area = 10._real32
    end if

    if(present(max_length)) then
       this%tolerance%maxlen = max_length
    else
       this%tolerance%maxlen = 20._real32
    end if

    if(present(max_area)) then
       this%tolerance%maxarea = max_area
    else
       this%tolerance%maxarea = 400._real32
    end if

    if(present(max_fit)) then
       this%tolerance%maxfit = max_fit
    else
       this%tolerance%maxfit = 5
    end if

    if(present(max_extension)) then
       this%tolerance%maxsize = max_extension
    else
       this%tolerance%maxsize = 5
    end if

    if(present(angle_weight)) then
       this%tolerance%ang_weight = angle_weight
    else
       this%tolerance%ang_weight = 1._real32
    end if

    if(present(area_weight)) then
       this%tolerance%area_weight = area_weight
    else
       this%tolerance%area_weight = 1._real32
    end if

  end subroutine set_tolerance
!###############################################################################


!###############################################################################
  subroutine generate_terminations( &
       this, basis, miller_plane, axis, num_layers, thickness &
  )
    !! Generate and prints terminations parallel to the supplied miller plane
    implicit none

    ! Arguments
    class(artemis_termination_generator_type), intent(inout) :: this
    !! Instance of artemis generator type
    type(basis_type), intent(in) :: basis
    !! Atomic structure data
    integer, dimension(3), intent(in) :: miller_plane
    !! Miller plane
    integer, intent(in) :: axis
    !! Axis along which to align the slab
    integer, intent(in), optional :: num_layers
    !! Number of layers in the slab
    real(real32), intent(in), optional :: thickness
    !! Thickness of the slab (in Å)

    type(basis_type), dimension(:), allocatable :: output
    !! Output structures

    ! Local variables
    integer :: itmp1, iterm, term_start, term_end, iterm_step, i
    !! Termination loop variables
    integer :: ncells, ntrans
    !! Number of cells in the slab
    integer :: num_structures
    !! Number of structures to be generated
    integer :: num_layers_
    !! Number of layers in the slab
    real(real32) :: height
    !! Height of the slab
    logical :: lcycle
    !! Boolean whether to cycle through the slab
    type(basis_type) :: tmp_bas1,tmp_bas2
    !! Temporary basis structures
    type(confine_type) :: confine
    !! Confine structure along the specified axis
    type(term_arr_type) :: term
    !! List of terminations
    real(real32), dimension(3,3) :: tfmat
    !! Transformation matrix

    character(len=256) :: warn_msg

    integer, allocatable, dimension(:,:,:) :: bas_map,t1bas_map
    real(real32), allocatable, dimension(:,:) :: trans


    !! copy lattice and basis for manipulating
    call tmp_bas1%copy(basis)
    allocate(bas_map(tmp_bas1%nspec,maxval(tmp_bas1%spec(:)%num,dim=1),2))
    bas_map = -1


    write(6,'(1X,"Using supplied plane...")')
    tfmat = planecutter(tmp_bas1%lat,real(miller_plane,real32))
    call transformer(tmp_bas1,tfmat,bas_map)
    !call err_abort_print_struc(bas,"check.vasp","stop")

    !---------------------------------------------------------------------------
    ! Finds smallest thickness of the slab and increases to ...
    ! ... user-defined thickness
    !---------------------------------------------------------------------------
    confine%l = .false.
    confine%axis = axis
    confine%laxis = .false.
    confine%laxis(axis) = .true.
    if(allocated(trans)) deallocate(trans)
    allocate(trans(minval(tmp_bas1%spec(:)%num+2),3))
    call gldfnd(confine, tmp_bas1, tmp_bas1, trans, ntrans)
    tfmat(:,:) = 0._real32
    tfmat(1,1) = 1._real32
    tfmat(2,2) = 1._real32
    if(ntrans.eq.0)then
       tfmat(3,3)=1._real32
    else
       itmp1=minloc(abs(trans(:ntrans,axis)),dim=1,&
            mask=abs(trans(:ntrans,axis)).gt.1.D-3/modu(tmp_bas1%lat(axis,:)))
       tfmat(3,:)=trans(itmp1,:)
    end if
    if(all(abs(tfmat(3,:)).lt.1.E-5_real32)) tfmat(3,3) = 1._real32
    call transformer(tmp_bas1,tfmat,bas_map)
    if(.not.compare_stoichiometry(tmp_bas1,basis))then
       write(0,'(1X,"ERROR: Internal error in generate_terminations")')
       write(0,'(2X,"The gldfnd subroutine could not reproduce a valid primitive cell for the material")')
       if(ierror.eq.1)then
          call err_abort_print_struc(tmp_bas1, "broken_primitive.vasp", &
           "Code exiting due to IPRINT = 1")
       end if
       write(0,'(2X,"Skipping this lattice match")')
       return
    end if

    ! get the terminations
    term = get_termination_info( &
         tmp_bas1, axis, &
         lprint = .true., layer_sep = this%layer_separation_cutoff, &
         break_on_fail = lbreak_on_no_term &
    )
    if(term%nterm .eq. 0)then
       write(warn_msg, '(A,I0,1X,I0,1X,I0,A)') &
            "No terminations found for Miller plane (",miller_plane,")"
       call print_warning(trim(warn_msg))
       return
    end if

    ! set thickness if provided by user
    if(present(num_layers))then
       num_layers_ = num_layers
    else
       num_layers_ = 1
    end if

    ! determine tolerance for layer separations (termination tolerance)
    ! ... this is different from layer_sep
    call set_layer_tol(term)

    ! determine required extension and perform that
    call set_slab_height(tmp_bas1,bas_map,term,lw_surf,&
         height,num_layers_, thickness, ncells,&
         term_start,term_end,iterm_step &
    )
    
    !---------------------------------------------------------------------------
    ! Normalise lattice
    !---------------------------------------------------------------------------
    if(lnorm_lat)then
       call reducer(tmp_bas1)
       tmp_bas1%lat = MATNORM(tmp_bas1%lat)
    end if
    

    !---------------------------------------------------------------------------
    ! loop over terminations and write them
    !---------------------------------------------------------------------------
    num_structures = ( term_end - term_start ) / iterm_step + 1
    allocate(output(num_structures))
    do iterm = term_start, term_end, iterm_step
       i = ( iterm - term_start ) / iterm_step + 1 
       call output(i)%copy(tmp_bas1)
       if(allocated(t1bas_map)) deallocate(t1bas_map)
       allocate(t1bas_map,source=bas_map)
       call build_slab(output(i),bas_map,term,[iterm,lw_surf(2)],&
            thickness, ncells, num_layers_, height,&
            "lw",lcycle,lortho,vacuum &
       )
    end do
    if(.not.allocated(this%structures))then
       call move_alloc(output,this%structures)
    else
       this%structures = [ this%structures, output ]
    end if

   end subroutine generate_terminations
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
  subroutine generate_intefaces_from_existing(this, basis)
    !! Generate interfaces for the given basis
    implicit none

    ! Arguments
    class(artemis_interface_generator_type), intent(inout) :: this
    !! Instance of artemis generator type
    type(basis_type), intent(in) :: basis
    !! Atomic structure data

    ! Local variables
    integer :: is,ia,js,ja
    !! Loop variables
    real(real32) :: dtmp1,min_bond,min_bond1,min_bond2
    !! Minimum bond length
    type(intf_info_type) :: intf
    !! Interface information
    real(real32), dimension(3) :: vtmp1
    !! Temporary vector


    min_bond1=huge(0._real32)
    min_bond2=huge(0._real32)
    if(any(udef_intf_loc.lt.0._real32))then
       if(ludef_axis)then
          intf=get_interface(basis%lat,basis,axis)
       else
          intf=get_interface(basis%lat,basis)
       end if
       intf%loc=intf%loc/modu(basis%lat(intf%axis,:))
       write(6,*) "interface axis:",intf%axis
       write(6,*) "interface loc:",intf%loc
       !! write interface location to a file for user to refer back to
       open(unit=10,file="interface_location.dat")
       write(10,'(1X,"AXIS = ",I0)') intf%axis
       write(10,'(1X,"INTF_LOC = ",2(2X,F9.6))') intf%loc
       close(10)
    else
       intf%axis = axis
       intf%loc = udef_intf_loc
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

    min_bond = ( min_bond1 + min_bond2 )/2._real32
    write(6,'(1X,"Avg min bulk bond: ",F0.3," Å")') min_bond
    write(6,'(1X,"Trans-interfacial scaling factor:",F0.3)') c_scale
    call gen_shifts_and_swaps(basis,intf%axis,intf%loc,min_bond,&
         ishift,nshift,&
         iswap,swap_den,nswap)


  end subroutine generate_intefaces_from_existing
!###############################################################################


!###############################################################################
  subroutine generate_interfaces(this, basis_lw, basis_up)
    !! Generate interfaces from two bulk structures
    implicit none

    ! Arguments
    class(artemis_interface_generator_type), intent(inout) :: this
    !! Instance of artemis generator type
    type(basis_type), intent(in) :: basis_lw
    !! Lower bulk structure
    type(basis_type), intent(in) :: basis_up
    !! Upper bulk structure

    ! Local variables
    real(real32) :: avg_min_bond
    !! Average minimum bond length

    type(basis_type) :: basis_lw_, basis_up_
    !! Temporary basis structures
    type(basis_type) :: supercell_lw, supercell_up
    !! Copy of the basis structures
    type(basis_type) :: slab_lw, slab_up
    !! Slab structures
    type(basis_type) :: interface
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



    integer :: ntrans,iunique,itmp1,old_intf
    integer :: lw_layered_axis,up_layered_axis
    real(real32) :: dtmp1,bondlength
    character(3) :: abc
    character(1024) :: pwd,intf_dir,dirpath,msg, filename
    logical :: ludef_lw_surf,ludef_up_surf,lcycle
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


!!!-----------------------------------------------------------------------------
!!! determines the primitive and niggli reduced cell for each bulk
!!!-----------------------------------------------------------------------------
    call basis_lw_%copy(basis_lw)
    call basis_up_%copy(basis_up)
    write(6,*)
    if(lw_use_pricel)then
       write(6,'(1X,"Using primitive cell for lower material")')
       call get_primitive_cell(basis_lw_)
    else
       write(6,'(1X,"Using supplied cell for lower material")')
       call reducer(basis_lw_)
       basis_lw_%lat=primitive_lat(basis_lw_%lat)
    end if
    if(up_use_pricel)then
       write(6,'(1X,"Using primitive cell for upper material")')
       call get_primitive_cell(basis_up_)
    else
       write(6,'(1X,"Using supplied cell for upper material")')
       call reducer(basis_up_)
       basis_up_%lat=primitive_lat(basis_up_%lat)
    end if
    write(6,*)


    ludef_lw_surf = .false.
    if(all(lw_surf.gt.0)) ludef_lw_surf = .true.
    ludef_up_surf = .false.
    if(all(up_surf.gt.0)) ludef_up_surf = .true.
    

    
!!!-----------------------------------------------------------------------------
!!! investigates individual bulks and their bondlengths
!!!-----------------------------------------------------------------------------
    avg_min_bond = &
         ( get_min_bulk_bond(basis_lw_) + get_min_bulk_bond(basis_up_) )/2._real32
    write(6,'(1X,"Avg min bulk bond: ",F0.3," Å")') avg_min_bond
    write(6,'(1X,"Trans-interfacial scaling factor: ",F0.3)') c_scale
    if(ishift.eq.-1) nshift=1
    

!!!-----------------------------------------------------------------------------
!!! gets bulk DONs, if ISHIFT = 4
!!!-----------------------------------------------------------------------------
    allocate(lw_map(basis_lw_%nspec,maxval(basis_lw_%spec(:)%num,dim=1),2))
    allocate(up_map(basis_up_%nspec,maxval(basis_up_%spec(:)%num,dim=1),2))    
    if(ishift.eq.4.or.ishift.eq.0)then
       lw_map=0
       bulk_DON(1)%spec=gen_DON(basis_lw_%lat,basis_lw_,&
            dist_max=max_bondlength,&
            scale_dist=.false.,&
            norm=.true.)
       do is = 1, inlw_bas%nspec
          if(all(abs(bulk_DON(1)%spec(is)%atom(:,:)).lt.1._real32))then
             bondlength = huge(0._real32)
             do ia = 1, inlw_bas%spec(is)%num
                dtmp1 = modu(get_min_bond(inlw_lat, inlw_bas, is, ia))
                if(dtmp1.lt.bondlength) bondlength = dtmp1
                if(dtmp1.gt.max_bondlength)then
                   write(filename,'("lw_DON_",I0,"_",I0,".dat")') is,ia
                   open(unit=13,file=filename)
                   do j=1,1000
                      write(13,*) &
                           (j-1)*max_bondlength/1000,&
                           bulk_DON(1)%spec(is)%atom(ia,j)
                   end do
                   close(13)
                  end if
             end do
             if(bondlength.gt.max_bondlength)then
                write(0,*) "Min bondlength for lower species ", &
                     is, " is ", bondlength
                write(0,*) "To account for this, increase MBOND_MAXLEN to at &
                     &least ",bondlength
             end if
             call err_abort("ISSUE WITH THE LOWER BULK DON!!!")
          end if
       end do
       up_map=0
       bulk_DON(2)%spec=gen_DON(basis_up_%lat,basis_up_,&
            dist_max=max_bondlength,&
            scale_dist=.false.,&
            norm=.true.)
       do is = 1, inup_bas%nspec
          if(all(abs(bulk_DON(2)%spec(is)%atom(:,:)).lt.1._real32))then
             bondlength = huge(0._real32)
             do ia = 1, inup_bas%spec(is)%num
                dtmp1 = modu(get_min_bond(inup_lat, inup_bas, is, ia))
                if(dtmp1.lt.bondlength) bondlength = dtmp1
                if(dtmp1.gt.max_bondlength)then
                   write(filename,'("up_DON_",I0,"_",I0,".dat")') is,ia
                   open(unit=13,file=filename)
                   do j=1,1000
                      write(13,*) &
                           (j-1)*max_bondlength/1000,&
                           bulk_DON(2)%spec(is)%atom(ia,j)
                   end do
                   close(13)
                  end if
             end do
             if(bondlength.gt.max_bondlength)then
                write(0,*) "Min bondlength for upper species ", &
                     is, " is ", bondlength
                write(0,*) "To account for this, increase MBOND_MAXLEN to at &
                     &least ",bondlength
             end if
             call err_abort("ISSUE WITH THE UPPER BULK DON!!!")
          end if
       end do
    else
       lw_map=-1
       up_map=-1       
    end if


!!!-----------------------------------------------------------------------------
!!! checks whether system appears layered
!!!-----------------------------------------------------------------------------
    lw_layered_axis=get_layered_axis(basis_lw_%lat,basis_lw_)
    if(.not.lw_layered.and.lw_layered_axis.gt.0)then
       ivtmp1=0
       ivtmp1(lw_layered_axis)=1
       if(ludef_lw_layered)then
          write(msg,'("Lower crystal appears layered along axis ",I0,"\n&
               &Partial layer terminations will be generated\n&
               &We suggest using LW_MILLER =",3(1X,I1))') lw_layered_axis,ivtmp1
          call print_warning(trim(msg))
       else
          write(msg,'("Lower crystal has been identified as layered\nalong",3(1X,I1),"\n&
               &Confining crystal to this plane and\nstoichiometric terminations.\n&
               &If you don''t want this, set\nLW_LAYERED = .FALSE.")') &
               ivtmp1
          call print_warning(trim(msg))
          lw_mplane=ivtmp1
          lw_layered=.true.
       end if
    elseif(lw_layered.and.lw_layered_axis.gt.0.and.all(lw_mplane.eq.0))then
       lw_mplane(lw_layered_axis)=1
    end if

    up_layered_axis=get_layered_axis(basis_up_%lat,basis_up_)
    if(.not.up_layered.and.up_layered_axis.gt.0)then
       ivtmp1=0
       ivtmp1(up_layered_axis)=1
       if(ludef_up_layered)then
          write(msg,'("Upper crystal appears layered along axis ",I0,"\n&
               &Partial layer terminations will be generated\n&
               &We suggest using UP_MILLER =",3(1X,I1))') up_layered_axis,ivtmp1
          call print_warning(trim(msg))
       else
          write(msg,'("Upper crystal has been identified as layered\nalong",3(1X,I1),"\n&
               &Confining crystal to this plane and\nstoichiometric terminations.\n&
               &If you don''t want this, set\nUP_LAYERED = .FALSE.")') &
               ivtmp1
          call print_warning(trim(msg))
          up_mplane=ivtmp1
          up_layered=.true.
       end if
    elseif(up_layered.and.up_layered_axis.gt.0.and.all(up_mplane.eq.0))then
       up_mplane(up_layered_axis)=1
    end if


!!!-----------------------------------------------------------------------------
!!! Finds and stores the best matches between the materials
!!!-----------------------------------------------------------------------------
    call getcwd(pwd)
    old_intf = -1
    intf=0
    abc="abc"
    if(imatch.ne.0.and.(any(lw_mplane.ne.0).or.any(up_mplane.ne.0)))then
       call err_abort( '&
            &Cannot use LW_MILLER or UP_MILLER with IMATCH>0\n&
            Exiting...', &
            fmtd=.true. &
       )
    elseif(imatch.ne.0)then
       write(msg,'("&
            &IMATCH /= 0 methods are experimental and may\n&
            &not work as expected.\n&
            &They are not intended to be thorough searches.\n&
            &This method is not recommended unless you\n&
            &are clear on its intended use and\n&
            &limitations.&
       &")')
       call print_warning(trim(msg))
    end if
    if(any(lw_mplane.ne.0))then
       if(imatch.ne.0)then
          abc="ab"
          tfmat=planecutter(basis_lw_%lat,real(lw_mplane,real32))
          call transformer(basis_lw_,tfmat,lw_map)
          SAV=get_best_match(&
               this%tolerance,&
               basis_lw_%lat,basis_up_%lat,&
               basis_lw_,basis_up_,&
               trim(abc),"abc",lprint_matches,ierror,imatch=imatch)
       elseif(any(up_mplane.ne.0))then
          SAV=get_best_match(&
               this%tolerance,&
               basis_lw_%lat,basis_up_%lat,&
               basis_lw_,basis_up_,&
               trim(abc),"abc",lprint_matches,ierror,imatch=imatch,&
               plane1=lw_mplane,plane2=up_mplane,nmiller=nmiller)
       else
          SAV=get_best_match(&
               this%tolerance,&
               basis_lw_%lat,basis_up_%lat,&
               basis_lw_,basis_up_,&
               trim(abc),"abc",lprint_matches,ierror,imatch=imatch,&
               plane1=lw_mplane,nmiller=nmiller)
       end if
    elseif(any(up_mplane.ne.0))then
       SAV=get_best_match(&
            this%tolerance,&
            basis_lw_%lat,basis_up_%lat,&
            basis_lw_,basis_up_,&
            trim(abc),"abc",lprint_matches,ierror,imatch=imatch,&
            plane2=up_mplane,nmiller=nmiller)
    else
       SAV=get_best_match(&
            this%tolerance,&
            basis_lw_%lat,basis_up_%lat,&
            basis_lw_,basis_up_,&
            trim(abc),"abc",lprint_matches,ierror,imatch=imatch,&
            nmiller=nmiller)
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
         nterm*nshift*this%tolerance%nstore
    if(.not.lgen_interfaces)then
       write(0,'(1X,"Told not to generate interfaces, just find matches.")')
       write(0,'("Exiting...")')
       call exit()
    end if

       
!!!-----------------------------------------------------------------------------
!!! Saves current directory and moves to new directory
!!!-----------------------------------------------------------------------------
    call system('mkdir -p '//trim(adjustl(dirname)))
    call chdir(dirname)
    call getcwd(intf_dir)

    if(iintf.gt.0)then
       intf_start=iintf
       intf_end=iintf
       write(6,'(1X,"Generating only interfaces for match ",I0)') iintf
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
       if(ishift.eq.4)then
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
               dist_max=max_bondlength,&
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
       confine%axis=axis
       confine%laxis=.false.
       confine%laxis(axis)=.true.
       if(allocated(trans)) deallocate(trans)
       allocate(trans(minval(supercell_lw%spec(:)%num+2),3))
       call gldfnd(confine,supercell_lw,supercell_lw,trans,ntrans)
       tfmat(:,:)=0._real32
       tfmat(1,1)=1._real32
       tfmat(2,2)=1._real32
       if(ntrans.eq.0)then
          tfmat(3,3)=1._real32
       else
          itmp1=minloc(abs(trans(:ntrans,axis)),dim=1,&
               mask=abs(trans(:ntrans,axis)).gt.1.D-3/modu(supercell_lw%lat(axis,:)))
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
            supercell_lw, axis, &
            lprint = lprint_terms, layer_sep = lw_layer_sep, &
            break_on_fail = lbreak_on_no_term &
       )
       if(lw_term%nterm .eq. 0)then
          write(0,'("WARNING: &
               &No terminations found for lower material Miller plane &
               &(",3(1X,I0)," )")' &
          ) SAV%tf1(ifit,3,1:3)
          cycle intf_loop
       end if
       if(any(lw_surf.gt.lw_term%nterm))then
          write(msg, '("LW_SURFACE VALUES INVALID!\nOne or more value &
               &exceeds the maximum number of terminations in the &
               structure.\n&
               &  Supplied values: ",I0,1X,I0,"\n&
               &  Maximum allowed: ",I0)') lw_surf, lw_term%nterm
          call err_abort(trim(msg),fmtd=.true.)
       end if


       !!-----------------------------------------------------------------------
       !! Sort out ladder rungs (checks whether the material is centrosymmetric)
       !!-----------------------------------------------------------------------
       !call setup_ladder(supercell_lw%lat,supercell_lw,axis,lw_term)
       if(sum(lw_term%arr(:)%natom)*lw_term%nstep.ne.supercell_lw%natom)then
          write(msg, '("ERROR: Number of atoms in lower layers not correct: "&
               &I0,2X,I0)') sum(lw_term%arr(:)%natom)*lw_term%nstep,supercell_lw%natom
          call err_abort(trim(msg),fmtd=.true.)
       end if
       call set_layer_tol(lw_term)


       !!-----------------------------------------------------------------------
       !! Defines height of lower slab from user-defined values
       !!-----------------------------------------------------------------------
       call set_slab_height(supercell_lw,t1lw_map,lw_term,lw_surf,&
            height_lw,lw_num_layers, lw_thickness,ncells_lw,&
            term_lw_start_idx,term_lw_end_idx,term_lw_step &
       )
       if(term_lw_end_idx.gt.this%max_num_term) term_lw_end_idx = this%max_num_term


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
          itmp1=minloc(abs(trans(:ntrans,axis)),dim=1,&
               mask=abs(trans(:ntrans,axis)).gt.1.D-3/modu(supercell_lw%lat(axis,:)))
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
            supercell_up, axis, &
            lprint = lprint_terms, layer_sep = up_layer_sep, &
            break_on_fail = lbreak_on_no_term &
       )
       if(up_term%nterm .eq. 0)then
          write(0,'("WARNING: &
               &No terminations found for upper material Miller plane &
               &(",3(1X,I0)," )")' &
          ) SAV%tf2(ifit,3,1:3)
          cycle intf_loop
       end if
       if(any(up_surf.gt.up_term%nterm))then
          write(msg, '("UP_SURFACE VALUES INVALID!\nOne or more value &
               &exceeds the maximum number of terminations in the &
               structure.\n&
               &  Supplied values: ",I0,1X,I0,"\n&
               &  Maximum allowed: ",I0)') up_surf, up_term%nterm
          call err_abort(trim(msg),fmtd=.true.)
       end if


       !!-----------------------------------------------------------------------
       !! Sort out ladder rungs (checks whether the material is centrosymmetric)
       !!-----------------------------------------------------------------------
       !call setup_ladder(supercell_up%lat,supercell_up,axis,up_term)
       if(sum(up_term%arr(:)%natom)*up_term%nstep.ne.supercell_up%natom)then
          write(msg, '("ERROR: Number of atoms in upper layers not correct: "&
               &I0,2X,I0)') sum(up_term%arr(:)%natom)*up_term%nstep,supercell_up%natom
          call err_abort(trim(msg),fmtd=.true.)
       end if
       call set_layer_tol(up_term)


       !!-----------------------------------------------------------------------
       !! Defines height of upper slab from user-defined values
       !!-----------------------------------------------------------------------
       call set_slab_height(supercell_up,t1up_map,up_term,up_surf,&
            height_up,up_num_layers, up_thickness, ncells_up,&
            term_up_start_idx,term_up_end_idx,term_up_step &
       )
       if(term_up_end_idx.gt.this%max_num_term) term_up_end_idx = this%max_num_term


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
          call build_slab(slab_lw,t2lw_map,lw_term,[iterm_lw,lw_surf(2)],&
          lw_thickness, ncells_lw, lw_num_layers, height_lw,&
               "lw",lcycle)
          if(lcycle) cycle lw_term_loop

          
          !!--------------------------------------------------------------------
          !! Cycles over terminations of upper material
          !!--------------------------------------------------------------------
          up_term_loop: do iterm_up = term_up_start_idx, term_up_end_idx, term_up_step
             call slab_up%copy(supercell_up)
             if(allocated(t2up_map)) deallocate(t2up_map)
             allocate(t2up_map,source=t1up_map)
             call build_slab(slab_up,t2up_map,up_term,[iterm_up,up_surf(2)],&
                  up_thickness, ncells_up, up_num_layers, height_up,&
                  "up",lcycle)
             if(lcycle) cycle up_term_loop

             
             !!-----------------------------------------------------------------
             !! Checks stoichiometry
             !!-----------------------------------------------------------------
             if(slab_lw%nspec.ne.basis_lw_%nspec.or.any(&
                  (basis_lw_%spec(1)%num*slab_lw%spec(:)%num)&
                  /slab_lw%spec(1)%num.ne.basis_lw_%spec(:)%num))then
                write(6,'("WARNING: This lower surface termination is not &
                     &stoichiometric")')
                if(lw_layered)then
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
                if(up_layered)then
                   write(6,'(2X,"As upper structure is layered, stoichiometric &
                        &surfaces are required.")')
                   write(6,'(2X,"Skipping this termination...")')
                   cycle up_term_loop
                end if
             end if


             !!-----------------------------------------------------------------
             !! Use the bulk moduli to determine the strain sharing
             !!-----------------------------------------------------------------
             if(lw_bulk_modulus.ne.0.E0.and.up_bulk_modulus.ne.0.E0)then
                call share_strain(slab_lw%lat,slab_up%lat,&
                     lw_bulk_modulus,up_bulk_modulus,lcompensate=.not.lc_fix)
             end if
             

             !!-----------------------------------------------------------------
             !! Merge the two bases and lattices and define the interface loc
             !!-----------------------------------------------------------------
             interface = basis_stack(&
                  basis1 = slab_lw, basis2 = slab_up, &
                  axis = axis, offset = init_offset(:), &
                  map1 = t2lw_map, map2 = t2up_map &
             )
             intf_loc(1) = ( modu(slab_lw%lat(axis,:)) + 0.5_real32*init_offset(axis) - &
                  tmp_vac)/modu(interface%lat(axis,:))
             intf_loc(2) = ( modu(slab_lw%lat(axis,:)) + modu(slab_up%lat(axis,:)) + &
                  1.5_real32*init_offset(axis) - 2._real32*tmp_vac )/modu(interface%lat(axis,:))
             if(ierror.ge.1)then
                write(0,*) "interface:",intf_loc
                if(ierror.eq.1.and.iunique.eq.icheck_intf-1)then
                   call chdir(intf_dir)
                   call err_abort_print_struc(slab_lw,"lw_term.vasp",&
                        "",.false.)
                   call err_abort_print_struc(slab_up,"up_term.vasp",&
                        "As IPRINT = 1 and ICHECK has been set, &
                        &code is now exiting...")
                elseif(ierror.eq.2.and.iunique.eq.icheck_intf-1)then
                   call chdir(intf_dir)
                   call err_abort_print_struc(interface,"test_intf.vasp",&
                        "As IPRINT = 2 and ICHECK has been set, &
                        &code is now exiting...")
                end if
             end if


             !!-----------------------------------------------------------------
             !! Saves current directory and moves to new directory
             !!-----------------------------------------------------------------
             if(intf.gt.old_intf)then
                iunique=iunique+1
                if(ishift.gt.0.and.nshift.gt.1) &
                     write(6,'(1X,"Generating shifts for unique interface ",&
                     &I0,":")') iunique
                write(dirpath,'(A,I0.2)') trim(adjustl(subdir_prefix)),iunique
                call system('mkdir -p '//trim(adjustl(dirpath)))
             else
                write(dirpath,'(A,I0.2)') trim(adjustl(subdir_prefix)),iunique
             end if
             call chdir(dirpath)
             old_intf = intf

             
             !!-----------------------------------------------------------------
             !! Writes information of current match to file in save directory
             !!-----------------------------------------------------------------
             call  output_intf_data(SAV, ifit, lw_term, iterm_lw, up_term, iterm_up,&
                  lw_use_pricel,up_use_pricel)


             !!-----------------------------------------------------------------
             !! Generates shifts and swaps and prints the subsequent structures
             !!-----------------------------------------------------------------
             call gen_shifts_and_swaps(interface,axis,intf_loc,avg_min_bond,&
                  ishift,nshift,&
                  iswap,swap_den,nswap,t2lw_map)

             if(intf.ge.nintf) exit intf_loop
             !call chdir(dirname)
             call chdir(intf_dir)

             if(ludef_up_surf) exit up_term_loop
          end do up_term_loop
          if(ludef_lw_surf) exit lw_term_loop
       end do lw_term_loop
       !!-----------------------------------------------------------------------
       !! Returns to working directory
       !!-----------------------------------------------------------------------
       call chdir(intf_dir)

    end do intf_loop

    call chdir(pwd)


    return
  end subroutine generate_interfaces
!###############################################################################


!!!#############################################################################
!!! Takes input interface structure and generates a set of shifts and swaps.
!!! Prints these new structures to POSCARs.
!!!#############################################################################
!!! ISWAP METHOD NOT YET SET UP
  subroutine gen_shifts_and_swaps(basis,axis,intf_loc,bond,&
       ishift,nshift,&
       iswap,swap_den,nswap,&
       map)
    implicit none
    type(basis_type), intent(in) :: basis
    integer :: shift_unit=10
    integer :: ounit,iaxis,k,l
    integer :: ngen_swaps,nswaps_per_cell
    real(real32) :: dtmp1
    type(basis_type) :: tbas
    type(bond_type) :: min_bond
    character(1024) :: filename,dirpath,pwd1,pwd2,msg
    integer, dimension(3) :: abc
    real(real32), dimension(2) :: intf_loc
    real(real32), dimension(3) :: toffset
    real(real32), dimension(3,3) :: tlat
    type(basis_type), allocatable, dimension(:) :: bas_arr
    real(real32), allocatable, dimension(:,:) :: output_shifts

    integer, intent(in) :: axis
    integer, intent(in) :: nshift,nswap
    integer, intent(in) :: ishift,iswap
    real(real32), intent(in) :: bond,swap_den

    integer, dimension(:,:,:), optional, intent(in) :: map


!!!-----------------------------------------------------------------------------
!!! Sets up shift axis
!!!-----------------------------------------------------------------------------
    abc = [ 1, 2, 3 ]
    abc = cshift(abc,axis)


!!!-----------------------------------------------------------------------------
!!! Sets up and moves to appropriate directories
!!!-----------------------------------------------------------------------------
    call getcwd(pwd1)
    if(ishift.gt.0.or.nshift.gt.1)then
       call system('mkdir -p '//trim(adjustl(shiftdir)))
       call chdir(shiftdir)
    end if
    call getcwd(pwd2)
    open(unit=shift_unit,file="shift_vals.txt")
    write(shift_unit,&
         '("# interface_num    shift (a,b,c) units=(direct,direct,Å)")')


!!!-----------------------------------------------------------------------------
!!! Generates sets of shifts based on shift version
!!!-----------------------------------------------------------------------------
    if(ishift.eq.0.or.ishift.eq.1) allocate(output_shifts(nshift,3))
    select case(ishift)
    case(1)
       output_shifts(1,:3)=0._real32
       do k=2,nshift
          do iaxis=1,2
             call random_number(output_shifts(k,iaxis))
          end do
       end do
    case(2)
       output_shifts = get_fit_shifts(&
            lat=basis%lat,bas=basis,&
            bond=bond,&
            axis=axis,&
            intf_loc=intf_loc,&
            depth=intf_depth,&
            nstore=nshift)
    case(3)
       output_shifts = get_descriptive_shifts(&
            lat=basis%lat,bas=basis,&
            bond=bond,&
            axis=axis,&
            intf_loc=intf_loc,&
            depth=intf_depth,c_scale=c_scale,&
            nstore=nshift,lprint=lprint_shifts)
    case(4)
       if(present(map))then
          output_shifts = get_shifts_DON(&
               lat=basis%lat,bas=basis,&
               axis=axis,&
               intf_loc=intf_loc,&
               nstore=nshift,c_scale=c_scale,offset=offset(1,:3),&
               lprint=lprint_shifts,bulk_DON=bulk_DON,bulk_map=map,&
               max_bondlength=max_bondlength)
       else
          output_shifts = get_shifts_DON(&
               lat=basis%lat,bas=basis,&
               axis=axis,&
               intf_loc=intf_loc,&
               nstore=nshift,c_scale=c_scale,offset=offset(1,:3),&
               lprint=lprint_shifts,&
               max_bondlength=max_bondlength)
       end if
       if(size(output_shifts(:,1)).eq.0)then
          write(0,'(2X,"No shifts were identified with ISHIFT = 4 for this lattice match")')
          write(0,'(2X,"We suggest increasing MBOND_MAXLEN to find shifts")')
          write(0,'("Skipping interface...")')
          return
       end if
    case default
      ! nshift=1 !!! SORT THIS OUT !!! RESET NSHIFT DUE TO ISHIFT
       if(.not.allocated(output_shifts)) allocate(output_shifts(1,3))
       output_shifts(:,:) = offset
       do iaxis=1,2
          output_shifts(1,iaxis) = output_shifts(1,iaxis)!/modu(lat(iaxis,:))
       end do
    end select
    if(ishift.gt.0)then
       output_shifts(:,axis) = output_shifts(:,axis)*modu(basis%lat(axis,:))
    end if


!!!-----------------------------------------------------------------------------
!!! Prints number of shifts to terminal
!!!-----------------------------------------------------------------------------
    write(6,'(3X,"Number of unique shifts structures: ",I0)') nshift


!!!-----------------------------------------------------------------------------
!!! Determines number of swaps across the interface
!!!-----------------------------------------------------------------------------
    nswaps_per_cell=nint(swap_den*get_area([basis%lat(abc(1),:)],[basis%lat(abc(2),:)]))
    if(iswap.ne.0)then
       write(6,&
            '(" Generating ",I0," swaps per structure ")') nswaps_per_cell
    end if


!!!-----------------------------------------------------------------------------
!!! Prints each unique shift structure
!!!-----------------------------------------------------------------------------
    shift_loop: do k=1,nshift
       call tbas%copy(basis)
       toffset=output_shifts(k,:3)
       do iaxis=1,2
          call shift_region(tbas,axis,&
               intf_loc(1),intf_loc(2),&
               shift_axis=iaxis,shift=toffset(iaxis),renorm=.true.)
       end do
       dtmp1=modu(tlat(axis,:))
       call set_vacuum(&
            basis=tbas,&
            axis=axis,loc=maxval(intf_loc(:)),&
            vac=toffset(axis))
       dtmp1=minval(intf_loc(:))*dtmp1/modu(tlat(axis,:))
       call set_vacuum(&
            basis=tbas,&
            axis=axis,loc=dtmp1,&
            vac=toffset(axis))
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
       write(shift_unit,'(2X,I0.2,15X,"(",2(" ",F9.6,", ")," ",F9.6," )")') &
            k,toffset(:)


       !!-----------------------------------------------------------------------
       !! Merges lower and upper materials
       !! Writes interfaces to output directories
       !!-----------------------------------------------------------------------
       intf=intf+1
       ounit=100+intf
       if(ishift.gt.0.or.nshift.gt.1)then
          write(dirpath,'(A,I0.2)') trim(adjustl(subdir_prefix)),k
          call system('mkdir -p '//trim(adjustl(dirpath)))
          write(filename,'(A,"/",A)') trim(adjustl(dirpath)),trim(out_filename)
       else
          filename = trim(out_filename)
       end if
       write(6,'(2X,"Writing interface ",I0,"...")') intf
       open(unit=ounit,file=trim(adjustl(filename)))
       call geom_write(ounit,tbas)
       close(ounit)
       if(intf.ge.nintf) return


       !!-----------------------------------------------------------------------
       !! Performs swaps within the shifted structures if requested
       !!-----------------------------------------------------------------------
       if_swap: if(iswap.ne.0)then
          bas_arr = rand_swapper(tlat,tbas,axis,swap_depth,&
               nswaps_per_cell,nswap,intf_loc,iswap,seed,sigma=swap_sigma,&
               require_mirror=lswap_mirror)
          ngen_swaps = nswap
          LOOPswaps: do l=1,nswap
             if (bas_arr(l)%nspec.eq.0) then
                ngen_swaps = l - 1
                exit LOOPswaps
             end if
          end do LOOPswaps
          if(ngen_swaps.eq.0)then
             exit if_swap
          end if
          call chdir(dirpath)
          call system('mkdir -p '//trim(adjustl(swapdir)))
          call chdir(swapdir)
          write(6,'(3X,"Number of unique swap structures: ",I0)') ngen_swaps
          do l=1,ngen_swaps
             write(dirpath,'(A,I0.2)') trim(adjustl(subdir_prefix)),l
             call system('mkdir -p '//trim(adjustl(dirpath)))
             write(filename,'(A,"/",A)') &
                  trim(adjustl(dirpath)),trim(out_filename)
             ounit=100+l
             write(6,'(3X,"Writing swap ",I0,"...")') l
             open(unit=ounit,file=trim(adjustl(filename)))
             call geom_write(ounit,bas_arr(l))
             close(ounit)
          end do
          deallocate(bas_arr)
          call chdir(pwd2)
       end if if_swap


    end do shift_loop
    call chdir(pwd1)
    close(unit=shift_unit)


  end subroutine gen_shifts_and_swaps
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


end module artemis__generator
