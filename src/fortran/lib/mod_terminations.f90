module artemis__terminations
  !! Module for handling termination identification and generation
  use artemis__constants, only: real32, tolerance
  use artemis__geom_rw,   only: basis_type, geom_write
  use artemis__misc,      only: sort_col, to_lower, to_upper
  use artemis__io_utils,  only: err_abort
  use artemis__io_utils_extd, only: err_abort_print_struc
  use misc_linalg,        only: modu, cross, uvec, det
  use artemis__sym,       only: sym_type, check_sym, sym_setup, clone_grp, s_end
  use edit_geom,          only: shifter, transformer, ortho_axis, set_vacuum
  implicit none

  private

  public :: term_arr_type
  public :: get_termination_info
  public :: set_layer_tol
  public :: set_slab_height
  public :: build_slab


  type term_type
     !real(real32) :: add
     real(real32) :: hmin
     real(real32) :: hmax
     integer :: natom
     integer :: nstep
     real(real32), allocatable, dimension(:) :: ladder
  end type term_type

  type term_arr_type
     integer :: nterm = 0, axis, nstep
     real(real32) :: tol
     logical :: lmirror=.false.
     type(term_type), allocatable, dimension(:) :: arr
  end type term_arr_type

  type term_list_type
     integer :: term
     real(real32) :: loc
  end type term_list_type



contains

!###############################################################################
  function get_termination_info( &
       basis, axis, lprint, layer_sep, break_on_fail &
  ) result(term)
    !! Function to find the terminations of a material along a given axis
    implicit none

    ! Arguments
    type(basis_type), intent(in) :: basis
    !! Basis structure
    integer, intent(in) :: axis
    !! Axis to find terminations along (1,2,3)
    !! 1=a, 2=b, 3=c
    real(real32), intent(in), optional :: layer_sep
    !! Minimum separation between layers
    logical, intent(in), optional :: lprint
    !! Boolean whether to print terminations
    logical, intent(in), optional :: break_on_fail
    !! Boolean whether to break on failure to find terminations


    integer :: i,j,is,nterm,mterm,dim,ireject
    integer :: itmp1,itmp2,init,min_loc
    logical :: ludef_print,lunique,ltmp1,lmirror, break_on_fail_
    real(real32) :: rtmp1,tol,height,max_sep,c_along,centre
    type(sym_type) :: grp1,grp_store
    type(term_arr_type) :: term
    integer, dimension(3) :: abc=(/1,2,3/)
    real(real32), dimension(3) :: vec_compare
    real(real32), dimension(3,3) :: inv_mat,ident
    type(basis_type),allocatable, dimension(:) :: basis_arr,basis_arr_reject
    type(term_type), allocatable, dimension(:) :: term_arr,term_arr_uniq
    integer, allocatable, dimension(:) :: success,tmpop
    integer, allocatable, dimension(:,:) :: reject_match
    real(real32), allocatable, dimension(:,:) :: basis_list
    real(real32), allocatable, dimension(:,:,:) :: tmpsym
    character(len=256) :: err_msg




!!!APPLY TRANSFORMATION MATRIX TO FIND TERMINATIONS ALONG OTHER PLANES
!!! E.G. (1 0 1)
    
    term%nterm = 0
    s_end=0
    grp_store%confine%l=.false.
    grp_store%confine%axis=axis
    grp_store%confine%laxis=.false.
!!!-----------------------------------------------------------------------------
!!! Sets printing option
!!!-----------------------------------------------------------------------------
    if(present(lprint))then
       ludef_print = lprint
    else
       ludef_print = .false.
    end if
    break_on_fail_ = .false.
    if(present(break_on_fail)) break_on_fail_ = break_on_fail


!!!-----------------------------------------------------------------------------
!!! Sets the surface identification tolerance
!!!-----------------------------------------------------------------------------
    if(present(layer_sep))then
       tol = layer_sep
    else
       tol = 1._real32  !!!tolerance of 1 Å for defining a layer
    end if

    abc=cshift(abc,3-axis)
    c_along = abs(dot_product(basis%lat(axis,:),&
         uvec(cross([basis%lat(abc(1),:)],[basis%lat(abc(2),:)]))))
    tol = tol / c_along
    lmirror=.false.


!!!-----------------------------------------------------------------------------
!!! Set up basis list that will order them wrt distance along 'axis'
!!!-----------------------------------------------------------------------------
    allocate(basis_list(basis%natom,3))
    init = 1
    do is=1,basis%nspec
       basis_list(init:init+basis%spec(is)%num-1,:3) = basis%spec(is)%atom(:,:3)
       init = init + basis%spec(is)%num
    end do
    call sort_col(basis_list,col=axis)


!!!-----------------------------------------------------------------------------
!!! Find largest separation between atoms
!!!-----------------------------------------------------------------------------
    max_sep = basis_list(1,axis) - (basis_list(basis%natom,axis)-1._real32)
    height = ( basis_list(1,axis) + (basis_list(basis%natom,axis)-1._real32) )/2._real32
    do i=1,basis%natom-1
       rtmp1 = basis_list(i+1,axis) - basis_list(i,axis)
       if(rtmp1.gt.max_sep)then
          max_sep = rtmp1
          height = ( basis_list(i+1,axis) + basis_list(i,axis) )/2._real32
       end if
    end do
    if(max_sep.lt.tol)then
       if(break_on_fail_)then
          write(0,'("ERROR: Error in artemis__sym.f90")')
       else
          write(0,'("WARNING:")')
       end if
       write(0,'(2X,"get_terminations subroutine unable to find a separation &
            &in the material that is greater than LAYER_SEP")')
       write(0,'(2X,"Writing material to ''unlayerable.vasp''")')
       open(13,file="unlayerable.vasp")
       call geom_write(13,basis)
       close(13)
       write(0,'(2X,"We suggest reducing LAYER_SEP to less than ",F6.4)') &
            max_sep
       write(0,'(2X,"NOTE: If LAYER_SEP < 0.7, the material likely does not &
            &support the Miller plane")')
       write(0,'(2X,"Please inform the developers of this and give details &
            &of what structure caused this")')
       if(break_on_fail_)then
          write( 0, &
               '("To allow the program to continue, set &
               &LBREAK_ON_NO_TERM = F")' &
          )
          write(0,'("Stopping...")')
          call exit()
       else
          return
       end if
    end if
    basis_list(:,axis) = basis_list(:,axis) - height
    basis_list(:,axis) = basis_list(:,axis) - floor(basis_list(:,axis))
    call sort_col(basis_list,col=axis)


!!!-----------------------------------------------------------------------------
!!! Finds number of non-unique terminations
!!!-----------------------------------------------------------------------------
    nterm=1
    allocate(term_arr(basis%natom))
    term_arr(:)%natom=0
    term_arr(:)%hmin=0
    term_arr(:)%hmax=0
    term_arr(1)%hmin=basis_list(1,axis)
    term_arr(1)%hmax=basis_list(1,axis)
    min_loc = 1
    itmp1 = 1
    term_loop1: do

       !! get the atom at that height.
       !vtmp1 = get_min_dist(basis%lat,basis,basis_list(itmp1,:3),.true.,axis,.true.,.false.)
       !itmp1 = minloc(basis_list(:,axis) - vtmp1(axis), dim=1, &
       !     mask = abs(basis_list(:,axis) - (basis_list(itmp1,axis) + vtmp1(axis))&
       !     ).lt.tolerance)
       
       itmp1 = minloc(basis_list(:,axis) - term_arr(nterm)%hmax, dim=1, &
            mask = basis_list(:,axis) - term_arr(nterm)%hmax.gt.0._real32)
       if(itmp1.gt.basis%natom.or.itmp1.le.0)then
          term_arr(nterm)%natom = basis%natom - min_loc + 1
          exit term_loop1
       end if

       !rtmp1 = modu(matmul(vtmp1,basis%lat))
       rtmp1 = basis_list(itmp1,axis) - term_arr(nterm)%hmax
       if(rtmp1.le.tol)then
          term_arr(nterm)%hmax = basis_list(itmp1,axis)
       else
          term_arr(nterm)%natom = itmp1 - min_loc
          min_loc = itmp1
          nterm = nterm + 1
          term_arr(nterm)%hmin = basis_list(itmp1,axis)
          term_arr(nterm)%hmax = basis_list(itmp1,axis)
       end if
       
    end do term_loop1
    term_arr(:nterm)%hmin = term_arr(:nterm)%hmin + height
    term_arr(:nterm)%hmax = term_arr(:nterm)%hmax + height


!!!-----------------------------------------------------------------------------
!!! Set up system symmetries
!!!-----------------------------------------------------------------------------
    allocate(basis_arr(2*nterm))
    allocate(basis_arr_reject(2*nterm))
    dim = size(basis%spec(1)%atom(1,:))
    do i=1,2*nterm
       allocate(basis_arr(i)%spec(basis%nspec))
       allocate(basis_arr_reject(i)%spec(basis%nspec))
       do is=1,basis%nspec
          allocate(basis_arr(i)%spec(is)%atom(&
               basis%spec(is)%num,dim))
          allocate(basis_arr_reject(i)%spec(is)%atom(&
               basis%spec(is)%num,dim))
       end do
    end do


!!!-----------------------------------------------------------------------------
!!! Print location of unique terminations
!!!-----------------------------------------------------------------------------
    mterm = 0
    ireject = 0
    grp_store%lspace = .true.
    grp_store%confine%l = .true.
    grp_store%confine%laxis(axis) = .true.
    call sym_setup(grp_store,basis%lat,predefined=.false.,new_start=.true.)



    !!--------------------------------------------------------------------------
    !! Handle inversion matrix (centre of inversion must be accounted for)
    !!--------------------------------------------------------------------------
    !! change symmetry constraints after setting up symmetries
    !! this is done to constrain the matching of two basises in certain directions
    grp_store%confine%l = .false.
    grp_store%confine%laxis(axis) = .false.
    call check_sym(grp_store,bas1=basis,iperm=-1,lsave=.true.)
    inv_mat = 0._real32
    do i=1,3
       inv_mat(i,i) = -1._real32
    end do
    itmp1 = 0
    do i=1,grp_store%nsym
       if(all(abs(grp_store%sym(i,:3,:3)-inv_mat).lt.tolerance))then
          itmp1 = i
          exit
       end if
    end do
    if(itmp1.eq.0)then
       write(err_msg,*) "No inversion symmetry found!"
       call err_abort(err_msg)
    end if
    do i=1,grp_store%nsymop
       if(all(abs(grp_store%sym(i,:3,:3)-inv_mat).lt.tolerance)) &
            grp_store%sym(itmp1,4,:3) = grp_store%sym(i,4,:3)
    end do



    !!--------------------------------------------------------------------------
    !! Determine unique surface terminations
    !!--------------------------------------------------------------------------
    grp_store%confine%l = .true.
    grp_store%confine%laxis(axis) = .true.
    allocate(term_arr_uniq(2*nterm))
    allocate(reject_match(nterm,2))
    shift_loop1:do i=1,nterm
       mterm = mterm + 1

       basis_arr(mterm) = basis
       centre = term_arr(i)%hmin + (term_arr(i)%hmax - term_arr(i)%hmin)/2._real32
       call shifter(basis_arr(mterm),axis,1-centre,.true.)
       !if(ludef_print) write(6,'(1X,I3,8X,F7.5,9X,F7.5,8X,I3)') &
       !     i,term_arr(i)%hmin,term_arr(i)%hmax,term_arr(i)%natom
       sym_if: if(i.ne.1)then
          sym_loop1:do j=1,mterm-1
             if(abs(abs(term_arr(i)%hmax-term_arr(i)%hmin) - &
                  abs(term_arr_uniq(j)%hmax-term_arr_uniq(j)%hmin)).gt.tolerance) &
                  cycle sym_loop1
             call clone_grp(grp_store,grp1)
             call check_sym(grp1,bas1=basis_arr(mterm),&
                  iperm=-1,tmpbas2=basis_arr(j),lsave=.true.)
             if(grp1%nsymop.ne.0)then
                !write(0,*) "we have a possible reject"
                !if(any(grp1%sym_save(:grp1%nsymop,axis,axis).eq.-1._real32))then
                if(grp1%sym_save(1,axis,axis).eq.-1._real32)then
                   ireject = ireject + 1
                   reject_match(ireject,:) = [ i, j ]
                   basis_arr_reject(ireject) = basis_arr(mterm)
                   lmirror=.true.
                else
                   term_arr_uniq(j)%nstep = term_arr_uniq(j)%nstep + 1
                   term_arr_uniq(j)%ladder(term_arr_uniq(j)%nstep) = &
                        term_arr(i)%hmin - term_arr_uniq(j)%hmin
                end if
                mterm = mterm - 1
                cycle shift_loop1
             end if
          end do sym_loop1
       end if sym_if
       term_arr_uniq(mterm) = term_arr(i)
       term_arr_uniq(mterm)%nstep = 1
       allocate(term_arr_uniq(mterm)%ladder(nterm))
       term_arr_uniq(mterm)%ladder(:) = 0._real32
       !open(100+mterm)
       !call geom_write(100+mterm,basis_arr(mterm))
       !close(100+mterm)
    end do shift_loop1


    !!--------------------------------------------------------------------------
    !! Set up mirror/inversion symmetries of the matrix
    !!--------------------------------------------------------------------------
    call sym_setup(grp_store,basis%lat,predefined=.false.,new_start=.true.)
    allocate(tmpsym(count(grp_store%sym(:,3,3).eq.-1._real32),4,4))
    allocate(tmpop(count(grp_store%sym(:,3,3).eq.-1._real32)))
    itmp1 = 0
    do i=1,grp_store%nsym
       if(grp_store%sym(i,3,3).eq.-1._real32)then
          itmp1=itmp1+1
          tmpsym(itmp1,:,:) = grp_store%sym(i,:,:)
          tmpop(itmp1) = i
       end if
    end do
    grp_store%nsym = itmp1
    grp_store%nlatsym = itmp1
    call move_alloc(tmpsym,grp_store%sym)
    allocate(grp_store%op(itmp1))
    grp_store%op(:) = tmpop(:itmp1)
    s_end = grp_store%nsym


    !!--------------------------------------------------------------------------
    !! Check rejects for inverse surface termination of saved
    !!--------------------------------------------------------------------------
    ident = 0._real32
    do i=1,3
       ident(i,i) = 1._real32
    end do
    vec_compare = 0._real32
    vec_compare(axis) = -1._real32
    allocate(success(ireject))
    success=0
    reject_loop1: do i=1,ireject
       lunique=.true.
       itmp1=reject_match(i,1)
       itmp2=reject_match(i,2)
       !! Check if comparison termination has already been compared successfully
       prior_check: if(any(success(1:i-1).eq.itmp2))then
          lunique=.false.
       else
          call clone_grp(grp_store,grp1)
          call check_sym(grp1,basis_arr(itmp2),&
               iperm=-1,lsave=.true.,lcheck_all=.true.)
          ltmp1=.false.

!!!HERE
          !! Check if pure translations are present in comparison termination?
          !do j=1,grp1%nsymop
          !   if(all(abs(grp1%sym_save(j,:3,:3)-ident).le.tolerance))then
          !      write(0,*) "FOUND TRANSLATION"
          !      cycle reject_loop1
          !   end if
          !end do
          !! Check if inversions are present in comparison termination
          do j=1,grp1%nsymop
             if(abs(det(grp1%sym_save(j,:3,:3))+1._real32).le.tolerance) ltmp1=.true.
          end do
          !! If they are not, then no point comparing. It is a new termination
          if(.not.ltmp1) exit prior_check 

          call clone_grp(grp_store,grp1)
          call check_sym(grp1,basis_arr(itmp2),&
               tmpbas2=basis_arr_reject(i),iperm=-1,lsave=.true.,lcheck_all=.true.)

          !! Check det of all symmetry operations. If any are 1, move on
          !! This is because they are just rotations as can be captured ...
          !! ... through lattice matches.
          !! Solely inversions are unique and must be captured.
          do j=1,grp1%nsymop
             !write(0,'(4(2X,F9.4))') grp1%sym_save(j,:,:)
             !write(0,*) det(grp1%sym_save(j,:3,:3))
             if(abs(det(grp1%sym_save(j,:3,:3))-1._real32).le.tolerance) lunique=.false.
          end do
          if(grp1%sym_save(1,4,axis).eq.&
               2._real32*min(term_arr_uniq(itmp2)%hmin,0.5_real32-term_arr_uniq(itmp2)%hmin))then
             lunique=.false.
          end if

          if(.not.(all(grp1%sym_save(1,axis,:3).eq.vec_compare(:)).and.&
               all(grp1%sym_save(1,:3,axis).eq.vec_compare(:)))) lunique=.false.
          
       end if prior_check

       if(lunique)then
          mterm=mterm+1
          success(i)=itmp2
          term_arr_uniq(mterm)=term_arr(reject_match(i,1))
          !if(ludef_print) write(6,'(1X,I3,8X,F7.5,9X,F7.5,8X,I3)') &
          !     mterm,&
          !     term_arr_uniq(mterm)%hmin,&
          !     term_arr_uniq(mterm)%hmax,term_arr_uniq(mterm)%natom
          reject_match(i,2)=0
          term_arr_uniq(mterm)%nstep = 1
          allocate(term_arr_uniq(mterm)%ladder(ireject+1))
          term_arr_uniq(mterm)%ladder(1) = 0._real32
       else
          term_arr_uniq(itmp2)%nstep = term_arr_uniq(itmp2)%nstep + 1
          term_arr_uniq(itmp2)%ladder(term_arr_uniq(itmp2)%nstep) = &
               term_arr(itmp1)%hmin - term_arr_uniq(itmp2)%hmin
       end if
    end do reject_loop1


    !!--------------------------------------------------------------------------
    !! Populate termination output
    !!--------------------------------------------------------------------------
    allocate(term%arr(mterm))
    term%tol=tol
    term%axis=axis
    term%nterm=mterm
    term%lmirror = lmirror
    if(ludef_print)&
         write(6,'(1X,"Term.",3X,"Min layer loc",3X,"Max layer loc",3X,"no. atoms")')
    rtmp1 = term_arr_uniq(1)%hmin-1.E-6_real32
    itmp1 = 1
    do i=1,mterm
       allocate(term%arr(i)%ladder(term_arr_uniq(i)%nstep))
       term%arr(i)%hmin = term_arr_uniq(itmp1)%hmin
       term%arr(i)%hmax = term_arr_uniq(itmp1)%hmax
       term%arr(i)%natom = term_arr_uniq(itmp1)%natom
       term%arr(i)%nstep = term_arr_uniq(itmp1)%nstep
       term%arr(i)%ladder(:term%arr(i)%nstep) = term_arr_uniq(i)%ladder(:term%arr(i)%nstep)
       if(ludef_print) write(6,'(1X,I3,8X,F7.5,9X,F7.5,8X,I3)') &
            i,term%arr(i)%hmin,term%arr(i)%hmax,term%arr(i)%natom
       itmp1 = minloc(term_arr_uniq(:)%hmin,&
            mask=term_arr_uniq(:)%hmin.gt.rtmp1+tol,dim=1)
       if(itmp1.eq.0) then
          itmp1 = minloc(term_arr_uniq(:)%hmin,&
               mask=term_arr_uniq(:)%hmin.gt.rtmp1+tol-1._real32,dim=1)
       end if
       rtmp1 = term_arr_uniq(itmp1)%hmin
    end do
    term%nstep = maxval(term%arr(:)%nstep)


    !!--------------------------------------------------------------------------
    !! Check to ensure equivalent number of steps for each termination
    !!--------------------------------------------------------------------------
    !! Not yet certain whether each termination should have samve number ...
    !! ... of ladder rungs. That's why this check is here.
    if(all(term%arr(:)%nstep.ne.term%nstep))then
       write(0,'("ERROR: Number of rungs in terminations no equivalent for &
            &every termination! Please report this to developers.\n&
            &Exiting...")')
       call exit()
    end if


  end function get_termination_info
!###############################################################################


!###############################################################################
  function get_term_list(term) result(list)
    !! Function to get a list of all terminations in the system
    implicit none

    ! Arguments
    type(term_arr_type), intent(in) :: term
    !! Termination info
    type(term_list_type), allocatable, dimension(:) :: list
    !! List of terminations

    ! Local variables
    integer :: i, j
    !! Loop indices
    integer :: itmp1, nlist, loc
    !! Temporary indices
    type(term_list_type) :: tmp_element
    !! Temporary element for swapping


    if(.not.allocated(term%arr(1)%ladder))then
       nlist = term%nterm
       allocate(list(nlist))
       list(:)%loc = term%arr(:)%hmin
       do i = 1, term%nterm
          list(i)%term = i
       end do
    else
       nlist = term%nstep*term%nterm
       allocate(list(nlist))
       itmp1 = 0
       do i = 1, term%nterm
          do j = 1, term%nstep
             itmp1=itmp1+1
             list(itmp1)%loc = term%arr(i)%hmin+term%arr(i)%ladder(j)
             list(itmp1)%loc = list(itmp1)%loc - &
                  ceiling( list(itmp1)%loc - 1._real32 )
             list(itmp1)%term = i
          end do
       end do
    end if

    !! sort the list now
    do i = 1, nlist
       loc = minloc(list(i:nlist)%loc,dim=1) + i - 1
       tmp_element = list(i)
       list(i)     = list(loc)
       list(loc)   = tmp_element
    end do

  end function get_term_list
!###############################################################################


!###############################################################################
  subroutine set_layer_tol(term)
    !! Set the tolerance for the layer definitions
    implicit none

    ! Arguments
    type(term_arr_type), intent(inout) :: term
    !! Termination info

    ! Local variables
    integer :: i
    !! Loop index
    real(real32) :: rtmp1
    !! Temporary variable for tolerance


    do i = 1, term%nterm
       if(i.eq.1)then
          rtmp1 = abs(term%arr(i)%hmin - &
               (term%arr(term%nterm)%hmax+term%arr(i)%ladder(term%nstep)-1._real32)&
               )/4._real32
       else
          rtmp1 = abs(term%arr(i)%hmin-term%arr(i-1)%hmax)/4._real32
       end if
       if(rtmp1.lt.term%tol)then
          term%tol = rtmp1
       end if
    end do

    !! add the tolerances to the edges of the layers
    !! this ensures that all atoms in the layers are captured
    term%arr(:)%hmin = term%arr(:)%hmin - term%tol
    term%arr(:)%hmax = term%arr(:)%hmax + term%tol

  end subroutine set_layer_tol
!###############################################################################


!###############################################################################
  subroutine set_slab_height( basis, map, term, surf,&
       height, num_layers, thickness, ncells,&
       term_start, term_end, term_step &
  )
    !! Extend the basis to the maximum required height for all terminations
    !!
    !! This procedure extends the basis to form a supercell of the required
    !! integer extension along the surface normal vector. This supercell is
    !! sufficiently large to be able to be cut down to all required
    !! terminations.
    implicit none

    ! Arguments
    type(basis_type), intent(inout) :: basis
    !! Basis to be extended
    integer, allocatable, dimension(:,:,:), intent(inout) :: map
    !! Map from the original basis to the extended basis
    type(term_arr_type), intent(inout) :: term
    !! List of termination information
    integer, dimension(2), intent(in) :: surf
    !! Surface termination indices (for a single slab with both surface indices)
    integer, intent(in) :: num_layers
    integer, intent(out) :: term_start, term_end, ncells
    integer, intent(out) :: term_step
    real(real32), intent(in) :: thickness
    real(real32), intent(out) :: height


    integer :: i,itmp1
    real(real32) :: rtmp1, slab_thickness, largest_sep
    character(1024) :: msg
    real(real32), dimension(3,3) :: tfmat
    real(real32), allocatable, dimension(:) :: vtmp1
    type(term_list_type), allocatable, dimension(:) :: list


    integer :: icell, istep, iterm
    real(real32) :: layer_thickness
    logical :: success
    logical :: ludef_surf
    

    !!--------------------------------------------------------------------
    !! Initialise variables
    !!--------------------------------------------------------------------
    height = 0._real32


    !!-----------------------------------------------------------------------
    !! Defines height of slab from user-defined values
    !!-----------------------------------------------------------------------
    ludef_surf = .false.
    term_start = 1
    term_end = term%nterm
    if(all(surf.ne.0))then
       if(any(surf.gt.term%nterm))then
          write(msg, '("INVALID SURFACE VALUES!\nOne or more value &
               &exceeds the maximum number of terminations in the &
               structure.\n&
               &  Supplied values: ",I0,1X,I0,"\n&
               &  Maximum allowed: ",I0)') surf, term%nterm
          call err_abort(trim(msg),fmtd=.true.)
       end if
       ludef_surf = .true.
       list = get_term_list(term)
       !! set term_start to first surface value
       term_start = surf(1)
       !! set term_end to first surface value as a user-defined surface ...
       !! ... should not be cycled over.
       !! it is just one, potentially assymetric, slab to be explored.
       term_end = surf(1)

       !! determines the maximum number of cells required
       allocate(vtmp1(size(list)))
       height = term%arr(term_start)%hmin
       do i=num_layers,2,-1
          vtmp1 = list(:)%loc - height
          vtmp1 = vtmp1 - ceiling( vtmp1 - 1._real32 )
          itmp1 = minloc( vtmp1(:), dim=1,&
               mask=&
               vtmp1(:).gt.0.and.&
               list(:)%term.eq.surf(1))
          height = height + vtmp1(itmp1)
       end do
       vtmp1 = list(:)%loc - height
       !vtmp1 = vtmp1 - ceiling( vtmp1 - 1._real32 )
       where(vtmp1.lt.-1.D-5)
          vtmp1 = vtmp1 - ceiling( vtmp1 + 1.D-5 - 1._real32 )
       end where
       itmp1 = minloc( vtmp1(:), dim=1,&
            mask=&
            vtmp1(:).ge.-1.D-5.and.&
            list(:)%term.eq.surf(2))
       height = height + vtmp1(itmp1) - term%arr(term_start)%hmin

       !if(.not.term%lmirror)then
          ! get thickness of top/surface layer
          rtmp1 = term%arr(surf(2))%hmax - term%arr(surf(2))%hmin
          if(rtmp1.lt.-1.D-5) rtmp1 = rtmp1 + 1._real32
          height = height + rtmp1 !(1._real32 - rtmp1)
       !end if

       ncells = ceiling(height)
       height = height/real(ncells,real32)
    end if

    
    !!-----------------------------------------------------------------------
    !! Define termination iteration counter
    !!-----------------------------------------------------------------------
    if(term_end.lt.term_start)then
       term_step = -1
    else
       term_step = 1
    end if

    
    !!-----------------------------------------------------------------------
    !! Extend slab to user-defined thickness
    !!-----------------------------------------------------------------------
    if(.not.ludef_surf) ncells = int((num_layers-1)/term%nstep)+1
    !! convert thickness, in angstroms to number of cells
    if(thickness.gt.0._real32)then
       select case(term%axis)
       case(1)
          slab_thickness = dot_product(uvec(cross(basis%lat(2,:),basis%lat(3,:))), basis%lat(1,:))
       case(2)
          slab_thickness = dot_product(uvec(cross(basis%lat(1,:),basis%lat(3,:))), basis%lat(2,:))
       case(3)
          slab_thickness = dot_product(uvec(cross(basis%lat(1,:),basis%lat(2,:))), basis%lat(3,:))
       end select
       ! get the largest separation between two terminations
       if(ludef_surf)then

          height = 0.E0
          largest_sep = abs( term%arr(surf(1))%hmin - &
               term%arr(surf(2))%ladder(term%nstep) - &
               term%arr(surf(2))%hmax + 1._real32 )
          if(largest_sep.lt.0._real32) largest_sep = 1._real32 + largest_sep
          ! check for all terminations that a certain step is sufficiently large to reproduce thickness
          cell_loop1: do icell = 0, ceiling(thickness/slab_thickness), 1
             layer_thickness = term%arr(surf(2))%hmax - term%arr(surf(1))%hmin - 2.E0 * term%tol
             success = .false.
             step_loop1: do istep = 1, term%nstep, 1
                if(surf(2).lt.surf(1))then
                   if(istep.eq.term%nstep)then
                      layer_thickness = &
                           term%arr(surf(2))%hmax - term%arr(surf(1))%hmin - &
                           2.E0 * term%tol + ( &
                                1.E0 + term%arr(surf(2))%ladder(1) - &
                                term%arr(surf(1))%ladder(term%nstep) &
                           )
                   else
                      layer_thickness = &
                           term%arr(surf(2))%hmax - term%arr(surf(1))%hmin - &
                           2.E0 * term%tol + ( &
                                term%arr(surf(2))%ladder(istep+1) - &
                                term%arr(surf(1))%ladder(istep) &
                           )
                   end if
                end if
                rtmp1 = &
                     ( &
                          icell + layer_thickness + &
                          term%arr(surf(2))%ladder(istep) - &
                          term%arr(surf(1))%ladder(1) &
                     ) * slab_thickness
                if(rtmp1.ge.thickness)then
                   success = .true.
                   height = rtmp1 + 2.E0 * term%tol * slab_thickness
                   exit step_loop1
                end if
             end do step_loop1
             if(.not.success) cycle cell_loop1
             ncells = icell + 1
             exit cell_loop1
          end do cell_loop1
          
       else
          largest_sep = abs( term%arr(1)%hmin - &
               term%arr(1)%ladder(term%nstep) - &
               term%arr(1)%hmax + 1._real32 )
          if(largest_sep.lt.0._real32) largest_sep = 1._real32 + largest_sep
          ! check for all terminations that a certain step is sufficiently large to reproduce thickness
          cell_loop2: do icell = 0, ceiling(thickness/slab_thickness), 1
             term_loop: do iterm = 1, term%nterm, 1
                layer_thickness = term%arr(iterm)%hmax - term%arr(iterm)%hmin - 2.E0 * term%tol
                success = .false.
                step_loop: do istep = 1, term%nstep, 1
                   rtmp1 = ( icell + layer_thickness + term%arr(iterm)%ladder(istep) ) * slab_thickness
                   if(rtmp1.ge.thickness)then
                      success = .true.
                      exit step_loop
                   end if
                end do step_loop
                if(.not.success) cycle cell_loop2
             end do term_loop
             ncells = icell + 1
             exit cell_loop2
          end do cell_loop2

       end if
       height = height/real(ncells * slab_thickness,real32)
    end if
    tfmat(:,:) = 0._real32
    tfmat(1,1) = 1._real32
    tfmat(2,2) = 1._real32
    tfmat(3,3) = ncells
    call transformer(basis,tfmat,map)

    
    !!-----------------------------------------------------------------------
    !! Readjust termination plane locations
    !! ... i.e. divide all termination values by the number of cells
    !!-----------------------------------------------------------------------
    term%arr(:)%hmin = term%arr(:)%hmin/real(ncells,real32)
    term%arr(:)%hmax = term%arr(:)%hmax/real(ncells,real32)
    term%tol = term%tol/real(ncells,real32)
    

  end subroutine set_slab_height
!###############################################################################



!###############################################################################
  subroutine build_slab( &
       basis, map, term, surf, thickness, ncells, num_layers, &
       height, lwup_in, lcycle, &
       orthogonalise, vacuum &
  )
    !! Build a slab of the specified terminations
    !!
    !! This procedure builds a slab of the specified terminations from a
    !! supplied supercell. The supercell must be large enough to be able to
    !! be cut down to the required slab size. The supercell is built by
    !! set_slab_height.
    implicit none

    ! Arguments
    type(basis_type), intent(inout) :: basis
    !! Basis to be extended
    integer, allocatable, dimension(:,:,:), intent(inout) :: map
    !! Map from the original basis to the extended basis
    type(term_arr_type), intent(in) :: term
    !! Termination info
    integer, dimension(2), intent(in) :: surf
    !! Surface termination indices (for a single slab with both surface indices)

    integer, intent(in) :: num_layers, ncells
    real(real32), intent(in) :: height, thickness
    character(2), intent(in) :: lwup_in
    logical, intent(out) :: lcycle

    logical, optional, intent(in) :: orthogonalise
    real(real32), optional, intent(in) :: vacuum


    ! Local variables
    integer :: term_btm_idx, term_top_idx
    !! Indices of the bottom and top terminations
    logical :: equivalent_surfaces
    !! Boolean whether the two surfaces are equivalent
    real(real32) :: vacuum_
    integer :: j, j_start, istep, natom_check
    real(real32) :: rtmp1, slab_thickness, shift_val
    character(2) :: lwup
    character(5) :: lowerupper
    character(1024) :: msg
    logical :: orthogonalise_
    integer, dimension(3) :: abc=(/1,2,3/)
    real(real32), dimension(3) :: surface_normal_vec
    real(real32), dimension(3,3) :: tfmat
    integer, allocatable, dimension(:) :: iterm_list


    integer :: icell, num_cells
    real(real32) :: layer_thickness, ladder_adjust

    !!--------------------------------------------------------------------
    !! Initialise variables
    !!--------------------------------------------------------------------
    lwup=to_lower(lwup_in)
    if(lwup.eq."lw") lowerupper="LOWER"
    if(lwup.eq."up") lowerupper="UPPER"
    lcycle = .false.
    rtmp1=0._real32
    tfmat=0._real32
    term_btm_idx = surf(1)
    if(surf(2).gt.0)then
       term_top_idx = surf(2)
    else
       term_top_idx = surf(1)
    end if
    equivalent_surfaces = .false.
    if(term_btm_idx.eq.term_top_idx) equivalent_surfaces = .true.
    select case(term%axis)
    case(1)
       surface_normal_vec = uvec(cross( [ basis%lat(2,:) ], [ basis%lat(3,:) ]))
       slab_thickness = abs( dot_product(surface_normal_vec, [ basis%lat(1,:) ]) )
    case(2)
       surface_normal_vec = uvec(cross( [ basis%lat(1,:) ], [ basis%lat(3,:) ]))
       slab_thickness = abs( dot_product(surface_normal_vec, [ basis%lat(2,:) ]) )
    case(3)
       surface_normal_vec = uvec(cross( [ basis%lat(1,:) ], [ basis%lat(2,:)] ))
       slab_thickness = abs( dot_product(surface_normal_vec, [ basis%lat(3,:) ]) )
    end select
    if(thickness.gt.0._real32)then
       rtmp1 = slab_thickness / ncells * ( ncells - 1 )
       istep = term%nstep
       num_cells = ncells - 1
       cell_loop: do icell = 0, ncells, 1
          layer_thickness = term%arr(term_top_idx)%hmax - term%arr(term_btm_idx)%hmin - 2.E0 * term%tol
          ladder_adjust = 0._real32
          step_loop: do j = 1, term%nstep
             if(term_top_idx.lt.term_btm_idx)then
                if(j.eq.term%nstep)then
                   layer_thickness = term%arr(term_top_idx)%hmax - term%arr(term_btm_idx)%hmin - 2.E0 * term%tol
                   ladder_adjust = 1.E0 + term%arr(term_top_idx)%ladder(1) - term%arr(term_btm_idx)%ladder(term%nstep)
                else
                   layer_thickness = term%arr(term_top_idx)%hmax - term%arr(term_btm_idx)%hmin - 2.E0 * term%tol
                   ladder_adjust = term%arr(term_top_idx)%ladder(j+1) - term%arr(term_btm_idx)%ladder(j)
                end if
             end if
             rtmp1 = ( icell / real(ncells,real32) + layer_thickness ) * slab_thickness + &
                  ( ladder_adjust + term%arr(term_top_idx)%ladder(j) - term%arr(term_btm_idx)%ladder(1) ) * slab_thickness / real(ncells,real32)
             if(rtmp1.ge.thickness)then
                istep = j
                num_cells = icell
                exit cell_loop
             end if
          end do step_loop
       end do cell_loop
    else
       istep = num_layers - (ncells-1)*term%nstep
       num_cells = ncells - 1
    end if
    natom_check = basis%natom

    orthogonalise_ = .true.
    if(present(orthogonalise)) orthogonalise_ = orthogonalise

    vacuum_ = 10._real32
    if(present(vacuum)) vacuum_ = vacuum


    !!--------------------------------------------------------------------
    !! Set up list for checking expected number of atoms
    !!--------------------------------------------------------------------
    allocate(iterm_list(term%nterm))
    do j = 1, term%nterm
       iterm_list(j) = j
    end do
    iterm_list = cshift( iterm_list, term_btm_idx - 1 )
    if(.not.equivalent_surfaces)then
       j_start = term_top_idx - term_btm_idx + 1
       if(j_start.le.0) j_start = j_start + term%nterm
       j_start = j_start + 1 !+ (istep-1)*term%nterm/term%nstep
    else
       !! handle ladder steps that are equivalent
       j_start = 2 !+ (istep-1)*term%nterm/term%nstep
    end if


    !!--------------------------------------------------------------------
    !! Shift lower material to specified termination
    !!--------------------------------------------------------------------
    call shifter(basis,term%axis,-term%arr(term_btm_idx)%hmin,.true.)


    !!--------------------------------------------------------------------
    !! Determine cell reduction to specified termination
    !!--------------------------------------------------------------------
    do j = 1, 3
       tfmat(j,j) = 1._real32
       if(j.eq.term%axis)then
          if(.not.equivalent_surfaces)then
             tfmat(j,j) = height
          else!if(term%lmirror)then
             if(istep.ne.0)then
                rtmp1 = num_cells + term%arr(term_btm_idx)%ladder(istep)
                rtmp1 = rtmp1/real(ncells, real32)
                tfmat(j,j) = rtmp1 + &
                     (term%arr(term_btm_idx)%hmax - term%arr(term_btm_idx)%hmin)
             end if
          !else
          !   tfmat(j,j) = tfmat(j,j) + (&
          !        term%arr(term_btm_idx)%hmax - &
          !        term%arr(term_btm_idx)%hmin) + term%tol*2._real32
          end if
       end if
    end do


    !!--------------------------------------------------------------------
    !! Check number of atoms is expected
    !!--------------------------------------------------------------------
    if(num_cells.ne.ncells-1)then
       do icell = num_cells + 2, ncells, 1
          natom_check = natom_check - nint( basis%natom / real(ncells) )
       end do
    end if


    !!--------------------------------------------------------------------
    !! Apply transformation and shift cell back to bottom of layer
    !! ... i.e. account for the tolerance that has been added to layer ...
    !! ... hmin and hmax
    !!--------------------------------------------------------------------
    shift_val = term%tol * slab_thickness / modu(basis%lat(term%axis,:))
    call transformer(basis,tfmat,map)
    call shifter(basis,term%axis,-shift_val/tfmat(term%axis,term%axis),.true.)


    !!--------------------------------------------------------------------
    !! Check number of atoms is expected
    !!--------------------------------------------------------------------
    if(term%nterm.gt.1.or.term%nstep.gt.1)then
       do j = 1, max(0,term%nstep-istep), 1
          natom_check = natom_check - sum(term%arr(:)%natom)
       end do
       do j = j_start, term%nterm, 1
          natom_check = natom_check - term%arr(iterm_list(j))%natom
       end do
    end if
    if(basis%natom.ne.natom_check)then
       write(msg, '("NUMBER OF ATOMS IN '//to_upper(lowerupper)//' SLAB! &
            &Expected ",I0," but generated ",I0," instead")') &
            natom_check,basis%natom
       if(tfmat(term%axis,term%axis).gt.1._real32)then
          write(0,'("THE TRANSFORMATION IS GREATER THAN ONE ",F0.9)') &
               tfmat(term%axis,term%axis)
       end if
       call err_abort_print_struc(basis,lwup//"_term.vasp",&
            trim(msg),.true.)
       lcycle = .true.
    end if


    !!--------------------------------------------------------------------
    !! Apply slab_cuber to orthogonalise lower material
    !!--------------------------------------------------------------------
    call set_vacuum(basis,term%axis,1._real32-term%tol/tfmat(term%axis,term%axis),vacuum_)
    !call err_abort_print_struc(basis,"check.vasp","stop")
    abc=cshift(abc,3-term%axis)
    if(orthogonalise_)then
       ortho_check: do j=1,2
          if(abs(dot_product(basis%lat(abc(j),:),basis%lat(term%axis,:))).gt.1.D-5)then
             call ortho_axis(basis%lat,basis,term%axis)
             exit ortho_check
          end if
       end do ortho_check
    end if
    call basis%normalise(ceil_val=0.9999_real32,floor_coords=.true.,zero_round=0._real32)


  end subroutine build_slab
!###############################################################################

end module artemis__terminations
