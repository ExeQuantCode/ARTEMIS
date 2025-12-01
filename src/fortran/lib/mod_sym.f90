!!!#############################################################################
!!! Code written by Ned Thaddeus Taylor and Francis Huw Davies
!!! Code part of the ARTEMIS group (Hepplestone research group).
!!! Think Hepplestone, think HRG.
!!!#############################################################################
!!!module contains symmetry-related functions and subroutines.
!!!module includes the following functions and subroutines:
!!! check_sym         (checks supplied symmetries against supplied basis or ...
!!!                    ... checks whether the two supplied bases match after ...
!!!                    ... applying symmetries)
!!! gldfnd            (output translations that maps two bases)
!!! mksym             (makes array of symmetries that apply to supplied lattice
!!! basis_map         (finds symmetry equivalent atoms in two bases based on ...
!!!                    ... the supplied transformation matrix)
!!!#############################################################################
module artemis__sym
  use artemis__constants,   only: real32, pi
  use misc_linalg,          only: modu, inverse_3x3, det, uvec
  use artemis__geom_rw,     only: basis_type
  implicit none


  private


  public :: tol_sym_default
  public :: sym_type
  public :: check_sym, gldfnd

  public :: confine_type

  public :: basis_map_type, basis_map



  real(real32) :: tol_sym_default = 1.E-6_real32
  integer, allocatable, dimension(:) :: symops_compare

  interface get_wyckoff_atoms
     procedure get_wyckoff_atoms_any,get_wyckoff_atoms_loc
  end interface get_wyckoff_atoms


  type spec_wyck_type
     integer :: num
     character(len=5) :: name
     integer, allocatable, dimension(:) :: atom
  end type spec_wyck_type
  type wyck_type
     integer :: nwyck
     type(spec_wyck_type), allocatable, dimension(:) :: spec
  end type wyck_type


  type spcmap_type
     integer, allocatable ,dimension(:) :: atom
  end type spcmap_type
  type basis_map_type
     type(spcmap_type), allocatable, dimension(:) :: spec
  end type basis_map_type

  type confine_type
     !! apply any confinement/constraints on symmetries
     logical :: l=.false.
     !! axis to confine
     integer :: axis=0
     !! states whether to consider mirrors in only one plane
     logical :: lmirror=.false.
     !! if l=.false. -> laxis defines which axes are free
     !! if l=.true.  -> laxis defines which axes  are confined
     logical, dimension(3) :: laxis=(/.false.,.false.,.false./)
  end type confine_type

  type sym_type
     integer :: nsym = 0
     integer :: nlatsym = 0
     integer :: nsymop = 0
     integer :: npntop = 0
     logical :: lspace = .true.
     logical :: lmolec = .false.
     integer :: start_idx = 1, end_idx  =0
     integer, allocatable, dimension(:) :: op
     real(real32), allocatable, dimension(:,:,:) :: sym
     type(confine_type) :: confine
     real(real32), allocatable, dimension(:,:,:) :: sym_save
   contains
     procedure, pass(this) :: init => initialise_sym_type
     procedure, pass(this) :: copy => copy_sym_type
  end type sym_type




contains

!###############################################################################
  subroutine initialise_sym_type(this,lat,predefined,new_start,tol_sym)
    !! Initialises the symmetry container
    implicit none

    ! Arguments
    class(sym_type), intent(inout) :: this
    real(real32), dimension(3,3), intent(in) :: lat
    logical, optional, intent(in) :: predefined
    logical, optional, intent(in) :: new_start
    real(real32), optional, intent(in) :: tol_sym


    real(real32) :: tol_sym_
    logical :: predefined_, new_start_


    tol_sym_ = tol_sym_default
    if(present(tol_sym)) tol_sym_ = tol_sym
    if(present(new_start))then
       if(new_start)then
          if(allocated(this%op)) deallocate(this%op)
          if(allocated(this%sym)) deallocate(this%sym)
       end if
    end if

    predefined_ = .true.
    if(present(predefined)) predefined_ = predefined
    if(predefined_)then
       call gen_fundam_sym_matrices(this, lat, tol_sym_)
    else
       call mksym(this, lat, tol_sym_)
    end if

    if(allocated(symops_compare)) deallocate(symops_compare)
    this%nsymop=0

    new_start_ = .true.
    if(present(new_start)) new_start_ = new_start
    if(new_start_.or.this%end_idx.eq.0)then
       this%end_idx = this%nsym
    end if

  end subroutine initialise_sym_type
!###############################################################################


!###############################################################################
  subroutine copy_sym_type(this, source)
    !! Copy symmetry container
    implicit none

    ! Arguments
    class(sym_type), intent(inout) :: this
    !! Destination symmetry group
    type(sym_type), intent(in) :: source
    !! Source symmetry group
   

    if(allocated(this%op)) deallocate(this%op)
    if(allocated(this%sym)) deallocate(this%sym)
    if(allocated(this%sym_save)) deallocate(this%sym_save)

    this%nsym = source%nsym
    this%nlatsym = source%nlatsym
    this%nsymop = source%nsymop
    this%npntop = source%npntop
    this%lspace = source%lspace
    this%lmolec = source%lmolec
    this%start_idx = source%start_idx
    this%end_idx = source%end_idx
    this%confine = source%confine

    if(allocated(source%op)) &
         allocate(this%op, source = source%op)
    if(allocated(source%sym)) &
         allocate(this%sym, source = source%sym)
    if(allocated(source%sym_save)) &
         allocate(this%sym_save, source = source%sym_save)

  end subroutine copy_sym_type
!###############################################################################


!!!#############################################################################
!!! builds an array of the symmetries that apply to the supplied lattice
!!!#############################################################################
!!! tfbas   : transformed basis
!!!#############################################################################
  subroutine check_sym( &
       grp, basis, iperm, tmpbas2, wyckoff, lsave, lat, loc, check_all_sym, &
       verbose, tol_sym &
  )
    implicit none
    type(basis_type), intent(in) :: basis
    type(sym_type), intent(inout) :: grp

    integer, optional, intent(in) :: iperm
    logical, optional, intent(in) :: lsave,check_all_sym
    type(basis_type), optional, intent(in) :: tmpbas2
    type(wyck_type), optional, intent(inout) :: wyckoff
    real(real32), dimension(3), optional, intent(in) :: loc
    real(real32), dimension(3,3), optional, intent(in) :: lat
    integer, optional, intent(in) :: verbose
    real(real32), optional, intent(in) :: tol_sym

    integer :: i,j,k,iatom,jatom,ispec,itmp1
    integer :: is,isym,jsym,count,ntrans
    integer :: samecount,oldnpntop
    logical :: lsave_,lwyckoff,ltransformed, is_a_symmetry
    integer :: verbose_
    logical :: check_all_sym_
    real(real32) :: tol_sym_
    type(basis_type) :: basis2, tfbas
    real(real32), dimension(3) :: diff
    real(real32), dimension(3,3) :: ident
    type(wyck_type), allocatable, dimension(:) :: wyck_check
    real(real32), allocatable, dimension(:,:) :: trans
    real(real32), allocatable, dimension(:,:,:) :: tmpsav


    verbose_ = 0
    tol_sym_ = tol_sym_default
    if(present(verbose)) verbose_ = verbose
    if(present(tol_sym)) tol_sym_ = tol_sym
204 format(4(F11.6),/,4(F11.6),/,4(F11.6),/,4(F11.6))

    ! check length of basis
    do is = 1, basis%nspec
       if(size(basis%spec(is)%atom,2).ne.4)then
          write(0,'("ERROR: error encountered in check_sym")')
          write(0,'(2X,"Internal error in subroutine check_sym in artemis__sym.f90")')
          write(0,'(2X,"size of basis is not 4")')
          return
       end if
    end do


!!!-----------------------------------------------------------------------------
!!! allocated grp%op
!!!-----------------------------------------------------------------------------
    if(allocated(grp%op)) deallocate(grp%op)
    allocate(grp%op(grp%nsym*minval(basis%spec(:)%num)))
    grp%op = 0

    if(present(lsave))then
       lsave_ = lsave
    else
       lsave_ = .false.
    end if


!!!-----------------------------------------------------------------------------
!!! checks for optional arguments and assigns values if not present
!!!-----------------------------------------------------------------------------
    check_all_sym_ = .true.
    if(present(tmpbas2)) then
       call basis2%copy(tmpbas2)
       if(present(check_all_sym)) check_all_sym_ = check_all_sym
    else
       call basis2%copy(basis)
    end if
    allocate(tmpsav(4,4,grp%nsym*minval(basis%spec(:)%num)))
    itmp1 = maxval(basis%spec(:)%num)


!!!-----------------------------------------------------------------------------
!!! initialises variables
!!!-----------------------------------------------------------------------------
    allocate(trans(minval(basis%spec(:)%num+2),3)); trans = 0._real32
    allocate(tfbas%spec(basis%nspec))
    itmp1 = size(basis%spec(1)%atom(1,:),dim=1)
    do is=1,basis%nspec
       allocate(tfbas%spec(is)%atom(basis%spec(is)%num,itmp1))
    end do
    grp%nsymop = 0
    grp%npntop = 0


!!!-----------------------------------------------------------------------------
!!! if present, initialises wyckoff arrays
!!!-----------------------------------------------------------------------------
    allocate(wyck_check(grp%nsym*minval(basis%spec(:)%num)))
    do isym=1,grp%nsym*minval(basis%spec(:)%num)
       allocate(wyck_check(isym)%spec(basis%nspec))
       do ispec=1,basis%nspec
          allocate(wyck_check(isym)%spec(ispec)%atom(basis%spec(ispec)%num))
          wyck_check(isym)%spec(ispec)%atom = 0
       end do
    end do
    if(present(wyckoff))then
       lwyckoff = .true.
       if(allocated(wyckoff%spec)) deallocate(wyckoff%spec)
       wyckoff%nwyck = 0
       allocate(wyckoff%spec(basis%nspec))
       do ispec=1,basis%nspec
          wyckoff%spec(ispec)%num = 0
          wyckoff%spec(ispec)%name = ""
          allocate(wyckoff%spec(ispec)%atom(basis%spec(ispec)%num))
          do iatom=1,basis%spec(ispec)%num
             wyckoff%spec(ispec)%atom(iatom) = iatom
          end do
       end do
    else
       lwyckoff = .false.
    end if


!!!-----------------------------------------------------------------------------
!!! set up identity matrix as reference
!!!-----------------------------------------------------------------------------
    ltransformed = .false.
    ident = 0._real32
    do i=1,3
       ident(i,i) = 1._real32
    end do


!!!-----------------------------------------------------------------------------
!!! applying symmetries to basis to see if the basis conforms to any of them
!!!-----------------------------------------------------------------------------
    itmp1 = 1
    symloop: do isym = grp%start_idx, grp%end_idx, 1
       if(verbose_.eq.2.or.verbose_.eq.3) write(*,204)  &
            grp%sym(1:4,1:4,isym)
       !------------------------------------------------------------------------
       ! apply symmetry operator to basis
       !------------------------------------------------------------------------
       do ispec = 1, basis%nspec, 1
          do iatom = 1, basis%spec(ispec)%num, 1
             tfbas%spec(ispec)%atom(iatom,1:3) = &
                  matmul(basis%spec(ispec)%atom(iatom,1:4),grp%sym(1:4,1:3,isym))
             do j=1,3
                tfbas%spec(ispec)%atom(iatom,j) = &
                     tfbas%spec(ispec)%atom(iatom,j) - &
                     ceiling(tfbas%spec(ispec)%atom(iatom,j)-0.5_real32)
             end do
          end do
       end do
       !------------------------------------------------------------------------
       ! check whether transformed basis matches original basis
       !------------------------------------------------------------------------
       count=0
       is_a_symmetry = .true.
       spcheck: do ispec = 1, basis%nspec, 1
          diff = 0._real32
          samecount = 0
          wyck_check(itmp1)%spec(ispec)%atom = 0
          atmcheck: do iatom = 1, basis%spec(ispec)%num, 1
             atmcyc: do jatom = 1, basis%spec(ispec)%num, 1
                !if(wyck_check(itmp1)%spec(ispec)%atom(jatom).ne.0) cycle atmcyc
                diff = tfbas%spec(ispec)%atom(iatom,1:3) - &
                     basis2%spec(ispec)%atom(jatom,1:3)
                diff(:) = diff(:) - floor(diff(:))
                where(abs(diff(:)-1._real32).lt.tol_sym_)
                   diff(:)=0._real32
                end where
                if(sqrt(dot_product(diff,diff)).lt.tol_sym_)then
                   samecount = samecount + 1
                   wyck_check(itmp1)%spec(ispec)%atom(iatom) = jatom
                end if
                if((iatom.eq.basis%spec(ispec)%num).and.&
                     (jatom.eq.basis%spec(ispec)%num))then
                   if (samecount.ne.basis%spec(ispec)%num)then
                     is_a_symmetry = .false.
                     exit spcheck
                   end if
                end if
             end do atmcyc
             count = count + samecount
          end do atmcheck
          if(samecount.ne.basis%spec(ispec)%num)then
             is_a_symmetry = .false.
             exit spcheck
          end if
       end do spcheck
       if(is_a_symmetry)then
          grp%npntop = grp%npntop + 1
          grp%nsymop = grp%nsymop + 1
          itmp1 = grp%nsymop + 1
          tmpsav(:,:,grp%nsymop) = grp%sym(:,:,isym)
          grp%op(grp%nsymop) = isym
          if(grp%nsymop.ne.0.and..not.check_all_sym_) exit symloop
       end if
       trans = 0._real32
       ntrans = 0
       !------------------------------------------------------------------------
       ! checks if translations are valid with the current symmetry operation
       !------------------------------------------------------------------------
       if(grp%lspace) then
          if(all(abs(grp%sym(1:3,1:3,isym)-ident).lt.tol_sym_))then
             ltransformed=.false.
          else
             ltransformed=.true.
          end if
          call gldfnd(grp%confine,&
               basis2,tfbas,&
               trans,ntrans,&
               tol_sym_,&
               transformed=ltransformed,&
               wyck_check=wyck_check(itmp1:))
          if(ntrans.gt.0) then
             if(.not.check_all_sym_.and..not.lsave_)then
                grp%nsymop = grp%nsymop + 1
                exit symloop
             end if
             transloop: do i = 1, ntrans, 1
                if(dot_product(trans(i,:),trans(i,:)).lt.tol_sym_) &
                     cycle transloop
                if(verbose_.eq.3) write(*,*) trans(i,:)
                if(isym.ne.1)then
                   do jsym=2,grp%nsymop
                      if(grp%op(jsym).eq.1) then
                         if(all(abs(trans(i,1:3)-tmpsav(4,1:3,jsym)).lt.&
                              tol_sym_)) cycle transloop
                         diff = trans(i,1:3) - tmpsav(4,1:3,jsym)
                         diff = diff - ceiling( diff - 0.5_real32 )
                         do k=1,i
                            if(all(abs(diff-trans(k,1:3)).lt.tol_sym_)) &
                                 cycle transloop
                         end do
                      end if
                   end do
                end if
                grp%nsymop = grp%nsymop + 1
                itmp1 = grp%nsymop + 1
                tmpsav(:,:,grp%nsymop) = grp%sym(:,:,isym)
                tmpsav(4,1:3,grp%nsymop) = trans(i,:)
                grp%op(grp%nsymop) = isym
             end do transloop
             if(.not.check_all_sym_) exit symloop
          end if
       end if
       oldnpntop = grp%npntop
    end do symloop


!!!-----------------------------------------------------------------------------
!!! allocates and saves the array sym_save if the first time submitted
!!!-----------------------------------------------------------------------------
    if(lsave_)then
       if(allocated(grp%sym_save)) deallocate(grp%sym_save)
       allocate(grp%sym_save(4,4,grp%nsymop))
       grp%sym_save=0._real32
       grp%sym_save(:,:,:grp%nsymop) = tmpsav(:,:,:grp%nsymop)
       grp%sym_save(4,4,:) = 1._real32
       deallocate(tmpsav)
    end if


    iperm_if: if(present(iperm))then
       select case(iperm)
       case(-1)
          return
       case(0)
          exit iperm_if
       case default
          if(.not.allocated(symops_compare))then
             write(0,'("ERROR: Internal error in check_sym")')
             write(0,'(2X,"check_sym in artemis__sym.f90 is trying to assign a &
                  &value to symops_compare, which hasn''t been allocated")')
             exit iperm_if
          end if
          symops_compare(iperm)=grp%nsymop
       end select
    end if iperm_if


    if(lsave_)then
       deallocate(grp%sym)
       call move_alloc(grp%sym_save, grp%sym)
       grp%nsym = grp%nsymop
    end if


!!!-----------------------------------------------------------------------------
!!! if wyckoff present, set up wyckoff atoms
!!!-----------------------------------------------------------------------------
    if(lwyckoff)then
       if(present(lat).and.present(loc))then
          wyckoff=get_wyckoff_atoms(wyck_check(:grp%nsymop),lat,basis,loc)
       else       
          wyckoff=get_wyckoff_atoms(wyck_check(:grp%nsymop))
       end if
    end if

  end subroutine check_sym
!!!#############################################################################


!!!#############################################################################
!!! supplies the glides (if any) that are required to match the two bases ...
!!! ... "basis1" and "basis2" onto one another
!!!#############################################################################
  subroutine gldfnd( &
       confine, basis1, basis2, &
       trans, ntrans, &
       tol_sym, &
       transformed, wyck_check &
  )
    implicit none
    type(confine_type), intent(in) :: confine
    type(basis_type), intent(in) :: basis1,basis2
    real(real32), dimension(:,:), intent(out) :: trans
    integer, intent(out) :: ntrans
    real(real32), intent(in) :: tol_sym

    logical, optional, intent(in) :: transformed

    type(wyck_type), dimension(:), optional, intent(inout) :: wyck_check

    integer :: i,j,ispec,iatom,jatom,katom,itmp1
    integer :: minspecloc,samecount
    logical :: lwyckoff
    real(real32), dimension(3) :: ttrans,tmpbas,diff
    real(real32), allocatable, dimension(:,:) :: sav_trans



!!!-----------------------------------------------------------------------------
!!! Allocate arrays and initialise variables
!!!-----------------------------------------------------------------------------
    ttrans=0._real32
    trans=0._real32
    samecount=0
    ntrans=0
    minspecloc=minloc(basis1%spec(:)%num,mask=basis1%spec(:)%num.ne.0,dim=1)

    if(present(transformed))then
       if(.not.transformed)then
          if(basis1%spec(minspecloc)%num.eq.1) return
       end if
    else
       if(basis1%spec(minspecloc)%num.eq.1) return
    end if
    allocate(sav_trans(basis1%natom,3))


!!!-----------------------------------------------------------------------------
!!! if present, initialises tmp_wyckoff arrays
!!!-----------------------------------------------------------------------------
    if(present(wyck_check))then
       lwyckoff=.true.
    else
       lwyckoff=.false.
    end if


!!!-----------------------------------------------------------------------------
!!! Cycles through each atom in transformed basis and finds translation ...
!!! ... vector that maps it back onto the 1st atom in the original, ...
!!! ... untransformed, basis.
!!! Then tests this translation vector on all other atoms to see if it works ...
!!! ... as a translation vector for the symmetry.
!!!-----------------------------------------------------------------------------
    trloop: do iatom = 1, basis1%spec(minspecloc)%num
       ttrans(:) = 0._real32
       ttrans(1:3) = basis1%spec(minspecloc)%atom(1,1:3)-&
            basis2%spec(minspecloc)%atom(iatom,1:3)
       if(all(abs(ttrans(1:3)-anint(ttrans(1:3))).lt.tol_sym)) cycle trloop
       if(confine%l)then
          if(confine%laxis(confine%axis).and.&
               abs(ttrans(confine%axis)-nint(ttrans(confine%axis)))&
               .gt.tol_sym) cycle trloop
       end if
       itmp1 = 0
       sav_trans = 0._real32
       if(lwyckoff.and.ntrans+1.gt.size(wyck_check))then
          write(0,'("ERROR: error encountered in gldfnd")')
          write(0,'(2X,"Internal error in subroutine gldfnd in artemis__sym.f90")')
          write(0,'(2X,"ntrans is greater than wyck_check")')
          write(0,'(2X,"EXITING SUBROUTINE")')
          return
       end if
       trcyc: do ispec = 1, basis1%nspec
          samecount=0
          if(lwyckoff) wyck_check(ntrans+1)%spec(ispec)%atom(:) = 0
          atmcyc2: do jatom=1,basis1%spec(ispec)%num
             itmp1 = itmp1 + 1
             tmpbas(1:3) = basis2%spec(ispec)%atom(jatom,1:3) + ttrans(1:3)
             tmpbas(:) = tmpbas(:) - ceiling(tmpbas(:)-0.5_real32)
             atmcyc3: do katom=1,basis1%spec(ispec)%num
                !if(lwyckoff.and.&
                !     wyck_check(ntrans+1)%spec(ispec)%atom(katom).ne.0) &
                !     cycle atmcyc3
                diff = tmpbas(1:3) - basis1%spec(ispec)%atom(katom,1:3)
                do j=1,3
                   diff(j) = mod((diff(j)+100._real32),1.0)
                   if((abs(diff(j)-1._real32)).lt.(tol_sym)) diff(j) = 0._real32
                end do
                if(sqrt(dot_product(diff,diff)).lt.tol_sym)then
                   samecount = samecount + 1
                   sav_trans(itmp1,:) = basis1%spec(ispec)%atom(katom,1:3) - &
                        basis2%spec(ispec)%atom(jatom,1:3)
                   sav_trans(itmp1,:) = sav_trans(itmp1,:) - &
                        ceiling(sav_trans(itmp1,:)-0.5_real32)
                   if(lwyckoff) &
                        wyck_check(ntrans+1)%spec(ispec)%atom(jatom) = katom
                   cycle atmcyc2
                end if
             end do atmcyc3
             !cycle trloop
          end do atmcyc2
          if (samecount.ne.basis1%spec(ispec)%num) cycle trloop
       end do trcyc
!!!-----------------------------------------------------------------------------
!!! Cleans up succeeded translation vector
!!!-----------------------------------------------------------------------------
       do j = 1, 3
          itmp1 = maxloc(abs(sav_trans(:,j)),dim=1)
          ttrans(j) = sav_trans(itmp1,j)
          ttrans(j) = ttrans(j) - ceiling(ttrans(j)-0.5_real32)
       end do
!!!-----------------------------------------------------------------------------
!!! If axis is confined, removes all symmetries not confined to the axis plane
!!!-----------------------------------------------------------------------------
       if(confine%l)then
          if(confine%laxis(confine%axis).and.&
               abs(ttrans(confine%axis)-nint(ttrans(confine%axis)))&
               .gt.tol_sym) cycle trloop
       else
          do i = 1, 3
             if(confine%laxis(i))then
                if(abs(ttrans(confine%axis)-floor(ttrans(confine%axis)))&
                     .lt.tol_sym) cycle trloop
             end if
          end do
       end if
!!!-----------------------------------------------------------------------------
!!! Checks whether this translation has already been saved
!!!-----------------------------------------------------------------------------
       do i = 1, ntrans
          if(all(abs(ttrans(:)-trans(i,:)).lt.tol_sym)) cycle trloop
       end do
       ntrans = ntrans + 1
       trans(ntrans,1:3) = ttrans(1:3)
       if(confine%l) return
    end do trloop


    return
  end subroutine gldfnd
!!!#############################################################################


!###############################################################################
  subroutine gen_fundam_sym_matrices(grp, lat, tol_sym)
    !! Generate fundamental symmetry matrices for the 3D space groups
    implicit none

    ! Arguments
    type(sym_type), intent(inout) :: grp
    !! Instance of the symmetry container
    real(real32), dimension(3,3), intent(in) :: lat
    !! The lattice matrix
    real(real32), intent(in) :: tol_sym
    !! Tolerance for symmetry operations

    ! Local variables
    integer :: i, count, old_count, jsym
    real(real32) :: cosPi3,sinPi3,mcosPi3,msinPi3
    real(real32), dimension(3,3) :: inversion,invlat,tmat1
    real(real32), dimension(3,3,64) :: fundam_mat
    real(real32), dimension(3,3,128) :: tmp_store


    cosPi3 = 0.5_real32
    sinPi3 = sin(pi/3._real32)
    mcosPi3 = -cosPi3
    msinPi3 = -sinPi3


    fundam_mat(1:3,1:3,1)=transpose(reshape((/&
         1._real32,  0._real32,  0._real32,  0._real32,  1._real32,  0._real32,  0._real32,  0._real32,  1._real32 /),&
         shape(inversion)))

    fundam_mat(1:3,1:3,2)=transpose(reshape((/&
         -1._real32,  0._real32,  0._real32,  0._real32, -1._real32,  0._real32,  0._real32, 0._real32,  1._real32 /),&
         shape(inversion)))

    fundam_mat(1:3,1:3,3)=transpose(reshape((/&
         -1._real32,  0._real32,  0._real32,  0._real32,  1._real32,  0._real32,  0._real32, 0._real32, -1._real32 /),&
         shape(inversion)))

    fundam_mat(1:3,1:3,4)=transpose(reshape((/&
         1._real32,  0._real32,  0._real32,  0._real32, -1._real32,  0._real32,  0._real32,  0._real32, -1._real32 /),&
         shape(inversion)))

    fundam_mat(1:3,1:3,5)=transpose(reshape((/&
         0._real32,  1._real32,  0._real32,  1._real32,  0._real32,  0._real32,  0._real32,  0._real32, -1._real32 /),&
         shape(inversion)))

    fundam_mat(1:3,1:3,6)=transpose(reshape((/&
         0._real32, -1._real32,  0._real32,  -1._real32,  0._real32,  0._real32,  0._real32, 0._real32, -1._real32 /),&
         shape(inversion)))

    fundam_mat(1:3,1:3,7)=transpose(reshape((/&
         0._real32, -1._real32,  0._real32,  1._real32,  0._real32,  0._real32,  0._real32,  0._real32,  1._real32 /),&
         shape(inversion)))

    fundam_mat(1:3,1:3,8)=transpose(reshape((/&
         0._real32,  1._real32,  0._real32,  -1._real32,  0._real32,  0._real32,  0._real32, 0._real32,  1._real32 /),&
         shape(inversion)))

    fundam_mat(1:3,1:3,9)=transpose(reshape((/&
         0._real32,  0._real32,  1._real32,  0._real32, -1._real32,  0._real32,  1._real32,  0._real32,  0._real32 /),&
         shape(inversion)))

    fundam_mat(1:3,1:3,10)=transpose(reshape((/&
         0._real32,  0._real32, -1._real32,  0._real32, -1._real32,  0._real32,  -1._real32, 0._real32,  0._real32 /),&
         shape(inversion)))

    fundam_mat(1:3,1:3,11)=transpose(reshape((/&
         0._real32,  0._real32, -1._real32,   0._real32,  1._real32,  0._real32,  1._real32, 0._real32,  0._real32 /),&
         shape(inversion)))

    fundam_mat(1:3,1:3,12)=transpose(reshape((/&
         0._real32,  0._real32,  1._real32,  0._real32,  1._real32,  0._real32,  -1._real32, 0._real32,  0._real32 /),&
         shape(inversion)))

    fundam_mat(1:3,1:3,13)=transpose(reshape((/&
         -1._real32,  0._real32,  0._real32,  0._real32,  0._real32,  1._real32,  0._real32, 1._real32,  0._real32 /),&
         shape(inversion)))

    fundam_mat(1:3,1:3,14)=transpose(reshape((/&
         -1._real32,  0._real32,  0._real32,  0._real32,  0._real32, -1._real32,  0._real32, -1._real32,  0._real32 /),&
         shape(inversion)))

    fundam_mat(1:3,1:3,15)=transpose(reshape((/&
         1._real32,  0._real32,  0._real32,  0._real32,  0._real32, -1._real32,  0._real32,  1._real32,  0._real32 /),&
         shape(inversion)))

    fundam_mat(1:3,1:3,16)=transpose(reshape((/&
         1._real32,  0._real32,  0._real32,  0._real32,  0._real32,  1._real32,  0._real32, -1._real32,  0._real32/),&
         shape(inversion)))

    fundam_mat(1:3,1:3,17)=transpose(reshape((/&
         0._real32,  0._real32,  1._real32,  1._real32,  0._real32,  0._real32,  0._real32,  1._real32,  0._real32 /),&
         shape(inversion)))

    fundam_mat(1:3,1:3,18)=transpose(reshape((/&
         0._real32,  0._real32, -1._real32, -1._real32,  0._real32,  0._real32,  0._real32,  1._real32,  0._real32 /),&
         shape(inversion)))

    fundam_mat(1:3,1:3,19)=transpose(reshape((/&
         0._real32,  0._real32, -1._real32,  1._real32,  0._real32,  0._real32,  0._real32, -1._real32,  0._real32 /),&
         shape(inversion)))

    fundam_mat(1:3,1:3,20)=transpose(reshape((/&
         0._real32,  0._real32,  1._real32, -1._real32,  0._real32,  0._real32,  0._real32, -1._real32,  0._real32 /),&
         shape(inversion)))

    fundam_mat(1:3,1:3,21)=transpose(reshape((/&
         0._real32,  1._real32,  0._real32,  0._real32,  0._real32,  1._real32,  1._real32,  0._real32,  0._real32 /),&
         shape(inversion)))

    fundam_mat(1:3,1:3,22)=transpose(reshape((/&
         0._real32, -1._real32,  0._real32,  0._real32,  0._real32, -1._real32,  1._real32,  0._real32,  0._real32 /),&
         shape(inversion)))

    fundam_mat(1:3,1:3,23)=transpose(reshape((/&
         0._real32, -1._real32,  0._real32,  0._real32,  0._real32,  1._real32, -1._real32,  0._real32,  0._real32 /),&
         shape(inversion)))

    fundam_mat(1:3,1:3,24)=transpose(reshape((/&
         0._real32,  1._real32,  0._real32,  0._real32,  0._real32, -1._real32, -1._real32,  0._real32,  0._real32 /),&
         shape(inversion)))

    fundam_mat(1:3,1:3,25)=transpose(reshape((/&
         cosPi3,  sinPi3, 0._real32, msinPi3,  cosPi3, 0._real32, 0._real32, 0._real32,  1._real32 /),&
         shape(inversion)))

    fundam_mat(1:3,1:3,26)=transpose(reshape((/&
         cosPi3, msinPi3, 0._real32,  sinPi3,  cosPi3, 0._real32, 0._real32, 0._real32,  1._real32 /),&
         shape(inversion)))

    fundam_mat(1:3,1:3,27)=transpose(reshape((/&
         mcosPi3,  sinPi3, 0._real32, msinPi3, mcosPi3, 0._real32, 0._real32, 0._real32, 1._real32 /),&
         shape(inversion)))

    fundam_mat(1:3,1:3,28)=transpose(reshape((/&
         mcosPi3, msinPi3, 0._real32,  sinPi3, mcosPi3, 0._real32, 0._real32, 0._real32, 1._real32 /),&
         shape(inversion)))

    fundam_mat(1:3,1:3,29)=transpose(reshape((/&
         cosPi3, msinPi3, 0._real32, msinPi3, mcosPi3, 0._real32, 0._real32, 0._real32, -1._real32 /),&
         shape(inversion)))

    fundam_mat(1:3,1:3,30)=transpose(reshape((/&
         cosPi3,  sinPi3, 0._real32,  sinPi3, mcosPi3, 0._real32, 0._real32, 0._real32, -1._real32 /),&
         shape(inversion)))

    fundam_mat(1:3,1:3,31)=transpose(reshape((/&
         mcosPi3, msinPi3, 0._real32, msinPi3,  cosPi3, 0._real32, 0._real32, 0._real32, -1._real32 /),&
         shape(inversion)))

    fundam_mat(1:3,1:3,32)=transpose(reshape((/&
         mcosPi3,  sinPi3, 0._real32,  sinPi3,  cosPi3, 0._real32, 0._real32, 0._real32, -1._real32 /),&
         shape(inversion)))

    inversion(:3,:3)=transpose(reshape((/&
         -1._real32,  0._real32,  0._real32,   0._real32,  -1._real32,  0._real32,   0._real32,  0._real32,  -1._real32 /),&
         shape(inversion)))


    do i=1,32
       fundam_mat(:3,:3,i+32) = matmul(inversion,fundam_mat(:3,:3,i))
    end do


    grp%nsym=0
    invlat=inverse_3x3(lat)
    old_count = 0
    count = 0
    do i = 1, 64, 1
       call add_sym(grp, fundam_mat(:3,:3,i), lat, invlat, tol_sym, tmp_store, count)
       if(old_count.ne.count) then
          same_check1: do jsym = 1, count-1, 1
             if(all(abs(tmp_store(:3,:3,count)-tmp_store(:3,:3,jsym)).lt.tol_sym))then
               count = count - 1
               exit same_check1
             end if
          end do same_check1
       end if
       old_count = count
       call add_sym_tf(grp, fundam_mat(:3,:3,i), lat, invlat, tol_sym, tmp_store, count)
       if(old_count.ne.count) then
          same_check2: do jsym = 1, count-1, 1
             if(all(abs(tmp_store(:3,:3,count)-tmp_store(:3,:3,jsym)).lt.tol_sym))then
               count = count - 1
               exit same_check2
             end if
          end do same_check2
       end if
       old_count = count
    end do


    grp%nsym = count
    allocate(grp%sym(4,4,grp%nsym), source = 0._real32)
    grp%sym(4,4,:) = 1._real32
    grp%sym(:3,:3,:grp%nsym) = tmp_store(:3,:3,:grp%nsym)
    grp%nlatsym=grp%nsym

  end subroutine gen_fundam_sym_matrices
!###############################################################################


!###############################################################################
  subroutine mksym(grp, lat, tol_sym)
    !! Generate the symmetry operations for a given lattice
    implicit none

    ! Arguments
    type(sym_type), intent(inout) :: grp
    !! Instance of the symmetry container
    real(real32), dimension(3,3), intent(in) :: lat
    !! Lattice matrix
    real(real32), intent(in) :: tol_sym
    !! Tolerance for symmetry operations

    ! Local variables
    integer :: amin,bmin,cmin
    integer :: i,j,ia,ib,ic,n,count,irot,nrot,isym,jsym, old_count
    real(real32) :: tht,a,b,c
    real(real32), dimension(3,3) :: rotmat,refmat,invlat,tmat1
    real(real32), allocatable, dimension(:,:,:) :: tsym1,tmp_store
    logical, dimension(3) :: laxis


    if(grp%confine%l)then
       laxis = grp%confine%laxis
    else
       laxis = .not.grp%confine%laxis
    end if


!!!-----------------------------------------------------------------------------
!!! initialise values and symmetry matrix
!!!-----------------------------------------------------------------------------
    allocate(tsym1(4,4,50000))
    tsym1 = 0._real32
    tsym1(4,4,:) = 1._real32
    count = 0


!!!-----------------------------------------------------------------------------
!!! rotation plane perp to z (1=E,2=C2,3=C3,4=C4,5=C5,6=C6)
!!!-----------------------------------------------------------------------------
    if(laxis(3))then
       mksyml: do n=1,10
          count=count+1
          if(n.gt.6)then
             tht = -2._real32*pi/real(n-4) !=2*pi/(n-4)
          else
             tht = 2._real32*pi/real(n) !=2*pi/n          
          end if
          tsym1(1:3,1:3,count)=transpose(reshape((/&
               cos(tht) ,  sin(tht),   0._real32,&
               -sin(tht),  cos(tht),   0._real32,&
               0._real32     ,      0._real32,   1._real32/), shape(rotmat)))
          do i=1,3
             do j=1,3
                if(abs(tsym1(i,j,count)).lt.tol_sym) tsym1(i,j,count)=0._real32
             end do
          end do
       end do mksyml
       nrot=count
    end if


!!!-----------------------------------------------------------------------------
!!! rotation plane perp to x
!!!-----------------------------------------------------------------------------
    if(laxis(1))then
       philoop: do n=1,10
          if(n.gt.6)then
             tht = -2._real32*pi/real(n-4) !=2*pi/n
          else
             tht = 2._real32*pi/real(n) !=2*pi/n
          end if
          rotmat = transpose(reshape((/&
               1._real32,      0._real32,      0._real32,  &
               0._real32,  cos(tht),  sin(tht),&
               0._real32, -sin(tht),  cos(tht)/), shape(rotmat)))
          rot2: do irot = 1, nrot
             count = count + 1
             tsym1(1:3,1:3,count) = matmul(rotmat(1:3,1:3),tsym1(1:3,1:3,irot))
          end do rot2
       end do philoop
       nrot=count
    end if


!!!-----------------------------------------------------------------------------
!!! rotation plane perp to y
!!!-----------------------------------------------------------------------------
    if(laxis(2))then
       psiloop: do n=1,10
          if(n.gt.6)then
             tht = -2._real32*pi/real(n-4) !=2*pi/n 
          else
             tht = 2._real32*pi/real(n) !=2*pi/n 
          end if
          rotmat = transpose(reshape((/&
               cos(tht) ,  0._real32,  sin(tht),&
               0._real32     ,  1._real32,      0._real32,    &
               -sin(tht),  0._real32,  cos(tht)/), shape(rotmat)))
          rot3: do irot=1,nrot
             count = count + 1
             tsym1(1:3,1:3,count) = matmul(rotmat(1:3,1:3),tsym1(1:3,1:3,irot))
             where (abs(tsym1(1:3,1:3,count)).lt.tol_sym)
                tsym1(1:3,1:3,count) = 0._real32
             end where
          end do rot3
       end do psiloop
       nrot=count
    end if


!!!-----------------------------------------------------------------------------
!!! inversion (i), x plane mirror (v), y plane mirror (v), z plane mirror (h)
!!!-----------------------------------------------------------------------------
    amin=1;bmin=1;cmin=1
    if(grp%confine%lmirror)then
       if(laxis(1)) amin=2
       if(laxis(2)) bmin=2
       if(laxis(3)) cmin=2
    end if
    aloop: do ia=amin,2
       a=(-1._real32)**ia
       bloop: do ib=bmin,2
          b=(-1._real32)**ib
          cloop: do ic=cmin,2
             c=(-1._real32)**ic
             !           if((a*b*c).ne.(-1._real32)) cycle cloop
             refmat(1:3,1:3) = transpose(reshape((/&
                  a,     0._real32,  0._real32,&
                  0._real32,  b   ,  0._real32,&
                  0._real32,  0._real32,     c/), shape(rotmat)))
             refloop: do irot = 1, nrot
                count = count + 1
                tsym1(1:3,1:3,count) = matmul(refmat(1:3,1:3),tsym1(1:3,1:3,irot))
             end do refloop
          end do cloop
       end do bloop
    end do aloop
    grp%nsym = count


    if(grp%lmolec)then
       allocate(grp%sym(4,4,grp%nsym))
       grp%sym(:,:,:grp%nsym)=tsym1(:,:,:grp%nsym)
       deallocate(tsym1)
       return
    else
       invlat = inverse_3x3(lat)
    end if
!!!-----------------------------------------------------------------------------
!!! checks all made symmetries to see if they apply to the supplied lattice
!!!-----------------------------------------------------------------------------
    allocate(tmp_store(3,3,grp%nsym))
    count = 0
    do i = 1, grp%nsym, 1
       call add_sym_tf(grp, tsym1(:3,:3,i), lat, invlat, tol_sym, tmp_store, count)
       if(old_count.ne.count) then
          same_check2: do jsym = 1, count-1, 1
             if(all(abs(tmp_store(:3,:3,count)-tmp_store(:3,:3,jsym)).lt.tol_sym))then
               count = count - 1
               exit same_check2
             end if
          end do same_check2
       end if
       old_count = count
    end do

    grp%nsym = count
    deallocate(tsym1)
    allocate(grp%sym(4,4,grp%nsym), source = 0._real32)
    grp%sym(4,4,:) = 1._real32
    grp%sym(:3,:3,:grp%nsym)=tmp_store(:3,:3,:grp%nsym)
    deallocate(tmp_store)

    grp%nlatsym = grp%nsym
    
  end subroutine mksym
!###############################################################################


!###############################################################################
  subroutine generate_all_symmetries(grp, lat, tol_sym)
    !! Generate all possible symmetry operations for a given lattice
    implicit none

    ! Arguments
    type(sym_type), intent(inout) :: grp
    !! Instance of the symmetry container
    real(real32), dimension(3,3), intent(in) :: lat
    !! Lattice matrix
    real(real32), intent(in) :: tol_sym
    !! Tolerance for symmetry operations

    ! Local variables
    integer :: i, j, count, n
    !! Counters
    real(real32) :: angle
    !! Angle for rotation
    real(real32), dimension(3,3) :: invlat
    !! Inverse lattice matrix
    real(real32), dimension(3,3) :: smat, mirror, ident
    !! Symmetry matrices
    real(real32), allocatable :: symm_matrices(:,:,:)
    !! Symmetry matrices array
    real(real32), dimension(3) :: axis
    !! Axis of rotation


    allocate(symm_matrices(3,3,20000))
    count = 0

    ident = 0._real32
    ident(1,1) = 1._real32; ident(2,2) = 1._real32; ident(3,3) = 1._real32
    invlat = inverse_3x3(lat)

    ! Off-axis mirrors (diagonal planes)
    smat = ident
    smat(1,1) = 0._real32; smat(1,2) = 1._real32;
    smat(2,1) = 1._real32; smat(2,2) = 0._real32
    smat(3,3) = 1._real32
    call add_sym(grp, smat, lat, invlat, tol_sym, symm_matrices, count)
    smat = ident
    smat(1,1) = 0._real32; smat(1,2) = -1._real32
    smat(2,1) = -1._real32; smat(2,2) = 0._real32
    smat(3,3) = 1._real32
    call add_sym(grp, smat, lat, invlat, tol_sym, symm_matrices, count)

    ! Rotations around x, y, z axes (common n-fold: 2, 3, 4, 6)
    do j = 1, 8
       mirror = ident
       if(j.gt.5) then
          mirror = -ident
          mirror(j-5,j-5) = 1._real32
       elseif(j.gt.2) then
          mirror(j-2,j-2) = -1._real32
       elseif(j.eq.2)then
          mirror = -ident
       end if

      do i = 1, 3
         do n = 1, 10
            if(n.gt.6)then
               angle = -2._real32*pi/real(n-4, real32) !=2*pi/(n-4)
           else
               angle = 2._real32*pi/real(n, real32) !=2*pi/n          
           end if
               smat = rotation_matrix(i, angle)
               smat = matmul(smat, mirror)
               call add_sym(grp, smat, lat, invlat, tol_sym, symm_matrices, count)
         end do
      end do

      ! Rotations around body diagonals (e.g. [111], [110])
      do i = 1, 11
         axis = uvec([1._real32, 1._real32, 1._real32])
         if(i.eq.11)then
            axis = -axis
         elseif(i.gt.7) then
            axis = -axis
            axis(i-7) = 1._real32
         elseif(i.gt.4)then
            axis(i-4) = -1._real32
         elseif(i.gt.1) then
            axis(i-1) = 0._real32
         end if
         do n = 1, 10
            if(n.gt.6)then
               angle = -2._real32*pi/real(n-4, real32) !=2*pi/(n-4)
         else
               angle = 2._real32*pi/real(n, real32) !=2*pi/n          
         end if
            smat = rotate_about_axis(axis, angle)
            smat = matmul(smat, mirror)
            call add_sym(grp, smat, lat, invlat, tol_sym, symm_matrices, count)
         end do
      end do
   end do

    ! Trim to valid
    grp%nsym = count
    allocate(grp%sym(4,4,grp%nsym), source=0._real32)
    grp%sym(:3,:3,:) = symm_matrices(:3,:3,1:count)
    count = 0
    sym_check: do i = 1, grp%nsym
       do j = 1, count, 1
          if(all(abs(grp%sym(:3,:3,i)-symm_matrices(:3,:3,j)).lt.tol_sym)) then
             cycle sym_check
          end if
       end do
       count = count + 1
       symm_matrices(1:3,1:3,count) = grp%sym(:3,:3,i)
    end do sym_check
    grp%nsym = count
    deallocate(grp%sym)
    allocate(grp%sym(4,4,grp%nsym), source=0._real32)
    grp%sym(:3,:3,:) = symm_matrices(:3,:3,1:count)
    grp%sym(4,4,:) = 1._real32
    deallocate(symm_matrices)
    grp%nlatsym = grp%nsym

   contains

    function rotation_matrix(axis, angle) result(output)
      implicit none
      integer, intent(in) :: axis
      real(real32), intent(in) :: angle
      real(real32) :: output(3,3), c, s
      c = cos(angle); s = sin(angle)
      if (axis == 1) then
         output = reshape([1._real32,0._real32,0._real32, 0._real32,c,s, 0._real32,-s,c], [3,3])
      elseif (axis == 2) then
         output = reshape([c,0._real32,-s, 0._real32,1._real32,0._real32, s,0._real32,c], [3,3])
      else
         output = reshape([c,s,0._real32, -s,c,0._real32, 0._real32,0._real32,1._real32], [3,3])
      end if
    end function rotation_matrix

    function rotate_about_axis(ax, angle) result(output)
      implicit none
      real(real32), intent(in) :: ax(3), angle
      real(real32) :: output(3,3), c, s, v
      real(real32) :: x, y, z
      x = ax(1); y = ax(2); z = ax(3)
      c = cos(angle); s = sin(angle); v = 1 - c
      output(1,1) = x*x*v + c
      output(1,2) = x*y*v - z*s
      output(1,3) = x*z*v + y*s
      output(2,1) = y*x*v + z*s
      output(2,2) = y*y*v + c
      output(2,3) = y*z*v - x*s
      output(3,1) = z*x*v - y*s
      output(3,2) = z*y*v + x*s
      output(3,3) = z*z*v + c
    end function rotate_about_axis

  end subroutine generate_all_symmetries
!###############################################################################


!###############################################################################
  subroutine add_sym(grp, mat, lat, invlat, tol_sym, store, count)
    !! Add symmetry matrix to the store if valid
    implicit none

    ! Arguments
    type(sym_type), intent(in) :: grp
    !! Instance of symmetry container
    real(real32), dimension(3,3), intent(in) :: mat
    !! Symmetry matrix
    real(real32), dimension(3,3), intent(in) :: lat, invlat
    !! Lattice and inverse lattice matrices
    real(real32), intent(in) :: tol_sym
    !! Tolerance for symmetry check
    real(real32), intent(inout) :: store(:,:,:)
    !! Store for symmetry matrices
    integer, intent(inout) :: count
    !! Counter for number of valid symmetries

    if (is_valid_symmetry(grp, mat, tol_sym))then
       count = count + 1
       store(:3,:3,count) = mat
    end if
  end subroutine add_sym
!-------------------------------------------------------------------------------
  subroutine add_sym_tf(grp, mat, lat, invlat, tol_sym, store, count)
    !! Add the coordinate transformed symmetry matrix to the store if valid
    implicit none

    ! Arguments
    type(sym_type), intent(in) :: grp
    !! Instance of symmetry container
    real(real32), dimension(3,3), intent(in) :: mat
    !! Symmetry matrix
    real(real32), dimension(3,3), intent(in) :: lat, invlat
    !! Lattice and inverse lattice matrices
    real(real32), intent(in) :: tol_sym
    !! Tolerance for symmetry check
    real(real32), intent(inout) :: store(:,:,:)
    !! Store for symmetry matrices
    integer, intent(inout) :: count
    !! Counter for number of valid symmetries

    ! Local variables
    real(real32) :: t(3,3)
    !! Transformed symmetry operation

    ! t = matmul(invlat, matmul(mat, lat))
    t = matmul(lat, matmul(mat, invlat))
    if (is_valid_symmetry(grp, t, tol_sym))then
       count = count + 1
       store(:3,:3,count) = t
    end if
  end subroutine add_sym_tf
!-------------------------------------------------------------------------------
  function is_valid_symmetry(grp, mat, tol_sym) result(output)
    !! Check if the symmetry matrix is valid
    implicit none

    ! Arguments
    type(sym_type), intent(in) :: grp
    !! Instance of symmetry container
    real(real32), dimension(3,3), intent(in) :: mat
    !! Symmetry matrix
    real(real32), intent(in) :: tol_sym
    !! Tolerance for symmetry check
    logical :: output
    !! Result of the symmetry check

    ! Local variables
    integer :: i
    !! Loop index
    real(real32), dimension(3) :: compare_vec, input_vec
    !! Vectors for comparison

    output = &
         all(abs(mat - nint(mat)) .lt. tol_sym) .and. &
         abs(abs(det(mat)) - 1._real32) .lt. tol_sym

    if(grp%lmolec) then
       output = output .and. all(abs(mat) .lt. 1._real32 + tol_sym)
    end if
    do i = 1, 3
       if(grp%confine%lmirror)then
          input_vec = mat(i,:)
       else
          input_vec = abs(mat(:,i))
       end if
       if( ( grp%confine%l .and. grp%confine%laxis(i) ) .or. &
           ( &
                .not.grp%confine%l .and. &
                grp%confine%lmirror .and. &
                grp%confine%laxis(i) &
           ) &
       ) then
             compare_vec = 0._real32
             compare_vec(i) = 1._real32
             output = output .and. &
                  all(abs(input_vec - compare_vec) .lt. tol_sym)
       end if
    end do

  end function is_valid_symmetry
!###############################################################################
 

!!!#############################################################################
!!! returns the wyckoff atoms of a basis (closest to a defined location)
!!!#############################################################################
  function get_wyckoff_atoms_any(wyckoff) result(wyckoff_atoms)
    implicit none
    integer :: i,is,ia,isym,imin,itmp1
    integer :: nsym,nspec
    type(wyck_type) :: wyckoff_atoms
    integer, allocatable, dimension(:) :: ivtmp1

    type(wyck_type), dimension(:), intent(in) :: wyckoff


    nsym = size(wyckoff)
    nspec = size(wyckoff(1)%spec(:))
    allocate(wyckoff_atoms%spec(nspec))
    wyckoff_atoms%spec(:)%num = 0
    do is=1,nspec
       allocate(ivtmp1(size(wyckoff(1)%spec(is)%atom)))
       ivtmp1 = 0
       do ia=1,size(wyckoff(1)%spec(is)%atom)

          imin = wyckoff(1)%spec(is)%atom(ia)
          if(imin.eq.0)then
             write(0,'("ERROR: imin in get_wyckoff_atoms is zero!!!")')
             write(0,'("Exiting...")')
             stop
          end if
          sym_loop1: do isym=2,nsym
             if(wyckoff(isym)%spec(is)%atom(ia).eq.0) cycle sym_loop1
             if(wyckoff(isym)%spec(is)%atom(ia).lt.imin)&
                  imin = wyckoff(isym)%spec(is)%atom(ia)
          end do sym_loop1
          sym_loop2: do 
             itmp1 = minval( (/ (wyckoff(i)%spec(is)%atom(imin),i=1,nsym) /),&
                  mask=(/ (wyckoff(i)%spec(is)%atom(imin),i=1,nsym) /).gt.0 )
             if(itmp1.ne.imin)then
                imin=itmp1
             else
                exit sym_loop2
             end if
          end do sym_loop2

          if(.not.any(ivtmp1(:).eq.imin))then
             wyckoff_atoms%spec(is)%num = wyckoff_atoms%spec(is)%num+1
             ivtmp1(wyckoff_atoms%spec(is)%num) = imin
          end if

       end do
       allocate(wyckoff_atoms%spec(is)%atom(wyckoff_atoms%spec(is)%num))
       wyckoff_atoms%spec(is)%atom(:)=ivtmp1(:wyckoff_atoms%spec(is)%num)
       deallocate(ivtmp1)
    end do
    wyckoff_atoms%nwyck = sum(wyckoff_atoms%spec(:)%num)

    
  end function get_wyckoff_atoms_any
!!!-----------------------------------------------------------------------------
!!!-----------------------------------------------------------------------------
  function get_wyckoff_atoms_loc(wyckoff,lat,bas,loc) result(wyckoff_atoms)
    implicit none
    integer :: is,ia,isym,imin,itmp1
    integer :: nsym
    real(real32) :: dist
    logical :: lfound_closer
    type(wyck_type) :: wyckoff_atoms
    real(real32), dimension(3) :: diff
    real(real32), allocatable, dimension(:) :: dists
    integer, allocatable, dimension(:) :: ivtmp1

    type(basis_type), intent(in) :: bas
    real(real32), dimension(3), intent(in) :: loc
    type(wyck_type), dimension(:), intent(in) :: wyckoff
    real(real32), dimension(3,3), intent(in) :: lat


    nsym = size(wyckoff)
    allocate(wyckoff_atoms%spec(bas%nspec))
    wyckoff_atoms%spec(:)%num = 0
    do is=1,bas%nspec
       allocate(ivtmp1(size(wyckoff(1)%spec(is)%atom)))
       ivtmp1 = 0

       allocate(dists(bas%spec(is)%num))
       do ia=1,bas%spec(is)%num
          diff = loc - bas%spec(is)%atom(ia,:3)
          diff = diff - ceiling(diff - 0.5_real32)
          dists(ia) = modu(matmul(diff,lat))
       end do

       wyckoff_loop1: do ia=1,size(wyckoff(1)%spec(is)%atom)

          dist = huge(0._real32)
          imin = wyckoff(1)%spec(is)%atom(ia)
          sym_loop1: do isym=1,nsym
             if(wyckoff(isym)%spec(is)%atom(ia).eq.0) cycle sym_loop1
             
             if(dists(wyckoff(isym)%spec(is)%atom(ia)).lt.dist)then
                dist = dists(wyckoff(isym)%spec(is)%atom(ia))
                imin = wyckoff(isym)%spec(is)%atom(ia)
             end if
          end do sym_loop1
          if(any(ivtmp1(:).eq.imin)) cycle wyckoff_loop1

          sym_loop2: do
             lfound_closer = .false.
             sym_loop3: do isym=1,nsym
                if(wyckoff(isym)%spec(is)%atom(imin).eq.0) cycle sym_loop3
                if(wyckoff(isym)%spec(is)%atom(imin).eq.imin) cycle sym_loop3
                if(dists(wyckoff(isym)%spec(is)%atom(imin)).lt.dist)then
                   dist = dists(wyckoff(isym)%spec(is)%atom(imin))
                   itmp1 = wyckoff(isym)%spec(is)%atom(imin)
                   lfound_closer = .true.
                elseif(dists(wyckoff(isym)%spec(is)%atom(imin)).eq.dist)then
                   if(any(ivtmp1(:).eq.wyckoff(isym)%spec(is)%atom(imin)))then
                      dist = dists(wyckoff(isym)%spec(is)%atom(imin))
                      itmp1 = wyckoff(isym)%spec(is)%atom(imin)
                      lfound_closer = .true.
                   end if
                end if
             end do sym_loop3
             if(lfound_closer)then
                imin = itmp1
             else
                exit sym_loop2
             end if
          end do sym_loop2


          if(.not.any(ivtmp1(:).eq.imin))then
             wyckoff_atoms%spec(is)%num = wyckoff_atoms%spec(is)%num+1
             ivtmp1(wyckoff_atoms%spec(is)%num) = imin
          end if
          if(imin.eq.0)then
             write(0,'("ERROR: imin in get_wyckoff_atoms is zero!!!")')
             write(0,'("Exiting...")')
             stop
          end if

       end do wyckoff_loop1
       allocate(wyckoff_atoms%spec(is)%atom(wyckoff_atoms%spec(is)%num))
       wyckoff_atoms%spec(is)%atom(:)=ivtmp1(:wyckoff_atoms%spec(is)%num)
       deallocate(ivtmp1)
       deallocate(dists)
    end do
    wyckoff_atoms%nwyck = sum(wyckoff_atoms%spec(:)%num)

    
  end function get_wyckoff_atoms_loc
!!!#############################################################################


!!!#############################################################################
!!! find corresponding basis2 atoms that the supplied symmetry operation ...
!!! ... maps basis1 atoms onto.
!!! Basis2 is optional. If missing, it uses basis1 for the comparison
!!!#############################################################################
  function basis_map(sym,bas1,tmpbas2, tol_sym) result(bas_map)
    implicit none
    real(real32), dimension(4,4), intent(in) :: sym
    type(basis_type), intent(in) :: bas1
    type(basis_type), optional, intent(in) :: tmpbas2
    real(real32), intent(in), optional :: tol_sym

    integer :: j,ispec,iatom,jatom,dim
    type(basis_map_type) :: bas_map
    type(basis_type) :: bas2,tfbas
    real(real32), dimension(3) :: diff


!!!-----------------------------------------------------------------------------
!!! checks for optional arguments and assigns values if not present
!!!-----------------------------------------------------------------------------
    allocate(bas2%spec(bas1%nspec))
    dim=size(bas1%spec(1)%atom(1,:),dim=1)
    do ispec=1,bas1%nspec
       allocate(bas2%spec(ispec)%atom(bas1%spec(ispec)%num,dim))
    end do
    if(present(tmpbas2)) then
       bas2 = tmpbas2
    else
       bas2 = bas1
    end if


!!!-----------------------------------------------------------------------------
!!! sets up basis map
!!!-----------------------------------------------------------------------------
    allocate(bas_map%spec(bas1%nspec))
    do ispec=1,bas1%nspec
       allocate(bas_map%spec(ispec)%atom(bas1%spec(ispec)%num))
       bas_map%spec(ispec)%atom(:)=0
    end do
    allocate(tfbas%spec(bas1%nspec))
    do ispec=1,bas1%nspec
       allocate(tfbas%spec(ispec)%atom(bas1%spec(ispec)%num,4))
    end do


!!!-----------------------------------------------------------------------------
!!! apply symmetry operator to bas1
!!!-----------------------------------------------------------------------------
    do ispec=1,bas1%nspec
       do iatom=1,bas1%spec(ispec)%num
          tfbas%spec(ispec)%atom(iatom,1:3) = &
               matmul(bas1%spec(ispec)%atom(iatom,1:4),sym(1:4,1:3))
          do j=1,3
             tfbas%spec(ispec)%atom(iatom,j) = &
                  tfbas%spec(ispec)%atom(iatom,j) - &
                  ceiling(tfbas%spec(ispec)%atom(iatom,j) - 0.5_real32)
             bas2%spec(ispec)%atom(iatom,j) = &
                  bas2%spec(ispec)%atom(iatom,j) - &
                  ceiling(bas2%spec(ispec)%atom(iatom,j) - 0.5_real32)
          end do
       end do
    end do


!!!-----------------------------------------------------------------------------
!!! check whether transformed basis matches original basis
!!!-----------------------------------------------------------------------------
    spcheck2: do ispec=1,bas1%nspec
       diff=0._real32
       atmcheck2: do iatom=1,bas1%spec(ispec)%num
          atmcyc2: do jatom=1,bas1%spec(ispec)%num
             if(any(bas_map%spec(ispec)%atom(:).eq.jatom)) cycle atmcyc2
             diff = tfbas%spec(ispec)%atom(iatom,1:3) - &
                  bas2%spec(ispec)%atom(jatom,1:3)
             diff = diff - ceiling(diff - 0.5_real32)
             if(sqrt(dot_product(diff,diff)).lt.tol_sym)then
                bas_map%spec(ispec)%atom(iatom) = jatom
             end if
          end do atmcyc2
       end do atmcheck2
    end do spcheck2


    return
  end function basis_map
!!!#############################################################################

end module artemis__sym
