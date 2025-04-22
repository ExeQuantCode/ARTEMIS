!!!#############################################################################
!!! Code written by Ned Thaddeus Taylor and Francis Huw Davies
!!! Code part of the ARTEMIS group (Hepplestone research group).
!!! Think Hepplestone, think HRG.
!!!#############################################################################
!!!module contains symmetry-related functions and subroutines.
!!!module includes the following functions and subroutines:
!!! sym_setup         (calls mksym and allocates unallocated symmetry arrays)
!!! check_sym         (checks supplied symmetries against supplied basis or ...
!!!                    ... checks whether the two supplied bases match after ...
!!!                    ... applying symmetries)
!!! gldfnd            (output translations that maps two bases)
!!! mksym             (makes array of symmetries that apply to supplied lattice
!!! clone_grp         (clones ingrp to outgrp)
!!! symwrite          (output human-readable supplied transformation matrix)
!!! basis_map         (finds symmetry equivalent atoms in two bases based on ...
!!!                    ... the supplied transformation matrix)
!!! setup_ladder      (sets up rungs of the layer ladder)
!!!#############################################################################
module artemis__sym
  use artemis__constants,   only: real32, pi
  use artemis__misc,        only: sort2D
  use misc_linalg,          only: modu,inverse_3x3,det,gcd,gen_group,cross
  use artemis__geom_rw,     only: basis_type
  use artemis__geom_utils,            only: reducer, primitive_lat
  implicit none
  real(real32) :: tol_sym_default = 1.E-6_real32
  integer, allocatable, dimension(:) :: symops_compare

  interface get_wyckoff_atoms
     procedure get_wyckoff_atoms_any,get_wyckoff_atoms_loc
  end interface get_wyckoff_atoms


  private


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
  end type sym_type


  public :: sym_type
  public :: clone_grp
  public :: sym_setup,check_sym,gldfnd

  public :: get_primitive_cell
  
  public :: confine_type

  public :: basis_map_type, basis_map


!!!updated 2023/02/14


contains

!!!#############################################################################
!!! calls mksym and allocates symops and wyckoff arrays
!!!#############################################################################
  subroutine sym_setup(grp,lat,predefined,new_start,tol_sym)
    implicit none
    type(sym_type), intent(inout) :: grp
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
          if(allocated(grp%op)) deallocate(grp%op)
          if(allocated(grp%sym)) deallocate(grp%sym)
       end if
    end if

    predefined_ = .false.
    if(present(predefined)) predefined_ = predefined
    if(predefined_)then
       call gen_fundam_sym_matrices(grp, lat, tol_sym_)
    else
       call mksym(grp, lat, tol_sym_)
    end if

    if(allocated(symops_compare)) deallocate(symops_compare)
    grp%nsymop=0

    new_start_ = .true.
    if(present(new_start)) new_start_ = new_start
    if(new_start_.or.grp%end_idx.eq.0)then
       grp%end_idx = grp%nsym
    end if

  end subroutine sym_setup
!!!#############################################################################


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
    if(.not.allocated(grp%op))then
       allocate(grp%op(grp%nsym*minval(basis%spec(:)%num)))
       grp%op = 0
    end if
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
    allocate(tmpsav(grp%nsym*minval(basis%spec(:)%num),4,4))
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
            grp%sym(isym,1:4,1:4)
       !------------------------------------------------------------------------
       ! apply symmetry operator to basis
       !------------------------------------------------------------------------
       do ispec = 1, basis%nspec, 1
          do iatom = 1, basis%spec(ispec)%num, 1
             tfbas%spec(ispec)%atom(iatom,1:3) = &
                  matmul(basis%spec(ispec)%atom(iatom,1:4),grp%sym(isym,1:4,1:3))
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
          tmpsav(grp%nsymop,:,:) = grp%sym(isym,:,:)
          grp%op(grp%nsymop) = isym
          if(grp%nsymop.ne.0.and..not.check_all_sym_) exit symloop
       end if
       trans = 0._real32
       ntrans = 0
       !------------------------------------------------------------------------
       ! checks if translations are valid with the current symmetry operation
       !------------------------------------------------------------------------
       if(grp%lspace) then
          if(all(abs(grp%sym(isym,1:3,1:3)-ident).lt.tol_sym_))then
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
                         if(all(abs(trans(i,1:3)-tmpsav(jsym,4,1:3)).lt.&
                              tol_sym_)) cycle transloop
                         diff = trans(i,1:3) - tmpsav(jsym,4,1:3)
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
                tmpsav(grp%nsymop,:,:) = grp%sym(isym,:,:)
                tmpsav(grp%nsymop,4,1:3) = trans(i,:)
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
       allocate(grp%sym_save(grp%nsymop,4,4))
       grp%sym_save=0._real32
       grp%sym_save(:grp%nsymop,:,:)=tmpsav(:grp%nsymop,:,:)
       grp%sym_save(:,4,4)=1._real32
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


!!!#############################################################################
!!! builds an array of the symmetries that apply to the supplied lattice
!!!#############################################################################
  subroutine gen_fundam_sym_matrices(grp, lat, tol_sym)
    implicit none
    type(sym_type), intent(inout) :: grp
    real(real32), dimension(3,3), intent(in) :: lat
    real(real32), intent(in) :: tol_sym

    integer :: i
    real(real32) :: cosPi3,sinPi3,mcosPi3,msinPi3
    real(real32), dimension(3,3) :: inversion,invlat,tmat1
    real(real32), dimension(64,3,3) :: fundam_mat


    cosPi3 = 0.5_real32
    sinPi3 = sin(pi/3._real32)
    mcosPi3 = -cosPi3
    msinPi3 = -sinPi3


    fundam_mat(1,1:3,1:3)=transpose(reshape((/&
         1._real32,  0._real32,  0._real32,  0._real32,  1._real32,  0._real32,  0._real32,  0._real32,  1._real32 /),&
         shape(inversion)))

    fundam_mat(2,1:3,1:3)=transpose(reshape((/&
         -1._real32,  0._real32,  0._real32,  0._real32, -1._real32,  0._real32,  0._real32, 0._real32,  1._real32 /),&
         shape(inversion)))

    fundam_mat(3,1:3,1:3)=transpose(reshape((/&
         -1._real32,  0._real32,  0._real32,  0._real32,  1._real32,  0._real32,  0._real32, 0._real32, -1._real32 /),&
         shape(inversion)))

    fundam_mat(4,1:3,1:3)=transpose(reshape((/&
         1._real32,  0._real32,  0._real32,  0._real32, -1._real32,  0._real32,  0._real32,  0._real32, -1._real32 /),&
         shape(inversion)))

    fundam_mat(5,1:3,1:3)=transpose(reshape((/&
         0._real32,  1._real32,  0._real32,  1._real32,  0._real32,  0._real32,  0._real32,  0._real32, -1._real32 /),&
         shape(inversion)))

    fundam_mat(6,1:3,1:3)=transpose(reshape((/&
         0._real32, -1._real32,  0._real32,  -1._real32,  0._real32,  0._real32,  0._real32, 0._real32, -1._real32 /),&
         shape(inversion)))

    fundam_mat(7,1:3,1:3)=transpose(reshape((/&
         0._real32, -1._real32,  0._real32,  1._real32,  0._real32,  0._real32,  0._real32,  0._real32,  1._real32 /),&
         shape(inversion)))

    fundam_mat(8,1:3,1:3)=transpose(reshape((/&
         0._real32,  1._real32,  0._real32,  -1._real32,  0._real32,  0._real32,  0._real32, 0._real32,  1._real32 /),&
         shape(inversion)))

    fundam_mat(9,1:3,1:3)=transpose(reshape((/&
         0._real32,  0._real32,  1._real32,  0._real32, -1._real32,  0._real32,  1._real32,  0._real32,  0._real32 /),&
         shape(inversion)))

    fundam_mat(10,1:3,1:3)=transpose(reshape((/&
         0._real32,  0._real32, -1._real32,  0._real32, -1._real32,  0._real32,  -1._real32, 0._real32,  0._real32 /),&
         shape(inversion)))

    fundam_mat(11,1:3,1:3)=transpose(reshape((/&
         0._real32,  0._real32, -1._real32,   0._real32,  1._real32,  0._real32,  1._real32, 0._real32,  0._real32 /),&
         shape(inversion)))

    fundam_mat(12,1:3,1:3)=transpose(reshape((/&
         0._real32,  0._real32,  1._real32,  0._real32,  1._real32,  0._real32,  -1._real32, 0._real32,  0._real32 /),&
         shape(inversion)))

    fundam_mat(13,1:3,1:3)=transpose(reshape((/&
         -1._real32,  0._real32,  0._real32,  0._real32,  0._real32,  1._real32,  0._real32, 1._real32,  0._real32 /),&
         shape(inversion)))

    fundam_mat(14,1:3,1:3)=transpose(reshape((/&
         -1._real32,  0._real32,  0._real32,  0._real32,  0._real32, -1._real32,  0._real32, -1._real32,  0._real32 /),&
         shape(inversion)))

    fundam_mat(15,1:3,1:3)=transpose(reshape((/&
         1._real32,  0._real32,  0._real32,  0._real32,  0._real32, -1._real32,  0._real32,  1._real32,  0._real32 /),&
         shape(inversion)))

    fundam_mat(16,1:3,1:3)=transpose(reshape((/&
         1._real32,  0._real32,  0._real32,  0._real32,  0._real32,  1._real32,  0._real32, -1._real32,  0._real32/),&
         shape(inversion)))

    fundam_mat(17,1:3,1:3)=transpose(reshape((/&
         0._real32,  0._real32,  1._real32,  1._real32,  0._real32,  0._real32,  0._real32,  1._real32,  0._real32 /),&
         shape(inversion)))

    fundam_mat(18,1:3,1:3)=transpose(reshape((/&
         0._real32,  0._real32, -1._real32, -1._real32,  0._real32,  0._real32,  0._real32,  1._real32,  0._real32 /),&
         shape(inversion)))

    fundam_mat(19,1:3,1:3)=transpose(reshape((/&
         0._real32,  0._real32, -1._real32,  1._real32,  0._real32,  0._real32,  0._real32, -1._real32,  0._real32 /),&
         shape(inversion)))

    fundam_mat(20,1:3,1:3)=transpose(reshape((/&
         0._real32,  0._real32,  1._real32, -1._real32,  0._real32,  0._real32,  0._real32, -1._real32,  0._real32 /),&
         shape(inversion)))

    fundam_mat(21,1:3,1:3)=transpose(reshape((/&
         0._real32,  1._real32,  0._real32,  0._real32,  0._real32,  1._real32,  1._real32,  0._real32,  0._real32 /),&
         shape(inversion)))

    fundam_mat(22,1:3,1:3)=transpose(reshape((/&
         0._real32, -1._real32,  0._real32,  0._real32,  0._real32, -1._real32,  1._real32,  0._real32,  0._real32 /),&
         shape(inversion)))

    fundam_mat(23,1:3,1:3)=transpose(reshape((/&
         0._real32, -1._real32,  0._real32,  0._real32,  0._real32,  1._real32, -1._real32,  0._real32,  0._real32 /),&
         shape(inversion)))

    fundam_mat(24,1:3,1:3)=transpose(reshape((/&
         0._real32,  1._real32,  0._real32,  0._real32,  0._real32, -1._real32, -1._real32,  0._real32,  0._real32 /),&
         shape(inversion)))

    fundam_mat(25,1:3,1:3)=transpose(reshape((/&
         cosPi3,  sinPi3, 0._real32, msinPi3,  cosPi3, 0._real32, 0._real32, 0._real32,  1._real32 /),&
         shape(inversion)))

    fundam_mat(26,1:3,1:3)=transpose(reshape((/&
         cosPi3, msinPi3, 0._real32,  sinPi3,  cosPi3, 0._real32, 0._real32, 0._real32,  1._real32 /),&
         shape(inversion)))

    fundam_mat(27,1:3,1:3)=transpose(reshape((/&
         mcosPi3,  sinPi3, 0._real32, msinPi3, mcosPi3, 0._real32, 0._real32, 0._real32, 1._real32 /),&
         shape(inversion)))

    fundam_mat(28,1:3,1:3)=transpose(reshape((/&
         mcosPi3, msinPi3, 0._real32,  sinPi3, mcosPi3, 0._real32, 0._real32, 0._real32, 1._real32 /),&
         shape(inversion)))

    fundam_mat(29,1:3,1:3)=transpose(reshape((/&
         cosPi3, msinPi3, 0._real32, msinPi3, mcosPi3, 0._real32, 0._real32, 0._real32, -1._real32 /),&
         shape(inversion)))

    fundam_mat(30,1:3,1:3)=transpose(reshape((/&
         cosPi3,  sinPi3, 0._real32,  sinPi3, mcosPi3, 0._real32, 0._real32, 0._real32, -1._real32 /),&
         shape(inversion)))

    fundam_mat(31,1:3,1:3)=transpose(reshape((/&
         mcosPi3, msinPi3, 0._real32, msinPi3,  cosPi3, 0._real32, 0._real32, 0._real32, -1._real32 /),&
         shape(inversion)))

    fundam_mat(32,1:3,1:3)=transpose(reshape((/&
         mcosPi3,  sinPi3, 0._real32,  sinPi3,  cosPi3, 0._real32, 0._real32, 0._real32, -1._real32 /),&
         shape(inversion)))

    inversion(:3,:3)=transpose(reshape((/&
         -1._real32,  0._real32,  0._real32,   0._real32,  -1._real32,  0._real32,   0._real32,  0._real32,  -1._real32 /),&
         shape(inversion)))


    do i=1,32
       fundam_mat(i+32,:3,:3) = matmul(inversion,fundam_mat(i,:3,:3))
    end do


    grp%nsym=0
    invlat=inverse_3x3(lat)
    do i=1,64
       tmat1=matmul(lat,fundam_mat(i,:3,:3))
       tmat1=matmul(tmat1,(invlat))
       !! ensure that the matrix preserves size of 1
       !! this is likely redundant
       if(abs(abs(det(tmat1))-1._real32).gt.tol_sym) cycle
       if(all(abs(tmat1-nint(tmat1)).le.tol_sym))then
          grp%nsym=grp%nsym+1
          fundam_mat(grp%nsym,:,:)=fundam_mat(i,:,:)
       end if
    end do


    allocate(grp%sym(grp%nsym,4,4))
    grp%sym(:,:,:)=0._real32
    grp%sym(:,4,4)=1._real32
    grp%sym(:grp%nsym,:3,:3)=fundam_mat(:grp%nsym,:3,:3)
    grp%nlatsym=grp%nsym


    !! REDUCE THIS SET BY DOING LTL^-1 AND JUST CHECK IF ANY BECOME NON-ZERO
    !! IF ONE DOES, SCRAP IT
    !! IF ONE DOESN'T, SAVE THE ORIGINAL (NOT THE NEWLY CREATED ONE)


  end subroutine gen_fundam_sym_matrices
!!!#############################################################################


!!!#############################################################################
!!! builds an array of the symmetries that apply to the supplied lattice
!!!#############################################################################
  subroutine mksym(grp, inlat, tol_sym)
    implicit none
    type(sym_type), intent(inout) :: grp
    real(real32), dimension(3,3), intent(in) :: inlat
    real(real32), intent(in) :: tol_sym

    integer :: amin,bmin,cmin
    integer :: i,j,ia,ib,ic,n,count,irot,nrot,isym,jsym
    real(real32) :: tht,a,b,c
    real(real32), dimension(3,3) :: rotmat,refmat,lat,invlat,tmat1
    real(real32), allocatable, dimension(:,:,:) :: tsym1,tsym2
    logical, dimension(3) :: laxis


    if(grp%confine%l)then
       laxis = grp%confine%laxis
    else
       laxis = .not.grp%confine%laxis
    end if


!!!-----------------------------------------------------------------------------
!!! set up inverse lattice
!!!-----------------------------------------------------------------------------
    lat = inlat
    if(grp%lmolec)then
       invlat = 0._real32
       lat    = 0._real32
    else
       invlat = inverse_3x3(lat)
    end if


!!!-----------------------------------------------------------------------------
!!! initialise values and symmetry matrix
!!!-----------------------------------------------------------------------------
    allocate(tsym1(50000,4,4))
    tsym1 = 0._real32
    tsym1(:,4,4)=1._real32
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
          tsym1(count,1:3,1:3)=transpose(reshape((/&
               cos(tht) ,  sin(tht),   0._real32,&
               -sin(tht),  cos(tht),   0._real32,&
               0._real32     ,      0._real32,   1._real32/), shape(rotmat)))
          do i=1,3
             do j=1,3
                if(abs(tsym1(count,i,j)).lt.tol_sym) tsym1(count,i,j)=0._real32
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
          rotmat=transpose(reshape((/&
               1._real32,      0._real32,      0._real32,  &
               0._real32,  cos(tht),  sin(tht),&
               0._real32, -sin(tht),  cos(tht)/), shape(rotmat)))
          rot2: do irot=1,nrot
             count=count+1
             tsym1(count,1:3,1:3)=matmul(rotmat(1:3,1:3),tsym1(irot,1:3,1:3))
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
          rotmat=transpose(reshape((/&
               cos(tht) ,  0._real32,  sin(tht),&
               0._real32     ,  1._real32,      0._real32,    &
               -sin(tht),  0._real32,  cos(tht)/), shape(rotmat)))
          rot3: do irot=1,nrot
             count=count+1
             tsym1(count,1:3,1:3)=matmul(rotmat(1:3,1:3),tsym1(irot,1:3,1:3))
             do i=1,3
                do j=1,3
                   if(abs(tsym1(count,i,j)).lt.tol_sym) tsym1(count,i,j)=0._real32
                end do
             end do
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
             refmat(1:3,1:3)=transpose(reshape((/&
                  a,     0._real32,  0._real32,&
                  0._real32,  b   ,  0._real32,&
                  0._real32,  0._real32,     c/), shape(rotmat)))
             refloop: do irot=1,nrot
                count=count+1
                tsym1(count,1:3,1:3)=matmul(refmat(1:3,1:3),tsym1(irot,1:3,1:3))
             end do refloop
          end do cloop
       end do bloop
    end do aloop
    grp%nsym=count


    if(grp%lmolec)then
       allocate(grp%sym(grp%nsym,4,4))
       grp%sym(:grp%nsym,:,:)=tsym1(:grp%nsym,:,:)
       deallocate(tsym1)
       return
    end if
    !! best so far
    !     sym(isym,1:3,1:3)=matmul(transpose(lat),sym(isym,1:3,1:3))
    !     sym(isym,1:3,1:3)=matmul(sym(isym,1:3,1:3),(invlat))
!!!-----------------------------------------------------------------------------
!!! checks all made symmetries to see if they apply to the supplied lattice
!!!-----------------------------------------------------------------------------
    allocate(tsym2(grp%nsym,4,4))
    tsym2=0._real32
    tsym2(:,4,4)=1._real32
    count=0
    samecheck: do isym=1,grp%nsym
       tmat1 = matmul((invlat),tsym1(isym,:3,:3))
       tmat1 = matmul(tmat1,(lat))
       do i=1,3
          do j=1,3
             if(abs(tmat1(i,j)).lt.tol_sym) tmat1(i,j)=0._real32
             if(abs(1._real32-abs(tmat1(i,j))).lt.tol_sym) &
                  tmat1(i,j)=sign(1._real32,tmat1(i,j))
          end do
       end do
       !!-----------------------------------------------------------------------
       !! Precautionary measure
       if(all(abs(tmat1).lt.tol_sym)) cycle samecheck
       if(abs(abs(det(tmat1))-1._real32).gt.tol_sym) cycle samecheck
       !!-----------------------------------------------------------------------
       if(.not.all(abs(tmat1-nint(tmat1)).lt.tol_sym)) cycle samecheck
       do jsym = 1, count, 1
          if(all(abs(tmat1-tsym2(jsym,:3,:3)).lt.tol_sym)) cycle samecheck
       end do
       count = count + 1
       tsym2(count,:3,:3) = tmat1
    end do samecheck
    grp%nsym=count
    deallocate(tsym1)
    allocate(grp%sym(grp%nsym,4,4))
    grp%sym(:grp%nsym,:4,:4)=tsym2(:grp%nsym,:4,:4)
    deallocate(tsym2)

    grp%nlatsym=grp%nsym


    return
  end subroutine mksym
!!!#############################################################################


!!!#############################################################################
!!! clone ingrp to outgrp
!!!#############################################################################
  subroutine clone_grp(from, to)
    implicit none
    type(sym_type), intent(in) :: from
    type(sym_type), intent(out) :: to
    
    
    if(allocated(from%op)) allocate(to%op(size(from%op)))
    if(allocated(from%sym)) allocate(to%sym(size(from%sym,dim=1),4,4))
    if(allocated(from%sym_save)) allocate(to%sym_save(size(from%sym_save,dim=1),4,4))
    to = from

  end subroutine clone_grp
!!!#############################################################################
 

!!!#############################################################################
!!! returns the primitive cell from a supercell
!!!#############################################################################
  subroutine get_primitive_cell(basis, tol_sym)
    implicit none
    type(basis_type), intent(inout) :: basis
    real(real32), intent(in), optional :: tol_sym

    integer :: is,ia,ja,i,j,k,itmp1
    integer :: ntrans,len
    real(real32) :: scale,proj,dtmp1
    real(real32) :: tol_sym_
    type(confine_type) :: confine
    real(real32), dimension(3,3) :: dmat1,invlat
    real(real32), allocatable, dimension(:,:) :: trans,atom_store
    

    
    !!-----------------------------------------------------------------------
    !! Allocate and initialise
    !!-----------------------------------------------------------------------
    tol_sym_ = tol_sym_default
    if(present(tol_sym)) tol_sym_ = tol_sym
    ntrans = 0
    dmat1=0._real32
    allocate(trans(minval(basis%spec(:)%num+2),3)); trans=0._real32

    
    !!-----------------------------------------------------------------------
    !! Find the translation vectors in the cell
    !!-----------------------------------------------------------------------
    call gldfnd(confine,basis,basis,trans,ntrans,tol_sym,.false.)
    len=size(basis%spec(1)%atom,dim=2)

    
    !!-----------------------------------------------------------------------
    !! For each translation, reduce the basis
    !!-----------------------------------------------------------------------
    if(ntrans.ge.1)then
       do i=ntrans+1,ntrans+3
          trans(i,:)=0._real32
          trans(i,i-ntrans)=1._real32
       end do
       !  trans=matmul(trans(1:ntrans,1:3),basis%lat)
       call sort2D( [ trans(1:ntrans+3,:) ] ,ntrans+3)
       !! for each lattice vector, determine the shortest translation ...
       !! ... vector that has a non-zero projection along that lattice vector.
       do i=1,3
          proj=1.D2
          trans_loop: do j=1,ntrans+3
             dtmp1 = dot_product(trans(j,:),trans(ntrans+i,:))
             if(dtmp1.lt.tol_sym) cycle trans_loop

             do k=1,i-1,1
                if(modu(abs(cross( [ trans(j,:) ], [ dmat1(k,:) ]))).lt.1.E-8_real32) cycle trans_loop
             end do

             dtmp1 = modu( [ trans(j,:) ] )
             if(dtmp1.lt.proj)then
                proj=dtmp1
                dmat1(i,:) = trans(j,:)
                trans(j,:) = 0._real32
             end if
          end do trans_loop
       end do
       !dmat1=trans(1:3,1:3)
       scale=det(dmat1)
       dmat1=matmul(dmat1,basis%lat)
       invlat=inverse_3x3(dmat1)
       do is=1,basis%nspec
          itmp1=0
          allocate(atom_store(nint(scale*basis%spec(is)%num),len))
          atcheck: do ia=1,basis%spec(is)%num
             !!-----------------------------------------------------------------
             !! Reduce the basis
             !!-----------------------------------------------------------------
             basis%spec(is)%atom(ia,1:3)=&
                  matmul(basis%spec(is)%atom(ia,1:3),basis%lat(1:3,1:3))
             basis%spec(is)%atom(ia,1:3)=&
                  matmul(transpose(invlat(1:3,1:3)),basis%spec(is)%atom(ia,1:3))
             do j=1,3
                basis%spec(is)%atom(ia,j)=&
                     basis%spec(is)%atom(ia,j)-floor(basis%spec(is)%atom(ia,j))
                if(basis%spec(is)%atom(ia,j).gt.1._real32-tol_sym) &
                     basis%spec(is)%atom(ia,j)=0._real32
             end do
             !!-----------------------------------------------------------------
             !! Check for duplicates in the cell
             !!-----------------------------------------------------------------
             do ja=1, itmp1
                if(all(abs(basis%spec(is)%atom(ia,1:3)-atom_store(ja,1:3)).lt.&
                     (/tol_sym,tol_sym,tol_sym/))) cycle atcheck
             end do
             itmp1=itmp1+1
             atom_store(itmp1,:)=basis%spec(is)%atom(ia,:)
             !!-----------------------------------------------------------------
             !! Check to ensure correct number of atoms remain after reduction
             !!-----------------------------------------------------------------
             if(itmp1.gt.size(atom_store,dim=1))then
                write(0,*) "ERROR! Primitive cell subroutine retained too &
                     &many atoms from supercell!", itmp1, size(atom_store,dim=1)
                call exit()
             end if
             !!-----------------------------------------------------------------
          end do atcheck
          deallocate(basis%spec(is)%atom)
          call move_alloc(atom_store,basis%spec(is)%atom)
          basis%spec(is)%num=size(basis%spec(is)%atom,dim=1)
          !deallocate(atom_store)
       end do
       !!-----------------------------------------------------------------------
       !! Reduce the lattice
       !!-----------------------------------------------------------------------
       basis%natom=sum(basis%spec(:)%num)
       basis%lat=dmat1
    end if

    
    !!-----------------------------------------------------------------------
    !! Reduce the lattice to symmetry definition
    !!-----------------------------------------------------------------------
    call reducer(basis)
    !! next line necessary as FCC and BCC do not conform to Niggli reduced ...
    !! ... cell definitions.
    basis%lat = primitive_lat(basis%lat)


    
  end subroutine get_primitive_cell
!!!#############################################################################


!!!#############################################################################
!!! takes in transformation matrix and outputs its (x,y,z) definition
!!!#############################################################################
  subroutine symwrite (sym,symchar)
    implicit none
    integer :: i,j,nt,nr,div
    real(real32), dimension(4,4) :: sym
    character(1024) :: symchar
    character(2) :: rm,c
    character(1), dimension(3) :: xyz

    xyz(1)="x";xyz(2)="y";xyz(3)="z"
    symchar=""
    do i=1,3
       select case (nint(100*sym(4,i)))
       case(0)
       case default
          div=abs(gcd(nint(100*sym(4,i)),100))
          write(symchar,'(A,I0,"aa",I0)') trim(symchar),nint(100*sym(4,i))/div,100/div
       end select

       do j=1,3
          select case (int(sym(j,i)))
          case(0)
             cycle
          case(1)
             c=""
          case default
             write(c,"(I2)") int(sym(j,i))
          end select
          symchar=trim(symchar) //"+"//trim(adjustl(c(1:1)))//xyz(j)
       end do
       if(i.ne.3) symchar=trim(symchar) //","
    end do

    rm="+-"
    nt=len_trim(symchar) ; nr=len_trim(symchar)
    remove: do
       i=index(symchar,trim(adjustl(rm)))
       if(i.eq.0) exit remove
       symchar = symchar(:i-1) //symchar(i+1:nt)
    end do remove

    rm=",+"
    nt=len_trim(symchar) ; nr=len_trim(symchar)
    remove2: do
       i=index(symchar,trim(adjustl(rm)))
       if(i.eq.0) exit remove2
       symchar = symchar(:i) //symchar(i+2:nt)
    end do remove2
    if(symchar(:1).eq."+") symchar=symchar(2:)

    rm="aa"
    nt=len_trim(symchar) ; nr=len_trim(symchar)
    remove3: do
       i=index(symchar,trim(adjustl(rm)))
       if(i.eq.0) exit remove3
       symchar = symchar(:i-1) //"/"//symchar(i+2:nt)
    end do remove3


    symchar = "("//trim(adjustl(symchar))//")"
    write(77,*) trim(adjustl(symchar))

  end subroutine symwrite
!!!#############################################################################


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
