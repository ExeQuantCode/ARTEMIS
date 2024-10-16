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
!!! get_terminations  (finds all possible terminations along an axis)
!!! print_terminations (prints the terminations to individual files)
!!!#############################################################################
module mod_sym
  use constants,   only: pi
  use misc,        only: sort1D,sort2D,sort_col,set
  use misc_linalg, only: modu,inverse_3x3,det,gcd,gen_group,cross,uvec
  use rw_geom,     only: bas_type,geom_write
  use edit_geom,   only: transformer,vacuumer,set_vacuum,shifter,&
       clone_bas,get_closest_atom,ortho_axis,reducer,primitive_lat,get_min_dist
  implicit none
  integer :: ierror_sym=0
  integer :: s_start=1,s_end=0
  double precision :: tol_sym=5.D-5
  character(1) :: verb_sym="n"
  integer, allocatable, dimension(:) :: symops_compare
  double precision, allocatable, dimension(:,:,:) :: savsym

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
  type basmap_type
     type(spcmap_type), allocatable, dimension(:) :: spec
  end type basmap_type

  type term_type
     !double precision :: add
     double precision :: hmin
     double precision :: hmax
     integer :: natom
     integer :: nstep
     double precision, allocatable, dimension(:) :: ladder
  end type term_type

  type term_arr_type
     integer :: nterm,axis,nstep
     double precision :: tol
     logical :: lmirror=.false.
     type(term_type), allocatable, dimension(:) :: arr
  end type term_arr_type


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
     integer :: nsym,nlatsym,nsymop,npntop
     logical :: lspace=.true.
     logical :: lmolec=.false.
     integer, allocatable, dimension(:) :: op
     double precision, allocatable, dimension(:,:,:) :: sym
     type(confine_type) :: confine
  end type sym_type


  public :: set_symmetry_tolerance
  public :: ierror_sym,s_start,s_end
  public :: sym_type
  public :: clone_grp
  public :: sym_setup,check_sym,gldfnd

  public :: get_primitive_cell
  
  public :: term_arr_type,confine_type
  public :: get_terminations

  public :: basmap_type,basis_map

  public :: wyck_type
  public :: get_wyckoff_atoms

  public :: symops_compare


!!!updated 2023/02/14


contains
!!!#############################################################################
!!! redefines the symmetry tolerance/precision
!!!#############################################################################
  subroutine set_symmetry_tolerance(tolerance)
    implicit none
    double precision, optional, intent(in) :: tolerance

    if(present(tolerance))then
       tol_sym = tolerance
    else
       tol_sym = 1.D-6
    end if

  end subroutine set_symmetry_tolerance
!!!#############################################################################


!!!#############################################################################
!!! calls mksym and allocates symops and wyckoff arrays
!!!#############################################################################
  subroutine sym_setup(grp,lat,predefined,new_start,tolerance)
    implicit none
    logical :: lpresent

    type(sym_type) :: grp

    double precision, dimension(3,3), intent(in) :: lat
    double precision, optional, intent(in) :: tolerance
    logical, optional, intent(in) :: predefined,new_start


    if(present(tolerance)) call set_symmetry_tolerance(tolerance)
    if(present(new_start))then
       if(new_start)then
          if(allocated(grp%op)) deallocate(grp%op)
          if(allocated(grp%sym)) deallocate(grp%sym)
       end if
    end if

    if(present(predefined))then
       if(predefined)then
          call gen_fundam_sym_matrices(grp,lat)
          goto 10
       end if
    end if
    call mksym(grp,lat)

10  if(allocated(savsym)) deallocate(savsym)
    if(allocated(symops_compare)) deallocate(symops_compare)
    grp%nsymop=0

    lpresent=.false.
    if(present(new_start))then
       if(new_start) lpresent=.true.
    end if
    if(.not.present(new_start).or.lpresent.or.s_end.eq.0)then
       s_end=grp%nsym
    end if


    return
  end subroutine sym_setup
!!!#############################################################################


!!!#############################################################################
!!! builds an array of the symmetries that apply to the supplied lattice
!!!#############################################################################
!!! tfbas   : transformed basis
!!!#############################################################################
  subroutine check_sym(grp,bas1,iperm,tmpbas2,wyckoff,lsave,lat,loc,lcheck_all)
    implicit none
    integer :: i,j,k,iatom,jatom,ispec,itmp1
    integer :: is,isym,jsym,count,ntrans
    integer :: samecount,oldnpntop
    logical :: lpresent,lsaving,lwyckoff,ltransformed
    type(bas_type) :: bas2,tfbas
    double precision, dimension(3) :: diff
    double precision, dimension(3,3) :: ident
    type(wyck_type), allocatable, dimension(:) :: wyck_check
    double precision, allocatable, dimension(:,:) :: trans
    double precision, allocatable, dimension(:,:,:) :: tmpsav

    type(bas_type), intent(in) :: bas1
    type(sym_type), intent(inout) :: grp

    integer, optional, intent(in) :: iperm
    logical, optional, intent(in) :: lsave,lcheck_all
    type(bas_type), optional, intent(in) :: tmpbas2
    type(wyck_type), optional, intent(inout) :: wyckoff
    double precision, dimension(3), optional, intent(in) :: loc
    double precision, dimension(3,3), optional, intent(in) :: lat


204 format(4(F11.6),/,4(F11.6),/,4(F11.6),/,4(F11.6))


!!!-----------------------------------------------------------------------------
!!! allocated grp%op
!!!-----------------------------------------------------------------------------
    if(allocated(grp%op)) deallocate(grp%op)
    allocate(grp%op(grp%nsym*minval(bas1%spec(:)%num)))
    grp%op = 0
    if(present(lsave))then
       lsaving = lsave
    else
       lsaving = .false.
    end if


!!!-----------------------------------------------------------------------------
!!! checks for optional arguments and assigns values if not present
!!!-----------------------------------------------------------------------------
    if(present(tmpbas2)) then
       bas2 = tmpbas2
       if(present(lcheck_all))then
          lpresent = .not.lcheck_all
       else
          lpresent = .true.
       end if
    else
       bas2 = bas1
       lpresent = .false.
    end if
    allocate(tmpsav(grp%nsym*minval(bas1%spec(:)%num),4,4))
    itmp1 = maxval(bas1%spec(:)%num)


!!!-----------------------------------------------------------------------------
!!! initialises variables
!!!-----------------------------------------------------------------------------
    allocate(trans(minval(bas1%spec(:)%num+2),3)); trans = 0.D0
    allocate(tfbas%spec(bas1%nspec))
    itmp1 = size(bas1%spec(1)%atom(1,:),dim=1)
    do is=1,bas1%nspec
       allocate(tfbas%spec(is)%atom(bas1%spec(is)%num,itmp1))
    end do
    grp%nsymop = 0
    grp%npntop = 0


!!!-----------------------------------------------------------------------------
!!! if present, initialises wyckoff arrays
!!!-----------------------------------------------------------------------------
    allocate(wyck_check(grp%nsym*minval(bas1%spec(:)%num)))
    do isym=1,grp%nsym*minval(bas1%spec(:)%num)
       allocate(wyck_check(isym)%spec(bas1%nspec))
       do ispec=1,bas1%nspec
          allocate(wyck_check(isym)%spec(ispec)%atom(bas1%spec(ispec)%num))
          wyck_check(isym)%spec(ispec)%atom = 0
       end do
    end do
    if(present(wyckoff))then
       lwyckoff = .true.
       if(allocated(wyckoff%spec)) deallocate(wyckoff%spec)
       wyckoff%nwyck = 0
       allocate(wyckoff%spec(bas1%nspec))
       do ispec=1,bas1%nspec
          wyckoff%spec(ispec)%num = 0
          wyckoff%spec(ispec)%name = ""
          allocate(wyckoff%spec(ispec)%atom(bas1%spec(ispec)%num))
          do iatom=1,bas1%spec(ispec)%num
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
    ident = 0.D0
    do i=1,3
       ident(i,i) = 1.D0
    end do


!!!-----------------------------------------------------------------------------
!!! applying symmetries to basis to see if the basis conforms to any of them
!!!-----------------------------------------------------------------------------
    itmp1 = 1
    symloop: do isym=s_start,s_end
       if(verb_sym.eq.'d') write(77,*) isym !,a,b,c
       if(verb_sym.eq.'d') write(77,204) grp%sym(isym,1:4,1:4)
       if(ierror_sym.eq.2.or.ierror_sym.eq.3) write(77,204)  &
            grp%sym(isym,1:4,1:4)
       !------------------------------------------------------------------------
       ! apply symmetry operator to basis
       !------------------------------------------------------------------------
       do ispec=1,bas1%nspec
          do iatom=1,bas1%spec(ispec)%num
             tfbas%spec(ispec)%atom(iatom,1:3) = &
                  matmul(bas1%spec(ispec)%atom(iatom,1:4),grp%sym(isym,1:4,1:3))
             do j=1,3
                tfbas%spec(ispec)%atom(iatom,j) = &
                     tfbas%spec(ispec)%atom(iatom,j) - &
                     ceiling(tfbas%spec(ispec)%atom(iatom,j)-0.5D0)
             end do
          end do
       end do
       !------------------------------------------------------------------------
       ! check whether transformed basis matches original basis
       !------------------------------------------------------------------------
       count=0
       spcheck: do ispec=1,bas1%nspec
          diff = 0.D0
          samecount = 0
          wyck_check(itmp1)%spec(ispec)%atom = 0
          atmcheck: do iatom=1,bas1%spec(ispec)%num
             atmcyc: do jatom=1,bas1%spec(ispec)%num
                !if(wyck_check(itmp1)%spec(ispec)%atom(jatom).ne.0) cycle atmcyc
                diff = tfbas%spec(ispec)%atom(iatom,1:3) - &
                     bas2%spec(ispec)%atom(jatom,1:3)
                diff(:) = diff(:) - floor(diff(:))
                where((abs(diff(:)-1.D0)).lt.(tol_sym))
                   diff(:)=0.D0
                end where
                if(sqrt(dot_product(diff,diff)).lt.tol_sym)then
                   samecount = samecount + 1
                   wyck_check(itmp1)%spec(ispec)%atom(iatom) = jatom
                end if
                if((iatom.eq.bas1%spec(ispec)%num).and.&
                     (jatom.eq.bas1%spec(ispec)%num))then
                   if (samecount.ne.bas1%spec(ispec)%num) goto 10
                end if
             end do atmcyc
             count = count + samecount
          end do atmcheck
          if(samecount.ne.bas1%spec(ispec)%num) goto 10
       end do spcheck
       grp%npntop = grp%npntop + 1
       grp%nsymop = grp%nsymop + 1
       itmp1 = grp%nsymop + 1
       tmpsav(grp%nsymop,:,:) = grp%sym(isym,:,:)
       grp%op(grp%nsymop) = isym
       if(grp%nsymop.ne.0.and.lpresent) exit symloop
10     trans = 0.D0
       ntrans = 0
       !------------------------------------------------------------------------
       ! checks if translations are valid with the current symmetry operation
       !------------------------------------------------------------------------
       if(grp%lspace) then
          if(all(abs(grp%sym(isym,1:3,1:3)-ident).lt.tol_sym))then
             ltransformed=.false.
          else
             ltransformed=.true.
          end if
          call gldfnd(grp%confine,&
               bas2,tfbas,&
               trans,ntrans,&
               transformed=ltransformed,&
               wyck_check=wyck_check(itmp1:))
          if(ntrans.gt.0) then
             if(lpresent.and..not.lsaving)then
                grp%nsymop = grp%nsymop + 1
                exit symloop
             end if
             transloop: do i=1,ntrans
                if(dot_product(trans(i,:),trans(i,:)).lt.tol_sym) &
                     cycle transloop
                if(ierror_sym.eq.3) write(77,*) trans(i,:)
                if(isym.ne.1)then
                   do jsym=2,grp%nsymop
                      if(grp%op(jsym).eq.1) then
                         if(all(abs(trans(i,1:3)-tmpsav(jsym,4,1:3)).lt.&
                              tol_sym)) cycle transloop
                         diff = trans(i,1:3) - tmpsav(jsym,4,1:3)
                         do j=1,3
                            diff(j) = diff(j) - floor(diff(j))
                            if(diff(j).gt.0.5) diff(j) = diff(j) - 1.D0
                         end do
                         do k=1,i
                            if(all(abs(diff-trans(k,1:3)).lt.tol_sym)) &
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
             if(lpresent) exit symloop
          end if
       end if
       oldnpntop = grp%npntop
    end do symloop


!!!-----------------------------------------------------------------------------
!!! allocates and saves the array savsym if the first time submitted
!!!-----------------------------------------------------------------------------
    if(lsaving)then
       if(allocated(savsym)) deallocate(savsym)
       allocate(savsym(grp%nsymop,4,4))
       savsym=0.D0
       savsym(:grp%nsymop,:,:)=tmpsav(:grp%nsymop,:,:)
       savsym(:,4,4)=1.D0
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
             write(0,'(2X,"check_sym in mod_sym.f90 is trying to assign a &
                  &value to symops_compare, which hasn''t been allocated")')
             exit iperm_if
          end if
          symops_compare(iperm)=grp%nsymop
       end select
    end if iperm_if


    if(lsaving)then
       deallocate(grp%sym)
       call move_alloc(savsym,grp%sym)
       grp%nsym = grp%nsymop
    end if


!!!-----------------------------------------------------------------------------
!!! if wyckoff present, set up wyckoff atoms
!!!-----------------------------------------------------------------------------
    if(lwyckoff)then
       if(present(lat).and.present(loc))then
          wyckoff=get_wyckoff_atoms(wyck_check(:grp%nsymop),lat,bas1,loc)
       else       
          wyckoff=get_wyckoff_atoms(wyck_check(:grp%nsymop))
       end if
    end if



    return
  end subroutine check_sym
!!!#############################################################################


!!!#############################################################################
!!! supplies the glides (if any) that are required to match the two bases ...
!!! ... "bas" and "tfbas" onto one another
!!!#############################################################################
  subroutine gldfnd(confine,bas,tfbas,trans,ntrans,transformed,wyck_check)
    implicit none
    integer :: i,j,ispec,iatom,jatom,katom,itmp1
    integer :: minspecloc,samecount
    logical :: lwyckoff
    double precision, dimension(3) :: ttrans,tmpbas,diff
    double precision, allocatable, dimension(:,:) :: sav_trans

    integer, intent(out) :: ntrans
    type(bas_type), intent(in) :: bas,tfbas
    type(confine_type), intent(in) :: confine
    double precision, dimension(:,:), intent(out) :: trans

    logical, optional, intent(in) :: transformed

    type(wyck_type), dimension(:), optional, intent(inout) :: wyck_check


!!!-----------------------------------------------------------------------------
!!! Allocate arrays and initialise variables
!!!-----------------------------------------------------------------------------
    ttrans=0.D0
    trans=0.D0
    samecount=0
    ntrans=0
    minspecloc=minloc(bas%spec(:)%num,mask=bas%spec(:)%num.ne.0,dim=1)

    if(present(transformed))then
       if(.not.transformed)then
          if(bas%spec(minspecloc)%num.eq.1) return
       end if
    else
       if(bas%spec(minspecloc)%num.eq.1) return
    end if
    allocate(sav_trans(bas%natom,3))


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
    trloop: do iatom=1,bas%spec(minspecloc)%num
       ttrans(:)=0.D0
       ttrans(1:3)=bas%spec(minspecloc)%atom(1,1:3)-&
            tfbas%spec(minspecloc)%atom(iatom,1:3)
       if(all(abs(ttrans(1:3)-anint(ttrans(1:3))).lt.tol_sym)) cycle trloop
       if(confine%l)then
          if(confine%laxis(confine%axis).and.&
               abs(ttrans(confine%axis)-nint(ttrans(confine%axis)))&
               .gt.tol_sym) cycle trloop
       end if
       itmp1 = 0
       sav_trans = 0.D0
       if(lwyckoff.and.ntrans+1.gt.size(wyck_check))then
          write(0,'("ERROR: error encountered in gldfnd")')
          write(0,'(2X,"Internal error in subroutine gldfnd in mod_sym.f90")')
          write(0,'(2X,"ntrans is greater than wyck_check")')
          write(0,'(2X,"EXITING SUBROUTINE")')
          return
       end if
       trcyc: do ispec=1,bas%nspec
          samecount=0
          if(lwyckoff) wyck_check(ntrans+1)%spec(ispec)%atom(:) = 0
          atmcyc2: do jatom=1,bas%spec(ispec)%num
             itmp1 = itmp1 + 1
             tmpbas(1:3) = tfbas%spec(ispec)%atom(jatom,1:3) + ttrans(1:3)
             tmpbas(:) = tmpbas(:) - ceiling(tmpbas(:)-0.5D0)
             atmcyc3: do katom=1,bas%spec(ispec)%num
                !if(lwyckoff.and.&
                !     wyck_check(ntrans+1)%spec(ispec)%atom(katom).ne.0) &
                !     cycle atmcyc3
                diff = tmpbas(1:3) - bas%spec(ispec)%atom(katom,1:3)
                do j=1,3
                   diff(j) = mod((diff(j)+100.D0),1.0)
                   if((abs(diff(j)-1.D0)).lt.(tol_sym)) diff(j) = 0.D0
                end do
                if(sqrt(dot_product(diff,diff)).lt.tol_sym)then
                   samecount = samecount + 1
                   !sav_trans(itmp1,:)=bas%spec(ispec)%atom(jatom,1:3)-&
                   !     bas%spec(ispec)%atom(katom,1:3)
                   sav_trans(itmp1,:) = bas%spec(ispec)%atom(katom,1:3) - &
                        tfbas%spec(ispec)%atom(jatom,1:3)
                   sav_trans(itmp1,:) = sav_trans(itmp1,:) - &
                        ceiling(sav_trans(itmp1,:)-0.5D0)
                   if(lwyckoff) &
                        wyck_check(ntrans+1)%spec(ispec)%atom(jatom) = katom
                   cycle atmcyc2
                end if
             end do atmcyc3
             !cycle trloop
          end do atmcyc2
          if (samecount.ne.bas%spec(ispec)%num) cycle trloop
       end do trcyc
!!!-----------------------------------------------------------------------------
!!! Cleans up succeeded translation vector
!!!-----------------------------------------------------------------------------
       do j=1,3
          itmp1 = maxloc(abs(sav_trans(:,j)),dim=1)
          ttrans(j) = sav_trans(itmp1,j)
          ttrans(j) = ttrans(j) - ceiling(ttrans(j)-0.5D0)
       end do
!!!-----------------------------------------------------------------------------
!!! If axis is confined, removes all symmetries not confined to the axis plane
!!!-----------------------------------------------------------------------------
       if(confine%l)then
          if(confine%laxis(confine%axis).and.&
               abs(ttrans(confine%axis)-nint(ttrans(confine%axis)))&
               .gt.tol_sym) cycle trloop
       else
          do i=1,3
             if(confine%laxis(i))then
                if(abs(ttrans(confine%axis)-floor(ttrans(confine%axis)))&
                     .lt.tol_sym) cycle trloop
             end if
          end do
       end if
!!!-----------------------------------------------------------------------------
!!! Checks whether this translation has already been saved
!!!-----------------------------------------------------------------------------
       do i=1,ntrans
          if(all(ttrans(:).eq.trans(i,:))) cycle trloop
          !if(all(abs(ttrans(:)-trans(i,:)).lt.tol_sym)) cycle trloop
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
  subroutine gen_fundam_sym_matrices(grp,lat)
    implicit none
    integer :: i
    type(sym_type) :: grp
    double precision :: cosPi3,sinPi3,mcosPi3,msinPi3
    double precision, dimension(3,3) :: inversion,invlat,tmat1
    double precision, dimension(64,3,3) :: fundam_mat
    double precision, dimension(3,3), intent(in) :: lat


    cosPi3 = 0.5D0
    sinPi3 = sin(pi/3.D0)
    mcosPi3 = -cosPi3
    msinPi3 = -sinPi3


    fundam_mat(1,1:3,1:3)=transpose(reshape((/&
         1.D0,  0.D0,  0.D0,  0.D0,  1.D0,  0.D0,  0.D0,  0.D0,  1.D0 /),&
         shape(inversion)))

    fundam_mat(2,1:3,1:3)=transpose(reshape((/&
         -1.D0,  0.D0,  0.D0,  0.D0, -1.D0,  0.D0,  0.D0, 0.D0,  1.D0 /),&
         shape(inversion)))

    fundam_mat(3,1:3,1:3)=transpose(reshape((/&
         -1.D0,  0.D0,  0.D0,  0.D0,  1.D0,  0.D0,  0.D0, 0.D0, -1.D0 /),&
         shape(inversion)))

    fundam_mat(4,1:3,1:3)=transpose(reshape((/&
         1.D0,  0.D0,  0.D0,  0.D0, -1.D0,  0.D0,  0.D0,  0.D0, -1.D0 /),&
         shape(inversion)))

    fundam_mat(5,1:3,1:3)=transpose(reshape((/&
         0.D0,  1.D0,  0.D0,  1.D0,  0.D0,  0.D0,  0.D0,  0.D0, -1.D0 /),&
         shape(inversion)))

    fundam_mat(6,1:3,1:3)=transpose(reshape((/&
         0.D0, -1.D0,  0.D0,  -1.D0,  0.D0,  0.D0,  0.D0, 0.D0, -1.D0 /),&
         shape(inversion)))

    fundam_mat(7,1:3,1:3)=transpose(reshape((/&
         0.D0, -1.D0,  0.D0,  1.D0,  0.D0,  0.D0,  0.D0,  0.D0,  1.D0 /),&
         shape(inversion)))

    fundam_mat(8,1:3,1:3)=transpose(reshape((/&
         0.D0,  1.D0,  0.D0,  -1.D0,  0.D0,  0.D0,  0.D0, 0.D0,  1.D0 /),&
         shape(inversion)))

    fundam_mat(9,1:3,1:3)=transpose(reshape((/&
         0.D0,  0.D0,  1.D0,  0.D0, -1.D0,  0.D0,  1.D0,  0.D0,  0.D0 /),&
         shape(inversion)))

    fundam_mat(10,1:3,1:3)=transpose(reshape((/&
         0.D0,  0.D0, -1.D0,  0.D0, -1.D0,  0.D0,  -1.D0, 0.D0,  0.D0 /),&
         shape(inversion)))

    fundam_mat(11,1:3,1:3)=transpose(reshape((/&
         0.D0,  0.D0, -1.D0,   0.D0,  1.D0,  0.D0,  1.D0, 0.D0,  0.D0 /),&
         shape(inversion)))

    fundam_mat(12,1:3,1:3)=transpose(reshape((/&
         0.D0,  0.D0,  1.D0,  0.D0,  1.D0,  0.D0,  -1.D0, 0.D0,  0.D0 /),&
         shape(inversion)))

    fundam_mat(13,1:3,1:3)=transpose(reshape((/&
         -1.D0,  0.D0,  0.D0,  0.D0,  0.D0,  1.D0,  0.D0, 1.D0,  0.D0 /),&
         shape(inversion)))

    fundam_mat(14,1:3,1:3)=transpose(reshape((/&
         -1.D0,  0.D0,  0.D0,  0.D0,  0.D0, -1.D0,  0.D0, -1.D0,  0.D0 /),&
         shape(inversion)))

    fundam_mat(15,1:3,1:3)=transpose(reshape((/&
         1.D0,  0.D0,  0.D0,  0.D0,  0.D0, -1.D0,  0.D0,  1.D0,  0.D0 /),&
         shape(inversion)))

    fundam_mat(16,1:3,1:3)=transpose(reshape((/&
         1.D0,  0.D0,  0.D0,  0.D0,  0.D0,  1.D0,  0.D0, -1.D0,  0.D0/),&
         shape(inversion)))

    fundam_mat(17,1:3,1:3)=transpose(reshape((/&
         0.D0,  0.D0,  1.D0,  1.D0,  0.D0,  0.D0,  0.D0,  1.D0,  0.D0 /),&
         shape(inversion)))

    fundam_mat(18,1:3,1:3)=transpose(reshape((/&
         0.D0,  0.D0, -1.D0, -1.D0,  0.D0,  0.D0,  0.D0,  1.D0,  0.D0 /),&
         shape(inversion)))

    fundam_mat(19,1:3,1:3)=transpose(reshape((/&
         0.D0,  0.D0, -1.D0,  1.D0,  0.D0,  0.D0,  0.D0, -1.D0,  0.D0 /),&
         shape(inversion)))

    fundam_mat(20,1:3,1:3)=transpose(reshape((/&
         0.D0,  0.D0,  1.D0, -1.D0,  0.D0,  0.D0,  0.D0, -1.D0,  0.D0 /),&
         shape(inversion)))

    fundam_mat(21,1:3,1:3)=transpose(reshape((/&
         0.D0,  1.D0,  0.D0,  0.D0,  0.D0,  1.D0,  1.D0,  0.D0,  0.D0 /),&
         shape(inversion)))

    fundam_mat(22,1:3,1:3)=transpose(reshape((/&
         0.D0, -1.D0,  0.D0,  0.D0,  0.D0, -1.D0,  1.D0,  0.D0,  0.D0 /),&
         shape(inversion)))

    fundam_mat(23,1:3,1:3)=transpose(reshape((/&
         0.D0, -1.D0,  0.D0,  0.D0,  0.D0,  1.D0, -1.D0,  0.D0,  0.D0 /),&
         shape(inversion)))

    fundam_mat(24,1:3,1:3)=transpose(reshape((/&
         0.D0,  1.D0,  0.D0,  0.D0,  0.D0, -1.D0, -1.D0,  0.D0,  0.D0 /),&
         shape(inversion)))

    fundam_mat(25,1:3,1:3)=transpose(reshape((/&
         cosPi3,  sinPi3, 0.D0, msinPi3,  cosPi3, 0.D0, 0.D0, 0.D0,  1.D0 /),&
         shape(inversion)))

    fundam_mat(26,1:3,1:3)=transpose(reshape((/&
         cosPi3, msinPi3, 0.D0,  sinPi3,  cosPi3, 0.D0, 0.D0, 0.D0,  1.D0 /),&
         shape(inversion)))

    fundam_mat(27,1:3,1:3)=transpose(reshape((/&
         mcosPi3,  sinPi3, 0.D0, msinPi3, mcosPi3, 0.D0, 0.D0, 0.D0, 1.D0 /),&
         shape(inversion)))

    fundam_mat(28,1:3,1:3)=transpose(reshape((/&
         mcosPi3, msinPi3, 0.D0,  sinPi3, mcosPi3, 0.D0, 0.D0, 0.D0, 1.D0 /),&
         shape(inversion)))

    fundam_mat(29,1:3,1:3)=transpose(reshape((/&
         cosPi3, msinPi3, 0.D0, msinPi3, mcosPi3, 0.D0, 0.D0, 0.D0, -1.D0 /),&
         shape(inversion)))

    fundam_mat(30,1:3,1:3)=transpose(reshape((/&
         cosPi3,  sinPi3, 0.D0,  sinPi3, mcosPi3, 0.D0, 0.D0, 0.D0, -1.D0 /),&
         shape(inversion)))

    fundam_mat(31,1:3,1:3)=transpose(reshape((/&
         mcosPi3, msinPi3, 0.D0, msinPi3,  cosPi3, 0.D0, 0.D0, 0.D0, -1.D0 /),&
         shape(inversion)))

    fundam_mat(32,1:3,1:3)=transpose(reshape((/&
         mcosPi3,  sinPi3, 0.D0,  sinPi3,  cosPi3, 0.D0, 0.D0, 0.D0, -1.D0 /),&
         shape(inversion)))

    inversion(:3,:3)=transpose(reshape((/&
         -1.D0,  0.D0,  0.D0,   0.D0,  -1.D0,  0.D0,   0.D0,  0.D0,  -1.D0 /),&
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
       if(abs(abs(det(tmat1))-1.D0).gt.tol_sym) cycle
       if(all(abs(tmat1-nint(tmat1)).le.tol_sym))then
          grp%nsym=grp%nsym+1
          fundam_mat(grp%nsym,:,:)=fundam_mat(i,:,:)
       end if
    end do


    allocate(grp%sym(grp%nsym,4,4))
    grp%sym(:,:,:)=0.D0
    grp%sym(:,4,4)=1.D0
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
  subroutine mksym(grp,inlat)
    implicit none
    integer :: amin,bmin,cmin
    integer :: i,j,ia,ib,ic,n,count,irot,nrot,isym,jsym
    double precision :: tht,a,b,c
    type(sym_type) :: grp
    double precision, dimension(3,3) :: rotmat,refmat,inlat,lat,invlat,tmat1
    double precision, allocatable, dimension(:,:,:) :: tsym1,tsym2
    logical, dimension(3) :: laxis


    if(grp%confine%l)then
       laxis=grp%confine%laxis
    else
       laxis=.not.grp%confine%laxis
    end if


!!!-----------------------------------------------------------------------------
!!! set up inverse lattice
!!!-----------------------------------------------------------------------------
    lat=inlat
    if(grp%lmolec)then
       invlat=0.D0
       lat=0.D0
    else
       invlat=inverse_3x3(lat)
    end if


!!!-----------------------------------------------------------------------------
!!! initialise values and symmetry matrix
!!!-----------------------------------------------------------------------------
    allocate(tsym1(50000,4,4))
    tsym1=0.D0
    tsym1(:,4,4)=1.D0
    count=0


!!!-----------------------------------------------------------------------------
!!! rotation plane perp to z (1=E,2=C2,3=C3,4=C4,5=C5,6=C6)
!!!-----------------------------------------------------------------------------
    if(laxis(3))then
       mksyml: do n=1,10
          count=count+1
          if(n.gt.6)then
             tht = -2.D0*pi/real(n-4) !=2*pi/(n-4)
          else
             tht = 2.D0*pi/real(n) !=2*pi/n          
          end if
          tsym1(count,1:3,1:3)=transpose(reshape((/&
               cos(tht) ,  sin(tht),   0.D0,&
               -sin(tht),  cos(tht),   0.D0,&
               0.D0     ,      0.D0,   1.D0/), shape(rotmat)))
          do i=1,3
             do j=1,3
                if(abs(tsym1(count,i,j)).lt.tol_sym) tsym1(count,i,j)=0.D0
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
             tht = -2.D0*pi/real(n-4) !=2*pi/n
          else
             tht = 2.D0*pi/real(n) !=2*pi/n
          end if
          rotmat=transpose(reshape((/&
               1.D0,      0.D0,      0.D0,  &
               0.D0,  cos(tht),  sin(tht),&
               0.D0, -sin(tht),  cos(tht)/), shape(rotmat)))
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
             tht = -2.D0*pi/real(n-4) !=2*pi/n 
          else
             tht = 2.D0*pi/real(n) !=2*pi/n 
          end if
          rotmat=transpose(reshape((/&
               cos(tht) ,  0.D0,  sin(tht),&
               0.D0     ,  1.D0,      0.D0,    &
               -sin(tht),  0.D0,  cos(tht)/), shape(rotmat)))
          rot3: do irot=1,nrot
             count=count+1
             tsym1(count,1:3,1:3)=matmul(rotmat(1:3,1:3),tsym1(irot,1:3,1:3))
             do i=1,3
                do j=1,3
                   if(abs(tsym1(count,i,j)).lt.tol_sym) tsym1(count,i,j)=0.D0
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
       a=(-1.D0)**ia
       bloop: do ib=bmin,2
          b=(-1.D0)**ib
          cloop: do ic=cmin,2
             c=(-1.D0)**ic
             !           if((a*b*c).ne.(-1.D0)) cycle cloop
             refmat(1:3,1:3)=transpose(reshape((/&
                  a,     0.D0,  0.D0,&
                  0.D0,  b   ,  0.D0,&
                  0.D0,  0.D0,     c/), shape(rotmat)))
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
    tsym2=0.D0
    tsym2(:,4,4)=1.D0
    count=0
    samecheck: do isym=1,grp%nsym
       tmat1=matmul((lat),tsym1(isym,:3,:3))
       tmat1=matmul(tmat1,(invlat))
       do i=1,3
          do j=1,3
             if(abs(tmat1(i,j)).lt.tol_sym) tmat1(i,j)=0.D0
             if(abs(1.D0-abs(tmat1(i,j))).lt.tol_sym) &
                  tmat1(i,j)=sign(1.D0,tmat1(i,j))
          end do
       end do
       !!-----------------------------------------------------------------------
       !! Precautionary measure
       if(all(abs(tmat1).lt.tol_sym)) cycle samecheck
       if(abs(abs(det(tmat1))-1.D0).gt.tol_sym) cycle samecheck
       !!-----------------------------------------------------------------------
       if(.not.all(abs(tmat1-nint(tmat1)).lt.tol_sym)) cycle samecheck
       do jsym=1,count
          if(all(tmat1.eq.tsym2(jsym,:3,:3))) cycle samecheck
          !if(all(tsym1(isym,:3,:3).eq.tsym2(jsym,:3,:3))) cycle samecheck
       end do
       count=count+1
       tsym2(count,:3,:3)=tmat1
       !tsym2(count,:4,:4)=tsym1(isym,:4,:4)
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
  subroutine clone_grp(ingrp,outgrp)
    implicit none
    type(sym_type), intent(in) :: ingrp
    type(sym_type), intent(out) :: outgrp
    
    
    allocate(outgrp%op(size(ingrp%op)))
    allocate(outgrp%sym(size(ingrp%sym(:,1,1)),4,4))
    outgrp = ingrp

  end subroutine clone_grp
!!!#############################################################################
 

!!!#############################################################################
!!! returns the primitive cell from a supercell
!!!#############################################################################
  subroutine get_primitive_cell(lat,bas)
    implicit none
    integer :: is,ia,ja,i,j,k,itmp1
    integer :: ntrans,len
    double precision :: scale,proj,dtmp1
    type(confine_type) :: confine
    double precision, dimension(3,3) :: dmat1,invlat
    double precision, allocatable, dimension(:,:) :: trans,atom_store
    
    type(sym_type) :: grp
    type(bas_type) :: bas,pbas
    double precision, dimension(3,3) :: lat

    
    !!-----------------------------------------------------------------------
    !! Allocate and initialise
    !!-----------------------------------------------------------------------
    ntrans = 0
    dmat1=0.D0
    allocate(trans(minval(bas%spec(:)%num+2),3)); trans=0.D0

    
    !!-----------------------------------------------------------------------
    !! Find the translation vectors in the cell
    !!-----------------------------------------------------------------------
    call gldfnd(confine,bas,bas,trans,ntrans,.false.)
    len=size(bas%spec(1)%atom,dim=2)

    
    !!-----------------------------------------------------------------------
    !! For each translation, reduce the basis
    !!-----------------------------------------------------------------------
    if(ntrans.ge.1)then
       do i=ntrans+1,ntrans+3
          trans(i,:)=0.D0
          trans(i,i-ntrans)=1.D0
       end do
       !  trans=matmul(trans(1:ntrans,1:3),lat)
       call sort2D(trans(1:ntrans+3,:),ntrans+3)
       !! for each lattice vector, determine the shortest translation ...
       !! ... vector that has a non-zero projection along that lattice vector.
       do i=1,3
          proj=1.D2
          trans_loop: do j=1,ntrans+3
             dtmp1 = dot_product(trans(j,:),trans(ntrans+i,:))
             if(dtmp1.lt.tol_sym) cycle trans_loop

             do k=1,i-1,1
                if(modu(abs(cross(trans(j,:),dmat1(k,:)))).lt.1.D-8) cycle trans_loop
             end do

             dtmp1 = modu(trans(j,:))
             if(dtmp1.lt.proj)then
                proj=dtmp1
                dmat1(i,:) = trans(j,:)
                trans(j,:) = 0.D0
             end if
          end do trans_loop
       end do
       !dmat1=trans(1:3,1:3)
       scale=det(dmat1)
       dmat1=matmul(dmat1,lat)
       invlat=inverse_3x3(dmat1)
       do is=1,bas%nspec
          itmp1=0
          allocate(atom_store(nint(scale*bas%spec(is)%num),len))
          atcheck: do ia=1,bas%spec(is)%num
             !!-----------------------------------------------------------------
             !! Reduce the basis
             !!-----------------------------------------------------------------
             bas%spec(is)%atom(ia,1:3)=&
                  matmul(bas%spec(is)%atom(ia,1:3),lat(1:3,1:3))
             bas%spec(is)%atom(ia,1:3)=&
                  matmul(transpose(invlat(1:3,1:3)),bas%spec(is)%atom(ia,1:3))
             do j=1,3
                bas%spec(is)%atom(ia,j)=&
                     bas%spec(is)%atom(ia,j)-floor(bas%spec(is)%atom(ia,j))
                if(bas%spec(is)%atom(ia,j).gt.1.D0-tol_sym) &
                     bas%spec(is)%atom(ia,j)=0.D0
             end do
             !!-----------------------------------------------------------------
             !! Check for duplicates in the cell
             !!-----------------------------------------------------------------
             do ja=1,ia-1
                if(all(abs(bas%spec(is)%atom(ia,1:3)-atom_store(ja,1:3)).lt.&
                     (/tol_sym,tol_sym,tol_sym/))) cycle atcheck
             end do
             itmp1=itmp1+1
             atom_store(itmp1,:)=bas%spec(is)%atom(ia,:)
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
          deallocate(bas%spec(is)%atom)
          call move_alloc(atom_store,bas%spec(is)%atom)
          bas%spec(is)%num=size(bas%spec(is)%atom,dim=1)
          !deallocate(atom_store)
       end do
       !!-----------------------------------------------------------------------
       !! Reduce the lattice
       !!-----------------------------------------------------------------------
       bas%natom=sum(bas%spec(:)%num)
       lat=dmat1
    end if

    
    !!-----------------------------------------------------------------------
    !! Reduce the lattice to symmetry definition
    !!-----------------------------------------------------------------------
    call reducer(lat, bas)
    !! next line necessary as FCC and BCC do not conform to Niggli reduced ...
    !! ... cell definitions.
    lat = primitive_lat(lat)


    
  end subroutine get_primitive_cell
!!!#############################################################################


!!!#############################################################################
!!! takes in transformation matrix and outputs its (x,y,z) definition
!!!#############################################################################
  subroutine symwrite (sym,symchar)
    implicit none
    integer :: i,j,nt,nr,div
    double precision, dimension(4,4) :: sym
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
    integer :: i,is,ia,isym,imin,itmp1
    integer :: nsym
    double precision :: dist
    logical :: lfound_closer
    type(wyck_type) :: wyckoff_atoms
    double precision, dimension(3) :: diff
    double precision, allocatable, dimension(:) :: dists
    integer, allocatable, dimension(:) :: ivtmp1

    type(bas_type), intent(in) :: bas
    double precision, dimension(3), intent(in) :: loc
    type(wyck_type), dimension(:), intent(in) :: wyckoff
    double precision, dimension(3,3), intent(in) :: lat


    nsym = size(wyckoff)
    allocate(wyckoff_atoms%spec(bas%nspec))
    wyckoff_atoms%spec(:)%num = 0
    do is=1,bas%nspec
       allocate(ivtmp1(size(wyckoff(1)%spec(is)%atom)))
       ivtmp1 = 0

       allocate(dists(bas%spec(is)%num))
       do ia=1,bas%spec(is)%num
          diff = loc - bas%spec(is)%atom(ia,:3)
          diff = diff - ceiling(diff - 0.5D0)
          dists(ia) = modu(matmul(diff,lat))
       end do

       wyckoff_loop1: do ia=1,size(wyckoff(1)%spec(is)%atom)

          dist = huge(0.D0)
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
  function basis_map(sym,bas1,tmpbas2) result(bas_map)
    implicit none
    integer :: j,ispec,iatom,jatom,dim
    type(basmap_type) :: bas_map
    type(bas_type) :: bas2,tfbas
    double precision, dimension(3) :: diff
    type(bas_type), intent(in) :: bas1
    double precision, dimension(4,4), intent(in) :: sym
    type(bas_type), optional, intent(in) :: tmpbas2


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
                  ceiling(tfbas%spec(ispec)%atom(iatom,j) - 0.5D0)
             bas2%spec(ispec)%atom(iatom,j) = &
                  bas2%spec(ispec)%atom(iatom,j) - &
                  ceiling(bas2%spec(ispec)%atom(iatom,j) - 0.5D0)
          end do
       end do
    end do


!!!-----------------------------------------------------------------------------
!!! check whether transformed basis matches original basis
!!!-----------------------------------------------------------------------------
    spcheck2: do ispec=1,bas1%nspec
       diff=0.D0
       atmcheck2: do iatom=1,bas1%spec(ispec)%num
          atmcyc2: do jatom=1,bas1%spec(ispec)%num
             if(any(bas_map%spec(ispec)%atom(:).eq.jatom)) cycle atmcyc2
             diff = tfbas%spec(ispec)%atom(iatom,1:3) - &
                  bas2%spec(ispec)%atom(jatom,1:3)
             diff = diff - ceiling(diff - 0.5D0)
             if(sqrt(dot_product(diff,diff)).lt.tol_sym)then
                bas_map%spec(ispec)%atom(iatom) = jatom
             end if
          end do atmcyc2
       end do atmcheck2
    end do spcheck2


    return
  end function basis_map
!!!#############################################################################


!!!#############################################################################
!!! finds all possible terminations along an axis
!!!#############################################################################
  function get_terminations(lat,bas,axis,lprint,layer_sep) result(term)
    implicit none
    integer :: i,j,is,nterm,mterm,dim,ireject
    integer :: itmp1,itmp2,init,min_loc
    logical :: ludef_print,lunique,ltmp1,lmirror
    double precision :: dtmp1,tol,height,max_sep,c_along,centre
    type(sym_type) :: grp1,grp_store
    type(term_arr_type) :: term
    integer, dimension(3) :: abc=(/1,2,3/)
    double precision, dimension(3) :: vec_compare,vtmp1
    double precision, dimension(3,3) :: inv_mat,ident
    type(bas_type),allocatable, dimension(:) :: bas_arr,bas_arr_reject
    type(term_type), allocatable, dimension(:) :: term_arr,term_arr_uniq
    integer, allocatable, dimension(:) :: success,tmpop
    integer, allocatable, dimension(:,:) :: reject_match
    double precision, allocatable, dimension(:,:) :: bas_list
    double precision, allocatable, dimension(:,:,:) :: tmpsym

    integer, intent(in) :: axis
    type(bas_type), intent(in) :: bas
    double precision, dimension(3,3), intent(in) :: lat

    double precision, optional, intent(in) :: layer_sep
    logical, optional, intent(in) :: lprint



!!!APPLY TRANSFORMATION MATRIX TO FIND TERMINATIONS ALONG OTHER PLANES
!!! E.G. (1 0 1)
    
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


!!!-----------------------------------------------------------------------------
!!! Sets the surface identification tolerance
!!!-----------------------------------------------------------------------------
    if(present(layer_sep))then
       tol = layer_sep
    else
       tol = 1.D0  !!!tolerance of 1 Å for defining a layer
    end if

    abc=cshift(abc,3-axis)
    c_along = abs(dot_product(lat(axis,:),&
         uvec(cross(lat(abc(1),:),lat(abc(2),:)))))
    tol = tol / c_along
    !tol = tol/modu(lat(axis,1:3))
    lmirror=.false.


!!!-----------------------------------------------------------------------------
!!! Set up basis list that will order them wrt distance along 'axis'
!!!-----------------------------------------------------------------------------
    allocate(bas_list(bas%natom,3))
    init = 1
    do is=1,bas%nspec
       bas_list(init:init+bas%spec(is)%num-1,:3) = bas%spec(is)%atom(:,:3)
       init = init + bas%spec(is)%num
    end do
    call sort_col(bas_list,col=axis)


!!!-----------------------------------------------------------------------------
!!! Find largest separation between atoms
!!!-----------------------------------------------------------------------------
    max_sep = bas_list(1,axis) - (bas_list(bas%natom,axis)-1.D0)
    height = ( bas_list(1,axis) + (bas_list(bas%natom,axis)-1.D0) )/2.D0
    do i=1,bas%natom-1
       dtmp1 = bas_list(i+1,axis) - bas_list(i,axis)
       if(dtmp1.gt.max_sep)then
          max_sep = dtmp1
          height = ( bas_list(i+1,axis) + bas_list(i,axis) )/2.D0
       end if
    end do
    if(max_sep.lt.tol)then
       write(0,'("ERROR: Error in mod_sym.f90")')
       write(0,'(2X,"get_terminations subroutine unable to find a separation &
            &in the material that is greater than LAYER_SEP")')
       write(0,'(2X,"Writing material to ''unlayerable.vasp''")')
       open(13,file="unlayerable.vasp")
       call geom_write(13,lat,bas)
       close(13)
       write(0,'(2X,"We suggest reducing LAYER_SEP to less than ",F6.4)') &
            max_sep
       write(0,'(2X,"Please inform the developers of this and give details &
            &of what structure caused this")')
       write(0,'("Stopping...")')
       stop
    end if
    bas_list(:,axis) = bas_list(:,axis) - height
    bas_list(:,axis) = bas_list(:,axis) - floor(bas_list(:,axis))
    call sort_col(bas_list,col=axis)


!!!-----------------------------------------------------------------------------
!!! Finds number of non-unique terminations
!!!-----------------------------------------------------------------------------
    nterm=1
    allocate(term_arr(bas%natom))
    term_arr(:)%natom=0
    term_arr(:)%hmin=0
    term_arr(:)%hmax=0
    term_arr(1)%hmin=bas_list(1,axis)
    term_arr(1)%hmax=bas_list(1,axis)
    min_loc = 1
    itmp1 = 1
    term_loop1: do

       !! get the atom at that height.
       !vtmp1 = get_min_dist(lat,bas,bas_list(itmp1,:3),.true.,axis,.true.,.false.)
       !itmp1 = minloc(bas_list(:,axis) - vtmp1(axis), dim=1, &
       !     mask = abs(bas_list(:,axis) - (bas_list(itmp1,axis) + vtmp1(axis))&
       !     ).lt.tol_sym)
       
       itmp1 = minloc(bas_list(:,axis) - term_arr(nterm)%hmax, dim=1, &
            mask = bas_list(:,axis) - term_arr(nterm)%hmax.gt.0.D0)
       if(itmp1.gt.bas%natom.or.itmp1.le.0)then
          term_arr(nterm)%natom = bas%natom - min_loc + 1
          exit term_loop1
       end if

       !dtmp1 = modu(matmul(vtmp1,lat))
       dtmp1 = bas_list(itmp1,axis) - term_arr(nterm)%hmax
       if(dtmp1.le.tol)then
          term_arr(nterm)%hmax = bas_list(itmp1,axis)
       else
          term_arr(nterm)%natom = itmp1 - min_loc
          min_loc = itmp1
          nterm = nterm + 1
          term_arr(nterm)%hmin = bas_list(itmp1,axis)
          term_arr(nterm)%hmax = bas_list(itmp1,axis)
       end if
       
    end do term_loop1
    term_arr(:nterm)%hmin = term_arr(:nterm)%hmin + height
    term_arr(:nterm)%hmax = term_arr(:nterm)%hmax + height


!!!-----------------------------------------------------------------------------
!!! Set up system symmetries
!!!-----------------------------------------------------------------------------
    allocate(bas_arr(2*nterm))
    allocate(bas_arr_reject(2*nterm))
    dim = size(bas%spec(1)%atom(1,:))
    do i=1,2*nterm
       allocate(bas_arr(i)%spec(bas%nspec))
       allocate(bas_arr_reject(i)%spec(bas%nspec))
       do is=1,bas%nspec
          allocate(bas_arr(i)%spec(is)%atom(&
               bas%spec(is)%num,dim))
          allocate(bas_arr_reject(i)%spec(is)%atom(&
               bas%spec(is)%num,dim))
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
    call sym_setup(grp_store,lat,predefined=.false.,new_start=.true.)

    !!WRITE OUT THE STRUCTURES HERE AND COMPARE
    !do i=1,grp_store%nsym
    !   write(0,*) i
    !   write(0,'(4(2X,F6.2))') grp_store%sym(i,:4,:3)
    !   write(0,*) det(grp_store%sym(i,:3,:3))
    !   write(0,*)
    !end do


    !!--------------------------------------------------------------------------
    !! Handle inversion matrix (centre of inversion must be accounted for)
    !!--------------------------------------------------------------------------
    !! change symmetry constraints after setting up symmetries
    !! this is done to constrain the matching of two bases in certain directions
    grp_store%confine%l = .false.
    grp_store%confine%laxis(axis) = .false.
    call check_sym(grp_store,bas1=bas,iperm=-1,lsave=.true.)
    inv_mat = 0.D0
    do i=1,3
       inv_mat(i,i) = -1.D0
    end do
    do i=1,grp_store%nsym
       if(all(abs(grp_store%sym(i,:3,:3)-inv_mat).lt.tol_sym))then
          itmp1 = i
          exit
       end if
    end do
    do i=1,grp_store%nsymop
       if(all(abs(savsym(i,:3,:3)-inv_mat).lt.tol_sym)) &
            grp_store%sym(itmp1,4,:3) = savsym(i,4,:3)
    end do
    !do i=1,grp_store%nsymop
    !   write(0,*) i
    !   write(0,'(4(2X,F9.4))') grp_store%sym(i,:4,:3)
    !   write(0,*) det(grp_store%sym(i,:3,:3))
    !   write(0,*)
    !end do


    !!--------------------------------------------------------------------------
    !! Determine unique surface terminations
    !!--------------------------------------------------------------------------
    grp_store%confine%l = .true.
    grp_store%confine%laxis(axis) = .true.
    allocate(term_arr_uniq(2*nterm))
    allocate(reject_match(nterm,2))
    shift_loop1:do i=1,nterm
       mterm = mterm + 1

       bas_arr(mterm) = bas
       centre = term_arr(i)%hmin + (term_arr(i)%hmax - term_arr(i)%hmin)/2.D0
       call shifter(bas_arr(mterm),axis,1-centre,.true.)
       !if(ludef_print) write(6,'(1X,I3,8X,F7.5,9X,F7.5,8X,I3)') &
       !     i,term_arr(i)%hmin,term_arr(i)%hmax,term_arr(i)%natom
       sym_if: if(i.ne.1)then
          sym_loop1:do j=1,mterm-1
             if(abs(abs(term_arr(i)%hmax-term_arr(i)%hmin) - &
                  abs(term_arr_uniq(j)%hmax-term_arr_uniq(j)%hmin)).gt.tol_sym) &
                  cycle sym_loop1
             call clone_grp(grp_store,grp1)
             call check_sym(grp1,bas1=bas_arr(mterm),&
                  iperm=-1,tmpbas2=bas_arr(j),lsave=.true.)
             if(grp1%nsymop.ne.0)then
                !write(0,*) "we have a possible reject"
                !if(any(savsym(:grp1%nsymop,axis,axis).eq.-1.D0))then
                if(savsym(1,axis,axis).eq.-1.D0)then
                   ireject = ireject + 1
                   reject_match(ireject,:) = [ i, j ]
                   bas_arr_reject(ireject) = bas_arr(mterm)
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
       term_arr_uniq(mterm)%ladder(:) = 0.D0
    end do shift_loop1


    !!--------------------------------------------------------------------------
    !! Set up mirror/inversion symmetries of the matrix
    !!--------------------------------------------------------------------------
    call sym_setup(grp_store,lat,predefined=.false.,new_start=.true.)
    allocate(tmpsym(count(grp_store%sym(:,3,3).eq.-1.D0),4,4))
    allocate(tmpop(count(grp_store%sym(:,3,3).eq.-1.D0)))
    itmp1 = 0
    do i=1,grp_store%nsym
       if(grp_store%sym(i,3,3).eq.-1.D0)then
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
    ident = 0.D0
    do i=1,3
       ident(i,i) = 1.D0
    end do
    vec_compare = 0.D0
    vec_compare(axis) = -1.D0
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
          call check_sym(grp1,bas_arr(itmp2),&
               iperm=-1,lsave=.true.,lcheck_all=.true.)
          ltmp1=.false.

!!!HERE
          !! Check if pure translations are present in comparison termination?
          !do j=1,grp1%nsymop
          !   if(all(abs(savsym(j,:3,:3)-ident).le.tol_sym))then
          !      write(0,*) "FOUND TRANSLATION"
          !      cycle reject_loop1
          !   end if
          !end do
          !! Check if inversions are present in comparison termination
          do j=1,grp1%nsymop
             if(abs(det(savsym(j,:3,:3))+1.D0).le.tol_sym) ltmp1=.true.
          end do
          !! If they are not, then no point comparing. It is a new termination
          if(.not.ltmp1) exit prior_check 

          call clone_grp(grp_store,grp1)
          call check_sym(grp1,bas_arr(itmp2),&
               tmpbas2=bas_arr_reject(i),iperm=-1,lsave=.true.,lcheck_all=.true.)

          !! Check det of all symmetry operations. If any are 1, move on
          !! This is because they are just rotations as can be captured ...
          !! ... through lattice matches.
          !! Solely inversions are unique and must be captured.
          do j=1,grp1%nsymop
             if(abs(det(savsym(j,:3,:3))-1.D0).le.tol_sym) lunique=.false.
          end do
          if(savsym(1,4,axis).eq.&
               2.D0*min(term_arr_uniq(itmp2)%hmin,0.5D0-term_arr_uniq(itmp2)%hmin))then
             lunique=.false.
          end if

          if(.not.(all(savsym(1,axis,:3).eq.vec_compare(:)).and.&
               all(savsym(1,:3,axis).eq.vec_compare(:)))) lunique=.false.
          
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
          term_arr_uniq(mterm)%ladder(1) = 0.D0
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
    dtmp1 = term_arr_uniq(1)%hmin-1.D-6
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
            mask=term_arr_uniq(:)%hmin.gt.dtmp1+tol,dim=1)
       if(itmp1.eq.0) then
          itmp1 = minloc(term_arr_uniq(:)%hmin,&
               mask=term_arr_uniq(:)%hmin.gt.dtmp1+tol-1.D0,dim=1)
       end if
       dtmp1 = term_arr_uniq(itmp1)%hmin
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


  end function get_terminations
!!!#############################################################################

end module mod_sym
