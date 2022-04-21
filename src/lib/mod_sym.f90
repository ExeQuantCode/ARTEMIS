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
       clone_bas,get_closest_atom,ortho_axis,reducer,primitive_lat
  implicit none
  integer :: ierror_sym=0
  integer :: s_start=1,s_end=0
  double precision :: tol_sym=5.D-5
  character(1) :: verb_sym="n"
  integer, allocatable, dimension(:) :: symops_compare
  integer, allocatable, dimension(:,:,:) :: wyckoff,tmpwyckoff
  double precision, allocatable, dimension(:,:,:) :: savsym


  private


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
     double precision, allocatable, dimension(:) :: ladder
  end type term_type

  type term_arr_type
     integer :: nterm,axis,nstep
     double precision :: tol
     logical :: lmirror
     type(term_type), allocatable, dimension(:) :: arr
  end type term_arr_type


  type confine_type
     integer :: axis=3
     logical :: l=.false.
     logical :: lmirror=.false.
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
  public :: sym_setup,check_sym,gldfnd

  public :: get_primitive_cell
  
  public :: setup_ladder
  public :: term_arr_type,confine_type
  public :: get_terminations,print_terminations

  public :: basmap_type,basis_map


!!!updated 2022/04/04


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
    if(allocated(wyckoff)) deallocate(wyckoff)
    if(allocated(tmpwyckoff)) deallocate(tmpwyckoff)
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
  subroutine check_sym(grp,bas1,iperm,tmpbas2,lsave,lcheck_all)
    implicit none
    integer :: i,j,k,iatom,jatom,ispec,itmp1
    integer :: is,isym,jsym,count,ntrans
    integer :: samecount,oldnpntop
    logical :: lpresent,lsaving,ltransformed
    type(bas_type) :: bas1,bas2,tfbas
    type(sym_type) :: grp
    double precision, dimension(3) :: diff
    double precision, dimension(3,3) :: ident
    double precision, allocatable, dimension(:,:) :: trans
    double precision, allocatable, dimension(:,:,:) :: tmpsav
    integer, optional, intent(in) :: iperm
    type(bas_type), optional :: tmpbas2
    logical, optional, intent(in) :: lsave, lcheck_all


204 format(4(F11.6),/,4(F11.6),/,4(F11.6),/,4(F11.6))
    

!!!-----------------------------------------------------------------------------
!!! allocated grp%op
!!!-----------------------------------------------------------------------------
    if(.not.allocated(grp%op))then
       allocate(grp%op(grp%nsym*minval(bas1%spec(:)%num)))
       grp%op=0
    end if
    if(present(lsave))then
       lsaving=lsave
    else
       lsaving=.false.
    end if


!!!-----------------------------------------------------------------------------
!!! checks for optional arguments and assigns values if not present
!!!-----------------------------------------------------------------------------
    if(present(tmpbas2)) then
       bas2=tmpbas2
       if(present(lcheck_all))then
          lpresent=.not.lcheck_all
       else
          lpresent=.true.
       end if
    else
       bas2=bas1
       lpresent=.false.
    end if
    allocate(tmpsav(grp%nsym*minval(bas1%spec(:)%num),4,4))
    itmp1=maxval(bas1%spec(:)%num)
    if(.not.allocated(wyckoff)) &
         allocate(wyckoff(bas1%nspec,itmp1,itmp1))
    if(.not.allocated(tmpwyckoff)) &
         allocate(tmpwyckoff(bas1%nspec,itmp1,itmp1))


!!!-----------------------------------------------------------------------------
!!! initialises variables
!!!-----------------------------------------------------------------------------
    allocate(trans(minval(bas1%spec(:)%num+2),3)); trans=0.D0
    allocate(tfbas%spec(bas1%nspec))
    itmp1=size(bas1%spec(1)%atom(1,:),dim=1)
    do is=1,bas1%nspec
       allocate(tfbas%spec(is)%atom(bas1%spec(is)%num,itmp1))
    end do
    grp%nsymop=0
    grp%npntop=0
    ! wyckoff section
    !##########################
    wyckoff=0
    tmpwyckoff=0
    do ispec=1,bas1%nspec
       do iatom=1,bas1%spec(ispec)%num
          wyckoff(ispec,iatom,1:bas1%spec(ispec)%num)=1
          tmpwyckoff(ispec,iatom,1:bas1%spec(ispec)%num)=1
       end do
    end do
    !##########################


!!!-----------------------------------------------------------------------------
!!! set up identity matrix as reference
!!!-----------------------------------------------------------------------------
    ltransformed=.false.
    ident = 0.D0
    do i=1,3
       ident(i,i) = 1.D0
    end do


!!!-----------------------------------------------------------------------------
!!! applying symmetries to basis to see if the basis conforms to any of them
!!!-----------------------------------------------------------------------------
    symloop: do isym=s_start,s_end
       tmpwyckoff=wyckoff
       if(verb_sym.eq.'d') write(77,*) isym !,a,b,c
       if(verb_sym.eq.'d') write(77,204) grp%sym(isym,1:4,1:4)
       if(ierror_sym.eq.2.or.ierror_sym.eq.3) write(77,204)  &
            grp%sym(isym,1:4,1:4)
       !------------------------------------------------------------------------
       ! apply symmetry operator to basis
       !------------------------------------------------------------------------
       do ispec=1,bas1%nspec
          do iatom=1,bas1%spec(ispec)%num
             tfbas%spec(ispec)%atom(iatom,1:3)=&
                  matmul(bas1%spec(ispec)%atom(iatom,1:4),grp%sym(isym,1:4,1:3))
             do j=1,3
                tfbas%spec(ispec)%atom(iatom,j)=&
                     tfbas%spec(ispec)%atom(iatom,j)-&
                     ceiling(tfbas%spec(ispec)%atom(iatom,j)-0.5D0)
             end do
          end do
       end do
       !------------------------------------------------------------------------
       ! check whether transformed basis matches original basis
       !------------------------------------------------------------------------
       count=0
       spcheck: do ispec=1,bas1%nspec
          diff=0.0
          samecount=0
          atmcheck: do iatom=1,bas1%spec(ispec)%num
             atmcyc: do jatom=1,bas1%spec(ispec)%num
                diff=tfbas%spec(ispec)%atom(iatom,1:3)-&
                     bas2%spec(ispec)%atom(jatom,1:3)
                do j=1,3
                   diff(j)=mod((diff(j)+100.D0),1.0)
                   diff(j)=diff(j)-floor(diff(j))
                   if((abs(diff(j)-1.D0)).lt.(tol_sym)) diff(j)=0.D0
                end do
                if(sqrt(dot_product(diff,diff)).lt.tol_sym)then
                   samecount=samecount+1
                   tmpwyckoff(ispec,iatom,jatom)=0
                   tmpwyckoff(ispec,jatom,iatom)=0
                end if
                if((iatom.eq.bas1%spec(ispec)%num).and.&
                     (jatom.eq.bas1%spec(ispec)%num))then
                   if (samecount.ne.bas1%spec(ispec)%num)then
                      goto 10
                   end if
                end if
             end do atmcyc
             count=count+samecount
          end do atmcheck
          if(samecount.ne.bas1%spec(ispec)%num) goto 10
       end do spcheck
       grp%npntop=grp%npntop+1
       grp%nsymop=grp%nsymop+1
       wyckoff=tmpwyckoff
       tmpsav(grp%nsymop,:,:)=grp%sym(isym,:,:)
       grp%op(grp%nsymop)=isym
       if(grp%nsymop.ne.0.and.lpresent) exit symloop
10     trans=0.D0
       ntrans=0
       tmpwyckoff=wyckoff
       !------------------------------------------------------------------------
       ! checks if translations are valid with the current symmetry operation
       !------------------------------------------------------------------------
       if(grp%lspace) then
          if(all(abs(grp%sym(isym,1:3,1:3)-ident).lt.tol_sym))then
             ltransformed=.false.
          else
             ltransformed=.true.
          end if
          call gldfnd(grp%confine,bas2,tfbas,trans,ntrans,transformed=ltransformed)
          if(ntrans.gt.0) then
             if(lpresent.and..not.lsaving)then
                grp%nsymop=grp%nsymop+1
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
                         diff=trans(i,1:3)-tmpsav(jsym,4,1:3)
                         do j=1,3
                            diff(j)=diff(j)-floor(diff(j))
                            if(diff(j).gt.0.5) diff(j)=diff(j)-1.D0
                         end do
                         do k=1,i
                            if(all(abs(diff-trans(k,1:3)).lt.tol_sym)) &
                                 cycle transloop
                         end do
                      end if
                   end do
                end if
                grp%nsymop=grp%nsymop+1
                tmpsav(grp%nsymop,:,:)=grp%sym(isym,:,:)
                tmpsav(grp%nsymop,4,1:3)=trans(i,:)
                grp%op(grp%nsymop)=isym
             end do transloop
             if(lpresent) exit symloop
          end if
       end if
       oldnpntop=grp%npntop
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
       call move_alloc(savsym,grp%sym)
       grp%nsym=grp%nsymop
    end if


    return
  end subroutine check_sym
!!!#############################################################################


!!!#############################################################################
!!! supplies the glides (if any) that are required to match the two bases ...
!!! ... "bas" and "tfbas" onto one another
!!!#############################################################################
  subroutine gldfnd(confine,bas,tfbas,trans,ntrans,transformed)
    implicit none
    integer :: i,j,ispec,iatom,jatom,katom,itmp1
    integer :: minspecloc,samecount,ntrans
    type(bas_type) :: bas,tfbas
    type(confine_type) :: confine
    !    integer, allocatable, dimension(:,:,:) :: tmpwyck
    double precision, dimension(3) :: ttrans,tmpbas,diff
    double precision, dimension(:,:) :: trans
    double precision, allocatable, dimension(:,:) :: sav_trans

    logical, optional, intent(in) :: transformed


!!!-----------------------------------------------------------------------------
!!! Allocate arrays and initialise variables
!!!-----------------------------------------------------------------------------
    !    allocate(tmpwyck(bas%nspec,maxval(bas%spec(:)%num),maxval(bas%spec(:)%num)))
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
       !       tmpwyck=wyckoff
       itmp1=0
       sav_trans=0.D0
       trcyc: do ispec=1,bas%nspec
          samecount=0
          atmcyc2: do jatom=1,bas%spec(ispec)%num
             itmp1=itmp1+1
             tmpbas(1:3)=tfbas%spec(ispec)%atom(jatom,1:3)+ttrans(1:3)
             tmpbas(:)=tmpbas(:)-ceiling(tmpbas(:)-0.5D0)
             atmcyc3: do katom=1,bas%spec(ispec)%num
                diff=tmpbas(1:3)-bas%spec(ispec)%atom(katom,1:3)
                do j=1,3
                   diff(j)=mod((diff(j)+100.D0),1.0)
                   if((abs(diff(j)-1.D0)).lt.(tol_sym)) diff(j)=0.D0
                end do
                if(sqrt(dot_product(diff,diff)).lt.tol_sym)then
                   samecount=samecount+1
                   !sav_trans(itmp1,:)=bas%spec(ispec)%atom(jatom,1:3)-&
                   !     bas%spec(ispec)%atom(katom,1:3)
                   sav_trans(itmp1,:)=bas%spec(ispec)%atom(katom,1:3)-&
                        tfbas%spec(ispec)%atom(jatom,1:3)
                   sav_trans(itmp1,:)=sav_trans(itmp1,:)-&
                        ceiling(sav_trans(itmp1,:)-0.5D0)
                   !                   tmpwyck(ispec,jatom,katom)=0
                   !                   tmpwyck(ispec,katom,jatom)=0
                   cycle atmcyc2
                end if
             end do atmcyc3
             cycle trloop
          end do atmcyc2
          if (samecount.ne.bas%spec(ispec)%num)then
             cycle trloop
          end if
       end do trcyc
!!!-----------------------------------------------------------------------------
!!! Cleans up succeeded translation vector
!!!-----------------------------------------------------------------------------
       do j=1,3
          itmp1=maxloc(abs(sav_trans(:,j)),dim=1)
          ttrans(j)=sav_trans(itmp1,j)
          ttrans(j)=ttrans(j)-ceiling(ttrans(j)-0.5D0)
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
             if(confine%laxis(i).and.&
                  abs(ttrans(confine%axis)-floor(ttrans(confine%axis)))&
                  .lt.tol_sym) cycle trloop
          end do
       end if
!!!-----------------------------------------------------------------------------
!!! Checks whether this translation has already been saved
!!!-----------------------------------------------------------------------------
       do i=1,ntrans
          if(all(ttrans(:).eq.trans(i,:))) cycle trloop
          !if(all(abs(ttrans(:)-trans(i,:)).lt.tol_sym)) cycle trloop
       end do
       !       wyckoff=tmpwyck
       ntrans=ntrans+1
       trans(ntrans,1:3)=ttrans(1:3)
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
       if(all(abs(tmat1-nint(tmat1)).lt.tol_sym))then
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
             tht = -8.D0*atan(1.D0)/real(n-4) !=2*pi/(n-4)
          else
             tht = 8.D0*atan(1.D0)/real(n) !=2*pi/n          
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
             tht = -8.D0*atan(1.D0)/real(n-4) !=2*pi/n
          else
             tht = 8.D0*atan(1.D0)/real(n) !=2*pi/n
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
             tht = -8.D0*atan(1.D0)/real(n-4) !=2*pi/n 
          else
             tht = 8.D0*atan(1.D0)/real(n) !=2*pi/n 
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
    integer :: is,ia,ja,i,j,itmp1
    integer :: ntrans,len
    double precision :: scale
    type(confine_type) :: confine
    double precision, dimension(3,3) :: dmat1, invlat
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
       dmat1=trans(1:3,1:3)
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
                     &many atoms from supercell!"
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
!!! sets up the ladder
!!!#############################################################################
  subroutine setup_ladder(lat,bas,axis,term)
    implicit none
    integer :: i
    integer :: nmirror,ntrans
    double precision :: dtmp1
    type(sym_type) :: grp
    logical :: lexclude_translation
    logical, dimension(2,2) :: mask
    
    double precision, allocatable, dimension(:) :: ladder !, ladder2
    double precision, allocatable, dimension(:,:,:) :: subgroup,group
    double precision, allocatable, dimension(:,:,:) :: mirror_mat,trans_mat

    integer, intent(in) :: axis
    type(bas_type), intent(in) :: bas
    double precision, dimension(3,3), intent(in) :: lat
    type(term_arr_type), optional, intent(inout) :: term

    
    !!--------------------------------------------------------------------------
    !! Test if mirror/inversion or translation exists
    !!--------------------------------------------------------------------------
    call sym_setup(grp,lat,predefined=.true.,new_start=.true.)
    call check_sym(grp,bas,lsave=.true.)
    term%lmirror = .false.
    allocate(mirror_mat(count(grp%sym(:,3,3).eq.-1.D0),2,2))
    allocate(trans_mat(&
         count(grp%sym(:,3,3).eq.1.D0.and.abs(grp%sym(:,4,axis)).gt.tol_sym),2,2))
    mirror_mat(:,:,:) = 0.D0
    do i=1,size(mirror_mat(:,1,1))
       mirror_mat(i,:,2) = [ 0.D0, 1.D0 ]
    end do
    trans_mat(:,:,:) = 0.D0
    do i=1,size(trans_mat(:,1,1))
       trans_mat(i,:,2) = [ 0.D0, 1.D0 ]
    end do

    nmirror = 0
    ntrans = 0
    mirror_loop: do i=1,grp%nsym
       !write(0,*) i
       !write(0,'(4(2X,F9.4))') grp%sym(i,:,:)
       !write(0,*)
       if(grp%sym(i,3,3).eq.-1.D0)then
          term%lmirror = .true.
          if(all(mirror_mat(:nmirror,2,1).ne.grp%sym(i,4,axis)))then
             nmirror = nmirror + 1
             mirror_mat(nmirror,:,1) = grp%sym(i,3:4,axis)
          end if
       elseif(grp%sym(i,3,3).eq.1.D0 .and. abs(grp%sym(i,4,axis)).gt.tol_sym)then
          if(all(trans_mat(:ntrans,2,1).ne.grp%sym(i,4,axis)))then
             ntrans = ntrans + 1
             trans_mat(ntrans,:,1) = grp%sym(i,3:4,axis)
          end if
       end if
    end do mirror_loop


    !!--------------------------------------------------------------------------
    !! If termination not present, then return
    !!--------------------------------------------------------------------------
    if(.not.present(term)) return


    !!--------------------------------------------------------------------------
    !! Handle no symmetry situation
    !!--------------------------------------------------------------------------
    if(ntrans+nmirror.eq.0)then
       term%nstep=1
       return
    end if


    !!--------------------------------------------------------------------------
    !! Set up rungs
    !!--------------------------------------------------------------------------
    allocate(subgroup(ntrans+nmirror,2,2))
    if(size(trans_mat).gt.0) subgroup(:ntrans,:,:) = trans_mat(:ntrans,:,:)
    if(size(mirror_mat).gt.0) subgroup(ntrans+1:ntrans+nmirror,:,:) = mirror_mat(:nmirror,:,:)
    mask = .false.
    mask(2,1) = .true.
    group = gen_group(subgroup,mask,tol_sym)
    !write(0,*) "-----------------"
    !do i=1,size(group(:,1,1))
    !   write(0,*) i
    !   write(0,'(2(2X,F9.6))') group(i,:,:)
    !   write(0,*)
    !end do
    allocate(ladder(size(group(:,1,1))))
    ladder = 0.D0

    lexclude_translation = .false.
    if(term%lmirror.and.abs(term%arr(1)%hmax-term%arr(1)%hmin).gt.tol_sym)&
         lexclude_translation = .true.
    term%nstep = 0
    ! account for when there are no translational symmetries in the cell
    if(size(subgroup).lt.1)then
       term%nstep = 1
    end if
    group_loop: do i=1,size(group(:,1,1))
       if(any(abs(ladder(:term%nstep)-group(i,2,1)).lt.tol_sym)) cycle group_loop
       if(any(abs(ladder(:term%nstep)-&
            floor(ladder(:term%nstep)+tol_sym)-group(i,2,1)).lt.tol_sym)) cycle group_loop
       if(lexclude_translation.and.&
            abs(group(i,1,1)-1.D0).lt.tol_sym) cycle group_loop
       term%nstep = term%nstep + 1
       ladder(term%nstep) = group(i,2,1)
       !if(abs(ladder(term%nstep)).lt.tol_sym) &
       !     ladder(term%nstep) = 0.D0
       if(abs(ladder(term%nstep)-nint(ladder(term%nstep))).lt.tol_sym) &
            ladder(term%nstep) = 0.D0
    end do group_loop
    call sort1D(ladder(:term%nstep))
    !call set(ladder2, tol_sym)
    !allocate(ladder2(term%nstep))
    !ladder2(:) = ladder(:term%nstep)
    !term%nstep = size(ladder2)

    if(term%nstep.gt.1)then
       dtmp1 = ladder(1)
       do i=1,term%nstep-1,1
          ladder(i)=ladder(i+1)
       end do
       ladder(term%nstep) = dtmp1+1.D0
    end if
    do i=1,term%nterm
       allocate(term%arr(i)%ladder(term%nstep))
       term%arr(i)%ladder(:) = ladder(:term%nstep)
    end do



  end subroutine setup_ladder
!!!#############################################################################


!!!#############################################################################
!!! finds all possible terminations along an axis
!!!#############################################################################
  function get_terminations(lat,bas,axis,lprint,layer_sep) result(term)
    implicit none
    integer :: i,j,is,nterm,mterm,dim,ireject
    integer :: itmp1,init,min_loc
    logical :: ludef_print
    double precision :: dtmp1,tol,height,max_sep,c_along
    type(sym_type) :: grp1,grp_store
    type(term_arr_type) :: term
    integer, dimension(3) :: abc=(/1,2,3/)
    double precision, dimension(3) :: vec_compare
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
    grp_store%confine%laxis(axis)=.true.
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
    term_loop1: do

       itmp1 = minloc(bas_list(:,axis) - term_arr(nterm)%hmax, dim=1, &
            mask = bas_list(:,axis) - term_arr(nterm)%hmax.gt.0.D0)
       if(itmp1.gt.bas%natom.or.itmp1.le.0)then
          term_arr(nterm)%natom = bas%natom - min_loc + 1
          exit term_loop1
       end if

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
    dim=size(bas%spec(1)%atom(1,:))
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
    mterm=0
    ireject=0
    grp_store%lspace=.true.
    call sym_setup(grp_store,lat)
    call check_sym(grp_store,bas1=bas)
    allocate(term_arr_uniq(2*nterm))
    allocate(reject_match(nterm,2))
    if(ludef_print)&
         write(6,'(1X,"Term.",3X,"Min layer loc",3X,"Max layer loc",3X,"no. atoms")')
    shift_loop1:do i=1,nterm
       mterm = mterm + 1

       bas_arr(mterm) = bas
       call shifter(bas_arr(mterm),axis,1-term_arr(i)%hmax,.true.)
       sym_if: if(i.ne.1)then
          sym_loop1:do j=1,mterm-1
             call clone_grp(grp_store,grp1)
             call check_sym(grp1,bas1=bas_arr(mterm),&
                  iperm=-1,tmpbas2=bas_arr(j),lsave=.true.)
             if(grp1%nsymop.ne.0)then
                if(savsym(1,axis,axis).eq.-1.D0)then
                   ireject = ireject + 1
                   reject_match(ireject,:) = [ i, j ]
                   bas_arr_reject(ireject) = bas_arr(mterm)
                end if
                mterm = mterm - 1
                cycle shift_loop1
             end if
          end do sym_loop1
       end if sym_if
       if(ludef_print) write(6,'(1X,I3,8X,F7.5,9X,F7.5,8X,I3)') &
            mterm,term_arr(i)%hmin,term_arr(i)%hmax,term_arr(i)%natom
       term_arr_uniq(mterm) = term_arr(i)
       !open(100+mterm)
       !call geom_write(100+mterm,lat,bas_arr(mterm))
       !close(100+mterm)
    end do shift_loop1


    !!--------------------------------------------------------------------------
    !! Set up mirror/inversion symmetries of the matrix
    !!--------------------------------------------------------------------------
    call sym_setup(grp_store,lat,predefined=.true.,new_start=.true.)
    allocate(tmpsym(count(grp_store%sym(:,3,3).eq.-1.D0),4,4))
    allocate(tmpop(count(grp_store%sym(:,3,3).eq.-1.D0)))
    itmp1=0
    do i=1,grp_store%nsym
       if(grp_store%sym(i,3,3).eq.-1.D0)then
          itmp1=itmp1+1
          tmpsym(itmp1,:,:)=grp_store%sym(i,:,:)
          tmpop(itmp1) = i
       end if
    end do
    grp_store%nsym=itmp1
    grp_store%nlatsym=itmp1
    call move_alloc(tmpsym,grp_store%sym)
    allocate(grp_store%op(itmp1))
    grp_store%op(:) = tmpop(:itmp1)
    s_end = grp_store%nsym


    !!--------------------------------------------------------------------------
    !! Check rejects for inverse surface termination of saved
    !!--------------------------------------------------------------------------
    vec_compare = 0.D0
    vec_compare(axis) = -1.D0
    allocate(success(ireject))
    success=0
    reject_loop1: do i=1,ireject
       itmp1=reject_match(i,2)
       if(any(success(1:i-1).eq.itmp1)) cycle reject_loop1
       call clone_grp(grp_store,grp1)
       call check_sym(grp1,bas_arr(itmp1),&
            tmpbas2=bas_arr_reject(i),iperm=-1,lsave=.true.,lcheck_all=.true.)

       !! CHECK DETERMINANT OF ALL SYMMETRY OPERATIONS. IF THERE ARE ANY THAT ARE 1, THEN MOVE ON
       !! This is because they are just rotations as can be captured through lattice matches.
       !! Solely inversions are unique and must be captured.
       do j=1,grp1%nsymop
          !write(0,'(4(2X,F9.4))') savsym(j,:,:)
          !write(0,*) det(savsym(j,:3,:3))
          if(abs(det(savsym(j,:3,:3))-1.D0).le.tol_sym) cycle reject_loop1
       end do
       
       if(all(savsym(1,axis,:3).eq.vec_compare(:)).and.&
            all(savsym(1,:3,axis).eq.vec_compare(:)))then
          !write(0,*) savsym(:,4,axis)
          !write(0,*) "test0",savsym(1,4,axis),2.D0*min(term_arr_uniq(itmp1)%hmin,0.5D0-term_arr_uniq(itmp1)%hmin)
          !write(0,*) "test1",term_arr_uniq(itmp1)%hmin,0.5D0-term_arr_uniq(itmp1)%hmin
          !write(0,*) "test2",term_arr_uniq(itmp1)%hmin,term_arr_uniq(itmp1)%hmax
          !write(0,*) "test3",term_arr(reject_match(i,1))%hmin
          !if(savsym(1,4,axis).eq.0.D0.and.&
          !     abs(term_arr_uniq(itmp1)%hmin-term_arr_uniq(itmp1)%hmax).lt.tol_sym)then
          !   cycle reject_loop1
          if(savsym(1,4,axis).eq.&
               2.D0*min(term_arr_uniq(itmp1)%hmin,0.5D0-term_arr_uniq(itmp1)%hmin))then
             cycle reject_loop1
          end if

          !dtmp1 = term_arr(reject_match(i,1))%hmin
          !if(dtmp1.gt.0.5D0)then
          !   dtmp1 = 0.5D0 - dtmp1
          !end if
          !dtmp1 = dtmp1 + savsym(1,4,axis)
          !write(0,*) "here",dtmp1,term_arr_uniq(itmp1)%hmin
          !if(abs(dtmp1 - term_arr_uniq(itmp1)%hmin).lt.tol_sym)then
          !   cycle reject_loop1
          !end if
          !
          !
          !! REMOVE THE j SYM_LOOP2 LOOP !!!
          !! IT IS UNECESSARY AS ALL INVERSION SYMS WILL HAVE SAME TRANSLATION !!!
          !reject_loop2: do j=1,i-1!count(reject_match(:,2).eq.reject_match(i,2))
          !   if(itmp1.eq.reject_match(j,2))then
          !      itmp2=reject_match(j,1)
          !      write(0,*) i,j,itmp1
          !      write(0,*) "min",term_arr(itmp2)%hmin,0.5D0-term_arr(itmp2)%hmin,&
          !           term_arr_uniq(itmp1)%hmin,term_arr(reject_match(i,1))%hmin
          !      dtmp1 = term_arr(reject_match(i,1))%hmin
          !      if(dtmp1.gt.0.5D0)then
          !         dtmp1 = 0.5D0 - dtmp1
          !      end if
          !      dtmp1 = dtmp1 + savsym(1,4,axis)
          !      if(abs(dtmp1 - term_arr(itmp2)%hmin).lt.tol_sym)then
          !         cycle reject_loop1
          !      end if
!         !       if(savsym(1,4,axis).eq.&
!         !            2.D0*min(term_arr(itmp2)%hmin,0.5D0-term_arr(itmp2)%hmin))then
!         !          cycle reject_loop1
!         !       end if
          !   end if
          !end do reject_loop2
!!! CYCLES ONE NOW SUCCESSFUL !!!
!!! PROBABLY SAVE itmp1 AND STOP LOOKING AT THOSE !!!
          mterm=mterm+1
          success(i)=itmp1
          term_arr_uniq(mterm)=term_arr(reject_match(i,1))
          if(ludef_print) write(6,'(1X,I3,8X,F7.5,9X,F7.5,8X,I3)') &
               mterm,&
               term_arr_uniq(mterm)%hmin,&
               term_arr_uniq(mterm)%hmax,term_arr_uniq(mterm)%natom
          reject_match(i,2)=0
       end if
    end do reject_loop1

    allocate(term%arr(mterm))
    term%tol=tol
    term%axis=axis
    term%nterm=mterm
    term%arr(:mterm)=term_arr_uniq(:mterm)

    do i=1,term%nterm
       if(i.eq.1)then
          dtmp1 = abs(term%arr(i)%hmin-term%arr(term%nterm)%hmax)/4.D0
       else
          dtmp1 = abs(term%arr(i)%hmin-term%arr(i-1)%hmax)/4.D0
       end if
       if(dtmp1.lt.term%tol)then
          term%tol = dtmp1
       end if
    end do
!!! THERE IS NO CLEAR TERMINATION PLANE
!!! DO THIS IF ANY TERMINATION IS LARGER THAN A CERTAIN SIZE, 


  end function get_terminations
!!!#############################################################################


!!!#############################################################################
!!! prints the terminations to individual files
!!!#############################################################################
  subroutine print_terminations(term,inlat,inbas,dirname,&
       thickness,vacuum,lortho)
    implicit none
    integer :: unit,i,j,istep,ncells,udef_thick
    double precision :: vac,dtmp1
    character(1024) :: filename,pwd
    logical :: udef_lortho
    type(term_arr_type) :: term
    type(bas_type) :: tbas,bas
    double precision, dimension(3,3) :: tfmat
    double precision, dimension(3,3) :: tlat,lat

    type(bas_type), intent(in) :: inbas
    double precision, dimension(3,3), intent(in) :: inlat

    integer, optional, intent(in) :: thickness
    double precision, optional, intent(in) :: vacuum
    character(*), optional, intent(in) :: dirname
    logical, optional, intent(in) :: lortho


    !!--------------------------------------------------------------------------
    !! Handles optional parameters
    !!--------------------------------------------------------------------------
    if(present(lortho))then
       udef_lortho = lortho
    else
       udef_lortho = .true.
    end if

    if(present(thickness))then
       udef_thick = thickness
    else
       udef_thick = 2
    end if
    
    if(present(vacuum))then
       vac = vacuum
    else
       vac = 14.D0
    end if


    !!--------------------------------------------------------------------------
    !! Makes directory and enters
    !!--------------------------------------------------------------------------
    call clone_bas(inbas,bas)
    if(present(dirname))then
       call system('mkdir -p '//trim(adjustl(dirname)))
       call getcwd(pwd)
       call chdir(dirname)
    end if
    

    !!--------------------------------------------------------------------------
    !! Increases crystal to max number of required unit cells thick
    !!--------------------------------------------------------------------------
    ncells = int((udef_thick-1)/term%nstep)+1
    tfmat(:,:)=0.D0
    tfmat(1,1)=1.D0
    tfmat(2,2)=1.D0
    tfmat(3,3)=ncells
    tbas = inbas
    tlat = inlat
    call transformer(tlat,tbas,tfmat)


    term%arr(:)%hmin = term%arr(:)%hmin/dble(ncells)
    term%arr(:)%hmax = term%arr(:)%hmax/dble(ncells)
    term%tol = term%tol/dble(ncells)
    

    !!--------------------------------------------------------------------------
    !! Generate each termination
    !!--------------------------------------------------------------------------
    do i=1,term%nterm
       bas = tbas
       lat = tlat
       tfmat = 0.D0
       !!-----------------------------------------------------------------------
       !! Shifts material to specified termination
       !!-----------------------------------------------------------------------
       call shifter(bas,term%axis,-term%arr(i)%hmin,.true.)
       
       
       !!-----------------------------------------------------------------------
       !! Determines how much extension is required and performs extension
       !!-----------------------------------------------------------------------
       do j=1,3
          tfmat(j,j) = 1.D0
          if(j.eq.term%axis)then
             istep = udef_thick - (ncells-1)*term%nstep
             dtmp1 = (ncells-1) + term%arr(i)%ladder(istep)
             dtmp1 = dtmp1/(ncells)
             dtmp1 = dtmp1 + (term%arr(i)%hmax - term%arr(i)%hmin)
             tfmat(j,j) = dtmp1 + term%tol/8.D0
             if(.not.term%lmirror)then
                tfmat(j,j) = tfmat(j,j) + (term%arr(i)%hmax - term%arr(i)%hmin)
             end if
          end if
       end do
       call transformer(lat,bas,tfmat)
  

       !!-----------------------------------------------------------------------
       !! If requested, orthogonalises interface axis wrt the other two axes
       !!-----------------------------------------------------------------------
       if(udef_lortho)then
          ortho_check1: do j=1,2
             if(abs(dot_product(lat(j,:),lat(term%axis,:))).gt.tol_sym)then
                call ortho_axis(lat,bas,term%axis)
                exit ortho_check1
             end if
          end do ortho_check1
       end if


       !!-----------------------------------------------------------------------
       !! Adds vacuum
       !!-----------------------------------------------------------------------
       call set_vacuum(lat,bas,term%axis,1.D0,vac)


       !!-----------------------------------------------------------------------
       !! Prints structure
       !!-----------------------------------------------------------------------
       unit=20+i
       write(filename,'("POSCAR_term",I0)') i
       open(unit,file=filename)
       call geom_write(unit,lat,bas)
       close(unit)
    end do
    

    
    !!--------------------------------------------------------------------------
    !! Returns to original directory
    !!--------------------------------------------------------------------------
    if(present(dirname))then
       call chdir(pwd)
    end if


    return
  end subroutine print_terminations
!!!#############################################################################

end module mod_sym
