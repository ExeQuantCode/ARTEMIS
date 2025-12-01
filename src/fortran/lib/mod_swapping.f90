!!!#############################################################################
!!! Code written by Ned Thaddeus Taylor and Isiah Edward Mikel Rudkin 
!!! Code part of the ARTEMIS group (Hepplestone research group).
!!! Think Hepplestone, think HRG.
!!!#############################################################################
module swapping
  use artemis__constants, only: real32
  use artemis__misc, only: sort1D
  use misc_maths, only: gauss
  use misc_linalg, only: modu
  use artemis__geom_rw, only: basis_type
  use artemis__sym, only: check_sym,sym_type,basis_map_type,basis_map
  use artemis__io_utils, only: err_abort
  implicit none
  real(real32) :: tiny=5.E-5_real32
  logical :: lmirror
  type(basis_map_type) :: bas_map

  private
  
  public :: rand_swapper


!!!updated 2020/02/07


contains
!!!#############################################################################
!!! Main function to be called from ARTEMIS
!!!#############################################################################
  function rand_swapper(lat,bas,axis,width,nswaps_per_cell,nswap,intf_loc,&
       iswap,seed_arr,tol_sym, verbose, sigma,require_mirror) result(bas_arr)
    implicit none
    integer :: i,j,is,iout,itmp,count1
    integer :: axis,nswap
    integer :: nabove,nbelow,nswaps_per_cell,nfail !,nperm
    real(real32) :: udef_sigma,small_sigma
    real(real32) :: dintf,dist
    type(basis_type) :: tmpbas,store_bas
    type(sym_type) :: grp
    !real(real32), dimension(4,4) :: intf_sym
    integer, allocatable, dimension(:) :: spec_list
    integer, allocatable, dimension(:) :: lw_close_list,up_close_list
    real(real32), allocatable, dimension(:) :: lw_dist_list,up_dist_list
    real(real32), allocatable, dimension(:) :: lw_weight_list,up_weight_list

    integer, allocatable, dimension(:,:) :: pos_list,up_list,lw_list
    real(real32), allocatable, dimension(:,:) :: bas_list
    real(real32), dimension(4,4) :: intf_sym

    integer, intent(in) :: iswap
    real(real32), intent(in) :: width
    real(real32), optional, intent(in) :: sigma
    logical, optional, intent(in) :: require_mirror
    type(basis_type), intent(in) :: bas
    integer, dimension(:), intent(in) :: seed_arr
    real(real32), dimension(2), intent(in) :: intf_loc !USE 1
    type(basis_type), allocatable, dimension(:) :: bas_arr
    real(real32), dimension(3,3), intent(in) :: lat
    real(real32), intent(in) :: tol_sym
    integer, intent(in) :: verbose


!!!-----------------------------------------------------------------------------
!!! initialises variables
!!!-----------------------------------------------------------------------------
    grp%nsymop = 1
    nfail=50
    if(present(sigma))then
       if(sigma.lt.0._real32)then
          udef_sigma = 0.05
       else
          udef_sigma = sigma
       end if
    else
       udef_sigma = 0.05
    end if
    udef_sigma = udef_sigma/modu(lat(axis,:))
    small_sigma = 0.01/modu(lat(axis,:))
    call random_seed(put=seed_arr)


!!!-----------------------------------------------------------------------------
!!! set up basis and positions list
!!!-----------------------------------------------------------------------------
    allocate(bas_list(bas%natom,3))
    allocate(pos_list(bas%natom,2))
    itmp=1
    do is=1,bas%nspec
       bas_list(itmp:bas%spec(is)%num+itmp-1,1:3)=bas%spec(is)%atom(:,1:3)
       pos_list(itmp:bas%spec(is)%num+itmp-1,1)=is
       do i=itmp,bas%spec(is)%num+itmp-1
          pos_list(i,2)=i-itmp+1
       end do
       !pos_list(itmp:bas%spec(is)%num+itmp-1,2)=(/ 1:bas%spec(is)%num /)
       itmp=itmp+bas%spec(is)%num
    end do


!!!-----------------------------------------------------------------------------
!!! find number of atoms within range of interface
!!!-----------------------------------------------------------------------------
    dist=width/modu(lat(axis,:))


!!!!-----------------------------------------------------------------------------
!!!! set number of permutations
!!!!-----------------------------------------------------------------------------
    !    dtmp=lnsum(min(nabove,nbelow))
    !    dtmp=dtmp-lnsum(min(nabove,nbelow)-nswaps_per_cell)-lnsum(nswaps_per_cell)
    !    nperm=nint(exp(dtmp))
    !    if(nperm.le.0) nperm=10
    !
    !    if(nswap.gt.nperm)then
    !       write(*,'(1X,A)') "Number of possible permutations is less than requested value."
    !       write(*,'(1X,A,I0)') "Resetting number of output structures to ",nperm
    !       nswap=nperm
    !    end if
!!!-----------------------------------------------------------------------------
    
    
!!!-----------------------------------------------------------------------------
!!! set up symmetries
!!!-----------------------------------------------------------------------------
    call grp%init(lat, tol_sym = tol_sym)
    call tmpbas%copy(bas, length = 4)
    call store_bas%copy(tmpbas, length = 4)


!!!-----------------------------------------------------------------------------
!!! set up array of bases
!!!-----------------------------------------------------------------------------
    allocate(bas_arr(nswap))
    do i=1,nswap
       allocate(bas_arr(i)%spec(bas%nspec))
       do is=1,bas%nspec
          allocate(bas_arr(i)%spec(is)%atom(bas%spec(is)%num,4))
       end do
       bas_arr(i)%nspec = 0
    end do


!!!-----------------------------------------------------------------------------
!!! find symmetry that maps top interface onto bottom interface
!!!-----------------------------------------------------------------------------
!!! NOT NEEDED?
!!! To Replace with
    lmirror = .false.
    call check_sym(grp,tmpbas,lsave=.true., tol_sym=tol_sym)
    intf_sym_loop: do i = 1, grp%nsymop
       !if(symops(i).eq.1) cycle intf_sym_loop
       if(abs(grp%sym(4,axis,i)).lt.tiny) cycle intf_sym_loop
       if(abs(grp%sym(axis,axis,i)+1._real32).gt.tiny) cycle intf_sym_loop
       intf_sym(1:4,1:4) = grp%sym(1:4,1:4,i)
       bas_map = basis_map(intf_sym,tmpbas, tol_sym=tol_sym)
       lmirror = .true.
       exit intf_sym_loop
    end do intf_sym_loop

    if(lmirror)then
       do i=1,bas%nspec
          if(any(bas_map%spec(i)%atom.eq.0))then
             call err_abort("&
                  &ERROR: Internal error in rand_swapper\n&
                  & Error in rand_swapper subroutine in mod_swapper.f90\n&
                  & atom missing a mapping even though mirror symmetry was found\n&
                  &Exiting...",fmtd=.true.)
          end if
       end do
       if(verbose.ge.1)then
          write(*,*) "mirror found for swaps"
          write(*,'(4(2X,F9.4))') intf_sym(:,:)
          write(*,*)
       end if
    else
       write(0,*) "WARNING: No mirror identified in interface"
       if(present(require_mirror))then
          if(.not.require_mirror)then
             write(0,*) "Performing swaps over only one interface then"
             goto 10
          end if
       end if
       write(0,*) "As such, cannot generate equivalent swaps on both interfaces"
       write(0,*) "Skipping..."
       return
    end if

10  deallocate(grp%sym)
    call grp%init(lat,new_start=.true., tol_sym = tol_sym)
    call check_sym(grp,tmpbas, tol_sym=tol_sym)!,lsave=.true.)
    
    
    dintf=intf_loc(1)
!!!-----------------------------------------------------------------------------
!!! finds number of atoms below and above the interface and records them
!!!-----------------------------------------------------------------------------
    select case(iswap)
    case(1)
       call check_intf(lat,bas,dintf,dist,lw_list,up_list,nbelow,nabove,bas_list,pos_list,axis) 
    case(2)
       call check_intf_depth(lat,bas,axis,intf_loc,udef_sigma,&
            spec_list,&
            lw_list,up_list,&
            lw_dist_list,up_dist_list,&
            lw_close_list,up_close_list,&
            lw_weight_list,up_weight_list)
       nbelow = size(lw_close_list)
       nabove = size(up_close_list)
    end select
    if(nswaps_per_cell.gt.min(nabove,nbelow))then
       write(*,'(1X,A)') "Number of possible swaps is less than requested value."
       write(*,'(1X,A,I0)') "Resetting number of swaps to ",min(nabove,nbelow)
       nswaps_per_cell=min(nabove,nbelow)
    end if


!!!-----------------------------------------------------------------------------
!!! swap atoms
!!!-----------------------------------------------------------------------------
    !randomly swaps atoms from below to above the interface (and vice versa)
    select case(iswap)
    case(1)
       call rand_swap(store_bas,tmpbas,nabove,nbelow,nswaps_per_cell,up_list,lw_list) 
    case(2)
       call rand_swap_depth(store_bas,tmpbas,&
            nswaps_per_cell,udef_sigma,small_sigma,&
            spec_list,&
            lw_list,up_list,&
            lw_dist_list,up_dist_list,&
            lw_close_list,up_close_list,&
            lw_weight_list,up_weight_list, &
            verbose=verbose &
       )
    end select
    bas_arr(1) = tmpbas
    iout = 1
    count1 = 0
    symloop:do
       count1 = count1 + 1
       if(count1.gt.nfail) exit symloop
       itmp = iout + 1
       if(iout.eq.nswap) exit symloop
       tmpbas = store_bas
       select case(iswap)
       case(1)
          !! randomly swaps atoms from below to above the interface (and vice versa)
          call rand_swap(store_bas,tmpbas,nabove,nbelow,nswaps_per_cell,up_list,lw_list)
       case(2)
          call rand_swap_depth(store_bas,tmpbas,&
               nswaps_per_cell,udef_sigma,small_sigma,&
               spec_list,&
               lw_list,up_list,&
               lw_dist_list,up_dist_list,&
               lw_close_list,up_close_list,&
               lw_weight_list,up_weight_list, &
               verbose=verbose &
          )
       end select
       !call check_sym(tmpbas,itmp,tol_sym=tol_sym)
       !call loadbar(iout,10)

       do j=1,iout
          call check_sym(grp,basis=tmpbas,tmpbas2=bas_arr(j),tol_sym=tol_sym)!,itmp,bas_arr(j))
          if(grp%nsymop.ne.0) cycle symloop
       end do

       iout=iout+1
       bas_arr(iout)=tmpbas
       count1=0
       !     call loadbar(iout,10,'y')
    end do symloop



end function rand_swapper
!!!#############################################################################



!!!#############################################################################
!!!#############################################################################
!!! M E T H O D   1
!!!#############################################################################
!!!#############################################################################


!!!#############################################################################
!!! finds number of atoms below and above the interface and records them
!!!#############################################################################
  subroutine check_intf(lat,bas,dintf,width,lw_list,up_list,nbelow,nabove,bas_list,pos_list,axis)
    implicit none
    integer :: i,itmp1,itmp2
    integer :: nbelow,nabove,axis
    real(real32) :: dintf,width
    type(basis_type) :: bas
    real(real32), dimension(3,3) :: lat
    real(real32), dimension(:,:) :: bas_list
    integer, allocatable, dimension(:,:) :: lw_list,up_list,pos_list


    nbelow=count(dintf-bas_list(:,axis).le.width.and.dintf-bas_list(:,axis).ge.0)
    nabove=count(bas_list(:,axis)-dintf.le.width.and.bas_list(:,axis)-dintf.gt.0)
    if(min(nabove,nbelow).eq.0)then
       write(*,'(1X,"No atoms found within ",F0.2," Å of the interface.")') width*modu(lat(axis,:))
       write(*,'(1X,"Exiting code...")')
       call exit()
    end if

    if(allocated(lw_list)) deallocate(lw_list)
    if(allocated(up_list)) deallocate(up_list)
    allocate(lw_list(nbelow,2))
    allocate(up_list(nabove,2))
    itmp1=0
    itmp2=0
    do i=1,bas%natom
       if(dintf-bas_list(i,axis).le.width.and.dintf-bas_list(i,axis).ge.0)then
          itmp1=itmp1+1
          lw_list(itmp1,:)=pos_list(i,:)
       end if
       if(bas_list(i,axis)-dintf.le.width.and.bas_list(i,axis)-dintf.gt.0)then
          itmp2=itmp2+1
          up_list(itmp2,:)=pos_list(i,:)
        end if
    end do


  end subroutine check_intf
!!!#############################################################################



!!!#############################################################################
!!! randomly swaps atoms from below to above the interface (and vice versa)
!!!#############################################################################
  subroutine rand_swap(bas,swap_bas,nabove,nbelow,nswaps_per_cell,up_list,lw_list)
    implicit none
    integer :: i,nfail
    integer :: itmp1,itmp2,old_itmp1
    integer :: lw_mirror,up_mirror
    integer :: lw_remove,up_remove,nabove,nbelow,nswaps_per_cell
    real(real32) :: r_rand
    integer, allocatable, dimension(:,:) :: swap_list,up_list,lw_list
    type(basis_type) :: bas,swap_bas

!!!-----------------------------------------------------------------------------
!!! randomly select atoms above and below the interface
!!!-----------------------------------------------------------------------------
    allocate(swap_list(nswaps_per_cell,2))
    swap_list=0
    itmp1 = 1
    old_itmp1 = 0
    itmp2 = 0
    nfail = 100

    swap_loop: do
       if(itmp1.eq.old_itmp1)then
          itmp2 = itmp2 + 1
       else
          itmp2 = 0
       end if

       if(itmp2.ge.nfail)then
          call err_abort("&
               &ERROR: Internal error in rand_swap\n&
               & Error in rand_swap subroutine in mod_swapper.f90\n&
               & all atoms either side of the interface appear to be the same species\n&
               &Exiting...",fmtd=.true.)
       end if
       !r_rand=rand(0)
       lw_loop: do
          call random_number(r_rand)
          lw_remove=nint(r_rand*(nbelow-1))
          lw_remove=lw_remove+1
          if(any(swap_list(1:itmp1-1,1).eq.lw_remove)) cycle lw_loop
          swap_list(itmp1,1)=lw_remove
          exit lw_loop
       end do lw_loop

       up_loop: do
          call random_number(r_rand)
          up_remove=nint(r_rand*(nabove-1))
          up_remove=up_remove+1
          if(any(swap_list(1:itmp1-1,2).eq.up_remove)) cycle up_loop
          swap_list(itmp1,2)=up_remove
          exit up_loop
       end do up_loop

       old_itmp1 = itmp1
       if(lw_list(lw_remove,1).eq.up_list(up_remove,1)) cycle swap_loop
       if(itmp1.ge.nswaps_per_cell) exit swap_loop
       itmp1 = itmp1 + 1
    end do swap_loop


!!!-----------------------------------------------------------------------------
!!! swap chosen atoms with each other
!!!-----------------------------------------------------------------------------
    do i=1,nswaps_per_cell
       swap_bas%spec(lw_list(swap_list(i,1),1))%atom(lw_list(swap_list(i,1),2),:)=&
            bas%spec(up_list(swap_list(i,2),1))%atom(up_list(swap_list(i,2),2),:)
       swap_bas%spec(up_list(swap_list(i,2),1))%atom(up_list(swap_list(i,2),2),:)=&
            bas%spec(lw_list(swap_list(i,1),1))%atom(lw_list(swap_list(i,1),2),:)
       
       if(lmirror)then
          lw_mirror=bas_map%spec(lw_list(swap_list(i,1),1))%atom(lw_list(swap_list(i,1),2))
          up_mirror=bas_map%spec(up_list(swap_list(i,2),1))%atom(up_list(swap_list(i,2),2))

          swap_bas%spec(lw_list(swap_list(i,1),1))%atom(lw_mirror,:)=&
               bas%spec(up_list(swap_list(i,2),1))%atom(up_mirror,:)
          swap_bas%spec(up_list(swap_list(i,2),1))%atom(up_mirror,:)=&
               bas%spec(lw_list(swap_list(i,1),1))%atom(lw_mirror,:)
       end if
    end do


  end subroutine rand_swap
!!!#############################################################################



!!!#############################################################################
!!!#############################################################################
!!! M E T H O D   2
!!!#############################################################################
!!!#############################################################################



!!!#############################################################################
!!! 
!!!#############################################################################
  subroutine check_intf_depth(lat,bas,axis,intf_loc,sigma,&
       spec_list,&
       lw_list,up_list,&
       lw_dist_list,up_dist_list,&
       lw_close_list,up_close_list,&
       lw_weight_list,up_weight_list)
    implicit none
    integer :: i,is,ia
    integer :: nbelow,nabove
    real(real32) :: rtol
    real(real32), dimension(2) :: midpoint
    integer, allocatable, dimension(:) :: tmp_list1,tmp_list2
    real(real32), allocatable, dimension(:) :: tmp_dist_list1,tmp_dist_list2

    integer, intent(in) :: axis
    real(real32), intent(in) :: sigma
    type(basis_type), intent(in) :: bas
    real(real32), dimension(2), intent(in) :: intf_loc
    real(real32), dimension(3,3), intent(in) :: lat

    integer, allocatable, dimension(:), intent(out) :: spec_list
    integer, allocatable, dimension(:), intent(out) :: lw_close_list,up_close_list
    real(real32), allocatable, dimension(:), intent(out) :: lw_dist_list,up_dist_list
    real(real32), allocatable, dimension(:), intent(out) :: lw_weight_list,up_weight_list
    integer, allocatable, dimension(:,:), intent(out) :: lw_list,up_list


!!!-----------------------------------------------------------------------------
!!! Initialise tolerances and set up midpoints
!!!-----------------------------------------------------------------------------
    rtol = 0.1/modu(lat(axis,:))

    midpoint(1) = (intf_loc(1) + intf_loc(2))/2
    midpoint(2) = (1._real32 + intf_loc(1) + intf_loc(2))/2

    if(midpoint(1).lt.intf_loc(1)) &
         midpoint(1) = midpoint(1) + 1._real32
    if(midpoint(2).gt.intf_loc(1)) &
         midpoint(2) = midpoint(2) - 1._real32


!!!-----------------------------------------------------------------------------
!!! Set up the list of atoms either side of the interface
!!!-----------------------------------------------------------------------------
    allocate(tmp_list1(bas%natom))
    allocate(tmp_list2(bas%natom))
    allocate(tmp_dist_list1(bas%natom))
    allocate(tmp_dist_list2(bas%natom))
    allocate(spec_list(bas%nspec))
    spec_list = 0
    nbelow = 0
    nabove = 0
    do is=1,bas%nspec
       do ia=1,bas%spec(is)%num

          if(bas%spec(is)%atom(ia,axis).lt.intf_loc(1).and.&
               bas%spec(is)%atom(ia,axis).ge.midpoint(2))then

             nbelow = nbelow + 1
             tmp_dist_list1(nbelow) = intf_loc(1) - bas%spec(is)%atom(ia,axis)
             tmp_list1(nbelow) = spec_list(is) + ia
          elseif(bas%spec(is)%atom(ia,axis).ge.intf_loc(1).and.&
               bas%spec(is)%atom(ia,axis).lt.midpoint(1))then

             nabove = nabove + 1
             tmp_dist_list2(nabove) = intf_loc(1) - bas%spec(is)%atom(ia,axis)
             tmp_list2(nabove) = spec_list(is) + ia
          end if

       end do
       if(is.ne.bas%nspec)then
          spec_list(is+1) = spec_list(is) + bas%spec(is)%num
       end if
    end do


!!!-----------------------------------------------------------------------------
!!! Move from temp to permanent arrays and sort
!!!-----------------------------------------------------------------------------
    allocate(lw_dist_list(nbelow))
    lw_dist_list(:nbelow) = tmp_dist_list1(:nbelow)
    call sort1D(lw_dist_list,tmp_list1(:nbelow))

    allocate(up_dist_list(nabove))
    up_dist_list(:nabove) = tmp_dist_list2(:nabove)
    call sort1D(up_dist_list,tmp_list2(:nabove))


!!!-----------------------------------------------------------------------------
!!! Set up lower and upper atom lists for later reference
!!!-----------------------------------------------------------------------------
    allocate(lw_list(nbelow,2))
    do i=1,nbelow
       lw_list(i,1) = minloc(tmp_list1(i) - spec_list(:), dim = 1,&
            mask = tmp_list1(i) - spec_list(:).gt.0)
       lw_list(i,2) = tmp_list1(i) - spec_list(lw_list(i,1))
    end do


    allocate(up_list(nabove,2))
    do i=1,nabove
       up_list(i,1) = minloc(tmp_list2(i) - spec_list(:), dim = 1,&
            mask = tmp_list2(i) - spec_list(:).gt.0)
       up_list(i,2) = tmp_list2(i) - spec_list(up_list(i,1))
    end do


!!!-----------------------------------------------------------------------------
!!! SET UP AN ATOM MAPPING HERE!!!!!
!!!-----------------------------------------------------------------------------
    

!!!-----------------------------------------------------------------------------
!!! Set up the weightings and closeness lists 
!!!-----------------------------------------------------------------------------
    allocate(lw_weight_list(nbelow))
    allocate(lw_close_list(nbelow))
    lw_weight_list(1) = gauss(pos=lw_dist_list(1),centre=0._real32,sigma=sigma)
    lw_close_list(1) = count(abs(lw_dist_list(1) - lw_dist_list(:nbelow)).le.rtol)
    do i=2,nbelow

       lw_weight_list(i) = lw_weight_list(i-1) + gauss(pos=lw_dist_list(i),centre=0._real32,sigma=sigma)
       lw_close_list(i) = count(abs(lw_dist_list(i) - lw_dist_list(:nbelow)).le.rtol)
       
    end do
    lw_weight_list = lw_weight_list/lw_weight_list(nbelow)

    allocate(up_weight_list(nabove))
    allocate(up_close_list(nabove))
    up_weight_list(1) = gauss(pos=up_dist_list(1),centre=0._real32,sigma=sigma)
    up_close_list(1) = count(abs(up_dist_list(1) - up_dist_list(:nabove)).le.rtol)
    do i=2,nabove

       up_weight_list(i) = up_weight_list(i-1) + gauss(pos=up_dist_list(i),centre=0._real32,sigma=sigma)
       up_close_list(i) = count(abs(up_dist_list(i) - up_dist_list(:nabove)).le.rtol)
       
    end do
    up_weight_list = up_weight_list/up_weight_list(nabove)
    


  end subroutine check_intf_depth
!!!#############################################################################


!!!#############################################################################
!!! 
!!!#############################################################################
  subroutine rand_swap_depth(bas,swap_bas,&
       nswaps_per_cell,sigma,small_sigma,&
       spec_list,&
       lw_list,up_list,&
       lw_dist_list,up_dist_list,&
       lw_close_list,up_close_list,&
       lw_weight_list,up_weight_list, &
       verbose &
  )
    implicit none
    integer :: i,loc1,loc2
    integer :: nbelow,nabove
    integer :: lw_mirror,up_mirror
    real(real32) :: r_rand1,r_rand2
    integer, allocatable, dimension(:) :: lw_convert,up_convert
    integer, allocatable, dimension(:,:) :: swap_list
    real(real32), allocatable, dimension(:) :: tlw_weight_list,tup_weight_list

    integer, dimension(:), intent(in) :: spec_list
    integer, dimension(:), intent(in) :: lw_close_list,up_close_list
    real(real32), dimension(:), intent(in) :: lw_dist_list,up_dist_list
    real(real32), dimension(:), intent(in) :: lw_weight_list,up_weight_list
    integer, dimension(:,:), intent(in) :: lw_list,up_list

    real(real32), intent(in) :: sigma,small_sigma
    type(basis_type), intent(inout) :: swap_bas
    integer, intent(in) :: nswaps_per_cell
    type(basis_type), intent(in) :: bas
    integer, intent(in) :: verbose


! make a list of natoms long, with each location pointing to a specific atomic species and number
! order the list based on distance from the interface
! choose a random number from that list ( nint(randnum()*) )

! use 1/distance as a weighting.

! exp(distance


!!!-----------------------------------------------------------------------------
!!! Identifies nbelow and above
!!!-----------------------------------------------------------------------------
    nbelow = size(lw_close_list)
    nabove = size(up_close_list)


!!!-----------------------------------------------------------------------------
!!! Allocates and initialises arrays
!!!-----------------------------------------------------------------------------
    allocate(tlw_weight_list, source=lw_weight_list)
    allocate(tup_weight_list, source=up_weight_list)


!!!-----------------------------------------------------------------------------
!!! Set up the converter lists between the unswapped and swapped system
!!!-----------------------------------------------------------------------------
    allocate(lw_convert(nbelow))
    do i=1,nbelow
       lw_convert(i) = i
    end do

    allocate(up_convert(nabove))
    do i=1,nabove
       up_convert(i) = i
    end do


!!!-----------------------------------------------------------------------------
!!! Choose swapping sets
!!!-----------------------------------------------------------------------------
    !allocate(swap_list(nswaps_per_cell,2))
    allocate(swap_list(min(nswaps_per_cell,nabove,nbelow),2))
    i = 1
    swap_loop: do while(i.le.nswaps_per_cell)
       if(nbelow.eq.1)then
          r_rand1 = 1
       else
          call random_number(r_rand1)
       end if
       if(nabove.eq.1)then
          r_rand2 = 1
       else
          call random_number(r_rand2)
       end if

       if(all(tlw_weight_list(:nbelow).le.0).or.nbelow.eq.0)then
          write(0,*) "Reached max number of atoms to remove from below"
          exit swap_loop
       end if
       loc1 = minloc(tlw_weight_list(:nbelow)-r_rand1,dim=1,mask=tlw_weight_list(:nbelow)-r_rand1.ge.0.0)
       if(loc1.eq.0)then
          loc1 = maxloc(tlw_weight_list(:nbelow)-r_rand1,dim=1,mask=tlw_weight_list(:nbelow).ge.0.0)
       end if
       swap_list(i,1) = lw_convert(loc1)

       if(loc1.ne.nbelow)then
          tlw_weight_list(loc1:) = tlw_weight_list(loc1+1:)
          lw_convert(loc1:) = lw_convert(loc1+1:)
          !lw_dist_list(loc1:) = lw_dist_list(loc1+1:)
       end if
       tlw_weight_list(nbelow) = 0
       lw_convert(nbelow) = 0
       nbelow = nbelow - 1


       if(all(tup_weight_list(:nabove).le.0).or.nabove.eq.0)then
          write(0,*) "Reached max number of atoms to remove from above"
          exit swap_loop
       end if
       loc2 = minloc(tup_weight_list(:nabove)-r_rand2,dim=1,mask=tup_weight_list(:nabove)-r_rand2.ge.0.0)
       if(loc2.eq.0)then
          loc2 = maxloc(tup_weight_list(:nabove)-r_rand2,dim=1,mask=tup_weight_list(:nabove).ge.0.0)
       end if
       swap_list(i,2) = up_convert(loc2)

       if(loc2.ne.nabove)then
          tup_weight_list(loc2:) = tup_weight_list(loc2+1:)
          up_convert(loc2:) = up_convert(loc2+1:)
          !up_dist_list(loc2:) = up_dist_list(loc2+1:)
       end if
       tup_weight_list(nabove) = 0
       up_convert(nabove) = 0
       nabove = nabove - 1


       tlw_weight_list(:nbelow) = recalc_rand_distrib(&
            lw_dist_list,&
            lw_convert(:nbelow),&
            lw_close_list,&
            swap_list(:i,1),&
            sigma,small_sigma)

       tup_weight_list(:nabove) = recalc_rand_distrib(&
            up_dist_list,&
            up_convert(:nabove),&
            up_close_list,&
            swap_list(:i,2),&
            sigma,small_sigma)
       
       if(verbose.ge.1) &
            write(0,'(&
            &I0,"th swap is ",I0,&
            &" with ",I0," at distances ",F7.3," and ",F7.3)') &
            i,swap_list(i,:),&
            lw_dist_list(swap_list(i,1)),up_dist_list(swap_list(i,2))
       i = i + 1

    end do swap_loop


    do i=1,nswaps_per_cell
       swap_bas%spec(lw_list(swap_list(i,1),1))%atom(lw_list(swap_list(i,1),2),:)=&
            bas%spec(up_list(swap_list(i,2),1))%atom(up_list(swap_list(i,2),2),:)
       swap_bas%spec(up_list(swap_list(i,2),1))%atom(up_list(swap_list(i,2),2),:)=&
            bas%spec(lw_list(swap_list(i,1),1))%atom(lw_list(swap_list(i,1),2),:)

       if(lmirror)then
          lw_mirror=bas_map%spec(lw_list(swap_list(i,1),1))%atom(lw_list(swap_list(i,1),2))
          up_mirror=bas_map%spec(up_list(swap_list(i,2),1))%atom(up_list(swap_list(i,2),2))

          swap_bas%spec(lw_list(swap_list(i,1),1))%atom(lw_mirror,:)=&
               bas%spec(up_list(swap_list(i,2),1))%atom(up_mirror,:)
          swap_bas%spec(up_list(swap_list(i,2),1))%atom(up_mirror,:)=&
               bas%spec(lw_list(swap_list(i,1),1))%atom(lw_mirror,:)
       end if

    end do


       ! Reduce the chance of removing an atom within 1 Å of the one just removed.
       ! If 1/2 of all gone from that region, stop anymore from being removed.
       ! Learn the number of atoms in that region and use the gaussian normalised to gaussian val at that point divided by the total number in that 1 Å region.

       ! remove the selected one from the lw_dist_list and lw_weight_list
       ! remove its effect on the lw_close_list also?


       ! we will have an issue with how to subtract from the weighting of similar level ions
       ! so, don't recalc lw_weight_list every time?
       ! instead, just undo the renormalisation, remove the effects of others on each part and
       ! then renormalise.
       ! have lw_weight_list be a changeable length. Have a list that points to the atomic number of ...
       ! ... each location in weight list.
       ! have a lw_dist_list and a reduced_lw_dist_list.
       ! have a lw_convert list between the two.

  contains
!!!-----------------------------------------------------------------------------
!!! Internal functions
!!!-----------------------------------------------------------------------------
    function recalc_rand_distrib(dist_list,conversion,close_list,swap_list,sigma,small_sigma) result(new_list)
      implicit none
      integer :: i,j
      integer :: nswaps,num
      real(real32) :: small_sigma

      real(real32), intent(in) :: sigma
      integer, dimension(:),intent(in) :: close_list,swap_list,conversion
      real(real32), dimension(:),intent(in) :: dist_list
      real(real32), allocatable, dimension(:) :: new_list


      num = size(conversion)
      nswaps = size(swap_list)

      allocate(new_list(num))
      do i=1,num

         new_list(i) = gauss(pos=dist_list(conversion(i)),centre=0._real32,sigma=sigma)

         do j=1,nswaps

            new_list(i) = new_list(i) - &
                 gauss(&
                 pos=dist_list(i),&
                 centre=dist_list(swap_list(j)),&
                 sigma=small_sigma)/close_list(conversion(i))

         end do

         if(i.ne.1) new_list(i) = new_list(i-1) + new_list(i)

      end do
      new_list = new_list/new_list(num)


    end function recalc_rand_distrib
!!!-----------------------------------------------------------------------------

  end subroutine rand_swap_depth
!!!#############################################################################



end module swapping
