!!!#############################################################################
!!! Code written by Isiah Edward Mikel Rudkin and Ned Thaddeus Taylor
!!! Code part of the ARTEMIS group (Hepplestone research group).
!!! Think Hepplestone, think HRG.
!!!#############################################################################
module shifting
  use artemis__constants, only: real32, pi, INF
  use misc_maths, only: get_nth_plane
  use misc_linalg, only: modu
  use artemis__geom_rw, only: basis_type,geom_write
  use artemis__geom_utils, only: split_bas,get_centre_atom,set_vacuum,shifter
  use artemis__io_utils
  use artemis__io_utils_extd, only: err_abort_print_struc
  use artemis__interface_identifier
  implicit none

  real(real32) :: f_scale = 0.5_real32
  real(real32) :: g_scale = 8._real32/3._real32

  private

  type bulk_DON_type
     type(den_of_neigh_type), allocatable, dimension(:) :: spec
  end type bulk_DON_type

  type map_type   
     integer, allocatable, dimension(:,:,:) :: spec       
  end type map_type


  public :: get_fit_shifts,get_descriptive_shifts,get_shifts_DON,bulk_DON_type


!!!updated  2020/02/25

  
contains 
!!!#############################################################################
!!! generates the top and bot bases near the interface
!!!#############################################################################
  subroutine get_top_bot_basis(lat,bas,bas_top,bas_bot,axis,intf_loc,depth)
    implicit none
    integer :: i,is,ia,itop,ibot,axis,count1
    real(real32) :: centre,dist,dist_max
    real(real32), optional :: depth
    type(basis_type) :: bas,bas_top,bas_bot
    real(real32), dimension(:) :: intf_loc
    real(real32), dimension(3,3) :: lat
    integer, allocatable, dimension(:) :: vtmp1
    integer, allocatable, dimension(:,:) :: intf_list
    real(real32), allocatable, dimension(:,:) :: regions
    type(basis_type), allocatable, dimension(:) :: splitbas


!!!-----------------------------------------------------------------------------
!!! allocating basis top and bottom
!!!-----------------------------------------------------------------------------
    allocate(bas_top%spec(bas%nspec))
    allocate(bas_bot%spec(bas%nspec))
    bas_top%nspec = bas%nspec
    bas_bot%nspec = bas%nspec

    bas_top%spec(:)%num = 0 !setting the number of atoms in each species to 0.
    bas_bot%spec(:)%num = 0 ! ""


    if(present(depth))then
       centre=intf_loc(1)
       dist=depth/modu(lat(axis,:))
!!!-----------------------------------------------------------------------------
!!! Loop to find the number of each species in each plane
!!!-----------------------------------------------------------------------------
       LOOP101: do is=1,bas%nspec ! Looping over the species
          LOOP102: do ia=1,bas%spec(is)%num ! Looping over the atom
             IF101: if ( ( bas%spec(is)%atom(ia,axis).gt.centre ) .and. &
                  ( bas%spec(is)%atom(ia,axis).le.(centre+dist) ) ) then ! checking if atom in the top plane.
                bas_top%spec(is)%num=bas_top%spec(is)%num+1 !update the number of the given species in the top plane.
             else if ( ( bas%spec(is)%atom(ia,axis).lt.centre ).and. &
                  ( bas%spec(is)%atom(ia,axis).ge.(centre-dist) ) ) then ! checking if atom in the bottom plane.
                bas_bot%spec(is)%num=bas_bot%spec(is)%num+1 !update the number of the given species in the bottom plane.
             end if IF101
          end do LOOP102
       end do LOOP101
       bas_top%natom = sum(bas_top%spec(:)%num)
       bas_bot%natom = sum(bas_bot%spec(:)%num)


!!!-----------------------------------------------------------------------------
!!! Allocate the required space in each bas and species.
!!!-----------------------------------------------------------------------------
       LOOP105: do is=1,bas%nspec
          allocate(bas_top%spec(is)%atom(bas_top%spec(is)%num,3))
          allocate(bas_bot%spec(is)%atom(bas_bot%spec(is)%num,3))
          bas_top%spec(is)%atom(:,:)=0._real32
          bas_bot%spec(is)%atom(:,:)=0._real32
       end do LOOP105


!!!-----------------------------------------------------------------------------
!!! Loop to add all the info from bas into bas_top and bas_bot
!!!-----------------------------------------------------------------------------
       LOOP103: do is=1,bas%nspec ! Looping over the species
          itop = 0
          ibot = 0
          LOOP104: do ia=1,bas%spec(is)%num ! Looping over the atom
             IF102: if ( ( bas%spec(is)%atom(ia,axis).gt.centre ).and. &
                  ( bas%spec(is)%atom(ia,axis).le.(centre+dist) ) ) then ! checking if atom in the top plane.
                itop = itop + 1
                bas_top%spec(is)%atom(itop,:) = bas%spec(is)%atom(ia,:)
             else if ( ( bas%spec(is)%atom(ia,axis).lt.centre ).and. &
                  ( bas%spec(is)%atom(ia,axis).ge.(centre-dist) ) ) then ! checking if atom in the bottom plane.
                ibot = ibot + 1
                bas_bot%spec(is)%atom(ibot,:) = bas%spec(is)%atom(ia,:)
             end if IF102
          end do LOOP104
       end do LOOP103
    else
       dist_max=4._real32/modu(lat(axis,:))
       allocate(vtmp1(bas%nspec))
       allocate(regions(size(intf_loc,dim=1),2))
       regions(1,1:2)=intf_loc(1:2)
       regions(2,1:2)=intf_loc(2:1:-1)
       splitbas=split_bas(bas,regions,axis)


       !!-----------------------------------------------------------------------
       !! Finds lower interfacial atoms near interface defined by intf_loc(1)
       !!-----------------------------------------------------------------------
       intf_list=gen_DONsim(gen_DON(lat,splitbas(1)),cutoff=4._real32)
101    do is=1,bas%nspec
          bas_bot%sysname=splitbas(1)%sysname
          bas_bot%spec(is)%name=splitbas(1)%spec(is)%name
          bas_bot%spec(is)%num = 0
          countloop1: do i=1,size(intf_list(:,1))
             if(intf_list(i,1).ne.is) cycle countloop1
             if(abs(splitbas(1)%spec(is)%atom(intf_list(i,2),axis)-&
                  intf_loc(1)).le.dist_max)then
                bas_bot%spec(is)%num = bas_bot%spec(is)%num + 1
             end if
          end do countloop1
          allocate(bas_bot%spec(is)%atom(bas_bot%spec(is)%num,3))
          count1=0
          atomloop1: do i=1,size(intf_list(:,1))
             if(intf_list(i,1).ne.is) cycle atomloop1
             if(abs(splitbas(1)%spec(is)%atom(intf_list(i,2),axis)-&
                  intf_loc(1)).le.dist_max)then
                count1=count1+1
                bas_bot%spec(is)%atom(count1,:3) = &
                     splitbas(1)%spec(is)%atom(intf_list(i,2),:3)
             end if
          end do atomloop1
       end do
       bas_bot%natom=sum(bas_bot%spec(:)%num)
       !!-----------------------------------------------------------------------
       !! If no lw_interfacial atoms found, redoes DONsim using averaging ...
       !! ... method 2
       !!-----------------------------------------------------------------------
       if(bas_bot%natom.eq.0)then
          intf_list=gen_DONsim(gen_DON(lat,splitbas(1)),cutoff=4._real32,avg_mthd=2)
          do is=1,bas%nspec
             deallocate(bas_bot%spec(is)%atom)
          end do
          goto 101
       end if


       !!-----------------------------------------------------------------------
       !! Finds upper interfacial atoms near interface defined by intf_loc(1)
       !!-----------------------------------------------------------------------
       intf_list=gen_DONsim(gen_DON(lat,splitbas(2)),cutoff=4._real32)
102    do is=1,bas%nspec
          bas_top%sysname=splitbas(2)%sysname
          bas_top%spec(is)%name=splitbas(2)%spec(is)%name
          bas_top%spec(is)%num = 0
          countloop2: do i=1,size(intf_list(:,1))
             if(intf_list(i,1).ne.is) cycle countloop2
             if(abs(splitbas(2)%spec(is)%atom(intf_list(i,2),axis)-&
                  intf_loc(1)).le.dist_max)then
                bas_top%spec(is)%num = bas_top%spec(is)%num + 1
             end if
          end do countloop2
          allocate(bas_top%spec(is)%atom(bas_top%spec(is)%num,3))
          count1=0
          atomloop2: do i=1,size(intf_list(:,1))
             if(intf_list(i,1).ne.is) cycle atomloop2
             if(abs(splitbas(2)%spec(is)%atom(intf_list(i,2),axis)-&
                  intf_loc(1)).le.dist_max)then
                count1=count1+1
                bas_top%spec(is)%atom(count1,:3) = &
                     splitbas(2)%spec(is)%atom(intf_list(i,2),:3)
             end if
          end do atomloop2
       end do
       bas_top%natom=sum(bas_top%spec(:)%num)
       !!-----------------------------------------------------------------------
       !! If no up_interfacial atoms found, redoes DONsim using averaging ...
       !! ... method 2
       !!-----------------------------------------------------------------------
       if(bas_top%natom.eq.0)then
          intf_list=gen_DONsim(gen_DON(lat,splitbas(2)),cutoff=4._real32,avg_mthd=2)
          do is=1,bas%nspec
             deallocate(bas_top%spec(is)%atom)
          end do
          goto 102
       end if
    end if


    return
  end subroutine get_top_bot_basis
!!!#############################################################################


!!!#############################################################################
!!! Function that figures out the best shift for the planes given the ...
!!! ... required minimum bulk bond length.            
!!!#############################################################################
  function get_fit_shifts(lat,bas,bond,axis,intf_loc,depth,nstore,itmp1,itmp2) result(best_shifts)
    real(real32) :: depth,bond ! the depth into the material we are interested (physical size in the c direction).
    integer :: i
    type(basis_type) :: bas_bot,bas_top

    real(real32) :: depth_bascoord
    real(real32), dimension(:) :: intf_loc
    real(real32), allocatable, dimension(:,:) :: min_atom_sep
    real(real32), allocatable, dimension(:,:,:) :: avg_min_atom_sep

    integer :: axis
    integer :: num_steps,num_c_shifts !number of pieces to divide the unit cell into in a and b direction.
    real(real32), allocatable, dimension(:,:) :: best_shifts
    integer :: nstore ! The required output number of the best shifts.

    integer, optional :: itmp1,itmp2
    type(basis_type) :: bas !The basis input by interfaces.f90
    real(real32), dimension(3,3) :: lat !The lattice input by interfaces.f90


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Setting up variables from interfaces.f90 !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    num_steps = 100
    num_c_shifts = 10
    if(present(itmp1)) num_steps=itmp1
    if(present(itmp2)) num_c_shifts=itmp1
    allocate(min_atom_sep(num_steps,num_steps))
    allocate(best_shifts(nstore,4))


    if(depth.eq.0._real32)then
       call get_top_bot_basis(lat,bas,bas_top,bas_bot,axis,intf_loc)
    else
       call get_top_bot_basis(lat,bas,bas_top,bas_bot,axis,intf_loc,depth)
    end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Looping through different plane allignments (shifts) and finding the best one !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    allocate(avg_min_atom_sep(num_steps,num_steps,num_c_shifts))
    avg_min_atom_sep = avgminsep(lat,bas_top,bas_bot,num_steps,num_c_shifts,depth_bascoord)



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Finding the best match from our shifts.   !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    best_shifts = findbestfits(bond,avg_min_atom_sep,num_steps,num_c_shifts,nstore,depth_bascoord)
    best_shifts(:,axis)=best_shifts(:,axis)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!       - Output-        !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! do j=1,num_c_shifts !output to text files
    !    unit=100+j
    !    open(unit=unit,file="output1.txt")
    !    write(unit,'(100(F0.8,2X))') (avg_min_atom_sep(:,i,j),i=1,num_steps)
    !    close(unit)
    ! end do

    write(*,'(4(F0.5,2X))') (best_shifts(i,:),i=1,nstore)

  end function get_fit_shifts
!!!#########################################################################


!!!#########################################################################
!!! Function that finds the best match between the average interface ...
!!! ... minimum bond length and the bulk minimum bond length            
!!!#########################################################################
  function findbestfits(bulkbond,avg_min_sep,num_steps,num_c_shifts,num_best_shifts,depth) result(best_shifts)
    implicit none
    real(real32), dimension(:,:,:) :: avg_min_sep
    real(real32) :: bulkbond,current_difference,min_difference,depth
    integer :: i,ia,ib,ic,num_steps,num_c_shifts,num_best_shifts,c_shift_low,c_shift_high
    real(real32), allocatable, dimension(:,:) :: best_shifts
    integer, dimension(3) :: placeholder

    allocate(best_shifts(num_best_shifts,4))

    if (mod(num_c_shifts,2) .eq. 0) then
       c_shift_low = -nint(real(num_c_shifts)/2.0)+1
       c_shift_high = nint(real(num_c_shifts)/2.0)
    else 
       c_shift_low = -(floor(real(num_c_shifts)/2.0))
       c_shift_high = (floor(real(num_c_shifts)/2.0))
    end if
    shiftloop: do i=1,num_best_shifts 

       placeholder = -1
       min_difference = huge(0._real32)
       LOOP5A: do ia=0,num_steps-1 !loop through shifts in a
          LOOP5B: do ib=0,num_steps-1 !loop through shifts in b
             LOOP5C: do ic=c_shift_low,c_shift_high,1 !Loop through shifts of the top plane in c
                current_difference = abs(avg_min_sep(ia+1,ib+1,ic-c_shift_low+1) - bulkbond)
                if (current_difference.lt.min_difference) then
                   min_difference = current_difference
                   best_shifts(i,1) = real(ia,real32)/real(num_steps,real32)
                   best_shifts(i,2) = real(ib,real32)/real(num_steps,real32)
                   best_shifts(i,3) = real(ic,real32)*depth*2._real32/real(num_c_shifts,real32)
                   best_shifts(i,4) = min_difference
                   placeholder(1) = ia
                   placeholder(2) = ib
                   placeholder(3) = ic
                end if
             end do LOOP5C
          end do LOOP5B
       end do LOOP5A
       if(any(placeholder.eq.-1)) then
          write(0,*) "ERROR: No shifts found for the given interface"
          stop
       end if
       avg_min_sep(placeholder(1)+1,placeholder(2)+1,placeholder(3)-c_shift_low+1) = huge(0._real32)
    end do shiftloop

  end function findbestfits
!!!#############################################################################


!!!#############################################################################
!!! Subroutine that finds the average minimum atomic seperation ...
!!! ... between any atoms in the top and bottom planes
!!!#############################################################################
  function avgminsep(lat,plane_up,plane_dw,num_steps,num_c_shifts,depth) result(avg_min_sep)
    implicit none
    type(basis_type) :: plane_up,plane_dw,tplane_up,tplane_dw
    real(real32) :: avg_sep_up,avg_sep_dw,depth
    integer :: num_steps,num_c_shifts !number of pieces to divide the unit cell into in a and b direction.
    real(real32), allocatable, dimension(:,:,:) :: avg_min_sep
    integer :: ia,ib,ic,is_up,ia_up,c_shift_low,c_shift_high
    real(real32), dimension(3,3) :: lat


    allocate(avg_min_sep(num_steps,num_steps,num_c_shifts))
    call tplane_up%copy(plane_up)
    call tplane_dw%copy(plane_dw)
    if (mod(num_c_shifts,2) .eq. 0) then
       c_shift_low = -nint(real(num_c_shifts)/2.0)+1
       c_shift_high = nint(real(num_c_shifts)/2.0)
    else 
       c_shift_low = -(floor(real(num_c_shifts)/2.0))
       c_shift_high = (floor(real(num_c_shifts)/2.0))
    end if


    avg_min_sep = huge(0._real32)
    LOOP4C: do ic=c_shift_low,c_shift_high,1 !Loop through shifts of the top plane in c
       LOOP4A: do ia=0,num_steps-1 !loop through shifts in a
          LOOP4B: do ib=0,num_steps-1 !loop through shifts in b

             do is_up=1,plane_up%nspec
                do ia_up=1,plane_up%spec(is_up)%num
                   plane_up%spec(is_up)%atom(ia_up,:) = &
                        plane_up%spec(is_up)%atom(ia_up,:) + &
                        (/&
                        (real(ia,real32)/real(num_steps,real32)),&
                        (real(ib,real32)/real(num_steps,real32)),&
                        (real(ic,real32)*depth*2._real32/real(num_c_shifts,real32)) /)
                end do
             end do

             avg_sep_up = find_avg_min_sep(lat,tplane_up, plane_dw)
             avg_sep_dw = find_avg_min_sep(lat, plane_dw,tplane_up)

             avg_min_sep(ia+1,ib+1,ic-c_shift_low+1) = &
                  (avg_sep_up + avg_sep_dw)/2._real32

          end do LOOP4B
       end do LOOP4A
    end do LOOP4C


  end function avgminsep
!!!#############################################################################


!!!#############################################################################
!!! finds average minimum separation between two planes
!!!#############################################################################
  function find_avg_min_sep(lat,plane_1,plane_2) result(avg_min_sep)
    implicit none
    integer :: is_1,ia_1,is_2,ia_2,j
    real(real32) :: avg_min_sep,min_sep,cur_sep
    type(basis_type) :: plane_1,plane_2
    real(real32), dimension(3) :: dvtmp1
    real(real32), dimension(3,3) :: lat


    avg_min_sep=0._real32
    LOOP401: do is_1=1,plane_1%nspec ! Loop though 1st plane
       LOOP402: do ia_1=1,plane_1%spec(is_1)%num

          min_sep = huge(0._real32)
          LOOP403: do is_2=1,plane_2%nspec ! Loop through 2nd plane
             LOOP404: do ia_2=1,plane_2%spec(is_2)%num

                dvtmp1 = &
                     plane_2%spec(is_2)%atom(ia_2,:) - &
                     plane_1%spec(is_1)%atom(ia_1,:)

                do j=1,3
                   dvtmp1(j) = dvtmp1(j) - ceiling( dvtmp1(j) - 0.5_real32 )
                end do
                dvtmp1 = dvtmp1(1) * lat(1,:) &
                     + dvtmp1(2) * lat(2,:) &
                     + dvtmp1(3) * lat(3,:)

                cur_sep = modu( dvtmp1 )


                if (cur_sep.lt.min_sep) min_sep = cur_sep
             end do LOOP404
          end do LOOP403
          avg_min_sep = avg_min_sep + min_sep
       end do LOOP402
    end do LOOP401

    avg_min_sep = avg_min_sep/plane_1%natom


  end function find_avg_min_sep
!!!#############################################################################



!!!#############################################################################
!!!#############################################################################
!!! M E T H O D   2
!!!#############################################################################
!!!#############################################################################


!!!#############################################################################
!!! Finds best c axis separation, then finds the most descriptive set of ...
!!! ... shifts for that separation (i.e. the ones that fit the best and ...
!!! ... worst to that of the average bulk bond).
!!! Outputs best, then worst, then 2nd best, then 2nd worst, etc.
!!!#############################################################################
  function get_descriptive_shifts(lat,bas,bond,axis,intf_loc,depth,nstore,c_scale,lprint) result(res_shifts)
    implicit none
    integer :: is
    integer :: nstore,axis,num_steps
    real(real32) :: bond,depth,cur_vac,c_shift
    type(basis_type) :: bas,bas_bot,bas_top
    real(real32), dimension(3,3) :: lat
    real(real32), dimension(:) :: intf_loc
    real(real32), allocatable, dimension(:) :: specval_bot,specval_top
    real(real32), allocatable, dimension(:,:) :: res_shifts
    real(real32), optional :: c_scale
    logical, optional :: lprint


    num_steps = 50
    allocate(res_shifts(nstore,3))
!!!-----------------------------------------------------------------------------
!!! separates basis into atoms above and below interface within a depth window
!!!-----------------------------------------------------------------------------
    if(depth.eq.0._real32)then
       call get_top_bot_basis(lat,bas,bas_top,bas_bot,axis,intf_loc)
    else
       call get_top_bot_basis(lat,bas,bas_top,bas_bot,axis,intf_loc,depth=depth)
    end if


!!!-----------------------------------------------------------------------------
!!! finds the current vacuum separation at the interface
!!!-----------------------------------------------------------------------------
    allocate(specval_bot(bas%nspec))
    allocate(specval_top(bas%nspec))
    specval_bot=-huge(0._real32)
    specval_top=huge(0._real32)
    do is=1,bas%nspec
       if(bas_bot%spec(is)%num.ne.0)then
          specval_bot(is)=maxval(bas_bot%spec(is)%atom(:,axis))
       end if
       if(bas_top%spec(is)%num.ne.0)then
          specval_top(is)=minval(bas_top%spec(is)%atom(:,axis))
       end if
    end do
    cur_vac=(minval(specval_top)-maxval(specval_bot))*modu(lat(axis,:))


!!!-----------------------------------------------------------------------------
!!! finds optimal separation for interface (based on average min bulk bond idea)
!!!-----------------------------------------------------------------------------
    do is=1,bas_top%nspec
       bas_top%spec(is)%atom(:,axis) = &
            bas_top%spec(is)%atom(:,axis) + (bond - cur_vac)/modu(lat(axis,:))
    end do
    c_shift = get_c_shift(lat,bas_top,bas_bot,bond,axis,num_steps)
    do is=1,bas_top%nspec
       bas_top%spec(is)%atom(:,axis) = bas_top%spec(is)%atom(:,axis) + c_shift
    end do


!!!-----------------------------------------------------------------------------
!!! finds descriptive set of shifts parallel to interface for supplied c shift
!!!-----------------------------------------------------------------------------
    !res_shifts(:,3) = c_shift + (bond - cur_vac)/modu(lat(axis,:))
    res_shifts(:,3) = c_shift + bond/modu(lat(axis,:))
    res_shifts(:,1:2) = get_descriptive_ab_shifts(lat,bas_top,bas_bot,bond,axis,nstore,num_steps)
    if(present(c_scale)) res_shifts(:,3) = res_shifts(:,3) * c_scale


    if(present(lprint))then
       if(lprint)then
          write(*,'(1X,"Shifts to be applied (Å)")')
          do is=1,nstore
             write(*,*) res_shifts(is,1),res_shifts(is,2), &
                  res_shifts(is,3)*modu(lat(axis,:))
          end do
       end if
    end if

!!! 1st, get it to find best c axis shift to find the best shift.
!!! start off by having the shift reduce the gap to the size of the bond and search from there by step sizes depending on the difference



  end function get_descriptive_shifts
!!!#############################################################################



!!!#############################################################################
!!! Subroutine that finds the average minimum atomic seperation ...
!!! ... between any atoms in the top and bottom planes
!!!#############################################################################
  function get_c_shift(lat,plane_up,plane_dw,bond,axis,num_steps) result(c_shift)
    implicit none
    integer :: num_steps,count1
    integer :: ia,ib,is_up,ia_up,axis
    real(real32) :: avg_sep_up,avg_sep_dw,bond,tol
    real(real32) :: c_shift,prev_c_shift,new_c_shift
    real(real32) :: prev_min_bond,min_bond
    type(basis_type) :: plane_up,plane_dw,tplane_up
    real(real32), allocatable, dimension(:,:) :: avg_min_sep
    real(real32), dimension(3,3) :: lat


!!!-----------------------------------------------------------------------------
!!! Clone upper basis for editing
!!!-----------------------------------------------------------------------------
    call tplane_up%copy(plane_up)
    allocate(avg_min_sep(num_steps,num_steps))


!!!-----------------------------------------------------------------------------
!!! Initialise variables
!!!-----------------------------------------------------------------------------
    tol=1.E-2_real32/modu(lat(axis,:))
    count1=0
    prev_min_bond=0._real32
    prev_c_shift=0._real32
    c_shift=0._real32
    avg_min_sep = huge(0._real32)


!!!-----------------------------------------------------------------------------
!!! Loop to change c_shift in order to find optimal c_shift
!!!-----------------------------------------------------------------------------
    LOOP5C: do !!! FIX WHAT avg_min_sep IS ALLOCATED AS
       count1=count1+1

       LOOP5A: do ia=0,num_steps-1 !loop through shifts in a
          LOOP5B: do ib=0,num_steps-1 !loop through shifts in b

             do is_up=1,plane_up%nspec
                do ia_up=1,plane_up%spec(is_up)%num
                   tplane_up%spec(is_up)%atom(ia_up,:) = plane_up%spec(is_up)%atom(ia_up,:) + &
                        (/&
                        (real(ia,real32)/real(num_steps,real32)),&
                        (real(ib,real32)/real(num_steps,real32)),&
                        c_shift /)
                end do
             end do

             avg_sep_up = find_avg_min_sep(lat,tplane_up, plane_dw)
             avg_sep_dw = find_avg_min_sep(lat, plane_dw,tplane_up)
             avg_min_sep(ia+1,ib+1) = (avg_sep_up + avg_sep_dw)/2._real32

          end do LOOP5B
       end do LOOP5A


       min_bond =  minval(avg_min_sep(:,:))
       !!-----------------------------------------------------------------------
       !! Checks if Infinity has been encountered
       !!-----------------------------------------------------------------------
       if(min_bond.gt.INF)then
          write(0,*)  "ERROR: &
               &Encountered Infinity when attempting to find optimal c_shift"
          stop
       end if
       !!-----------------------------------------------------------------------
       !! Using previous and current min_bond values, it predicts best min_bond
       !!-----------------------------------------------------------------------
       if(count1.ne.1)then
          if(abs(prev_min_bond - min_bond).lt.1.D-4)then
             exit
          else
             new_c_shift = prev_c_shift - &
                  (prev_c_shift - c_shift)*( prev_min_bond - bond )/( prev_min_bond - min_bond )
          end if
       else
          new_c_shift = 0.5_real32/modu(lat(axis,:))
       end if
       !!-----------------------------------------------------------------------
       !! Breaks afer 50 failed steps
       !!-----------------------------------------------------------------------
       if(count1.ge.50)then
          call err_abort('ERROR: Internal error in get_c_shift\n&
               &get_c_shift in mod_shifting.f90 hasn'' worked to&
               &find a good shift.\nSuggest using ISHIFT≠3̄')
       end if
       !write(0,*) c_shift,min_bond

       prev_c_shift = c_shift
       prev_min_bond = min_bond
       c_shift = new_c_shift

    end do LOOP5C



  end function get_c_shift
!!!#############################################################################



!!!#############################################################################
!!! Subroutine that finds the average minimum atomic seperation ...
!!! ... between any atoms in the top and bottom planes
!!!#############################################################################
  function get_descriptive_ab_shifts(lat,plane_up,plane_dw,bond,axis,nstore,num_steps) result(ab_shifts)
    implicit none
    integer :: nstore,num_steps,count1
    integer :: ia,ib,is_up,ia_up,axis,iden,inum
    real(real32) :: avg_sep_up,avg_sep_dw,bond
    real(real32) :: min_sep,max_sep
    type(basis_type) :: plane_up,plane_dw,tplane_up,tplane_dw
    real(real32), allocatable, dimension(:,:) :: ab_shifts,avg_min_sep
    real(real32), dimension(3,3) :: lat


    call tplane_up%copy(plane_up)
    call tplane_dw%copy(plane_dw)
    allocate(avg_min_sep(num_steps,num_steps))
    allocate(ab_shifts(nstore,2))


    count1=0
    avg_min_sep = huge(0._real32)
    LOOP5A: do ia=0,num_steps-1 !loop through shifts in a
       LOOP5B: do ib=0,num_steps-1 !loop through shifts in b

          do is_up=1,plane_up%nspec
             do ia_up=1,plane_up%spec(is_up)%num
                tplane_up%spec(is_up)%atom(ia_up,:) = plane_up%spec(is_up)%atom(ia_up,:) + &
                     (/&
                     (real(ia,real32)/real(num_steps,real32)),&
                     (real(ib,real32)/real(num_steps,real32)),&
                     0._real32 /)
             end do
          end do

          avg_sep_up = find_avg_min_sep(lat,tplane_up, plane_dw)
          avg_sep_dw = find_avg_min_sep(lat, plane_dw,tplane_up)


          avg_min_sep(ia+1,ib+1) = (avg_sep_up + avg_sep_dw)/2._real32
       end do LOOP5B
    end do LOOP5A

    min_sep=minval(avg_min_sep)
    max_sep=maxval(avg_min_sep)
    !write(0,*) "min:",minloc(avg_min_sep),min_sep
    !write(0,*) "max:",maxloc(avg_min_sep),max_sep



    ab_shifts(1,:)=real((/minloc(avg_min_sep)/),real32)/real(num_steps,real32)
    ab_shifts(2,:)=real((/maxloc(avg_min_sep)/),real32)/real(num_steps,real32)
    iden=1
    count1=2
    denom_loop: do
       iden=iden*2
       if(iden.gt.nstore) exit denom_loop
       do inum=1,iden,2
          count1=count1+1
          if(count1.gt.nstore) exit denom_loop
          ab_shifts(count1,:) = real((/ &
               minloc( &
               abs( avg_min_sep - ( min_sep + (max_sep-min_sep)*real(inum,real32)/real(iden,real32) ) ) )&
               /),real32)/real(num_steps,real32)

       end do
    end do denom_loop

    !! REMOVE DUPLICATES

    !! HAVE IT CYCLE OVER VARYING INCREMENTS OF X AND Y, SOMEHOW
    !ab_shifts(3,1)=ab_shifts(1,1)
    !ab_shifts(3,2)=maxloc(avg_min_sep(ab_shifts(1,1),:))

    !ab_shifts(4,2)=ab_shifts(1,2)
    !ab_shifts(4,2)=maxloc(avg_min_sep(:,ab_shifts(1,2)))



  end function get_descriptive_ab_shifts
!!!#############################################################################



!!!#############################################################################
!!!#############################################################################
!!! M E T H O D   4
!!!#############################################################################
!!!#############################################################################



!!!#############################################################################
!!! generate shifts by filling missing neighours for surface atoms
!!!#############################################################################
  function get_shifts_DON(bas,axis,intf_loc,nstore,tol_sym,c_scale,offset,&
       bulk_DON,bulk_map,verbose,max_bondlength) result(res_shifts)
    use artemis__sym, only: gldfnd,confine_type
    use artemis__geom_utils, only: get_bulk,wyck_spec_type,get_wyckoff
    use artemis__interface_identifier, only: gen_single_DON,nstep_default,den_of_neigh_type
    implicit none
    type(basis_type), intent(in) :: bas
    !! Interface structure
    integer, intent(in) :: axis
    !! Axis of the interface
    real(real32), dimension(:), intent(in) :: intf_loc
    !! Location of the interfaces
    integer, intent(in) :: nstore
    !! Number of shifts to be generated
    real(real32), intent(in) :: tol_sym
    !! Tolerance for symmetry
    real(real32), intent(in), optional :: c_scale
    !! Scaling factor for the interface separation
    real(real32), dimension(3), optional, intent(in) :: offset
    !! Input offset of the two interface substructures
    integer, intent(in), optional :: verbose
    !! Boolean whether to print the shifts
    type(bulk_DON_type), dimension(:), optional, intent(in) :: bulk_DON
    !! Bulk DONs to be used for the interface
    integer, dimension(:,:,:), optional, intent(in) :: bulk_map
    !! Mapping of bulk atoms to the interface atoms
    real(real32), intent(in), optional :: max_bondlength
    !! Cutoff bondlength to consider first neighbours


    integer :: i,j,k,l,is,ia,ja,jb,jc,count1,itmp1
    integer :: ntrans,iatom,nneigh,ncheck
    integer :: verbose_
    real(real32) :: stepsize,max_sep,dist_max
    real(real32) :: rtmp1,rtmp2,rtmp3
    real(real32) :: val,dtmp1,dtmp2
    logical :: lbulk, lpresent
    type(confine_type) :: confine
    integer, dimension(2) :: plane_loc
    integer, dimension(3) :: ngrid,nstep,ivtmp1
    real(real32), dimension(2) :: min_trans,lowest_atom,highest_atom
    real(real32), dimension(3) :: pos,vtmp1,vtmp2,vtmp3,gridsize,add
    logical, dimension(2) :: lwyckoff
    type(map_type), dimension(2) :: map
    type(wyck_spec_type), dimension(2) :: wyckoff
    real(real32), allocatable, dimension(:) :: fit_store,tmp_neigh
    type(basis_type), allocatable, dimension(:) :: splitbas
    type(den_of_neigh_type), allocatable, dimension(:,:) :: DON_missing
    integer, allocatable, dimension(:,:) :: shift_store
    real(real32), allocatable, dimension(:,:) :: res_shifts,trans,regions



    !integer :: OMP_GET_NUM_THREADS,OMP_GET_MAX_THREADS,OMP_GET_THREAD_NUM,CHUNK
    !integer :: nthreads


    type neighbour_type
       integer :: num
       real(real32) :: bond
       real(real32), dimension(3) :: pos
    end type neighbour_type
    type(neighbour_type), allocatable, dimension(:,:) :: neighbour
    type intf_type
       type(neighbour_type), allocatable, dimension(:) :: neigh
    end type intf_type
    type(intf_type), dimension(2) :: intf


    type grid_type
       real(real32), allocatable, dimension(:) :: neigh
    end type grid_type
    type(grid_type), allocatable, dimension(:,:,:,:) :: course_grid



    verbose_ = 0
    if(present(verbose)) verbose_ = verbose
!!!-----------------------------------------------------------------------------
!!! check if bulk DONs supplied
!!!-----------------------------------------------------------------------------
    if(present(bulk_DON).and.present(bulk_map))then
       lbulk=.true.
       allocate(map(1)%spec,source=bulk_map)
    else
       lbulk=.false.
    end if


!!!-----------------------------------------------------------------------------
!!! sets up step size
!!!-----------------------------------------------------------------------------
    allocate(res_shifts(nstore,3))
    res_shifts=0._real32


!!!-----------------------------------------------------------------------------
!!! separates basis into atoms above and below interface within a depth window
!!!-----------------------------------------------------------------------------
    confine%l=.true.
    confine%axis=axis
    allocate(trans(minval(bas%spec(:)%num,dim=1),3))
    allocate(regions(size(intf_loc,dim=1),2))
    if(intf_loc(1).lt.intf_loc(2))then
       regions(2,1:2)=intf_loc(1:2)
       regions(1,1:2)=intf_loc(2:1:-1)
    else
       regions(1,1:2)=intf_loc(1:2)
       regions(2,1:2)=intf_loc(2:1:-1)
    end if
    if(lbulk)then
       splitbas=split_bas(bas,regions,axis,lall_same_nspec=.false.,&
            map1=map(1)%spec,map2=map(2)%spec)
    else
       splitbas=split_bas(bas,regions,axis)!,lall_same_nspec=.false.)
    end if


!!!-----------------------------------------------------------------------------
!!! determines repeated translations within the cell (reduces shift by that)
!!!-----------------------------------------------------------------------------
    min_trans=1._real32
    do i=1,2
       call gldfnd(confine, splitbas(i), splitbas(i), trans, ntrans, tol_sym)
       if(ntrans.eq.0) cycle
       do j=1,ntrans
          do k=1,2
             if(abs(trans(j,k)).lt.abs(min_trans(k))) min_trans(k)=trans(j,k)
          end do
       end do
    end do
    min_trans=abs(min_trans)
    where(abs(min_trans).lt.1.E-5_real32)
       min_trans=1._real32
    end where
    if(verbose_.eq.1) write(*,*) "repeated_trans:",min_trans


!!!-----------------------------------------------------------------------------
!!! If given bulk, then use the DOS' given
!!! Else, work out atom in slab that is the same
!!!-----------------------------------------------------------------------------
    if(.not.lbulk)then
       lwyckoff=.true.
       do i=1,2
          wyckoff(i)=get_wyckoff(splitbas(i),axis)
          if(.not.allocated(wyckoff(i)%spec))then
             write(*,'(1X,"Using centre atoms as bulk representation for parent slab", I0)') i
             lwyckoff(i)=.false.
          else
             write(*,'(1X,"Using Wyckoff atoms as bulk representation for parent slab", I0)') i
          end if
       end do
    else
       lwyckoff=.false.
    end if


!!!-----------------------------------------------------------------------------
!!! Evaluates DON for each atom.
!!! Determines whether it is an atom to consider by calculating its ...
!!! ... dissimilarity to that of a same-species atom in the centre of the slab.
!!! For dissimilar atoms, number of nearest "missing" bonds is stored.
!!!-----------------------------------------------------------------------------
    allocate(neighbour(2,bas%natom))
    if(present(max_bondlength))then
       dist_max = max_bondlength
    else
       dist_max = 4.0
    end if
    allocate(DON_missing(2,bas%nspec))
    if(verbose_.ge.1) write(*,*)
    region_loop: do i=1,2
       if(verbose_.ge.1) write(*,'&
            &(2X,"is",2X,"ia",4X,"nmissing",4X,"bond size (Å)",8X,"position")')

       count1 = 0
       DON_missing(i,:) = &
            gen_DON(bas%lat,splitbas(i),dist_max,scale_dist=.false.,norm=.true.)
       !!-----------------------------------------------------------------------
       !! Loops through the basis and finds the missing bonds of surface atoms.
       !! Does this by minusing the DON of the wyckoff atom of the surface ...
       !! ... from the surface atom's DON.
       !! Using that, it can find the missing DON peaks and select only the ...
       !! ... missing 1st nearest neighbour
       !!-----------------------------------------------------------------------
       spec_loop: do is=1,splitbas(i)%nspec
          if(splitbas(i)%spec(is)%num.le.0) cycle spec_loop
          if(.not.lbulk.and..not.lwyckoff(i))then
             iatom = get_centre_atom(&
                  splitbas(i),is,axis,lw=regions(i,1),up=regions(i,2))
             if(iatom.eq.0)&
                  call err_abort("ERROR: Internal error in get_shifts_DON\n&
                  &  No centre atom found in get_shifts_DON.",.true.)
          end if
          if(lbulk)then
             if(any(map(i)%spec(is,:splitbas(i)%spec(is)%num,:).le.0))then
                write(0,'("parent  species  atom")')
                write(0,'(2X,I2,6X,I2,4X,I4)') i,is,ia
                call err_abort("ERROR: Internal error in get_shifts_DON\n&
                     &  Mapping of bulk missing",.true.)
             end if
          end if
          atom_loop1: do ia=1,splitbas(i)%spec(is)%num
             !!-----------------------------------------------------------------
             !! check if atom is closer to this interface or the other
             !!-----------------------------------------------------------------
             dtmp1 = splitbas(i)%spec(is)%atom(ia,axis) - intf_loc(1)
             dtmp2 = splitbas(i)%spec(is)%atom(ia,axis) - intf_loc(2)
             if( abs(dtmp1 - ceiling(dtmp1 - 0.5_real32)) .gt. &
                  abs(dtmp2 - ceiling(dtmp2 - 0.5_real32)))then
                cycle atom_loop1
             end if


             !!-----------------------------------------------------------------
             !! evaluates the missing components of the atom's DON
             !!-----------------------------------------------------------------
             if(lbulk)then
                DON_missing(i,is)%atom(ia,:) = &
                     bulk_DON(i)%spec(map(i)%spec(is,ia,1))%atom(map(i)%spec(is,ia,2),:) - &
                     DON_missing(i,is)%atom(ia,:)
             elseif(lwyckoff(i))then
                DON_missing(i,is)%atom(ia,:) = &
                     DON_missing(i,is)%atom(wyckoff(i)%spec(is)%atom(ia),:) - &
                     DON_missing(i,is)%atom(ia,:)
             else
                DON_missing(i,is)%atom(ia,:) = &
                     DON_missing(i,is)%atom(iatom,:) - &
                     DON_missing(i,is)%atom(ia,:)
             end if
             !where(DON_missing(i,is)%atom(ia,:).lt.0._real32)
             !   DON_missing(i,is)%atom(ia,:)=0._real32
             !end where
             if(all(abs(DON_missing(i,is)%atom(ia,:)).lt.1.E-2_real32))&
                  cycle atom_loop1


             !!-----------------------------------------------------------------
             !! checks only 1st missing bond
             !!-----------------------------------------------------------------
             plane_loc(:)=&
                  get_nth_plane(invec=real(DON_missing(i,is)%atom(ia,:),real32),&
                  nth=2,window=20,is_periodic=.false.) !! WINDOW WAS 10, NOW 20
             itmp1=nint( &
                  sum(DON_missing(i,is)%atom(ia,:plane_loc(1)))*&
                  (dist_max/nstep_default) )
             if(itmp1.gt.0)then
                count1 = count1 +1
                neighbour(i,count1)%pos = splitbas(i)%spec(is)%atom(ia,:3)
                !neighbour(i,count1)%pos = neighbour(i,count1)%pos - &
                !     ceiling(neighbour(i,count1)%pos - 1._real32)
                neighbour(i,count1)%bond = &
                     ( maxloc(DON_missing(i,is)%atom(ia,:plane_loc(1)),dim=1) &
                     - 1 ) * dist_max/nstep_default
                neighbour(i,count1)%num = itmp1
                if(verbose_.ge.1)&
                     write(*,'(2X,I2,3X,I3,7X,I2,9X,F0.3,8X,3(1X,F5.2))') &
                     is,ia,&
                     neighbour(i,count1)%num,&
                     neighbour(i,count1)%bond,&
                     neighbour(i,count1)%pos
             elseif(itmp1.lt.0.and.lbulk)then
                write(0,'("parent  species  atom  nmissing")')
                write(0,'(2X,I2,6X,I2,4X,I4,4X,I4)') i,is,ia,itmp1
                write(0,'("species  atom")')
                write(0,'(2X,I2,4X,I4)') is, map(i)%spec(is,ia,2)
                write(0,*) "Writing failed DON to output file &
                     &'full_broken_DON.dat'"
                open(unit=13,file="full_broken_DON.dat")
                write(13,'("#nstep   nmissing   n_in_bulk")')
                do j=1,nstep_default
                   write(13,*) &
                        (j-1)*dist_max/nstep_default,&
                        DON_missing(i,is)%atom(ia,j),&
                        bulk_DON(i)%spec(map(i)%spec(is,ia,1))%atom(map(i)%spec(is,ia,2),j)
                end do
                close(13)
                write(0,*) "Writing failed DON to output file &
                     &'zoom_broken_DON.dat'"
                open(unit=14,file="zoom_broken_DON.dat")
                write(14,'("#nstep   nmissing   n_in_bulk")')
                do j=1,plane_loc(1)
                   write(14,*) &
                        (j-1)*dist_max/nstep_default,&
                        DON_missing(i,is)%atom(ia,j),&
                        bulk_DON(i)%spec(map(i)%spec(is,ia,1))%atom(map(i)%spec(is,ia,2),j)
                end do
                close(14)
                call err_abort_print_struc(splitbas(1),"lw_term.vasp",&
                     "",.false.)
                call err_abort_print_struc(splitbas(2),"up_term.vasp",&
                     "",.false.)
                call err_abort("ERROR: Internal error in get_shifts_DON\n&
                  &  More neighbours found in slab than in bulk.",.true.)
             end if
          end do atom_loop1
       end do spec_loop
       if(verbose_.ge.1)then
          write(*,*) "nneigh:",count1
          write(*,*)
       end if
       if(count1.le.0)then
          write(0,'("WARNING: No missing bonds identified for parent slab ",I0)') i
          deallocate(res_shifts)
          allocate(res_shifts(0,0))
          return
       end if
       allocate(intf(i)%neigh(count1))
       intf(i)%neigh(:)=neighbour(i,:count1)

    end do region_loop
    deallocate(neighbour)


!!!-----------------------------------------------------------------------------
!!! Zeroes interface atoms to the lowest atom on top slab
!!!-----------------------------------------------------------------------------
    highest_atom(1) = maxval(intf(1)%neigh(:)%pos(3),dim=1)
    lowest_atom(2) = minval(intf(2)%neigh(:)%pos(3),dim=1)
    intf(1)%neigh(:)%pos(3) = intf(1)%neigh(:)%pos(3) - highest_atom(1)
    intf(2)%neigh(:)%pos(3) = intf(2)%neigh(:)%pos(3) - lowest_atom(2)
    lowest_atom(1) = minval(intf(1)%neigh(:)%pos(3),dim=1)
    highest_atom(2) = maxval(intf(2)%neigh(:)%pos(3),dim=1)
    if(abs(verbose_).ge.1)then
       write(*,*) "lowest atom:",lowest_atom
       write(*,*) "highest atom:",highest_atom
    end if


!!!-----------------------------------------------------------------------------
!!! Defines grid size and equivalent step size
!!!-----------------------------------------------------------------------------
    lpresent=.false.
    if(present(offset))then
       if(offset(axis).ge.0._real32)then
          max_sep = max(abs(highest_atom(2)),abs(lowest_atom(1)))*modu(bas%lat(axis,:))
          lpresent=.true.
       end if
    end if
    if(.not.lpresent)then
       max_sep = max(abs(highest_atom(2)),abs(lowest_atom(1)))*modu(bas%lat(axis,:)) + 6._real32
       add = 0._real32
    end if

    stepsize=0.1
    ngrid(1)=nint(modu(bas%lat(1,:))/stepsize)
    ngrid(2)=nint(modu(bas%lat(2,:))/stepsize)
    ngrid(3)=ceiling(max_sep/stepsize)+1
    allocate(course_grid(2,ngrid(1),ngrid(2),ngrid(3)))
    allocate(tmp_neigh(max(size(intf(1)%neigh),size(intf(2)%neigh))*9))
    gridsize(1) = stepsize/modu(bas%lat(1,:))
    gridsize(2) = stepsize/modu(bas%lat(2,:))
    gridsize(3) = stepsize/modu(bas%lat(3,:))

    nstep(:2) = nint( min_trans(:2) * ngrid(:2) )
    nstep(3) = 0
    do jc=1,ngrid(3)
       pos(3) = real(jc-1,real32)*gridsize(3)
       if(pos(3)+highest_atom(2).gt.(ngrid(3)-1)*gridsize(3)) exit
       if(pos(3)-lowest_atom(1).gt.(ngrid(3)-1)*gridsize(3)) exit
       nstep(3) = nstep(3) + 1
    end do
    if(present(offset))then
       if(verbose_.ge.1) write(*,'(1X,"user-defined offset:",3(3X,F7.3))') offset
       add = -1.0
       do i=1,3
          if(offset(i).ge.0._real32)then
             nstep(i) = 1
             add(i) = offset(i)
          end if
       end do

       do i=1,3
          if(add(i).lt.0.0)then
             add(i) = 0.0
          end if
       end do
       add(axis) = add(axis)/modu(bas%lat(axis,:))
    end if

    !nthreads=8
    !call OMP_SET_NUM_THREADS(nthreads)
    !CHUNK = 2
!!!-----------------------------------------------------------------------------
!!! Determines neighbours for each grid point 
!!!-----------------------------------------------------------------------------
    if(abs(verbose_).ge.1)then
       write(*,'(1X,A,3(2X,F8.4))') &
            "lat:",modu(bas%lat(1,:)),modu(bas%lat(2,:)),modu(bas%lat(3,:))
       write(*,'(1X,A,3(2X,F8.4))') "gridsize:",gridsize
       write(*,*) "add:",add
       write(*,*) "nstep:",nstep
       write(*,*) "ngrid:",ngrid
       write(*,*) "max_sep:",max_sep
    end if

    if(any(nstep(:).le.0))then
       write(0,*) "ERROR: Internal error in get_shifts_DON"
       write(0,*) "nstep:",nstep
       write(0,*) "ngrid:",ngrid
       call err_abort_print_struc(splitbas(1),"lw_term.vasp",&
            "",.false.)
       call err_abort_print_struc(splitbas(2),"up_term.vasp",&
            "",.false.)
       call err_abort("ERROR: Internal error in get_shifts_DON",.true.)
    end if
!$OMP PARALLEL DO &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(is,ja,jb,jc,pos,vtmp1,vtmp2,vtmp3,count1,tmp_neigh) &
!$OMP SCHEDULE(DYNAMIC,2)
    do k=1,2
       nneigh = size(intf(k)%neigh,dim=1)

       do ja=1,ngrid(1)
          pos(1) = real(ja-1,real32)*gridsize(1) + add(1)
          do jb=1,ngrid(2)
             pos(2) = real(jb-1,real32)*gridsize(2) + add(2)
             do jc=1,ngrid(3)

                count1 = 0
                tmp_neigh = 0
                pos(3) = real(jc-1,real32)*gridsize(3) + add(3)
                do is=1,nneigh
                   vtmp1 = ( &
                        pos*(-1)**real(k-1,real32) - &
                        intf(k)%neigh(is)%pos(:3) )*(-1)**real(k-1,real32)
                   !vtmp1 = ( pos - intf(k)%neigh(is)%pos(:3) )!*(-1)**real(k,real32)
                   vtmp2(3) = vtmp1(3)
                   a_extend_loop: do i=-1,1,1
                      vtmp2(1) = vtmp1(1) + real(i,real32)
                      b_extend_loop: do j=-1,1,1
                         vtmp2(2) = vtmp1(2) + real(j,real32)
                         vtmp3 = matmul(vtmp2,bas%lat)
                         if(modu(vtmp3).gt.dist_max) cycle b_extend_loop
                         count1 = count1 + 1
                         tmp_neigh(count1) = modu(vtmp3)

                      end do b_extend_loop
                   end do a_extend_loop

                end do
                allocate(course_grid(k,ja,jb,jc)%neigh(count1))
                course_grid(k,ja,jb,jc)%neigh=tmp_neigh(:count1)

             end do
          end do
       end do
    end do
!$OMP END PARALLEL DO


!!!-----------------------------------------------------------------------------
!!! Shifts lower interface atoms within the grid and evaluates match
!!!-----------------------------------------------------------------------------
    allocate(fit_store(nstore))
    allocate(shift_store(nstore,3))
    fit_store=huge(0.0)
    shift_store=0._real32
! !$OMP PARALLEL DEFAULT(SHARED) NUM_NHREADS(nthreads) 
! !$OMP DO PRIVATE(ja,jb,jc,pos,val,l,is,nneigh,vtmp1,ivtmp1,ncheck,rtmp1,rtmp2,val) SCHEDULE(DYNAMIC,CHUNK)
    do ja=1,nstep(1)
       pos(1) = real(ja-1,real32)*gridsize(1)
       b_loop1: do jb=1,nstep(2)
          pos(2) = real(jb-1,real32)*gridsize(2)
          c_loop1: do jc=1,nstep(3)
             pos(3) = real(jc-1,real32)*gridsize(3)

             val = 0._real32
             do k=1,2
                l=minval([1,2],mask=[1,2].ne.k)
                nneigh = size(intf(l)%neigh,dim=1)
                do is=1,nneigh
                   vtmp1 = ( &
                        pos*(-1)**real(l,real32) + &
                        intf(l)%neigh(is)%pos  )*(-1)**real(l,real32)
                   vtmp1(:2) = vtmp1(:2) - floor( vtmp1(:2) )
                   ivtmp1 = nint(vtmp1/gridsize)
                   ivtmp1 = ivtmp1 + 1
                   !if(any(ivtmp1.gt.ngrid)) write(0,*) l,is,ivtmp1
                   where(ivtmp1(:2).gt.ngrid(:2))
                      ivtmp1(:2) = ivtmp1(:2) - ngrid(:2)
                   end where
                   ncheck = size(course_grid(k,ivtmp1(1),ivtmp1(2),ivtmp1(3))%neigh(:),dim=1)
                   !!-----------------------------------------------------------
                   !! Checks for bonds that match the missing bond set
                   !!-----------------------------------------------------------
                   rtmp2 = 0.0
                   grid_bond_loop: do i=1,ncheck
                      rtmp1 = abs( intf(l)%neigh(is)%bond - &
                           course_grid(k,ivtmp1(1),ivtmp1(2),ivtmp1(3))%neigh(i) )
                      if(rtmp1.gt.0.5) cycle grid_bond_loop
                      rtmp1 = ( 1.0 - tanh( 9.0*(rtmp1-0.25) ) )
                      rtmp2 = rtmp2 + rtmp1
                   end do grid_bond_loop
                   !val = val + abs(intf(l)%neigh(is)%num - rtmp2)!*(intf(l)%neigh(is)%bond)**2.5
                   !!-----------------------------------------------------------
                   !! Checks for atoms that are too close to the surface atoms
                   !!-----------------------------------------------------------
                   rtmp3 = 0.0
                   under_bond_loop: do i=1,ncheck
                      rtmp1 = course_grid(k,ivtmp1(1),ivtmp1(2),ivtmp1(3))%neigh(i) - &
                           intf(l)%neigh(is)%bond
                      if(rtmp1.ge.0.0) cycle under_bond_loop
                      !rtmp1 = rtmp1/course_grid(k,ivtmp1(1),ivtmp1(2),ivtmp1(3))%neigh(i)
                      rtmp1 = rtmp1/intf(l)%neigh(is)%bond
                      rtmp1 = abs(tan(pi*rtmp1/2.0))
                      rtmp3 = rtmp3 + rtmp1
                   end do under_bond_loop
                   
                   !val = val + abs(&
                   !     intf(l)%neigh(is)%num - rtmp2*f_scale - rtmp3*g_scale )
                   val = val + abs(intf(l)%neigh(is)%num - rtmp2*f_scale) + rtmp3*g_scale


                end do
                

             end do

             if(val.lt.fit_store(nstore))then
                fit_store(nstore) = val
                shift_store(nstore,:) = [ ja-1, jb-1, jc-1 ]
                call sort_shifts(fit_store,shift_store)
             end if

          end do c_loop1
       end do b_loop1
    end do
! !$OMP END DO
! !$OMP END PARALLEL


!!!-----------------------------------------------------------------------------
!!! Checks whether any shifts have been identified
!!!-----------------------------------------------------------------------------
    if(all(shift_store.eq.0))then
       call err_abort("ERROR: Internal error in get_shifts_DON\n&
               &  No shifts found.",.true.)
    end if


!!!-----------------------------------------------------------------------------
!!! Sets output of shifts
!!!-----------------------------------------------------------------------------
    if(verbose_.gt.0)then
       write(*,'("Determined shifts (gridsize:",3(2X,F6.4),")")') gridsize
       write(*,'(" num   fit_val   x    y    z")')
    end if
    do i = 1, nstore, 1
       res_shifts(i,:) = real(shift_store(i,:),real32)/real(ngrid(:)-1,real32)
       res_shifts(i,:2) = res_shifts(i,:2) + add(:2)
       if(verbose_.gt.0) &
            write(*,'(1X,I3,":",2X,F6.2,3(2X,I3))') i,fit_store(i),shift_store(i,:)
    end do
    res_shifts(:,axis) = (res_shifts(:,axis)*max_sep)/modu(bas%lat(axis,:)) + &
         add(axis)
    if(present(c_scale)) res_shifts(:,axis) = res_shifts(:,axis) * c_scale


    if(verbose_.gt.0)then
       write(*,'(1X,"Shifts to be applied (Å)")')
       do i = 1, nstore, 1
          write(*,'(I3,":",2X,3(2X,F7.4))') &
               i,res_shifts(i,:2),res_shifts(i,3)*modu(bas%lat(axis,:))
       end do
    end if
       

  end function get_shifts_DON
!!!#############################################################################


!!!#############################################################################
!!! sorts shifts by fit values
!!!#############################################################################
  subroutine sort_shifts(fits,shifts)
    implicit none
    integer :: i,loc,num
    real(real32) :: dbuff
    integer, dimension(3) :: ivtmp1
    integer, dimension(:,:), intent(inout) :: shifts
    real(real32), dimension(:), intent(inout) :: fits


    num = size(fits,dim=1)
    do i=1,num
       loc = minloc(fits(i:),dim=1) + i - 1
       dbuff = fits(i)
       fits(i) = fits(loc)
       fits(loc) = dbuff

       ivtmp1(:3) = shifts(i,:3)
       shifts(i,:3) = shifts(loc,:3)
       shifts(loc,:3) = ivtmp1(:3)
    end do


  end subroutine sort_shifts
!!!#############################################################################

end module shifting

