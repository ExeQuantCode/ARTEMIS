!!!#############################################################################
!!! Module compare two lattices to find potential matches
!!! Code written by Ned Thaddeus Taylor and Isiah Edward Mikel Rudkin
!!! Code part of the ARTEMIS group
!!!#############################################################################
!!! module contains various miscellaneous functions and subroutines.
!!! module includes the following functionsand subroutines:
!!! get_best_match
!!! pick_axis         (shifts lattices according to the static planes)
!!! cyc_lat1          (cycles through lat1 combinations)
!!! cyc_lat2          (finds best fit of lat2 to the designated lat1)
!!! tol_check
!!! lat_check
!!! convert_n_tf1!!! endcode
!!!#############################################################################
module lat_compare
  use artemis__constants, only: real32, pi, INF
  use artemis__misc_types, only: latmatch_type, tol_type
  use misc_linalg, only: cross,uvec,modu,get_area,find_tf,det,reduce_vec_gcd,&
       inverse_3x3,get_vec_multiple,get_frac_denom
  use artemis__geom_rw, only: basis_type
  use artemis__geom_utils, only: MATNORM,planecutter
  implicit none
  integer :: ierr_compare = 0
  logical :: lstop=.true.
  logical :: reduce=.false.


  
!!!updated  2021/11/19
 
  
contains

!!!#############################################################################
!!! cycles lattice 1
!!!#############################################################################
!!! Conditions used to stop transformation matrix for lattice 1 from ...
!!! ... cycling over ones already used:
!!! For n array:
!!!   - diagonals and lower off-diagonals are, at most, half the n_num
!!!   - upper off-diagonals are, at most, n_num
!!!   - after n(1,:) reaches (/ (n_num-1)/2, (n_num-1), (n_num-1) /), cycle ...
!!!     ... onto n_num = n_num + 1. This stops it from going over all ...
!!!     ... previously checked n arrays.
!!!   - if all n(1,:) or n(2,:) = 0, skip this configuration, 
!!! Means for transformation matrix:
!!!   - stops diagonal and lower off-diagonal elements from being negative. ...
!!!     ... These are already covered by the negative of that vector
!!!   - allows for negative values on the upper off-diagonal elements
!!!   - stops transformation matrix from checking over previous superlattices
!!!   - equivalent transformation matrix will have a determinant of zero
  subroutine cyc_lat1(SAV, tol, match_method, verbose)
    implicit none
    type(latmatch_type), intent(inout) :: SAV
    integer, intent(in) :: match_method
    integer, intent(in) :: verbose
    integer :: i,j,k
    integer :: n_num,count1
    logical :: l1change
    type(tol_type) :: tol
    integer, dimension(3,3) :: tf1,tf2
    integer, dimension(2,3) :: n
    real(real32), dimension(3,3) :: tlat1,tlat2
    real(real32), allocatable, dimension(:,:,:) :: match_tfs

    
!!!-----------------------------------------------------------------------------
!!! Initialised varaibles and allocates arrays
!!!-----------------------------------------------------------------------------
    allocate(match_tfs(tol%maxfit,3,3))
    match_tfs=0._real32
    SAV%nfit=0
    count1=0
    n=0
    n_num=0
    l1change=.false.


!!!-----------------------------------------------------------------------------
!!! Sets up the n array for the current value of n_num
!!!-----------------------------------------------------------------------------
101 n_num=n_num+1
    if(n_num.gt.tol%maxsize) return
    n(:,:)=0
    do i=1,2
       do j=1,SAV%axes(1)
          n(i,j)=n_num*(2**abs(floor((i-j)/2.0)))!CHECK WHY i-j AND NOT abs(i-j)
       end do
    end do
    

!!!-----------------------------------------------------------------------------
!!! Loops over the n array to check whether values are allowed
!!!-----------------------------------------------------------------------------
102 nloop: do
       chngloop: do i=2,1,-1
          do j=1,SAV%axes(1)
             if(n(i,j).lt.0)then
                if(j.ne.SAV%axes(1))then
                   n(i,j+1)=n(i,j+1)-1
                elseif(i.eq.2)then
                   n(1,1)=n(1,1)-1
                   l1change=.true.
                   cycle chngloop
                else
                   goto 101
                end if
                n(i,j)=n_num*(2**abs(floor((i-j)/2.0)))
             end if
          end do
       end do chngloop
       !!-----------------------------------------------------------------------
       !! Checks whether to cycle loop.
       !! Will do so if current n values are valid ...
       !! ... or if sequence has already been used
       !!-----------------------------------------------------------------------
       if(all(n(2,:SAV%axes(1)).eq.0))then
          n(2,1)=n(2,1)-1
          cycle nloop
       elseif(all(n(1,:SAV%axes(1)).eq.0).or.&
            all(n(1,:SAV%axes(1)).le.(n_num-1)*&
            (2**abs(floor((1- (/(i,i=1,SAV%axes(1),1)/) )/2.0)))))then
          goto 101
       end if
       if(l1change)then
          l1change=.false.
          if(any(mod(n(1,:),2).eq.0.and.n(1,:).ne.0))then
             do k=1,SAV%axes(1)
                n(2,k)=n_num*(2**abs(floor((2-k)/2.0)))
             end do
          else
             n(2,:)=n(1,:)
          end if
       end if
       

!!!-----------------------------------------------------------------------------
!!! Creates transformation matrix using n array
!!!-----------------------------------------------------------------------------
       tf1=convert_n_tf1(n,SAV%axes(1))
!       if(abs(nint(get_area(real(tf1(1,:),real32),real(tf1(2,:),real32)))).gt.tol%area)then
!          n(1,1)=n(1,1)-1
!          l1change=.true.
!          cycle
       if(abs(get_area(real(tf1(1,:),real32),real(tf1(2,:),real32))).lt.0.99_real32) goto 103


!       tf1(1,:)=(/1,0,0/)
!       tf1(2,:)=(/0,1,0/)
!       tf1(3,:)=(/0,0,1/)
!!!-----------------------------------------------------------------------------
!!! Creates superlattice from transformation matrix and original lattice
!!!-----------------------------------------------------------------------------
       tlat1=matmul(tf1,SAV%lat1)
       tlat1=MATNORM(tlat1)


!!!-----------------------------------------------------------------------------
!!! Compares superlattice to previously saved superlattices and cycles if same
!!!-----------------------------------------------------------------------------
       if(lat_check(SAV,tol,tlat1)) goto 103


!!!-----------------------------------------------------------------------------
!!! Generates corresponding superlattice of 2nd lattice that will best ...
!!! ... fit with current superlattice of 1st lattice.
!!!-----------------------------------------------------------------------------
       tlat2=cyc_lat2(SAV,tol,tlat1,tf1,tf2,match_tfs, match_method)
       count1=count1+1


!!!-----------------------------------------------------------------------------
!!! Compares superlattice to previously saved superlattices and cycles if same
!!!-----------------------------------------------------------------------------
       if(any(isnan(match_tfs(SAV%nfit+1,:,:)))) goto 103
       if(lat_check(SAV,tol,tlat1)) goto 103


!!!-----------------------------------------------------------------------------
!!! Checks whether similar transformation has already been saved
!!!-----------------------------------------------------------------------------
       do i=1,SAV%nfit
          if(all(abs(&
               match_tfs(SAV%nfit+1,:2,:3)-&
               match_tfs(i,:2,:3)).lt.1.E-5_real32)) goto 103
       end do


!!!-----------------------------------------------------------------------------
!!! Checks whether corresponding superlattices are within tolerance factors
!!!-----------------------------------------------------------------------------
       if(tol_check(SAV,tol,tlat1,tlat2,tf1,tf2,verbose))then
          !!--------------------------------------------------------------------
          !! Handles counters accordingly
          !!--------------------------------------------------------------------
          count1=0
          SAV%nfit=SAV%nfit+1
          if(ierr_compare.ge.1)then
             write(0,'(3(F0.2,1X))') (match_tfs(SAV%nfit,i,:),i=1,3)
             if(ierr_compare.eq.2) call exit()
          end if
       end if


!!!-----------------------------------------------------------------------------
!!! Checks whether any stop conditions are met
!!!-----------------------------------------------------------------------------
       if(SAV%nfit.eq.tol%maxfit) then
          if(verbose.gt.0) &
               write(*,'(/,"Number of fits reached maxfits ",I0)') SAV%nfit
          return
       end if
       if(lstop.and.count1.gt.100) then
          if(verbose.gt.0) &
               write(*,'(/,"Stopped as we reached ",I0," failed checks.")')&
               count1
          return
       end if


!!!-----------------------------------------------------------------------------
!!! Cycles n array to new values and cycles loop
!!!-----------------------------------------------------------------------------
103    n(2,1)=n(2,1)-1


    end do nloop


    return
  end subroutine cyc_lat1
!!!#############################################################################


!!!#############################################################################
!!! cycles lattice 2
!!!#############################################################################
  function cyc_lat2(SAV,tol,tlat1,tf1,tf2,match_tfs, match_method) result(tlat2)
    implicit none
    integer :: i,j
    type(tol_type) :: tol
    type(latmatch_type) :: SAV
    integer, dimension(3,3) :: tf1,tf2
    integer, dimension(3,3) :: it1_mat,it2_mat
    real(real32), dimension(3,3) :: t_mat,tlat1,tlat2
    real(real32), dimension(:,:,:) :: match_tfs
    integer, intent(in) :: match_method


    select case(match_method)
    case(1)
       tf2=get_lat2(SAV,tlat1)
    case(2)
       tf2=get_lat2_alt(SAV,tol,tlat1)
    end select
    if(all(tf2.eq.0)) goto 201


!!!-----------------------------------------------------------------------------
!!! Finds the transformation matrix between the two transformation matrices ...
!!! ... for lat1 and lat2 to their respective supercells.
!!! This can be used to make the simplest conversion from an identity ...
!!! ... transformation of lat1 and the corresponding transformation of lat2.
!!!-----------------------------------------------------------------------------
    SAV%reduced=.false.
    match_tfs(SAV%nfit+1,:,:)=find_tf((real(tf1,real32)),(real(tf2,real32)))
    if(any(isnan(match_tfs(SAV%nfit+1,:,:)))) goto 201
    t_mat(:,:)=match_tfs(SAV%nfit+1,:,:)
    if(ierr_compare.ge.1) then
       write(0,'("###############################")')
       write(0,'(3(I0,1X))') ((tf1(i,:)),i=1,3)
       write(0,*)
       write(0,'(3(I0,1X))') ((tf2(i,:)),i=1,3)
       write(0,*)
       write(0,'(3(F9.3,1X))') (match_tfs(SAV%nfit+1,i,:),i=1,3)
    end if


    reduce_if: if(SAV%reduce)then
       !t_mat=transpose(t_mat) !NOT SURE WHY TRANSPOSED
       it1_mat(:,:)=0
       it2_mat(:,:)=0
       it1_mat(3,3)=1
       it2_mat(3,:)=nint(t_mat(3,:))
       do i=1,2
          t_mat(i,:)=reduce_vec_gcd(t_mat(i,:))
          if(any(abs(t_mat(i,:)-nint(t_mat(i,:))).gt.1.E-5_real32)) exit reduce_if
          it2_mat(i,:)=nint(t_mat(i,:))
          do j=1,3
             if(match_tfs(SAV%nfit+1,j,i).ne.0._real32)then
                it1_mat(i,i)=nint(t_mat(i,j)/match_tfs(SAV%nfit+1,j,i))
                exit
             end if
          end do
          if(all(it2_mat(i,:).eq.0)) exit reduce_if
          if(it1_mat(i,i).eq.0) exit reduce_if
       end do
       if(abs(get_area(real(tf1(1,:),real32),real(tf1(2,:),real32))).lt.&
            abs(get_area(real(it1_mat(1,:),real32),real(it1_mat(2,:),real32))))&
            exit reduce_if
       SAV%reduced=.true.
       tf1=it1_mat
       tf2=it2_mat
       !tf2=matmul(tf2,it2_mat) !WHY WAS THIS USED?
       tlat1=matmul(tf1,SAV%lat1)
    end if reduce_if


!!!-----------------------------------------------------------------------------
!!! Writes the new superlattice of lat2
!!!-----------------------------------------------------------------------------
201 tlat2=matmul(tf2,SAV%lat2)


    return
  end function cyc_lat2
!!!#############################################################################


!!!#############################################################################
!!! finds lat 2 transformation matrix
!!!#############################################################################
  function get_lat2(SAV,tlat1) result(tf)
    implicit none
    integer :: i,kmax
    real(real32) :: dtmp,t_area,ang1,ang2,t_ang
    type(latmatch_type) :: SAV
    integer, dimension(3,3) :: tf,it_mat
    real(real32), dimension(3,3) :: t_mat,t_lat,tlat1,tlat2


    tf=0
!!!-----------------------------------------------------------------------------
!!! Finds the exact transformation matrix between superlattice of lat1 and ...
!!! ... basic lattice of lat2.
!!! Converts each matrix element to nint (to make it still maintain the same ...
!!! ... basis for lat2).
!!! Tests this for the three different faces of lat2 and saves smallest ...
!!! ... difference between superlattices for lat1 and lat2.
!!! NEED TO FIX THE LAST STATEMENT TO ONLY APPLY IT TO THE FACES UP TO MAXAXIS
!!!-----------------------------------------------------------------------------
    ang1=acos(dot_product(tlat1(1,:),tlat1(2,:))/&
         (modu(tlat1(1,:)*modu(tlat1(2,:)))))
    t_area=1000._real32
    t_ang=5._real32
    kmax=1
    if(SAV%axes(2).eq.3) kmax=3
    t_lat=SAV%lat2
    if(ierr_compare.ge.3) then
       write(0,*) tlat1
       write(0,*)
       write(0,*) t_lat
       write(0,*)
       write(0,*) find_tf(t_lat,tlat1)/det(find_tf(t_lat,tlat1))
       write(0,*)
       write(0,*) matmul(nint(find_tf(t_lat,tlat1)),t_lat)
    end if
    abc_loop: do i=1,kmax
       !! MAKE SURE THIS WORKS WITHIN THE CONFINES OF MAXAXIS
       !! REPLACE WITH SWAPPING VEC1 WITH VEC2 AND THEN VEC3

       t_lat(:,:)=cshift(SAV%lat2(:,:),shift=i-1,dim=1)
       t_lat(:,:)=cshift(t_lat(:,:),shift=i-1,dim=2)
       !t_lat(:,:)=cshift(SAV%lat2(:,:),shift=k-1,dim=2)
       !t_mat=find_tf(tlat1,t_lat)
       t_mat=find_tf(t_lat,tlat1)
       it_mat=nint(t_mat)

       tlat2=matmul(real(it_mat,real32),t_lat)
       ang2=acos(dot_product(tlat2(1,:),tlat2(2,:))/&
            (modu(tlat2(1,:))*modu(tlat2(2,:))))
       t_mat=tlat1-tlat2
       dtmp=get_area(t_mat(1,:),t_mat(2,:))
       t_area=1000._real32
       !! SORT OUT HANDLING OF AREA COMPARISON
       if(dtmp.le.t_area.and.&!-1.E-8_real32.and.&
            abs(ang1-ang2).lt.t_ang)then
          if(i.ne.1) it_mat(:,:)=cshift(it_mat(:,:),shift=1-i,dim=2)
          tf=it_mat
          t_area=dtmp
          t_ang=abs(ang1-ang2)
          if(ierr_compare.ge.1) write(0,*) i,ang1,ang2
       end if
    end do abc_loop


  end function get_lat2
!!!#############################################################################


!!!#############################################################################
!!! alternative method: finds lat 2 transformation matrix
!!!#############################################################################
  function get_lat2_alt(SAV,tol,tlat1) result(tf)
    implicit none
    integer :: i,j,k,m_max,m_num
    type(tol_type) :: tol
    type(latmatch_type) :: SAV
    integer, dimension(2,3) :: m
    integer, dimension(3,3) :: tf
    real(real32), dimension(3,3) :: tlat1
    real(real32), dimension(3,3) :: mA,mB,S,newlat
    logical :: lchange

    
!!! GET THE VOLUME OF tlat1, COMPARE TO HOW MUCH LARGER IT IS THAN SAV%lat2
!!! THAT WILL DEFINE THE MAX VALUE OF m (m_max)

!!! IF tf RETURNED AS ALL 0, THEN NO MATCH FOUND

    lchange = .false.
    m_num=0
    m_max=ceiling(&
         get_area(tlat1(1,:),tlat1(2,:))/get_area(SAV%lat2(1,:),SAV%lat2(2,:)))
    do i=1,3
       do j=1,3
          mA(i,j)=dot_product(tlat1(i,:),tlat1(j,:))
          mB(i,j)=dot_product(SAV%lat2(i,:),SAV%lat2(j,:))
       end do
    end do
    S=1._real32
    S(:,:)=sqrt(mA(1,1))*sqrt(mA(2,2))*cos(pi/2._real32-tol%ang)
    do i=1,3
       S(i,i)=(2._real32*tol%vec)*mA(i,i)
    end do


!!!-----------------------------------------------------------------------------
!!! Sets up the n array for the current value of n_num
!!!-----------------------------------------------------------------------------
301 m_num=m_num+1
    if(m_num.gt.m_max)then
       tf=0
       return
    end if
    m(:,:)=0
    do i=1,2
       do j=1,SAV%axes(2)
          m(i,j)=m_num*(2**abs(floor((i-j)/2.0)))!CHECK WHY i-j AND NOT abs(i-j)
       end do
    end do
    

!!!-----------------------------------------------------------------------------
!!! Loops over the n array to check whether values are allowed
!!!-----------------------------------------------------------------------------
    mloop: do
       chngloop2: do i=2,1,-1
          do j=1,SAV%axes(2)
             if(m(i,j).lt.0)then
                if(j.ne.SAV%axes(2))then
                   m(i,j+1)=m(i,j+1)-1
                elseif(i.eq.2)then
                   m(1,1)=m(1,1)-1
                   lchange=.true.
                   cycle chngloop2
                else
                   goto 301
                end if
                m(i,j)=m_num*(2**abs(floor((i-j)/2.0)))
             end if
          end do
       end do chngloop2
       !!-----------------------------------------------------------------------
       !! Checks whether to cycle loop.
       !! Will do so if current n values are valid ...
       !! ... or if sequence has already been used
       !!-----------------------------------------------------------------------
       if(all(m(2,:SAV%axes(2)).eq.0))then
          m(2,1)=m(2,1)-1
          cycle mloop
!       elseif(all(m(:,1).eq.0).or.all(m(:,2).eq.0))then
!          m(2,1)=m(2,1)-1
!          cycle mloop
       elseif(all(m(1,:SAV%axes(2)).eq.0))then
          !.or.all(m(1,:SAV%axes(2)).le.(m_num-1)*&
          !(2**abs(floor((1- (/(i,i=1,SAV%axes(2),1)/) )/2.0)))))then
          goto 301
       end if
       if(lchange)then
          lchange=.false.
          if(any(mod(m(1,:),2).eq.0.and.m(1,:).ne.0))then
             do k=1,SAV%axes(2)
                m(2,k)=m_num*(2**abs(floor((2-k)/2.0)))
             end do
          else
             m(2,:)=m(1,:)
          end if
       end if


       tf=convert_n_tf1(m,SAV%axes(2))
       newlat=matmul((tf),(mB))
       newlat=matmul(newlat,transpose(tf))


!!!using 1 Å as the tolerance
!!! probably want smaller off, diagonal differences
       if(all((abs(newlat(:2,:2)-mA(:2,:2))-S(:2,:2)).lt.0._real32))then
          if(ierr_compare.gt.1)then
             write(0,*) "success"
             write(0,'(3(I0,1X))') tf
             write(0,*)
             write(0,'(2(F0.2,1X))') mA(:2,:2)
             write(0,*)
             write(0,'(2(F0.2,1X))') newlat(:2,:2)-mA(:2,:2)
          end if
          exit
       end if


       m(2,1)=m(2,1)-1
    end do mloop

    
  end function get_lat2_alt
!!!#############################################################################


!!!#############################################################################
!!! Checks whether the supplied superlattices fit within the tolerances
!!!#############################################################################
  function tol_check(SAV,tol,tlat1,tlat2,tf1,tf2,verbose) result(lmatch)
    implicit none
    type(latmatch_type), intent(inout) :: SAV
    real(real32), dimension(3,3), intent(in) :: tlat1, tlat2
    integer, dimension(3,3), intent(inout) :: tf1, tf2
    integer, intent(in) :: verbose

    integer :: i,j
    real(real32) :: ang1,ang2,t_area1,t_area2,diff
    logical :: la1a2,la1b2,l12,lmatch
    type(tol_type) :: tol
    real(real32), dimension(2) :: mag_mat1,mag_mat2
    real(real32), dimension(3) :: tvec


    lmatch=.false.
!!!-----------------------------------------------------------------------------
!!! Generates the corresponding areas and vector lengths of both lattices
!!!-----------------------------------------------------------------------------
    t_area1=get_area(tlat1(1,:),tlat1(2,:))
    t_area2=get_area(tlat2(1,:),tlat2(2,:))
    mag_mat1(1)=modu(tlat1(1,:))
    mag_mat1(2)=modu(tlat1(2,:))
    mag_mat2(1)=modu(tlat2(1,:))
    mag_mat2(2)=modu(tlat2(2,:))


    if(ierr_compare.gt.1) write(0,*) "area:",t_area1,t_area2
!!!-----------------------------------------------------------------------------
!!! Compares lattices using tolerances to check for a potential match
!!!-----------------------------------------------------------------------------
    if(abs((t_area1-t_area2)/t_area1).gt.tol%area) then
       return
    elseif(abs((t_area1-t_area2)/t_area1).le.tol%area) then
       ang1=acos(dot_product(tlat1(1,:),tlat1(2,:))/(mag_mat1(1)*mag_mat1(2)))
       ang2=acos(dot_product(tlat2(1,:),tlat2(2,:))/(mag_mat2(1)*mag_mat2(2)))
       !!-----------------------------------------------------------------------
       !! Changed angles to all less than pi/2 to deal with negative vectors
       !!-----------------------------------------------------------------------
       if(ang1.gt.pi/2._real32) ang1=pi-ang1
       if(ang2.gt.pi/2._real32) ang2=pi-ang2
       if(ierr_compare.gt.1) write(0,*) ang1,ang2
       !!-----------------------------------------------------------------------
       la1a2=(abs((mag_mat1(1)-mag_mat2(1))/mag_mat1(1)).lt.tol%vec.and.&
            abs((mag_mat1(2)-mag_mat2(2))/mag_mat1(2)).lt.tol%vec)
       la1b2=(abs((mag_mat1(1)-mag_mat2(2))/mag_mat1(1)).lt.tol%vec.and.&
            abs((mag_mat1(2)-mag_mat2(1))/mag_mat1(2)).lt.tol%vec)
       !!needs .or. below to cover either a or b matching lat2 a or b
       l12=(la1a2.or.la1b2)
       if(abs(ang1-ang2)*180/pi.lt.tol%ang.and.l12)then
          lmatch=.true.
       else
          return
       end if
       if(la1a2)then
          diff=max(&
               abs((mag_mat1(1)-mag_mat2(1))/mag_mat1(1)),&
               abs((mag_mat1(2)-mag_mat2(2))/mag_mat1(2)))
       else
          tvec=tf2(1,:)
          tf2(1,:)=tf2(2,:)
          tf2(2,:)=-nint(tvec)
          diff=max(&
               abs((mag_mat1(1)-mag_mat2(2))/mag_mat1(1)),&
               abs((mag_mat1(2)-mag_mat2(1))/mag_mat1(2)))
       end if
       !!-----------------------------------------------------------------------
       !! Generating unit vector c axis for both superlattices ...
       !! ... perpendicular to the interface plane.
       !!-----------------------------------------------------------------------
       tf1(3,:)=nint(uvec(cross(real(tf1(1,:),real32),real(tf1(2,:),real32))))
       tf2(3,:)=nint(uvec(cross(real(tf2(1,:),real32),real(tf2(2,:),real32))))
       !!-----------------------------------------------------------------------
       !! Prints the mismatches for the current successful match
       !!-----------------------------------------------------------------------
       if(verbose.gt.0)then
          write(*,'(/,A,I0,2X,A,I0)') &
               "Fit number: ",SAV%nfit+1,&
               "Area increase: ",&
               nint(get_area(real(tf1(1,:),real32),real(tf1(2,:),real32)))
          write(*,'("   Transmat 1:    Transmat 2:")')
          write(*,'((/,1X,3(3X,A1),3X,3(3X,A1)))') SAV%abc,SAV%abc
          write(*,'(3(/,2X,3(I3," "),3X,3(I3," ")))') &
               tf1(1,1:3),tf2(1,1:3),&
               tf1(2,1:3),tf2(2,1:3),&
               tf1(3,1:3),tf2(3,1:3)
          write(*,'(" vector mismatch (%) = ",F0.9)') diff*100._real32
          write(*,'(" angle mismatch (°)  = ",F0.9)') abs(ang1-ang2)*180/pi
          write(*,'(" area mismatch (%)   = ",F0.9)') (&
               1-abs(t_area1/t_area2))*100._real32
          write(*,*) "reduced:",SAV%reduced
       end if
       !!-----------------------------------------------------------------------
       !! Checks if best mismatch and saves accordingly
       !!-----------------------------------------------------------------------
       best_check: do i=1,SAV%max_num_matches
          if(i.gt.SAV%nfit)then
             SAV%tol(i,1)=diff*100._real32
             SAV%tol(i,2)=abs(ang1-ang2)
             SAV%tol(i,3)=(1-abs(t_area1/t_area2))*100._real32
             SAV%tf1(i,:,:)=tf1(:,:)
             SAV%tf2(i,:,:)=tf2(:,:)
             exit best_check
          end if
          if(diff*100._real32.le.SAV%tol(i,1).and.&
               abs(ang1-ang2).le.SAV%tol(i,2).and.&
               (1-abs(t_area1/t_area2))*100._real32.le.SAV%tol(i,3)) then
             if(nint(get_area(real(tf1(1,:),real32),real(tf1(2,:),real32))).ge.&
                  nint(get_area(real(SAV%tf1(i,1,:),real32),real(SAV%tf1(i,2,:),real32))))&
                  cycle best_check
             do j=SAV%max_num_matches,i+1,-1
                SAV%tol(j,:)=SAV%tol(j-1,:)
                SAV%tf1(j,:,:)=SAV%tf1(j-1,:,:)
                SAV%tf2(j,:,:)=SAV%tf2(j-1,:,:)
             end do
             SAV%tol(i,1)=diff*100._real32
             SAV%tol(i,2)=abs(ang1-ang2)
             SAV%tol(i,3)=(1-abs(t_area1/t_area2))*100._real32
             SAV%tf1(i,:,:)=tf1(:,:)
             SAV%tf2(i,:,:)=tf2(:,:)
             exit best_check
          end if
       end do best_check
    end if

    
    return
  end function tol_check
!!!#############################################################################


!!!#############################################################################
!!! checks whether, after applying transmat to lat, if it is the same as a ...
!!! ... previously successful one
!!!#############################################################################
  function lat_check(SAV,tol,lat) result(lcheck)
    implicit none
    integer :: i
    real(real32) :: ang1,ang2,tiny
    logical :: lcheck,lmatch_aa,lmatch_ab
    type(tol_type) :: tol
    type(latmatch_type) :: SAV
    real(real32), dimension(3,3) :: lat,tlat


    tiny=1.E-6_real32
    lcheck=.false.
    lat_loop: do i=1,min(SAV%max_num_matches,SAV%nfit)
       tlat=matmul(SAV%tf1(i,:,:),SAV%lat1)
       ang1=acos(dot_product(lat(1,:),lat(2,:))/(&
            sqrt(dot_product(lat(1,:),lat(1,:)))*&
            sqrt(dot_product(lat(2,:),lat(2,:)))))
       ang2=acos(dot_product(tlat(1,:),tlat(2,:))/(&
            sqrt(dot_product(tlat(1,:),tlat(1,:)))*&
            sqrt(dot_product(tlat(2,:),tlat(2,:)))))
       if(ang1.gt.pi/2._real32) ang1=pi-ang1
       if(ang2.gt.pi/2._real32) ang2=pi-ang2
       if(abs(ang1-ang2).lt.tiny)then
          lmatch_aa=&
               (abs(dot_product(lat(1,:),lat(1,:))-&
               dot_product(tlat(1,:),tlat(1,:))).lt.tiny).and.&
               (abs(dot_product(lat(2,:),lat(2,:))-&
               dot_product(tlat(2,:),tlat(2,:))).lt.tiny)

          lmatch_ab=&
               (abs(dot_product(lat(1,:),lat(1,:))-&
               dot_product(tlat(2,:),tlat(2,:))).lt.tiny).and.&
               (abs(dot_product(lat(2,:),lat(2,:))-&
               dot_product(tlat(1,:),tlat(1,:))).lt.tiny)
          if(lmatch_aa.or.lmatch_ab)then
             lcheck=.true.
             return
          end if
       end if
    end do lat_loop


  end function lat_check
!!!#############################################################################


!!!#############################################################################
!!! converts n array to tf1
!!!#############################################################################
  function convert_n_tf1(n,maxa) result(converted)
    implicit none
    integer :: i,j
    integer :: maxa
    integer, dimension(2,3) :: n
    integer, dimension(3,3) :: converted
  

    converted=0
    do i=1,2
       do j=1,maxa
          if(j.le.i)then
             converted(i,j)=n(i,j)
          else
             converted(i,j)=&
                  nint(floor(abs(1.0+(n(i,j)-1.0)/2.0))*((-1.0)**(n(i,j)-1.0)))
          end if
       end do
    end do
    converted(3,:)=(/0,0,1/)


    return    
  end function convert_n_tf1
!!!#############################################################################


!!!#############################################################################
!!! ends the code
!!!#############################################################################
  subroutine endcode(SAV)
    implicit none
    type(latmatch_type) :: SAV


    write(*,*)
    if(SAV%nfit.eq.0)then
       write(*,'(" No matches were found within the tolerances supplied.")')
       write(*,*)
       call exit(1)
    end if

    write(*,'(1X,"BEST MATCH      Area increase: ",I0)') &
         nint(get_area(real(SAV%tf1(1,1,:),real32),real(SAV%tf1(1,2,:),real32)))
    write(*,'("   Transmat 1:    Transmat 2:")')
    write(*,'((/,1X,3(3X,A1),3X,3(3X,A1)),3(/,2X,3(I3," "),3X,3(I3," ")))') &
         SAV%abc,SAV%abc,&
         SAV%tf1(1,1,1:3),SAV%tf2(1,1,1:3),&
         SAV%tf1(1,2,1:3),SAV%tf2(1,2,1:3),&
         SAV%tf1(1,3,1:3),SAV%tf2(1,3,1:3)
    write(*,'(" vector mismatch (%) = ",F0.9)') SAV%tol(1,1)
    write(*,'(" angle mismatch (°)  = ",F0.9)') SAV%tol(1,2)*180/pi
    write(*,'(" area mismatch (%)   = ",F0.9)') SAV%tol(1,3)
    write(*,*)

    write(*,'(A)') "EXITING"


    return
  end subroutine endcode
!!!#############################################################################


!!!#############################################################################
!!! Steve lattice match
!!!#############################################################################
  function vec_comp(S1,S1p,S2p,delta) result(match)
    implicit none
    real(real32) :: ct,cp,cv,th,ph,va
    real(real32) :: beta,pm1,alpha,pm2
    real(real32) :: mS1,mS1p,mS2p,tiny,md
    real(real32), dimension(2) :: match
    real(real32), dimension(3) :: S1,S1p,S2p,delta


    match=0._real32
    tiny=1.E-8_real32
    mS1=modu(S1)
    mS1p=modu(S1p)
    mS2p=modu(S2p)
    md=modu(delta)

    ct=dot_product(S1p,S2p)/(mS1p*mS2p)
    cp=dot_product(S1,S1p) /(mS1* mS1p)
    cv=dot_product(S1,S2p) /(mS1* mS2p)
    th=acos(dot_product(S1p,S2p)/(mS1p*mS2p))
    ph=acos(dot_product(S1,S1p) /(mS1* mS1p))
    va=acos(dot_product(S1,S2p) /(mS1* mS2p))

    beta=mS1*(cv-ct*cp)/(mS2p*sin(acos(ct))**2.0)
    pm1=(cv-ct*cp)**2.0 - (sin(th)*sin(ph))**2.0 !- md*(sin(th)/mS1)**2._real32
    if(abs(pm1).lt.tiny.or.pm1+(md*sin(th)/mS1)**2._real32.gt.0._real32) pm1=0._real32
    pm1=mS1*sqrt(pm1)/( mS2p*(1-ct**2._real32) )
    if(abs(beta+pm1-nint(beta+pm1)).lt.&
         abs(beta-pm1-nint(beta-pm1)))then
       match(1)=beta+pm1
    else
       match(1)=beta-pm1
    end if
    beta=match(1)

    !t_beta=beta+pm1*(-1.0)**real(i,real32)
    alpha=-( beta*mS2p*ct - mS1*cp )/mS1p
    pm2=-(beta*mS2p*sin(th))**2._real32 &
         -(mS1*sin(ph))**2._real32 + &
         2._real32*beta*mS1*mS2p*(cv - ct*cp) !- md
    if(abs(pm2).lt.tiny.or.pm2+md**2._real32.gt.0._real32) pm2=0._real32
    pm2=sqrt(pm2)/mS1p


    if(abs(alpha+pm2-nint(alpha+pm2)).lt.&
         abs(alpha-pm2-nint(alpha-pm2)))then
       match(2)=alpha+pm2
    else
       match(2)=alpha-pm2
    end if


  end function vec_comp
!!!#############################################################################


!!!#############################################################################
!!! Isiah lattice match
!!! Program to match lattices of two position cards.
!!!#############################################################################
  subroutine lattice_matching( &
       SAV, tol, structure_lw, structure_up, &
       miller_lw, miller_up, max_num_planes, &
       verbose, tol_sym &
  )
    use artemis__sym
    use plane_matching
    implicit none

    type(latmatch_type), intent(inout) :: SAV
    type(basis_type), intent(in) :: structure_lw,structure_up
    integer, dimension(3), intent(in) :: miller_lw,miller_up
    integer, intent(in) :: max_num_planes
    integer, intent(in) :: verbose
    real(real32), intent(in) :: tol_sym
    
    type(sym_type) :: grp1,grp2
    type(tol_type) :: tol
    type(tol_type) :: pm_tol
    real(real32), dimension(3,3) :: tf
    real(real32), dimension(3,3) :: lat1,lat2 !original lattices.
    real(real32), dimension(3,3) :: templat1,templat2 !tmp lattices to feed into plane matching.
    integer :: itmp1,nsym1,nsym2
    integer :: m1,m2,m3,i1,i2,i3,loc
    integer :: loopsize !size of the main loops
    integer :: i,j,num_of_transforms ! n = number of output transforms
    real(real32) :: dtmp1
    logical, allocatable, dimension(:) :: lvec1

    integer, dimension(3,3) :: tmat1,tmat2
    integer, dimension(3,3) :: transform1,transform2 !The transformations output by planecutter.

    real(real32), dimension(3) :: rvec1, rvec2
    real(real32), dimension(3,3) :: rmat1
    
    real(real32), allocatable, dimension(:,:,:) :: tmpsym1,tmpsym2,tmpsym
    real(real32), allocatable, dimension(:,:,:) :: transform1_saved,transform2_saved !The transformations output by plane cutter

    integer, allocatable, dimension(:,:,:) :: Tcellmatch_1,Tcellmatch_2 !The transformation matrices output from the cell_match program for lattices 1 and 2.
    real(real32), allocatable, dimension(:,:,:) :: Tsaved_1,Tsaved_2
    real(real32), allocatable, dimension(:,:,:) :: big_T_1,big_T_2 ! 3x3 versions of the matrices output by cell_match
    real(real32), dimension(3,3) :: dummy_mat1,dummy_mat2 ! temporary matrices used when the info is stored in a tensor.
    real(real32), dimension(2,2) :: temp_mat1,temp_mat2 ! temporary matrices used when the info is stored in a tensor.
    real(real32), allocatable, dimension(:,:,:) :: comb_trans_1,comb_trans_2 !The combined transformations (planecutter output)x(cellmatch output).

    real(real32), allocatable, dimension(:,:) :: tolerances,saved_tolerances
    integer, allocatable, dimension(:,:) :: ivtmp1,miller1,miller2
    integer, dimension(3) :: ivtmp2



    !!--------------------------------------------------------------------------
    !! sets initial variables
    !!--------------------------------------------------------------------------
    SAV%nfit = 0
    allocate(transform1_saved(SAV%max_num_matches,3,3))
    allocate(transform2_saved(SAV%max_num_matches,3,3))
    allocate(Tsaved_1(SAV%max_num_matches,2,2))
    allocate(Tsaved_2(SAV%max_num_matches,2,2))
    transform1_saved = 0._real32
    transform2_saved = 0._real32
    Tsaved_1 = 0._real32
    Tsaved_2 = 0._real32
    allocate(tolerances(SAV%max_num_matches,3))
    allocate(saved_tolerances(SAV%max_num_matches,3))
    saved_tolerances = INF
    lat1 = SAV%lat1
    lat2 = SAV%lat2
    pm_tol%maxsize=tol%maxsize
    pm_tol%maxfit=tol%maxfit
    pm_tol%vec=tol%vec
    pm_tol%ang=tol%ang
    pm_tol%area=tol%area
    pm_tol%ang_weight=tol%ang_weight
    pm_tol%area_weight=tol%area_weight
    pm_tol%maxlen=tol%maxlen
    pm_tol%maxarea=tol%maxarea


    !!--------------------------------------------------------------------------
    !! finds and stores symmetry operations for each lattice
    !!--------------------------------------------------------------------------
    call grp1%init(lat1, tol_sym=tol_sym,new_start=.true.)
    call check_sym(grp1,structure_lw,lsave=.true.,tol_sym=tol_sym)
    allocate(tmpsym1(3,3,grp1%nsym))
    

    call grp2%init(lat2, tol_sym=tol_sym,new_start=.true.)
    call check_sym(grp2,structure_up,lsave=.true.,tol_sym=tol_sym)
    allocate(tmpsym2(3,3,grp2%nsym))


    !!--------------------------------------------------------------------------
    !! initialises temporary Miller plane storage
    !!--------------------------------------------------------------------------
    loopsize=10
    allocate(ivtmp1((2*loopsize+1)**3,3))
    !allocate(ivtmp1(max_num_planes,3))


    !!--------------------------------------------------------------------------
    !! generate all unique planes for lattice 1
    !!--------------------------------------------------------------------------
    ivtmp1=0
    itmp1=0
    if(any(miller_lw.ne.0))then
       allocate(miller1(1,size(miller_lw)))
       miller1(1,:3)=miller_lw(:3)
    else
       mloop1: do i1=1,loopsize
          m1=floor((i1)/2.0)*(-1)**i1
          mloop2: do i2=1,loopsize
             m2=floor((i2)/2.0)*(-1)**i2
             mloop3: do i3=1,loopsize
                m3=floor((i3)/2.0)*(-1)**i3
                if ( .not.is_unique( [ m1, m2, m3 ], grp1%sym(:3,:3,:) ) ) &
                     cycle mloop3
                itmp1 = itmp1 + 1
                ivtmp1(itmp1,:) = [ m1, m2, m3 ]
                !if(itmp1.eq.max_num_planes) exit mloop1
             end do mloop3
          end do mloop2
       end do mloop1
       do i=1,itmp1
          loc = minloc(&
               abs(ivtmp1(i:itmp1,1)) + &
               abs(ivtmp1(i:itmp1,2)) + &
               abs(ivtmp1(i:itmp1,3)),dim=1) + i - 1
          ivtmp2(:) = ivtmp1(i,:)
          ivtmp1(i,:) = ivtmp1(loc,:)
          ivtmp1(loc,:) = ivtmp2(:)
       end do
       itmp1 = min(itmp1,max_num_planes)
       allocate(miller1(itmp1,3))
       miller1(:,:) = ivtmp1(:itmp1,:)
    end if


    !!--------------------------------------------------------------------------
    !! generate all unique planes for lattice 2
    !!--------------------------------------------------------------------------
    itmp1 = 0
    ivtmp1 = 0
    if(any(miller_up.ne.0))then
       allocate(miller2(1,size(miller_up)))
       miller2(1,:3)=miller_up(:3)
    else
       mloop4: do i1=1,loopsize
          m1=floor((i1)/2.0)*(-1)**i1
          mloop5: do i2=1,loopsize
             m2=floor((i2)/2.0)*(-1)**i2
             mloop6: do i3=1,loopsize
                m3=floor((i3)/2.0)*(-1)**i3
                if ( .not.is_unique( (/m1,m2,m3/), grp2%sym(:3,:3,:) ) ) &
                     cycle mloop6
                itmp1=itmp1+1
                ivtmp1(itmp1,:)=(/m1,m2,m3/)
                !if(itmp1.eq.max_num_planes) exit mloop4
             end do mloop6
          end do mloop5
       end do mloop4
       do i=1,itmp1
          loc = minloc(&
               abs(ivtmp1(i:itmp1,1)) + &
               abs(ivtmp1(i:itmp1,2)) + &
               abs(ivtmp1(i:itmp1,3)),dim=1) + i - 1
          ivtmp2(:) = ivtmp1(i,:)
          ivtmp1(i,:) = ivtmp1(loc,:)
          ivtmp1(loc,:) = ivtmp2(:)
       end do
       itmp1 = min(itmp1,max_num_planes)
       allocate(miller2(itmp1,3))
       miller2(:,:) = ivtmp1(:itmp1,:)
    end if
    if(verbose.gt.0)then
       write(*,*)
       write(*,'(1X,"Miller planes considered for lower material: ",I0)') &
            size(miller1(:,1))
       do i=1,size(miller1(:,1))
          write(*,'(2X,I2,")",3X,3(3X,I0))') i,miller1(i,:)
       end do
       write(*,*)
       write(*,'(1X,"Miller planes considered for upper material: ",I0)') &
            size(miller2(:,1))
       do i=1,size(miller2(:,1))
          write(*,'(2X,I2,")",3X,3(3X,I0))') i,miller2(i,:)
       end do
       write(*,*)
    end if


    !!--------------------------------------------------------------------------
    !! cycles through the unique miller planes to find matches
    !!--------------------------------------------------------------------------
    allocate(tmpsym(3,3,max(grp1%nsym,grp2%nsym)))
    MAINLOOP1: do m1 = 1, size( miller1, dim = 1 )
       transform1 = nint(planecutter(lat1,real(miller1(m1,:),real32)))
       if (all(transform1 .eq. 0)) cycle MAINLOOP1
       templat1 = matmul(transform1,lat1)
       tmpsym = 0._real32
       do i=1,grp1%nsym
          tmpsym(:3,:3,i) = &
               matmul(grp1%sym(:3,:3,i),inverse_3x3(real(transform1,real32)))
          ! next step required to transform properly into the space?
          tmpsym(:3,:3,i) = &
               matmul(real(transform1,real32),tmpsym(:3,:3,i))
       end do

       nsym1=0
       tmpsym1=0._real32
!!! IS THIS REASONABLE TO DO IT THIS WAY? OR DO WE NEED TO CHANGE sym TO BE IN THE NEW LAT?
!!! Wait, should it be instead that the cross product of the a-b plane is always consistent?
       rvec1=cross([templat1(1,:)],[templat1(2,:)])
       do i = 1, grp1%nsym, 1
          rmat1=matmul(tmpsym(:3,:3,i),templat1(:,:))
          rvec2=cross([rmat1(1,:)],[rmat1(2,:)])
          if(all(abs( rvec1(:) - rvec2(:) ).lt.1.E-8_real32).or.&
               all(abs( rvec1(:) + rvec2(:) ).lt.1.E-8_real32))then
             nsym1 = nsym1 + 1
             tmpsym1(:3,:3,nsym1) = tmpsym(:3,:3,i)
          else
             cycle
          end if
          ! redundant if a-b plane works instead.
          !if(all(&
          !     abs( templat1(3,:) - matmul(templat1(3,:),tmpsym(:3,:3,i)) )&
          !     .lt.1.E-8_real32).or.&
          !     all(&
          !     abs( templat1(3,:) + matmul(templat1(3,:),tmpsym(:3,:3,i)) )&
          !     .lt.1.E-8_real32))then
          !   nsym1=nsym1+1
          !   tmpsym1(:3,:3,nsym1) = tmpsym(:3,:3,i)
          !end if
          !write(0,*) "################################"
          !write(0,*) i
          !write(0,'(3(2X,F7.2))') tmpsym(:3,:3,i)
          !write(0,*)
          !write(0,'(3(2X,F7.2))') rvec1!(templat1(j,:),j=1,3)!(grp1%sym(j,:3,i),j=1,3) !tmpsym(:3,:3,i)
          !write(0,*)
          !write(0,'(3(2X,F7.2))') rvec2!matmul(templat1(3,:),tmpsym(:3,:3,i))!(tmpsym(j,:3,i),j=1,3)
       end do
       !stop


       MAINLOOP2: do m2=1,size(miller2(:,1),dim=1)
          transform2 = nint(planecutter(lat2,real(miller2(m2,:),real32)))
          if (all(transform2 .eq. 0)) cycle MAINLOOP2
          templat2 = matmul(transform2,lat2)
          
          tmpsym=0._real32
          do i = 1, grp2%nsym, 1
             tmpsym(:3,:3,i) = &
                  matmul(grp2%sym(:3,:3,i),inverse_3x3(real(transform2,real32)))
             ! next step required to transform properly into the space?
             tmpsym(:3,:3,i) = &
                  matmul(real(transform2,real32),tmpsym(:3,:3,i))
          end do
          nsym2=0
          tmpsym2 = 0._real32
          do i = 1, grp2%nsym, 1
             !write(0,*) "################################"
             !write(0,*) i
             !write(0,'(3(2X,F7.2))') (grp2%sym(j,:3,i),j=1,3) !tmpsym(:3,:3,i)
             !write(0,*)
             !write(0,'(3(2X,F7.2))') (tmpsym(j,:3,i),j=1,3)
             if(all(&
                  abs( templat2(3,:) - matmul(templat2(3,:),tmpsym(:3,:3,i)) )&
                  .lt.1.E-8_real32).or.&
                  all(&
                  abs( templat2(3,:) + matmul(templat2(3,:),tmpsym(:3,:3,i)) )&
                  .lt.1.E-8_real32))then
                nsym2 = nsym2 + 1
                tmpsym2(:3,:3,nsym2) = tmpsym(:3,:3,i)
             end if
          end do


          !!--------------------------------------------------------------------
          !! Calls the function cell_match which matches the ab plane for ...
          !! ... two input lattices
          !!--------------------------------------------------------------------
          call cell_match(pm_tol,&
               lat1=templat1,lat2=templat2,&
               transforms1=Tcellmatch_1,&
               transforms2=Tcellmatch_2,&
               ntransforms=num_of_transforms,&
               matched_tols=tolerances,&
               sym1=tmpsym1(:,:,:nsym1),sym2=tmpsym2(:,:,:nsym2))


          !!--------------------------------------------------------------------
          !! Find the (SAV%max_num_matches) best matches overall
          !!--------------------------------------------------------------------
          loop110: do i=1,num_of_transforms
             IF101: if ( dot_product(tolerances(i,:),vaa_weighting).le.&
                  dot_product(saved_tolerances(SAV%max_num_matches,:),vaa_weighting) )then
                temp_mat1(:,:) = real(Tcellmatch_1(i,:,:),real32)
                temp_mat2(:,:) = real(Tcellmatch_2(i,:,:),real32)
                IF102: if (.not.is_duplicate(&
                     (Tsaved_1),(Tsaved_2),&
                     (temp_mat1),(temp_mat2),&
                     tmpsym1,tmpsym2) ) then
                   saved_tolerances(SAV%max_num_matches,:) = tolerances(i,:)
                   Tsaved_1(SAV%max_num_matches,:,:) = temp_mat1(:,:)
                   Tsaved_2(SAV%max_num_matches,:,:) = temp_mat2(:,:)
                   transform1_saved(SAV%max_num_matches,:,:) = real(transform1(:,:),real32)
                   transform2_saved(SAV%max_num_matches,:,:) = real(transform2(:,:),real32)


                   if(SAV%nfit.lt.SAV%max_num_matches) SAV%nfit = SAV%nfit + 1
                   call datasortmain_tols(saved_tolerances,&
                        Tsaved_1,Tsaved_2,transform1_saved,transform2_saved)
                end if IF102
             end if IF101
          end do loop110


       end do MAINLOOP2
    end do MAINLOOP1


!!!-----------------------------------------------------------------------------
!!! Convert the 2x2 transformations to 3x3 matrices
!!!-----------------------------------------------------------------------------
    allocate(big_T_1(SAV%max_num_matches,3,3))
    allocate(big_T_2(SAV%max_num_matches,3,3))
    big_T_1(:,:,:) = 0
    big_T_2(:,:,:) = 0
    loop101: do i=1,SAV%max_num_matches
       big_T_1(i,3,3) = 1
       big_T_2(i,3,3) = 1
    end do loop101
    loop103: do i=1,SAV%max_num_matches
       big_T_1(i,:2,:2) = (Tsaved_1(i,:,:))
       big_T_2(i,:2,:2) = (Tsaved_2(i,:,:))
    end do loop103


!!!-----------------------------------------------------------------------------
!!! Combine 3x3 planecutter matrix with 3x3 plane matching matrix 
!!!-----------------------------------------------------------------------------
    allocate(comb_trans_1(SAV%max_num_matches,3,3))
    allocate(comb_trans_2(SAV%max_num_matches,3,3))
    loop104: do i=1,SAV%max_num_matches
       dummy_mat1(:,:) = big_T_1(i,:,:)
       dummy_mat2(:,:) = transform1_saved(i,:,:)
       comb_trans_1(i,:,:) = matmul((dummy_mat1),(dummy_mat2))

       dummy_mat1(:,:) = big_T_2(i,:,:)
       dummy_mat2(:,:) = transform2_saved(i,:,:)
       comb_trans_2(i,:,:) = matmul((dummy_mat1),(dummy_mat2))
    end do loop104


!!!-----------------------------------------------------------------------------
!!! Reduce transformation matrices if necessary
!!!-----------------------------------------------------------------------------
    write(*,*) "Performing lattice match reduction"
    allocate(lvec1(SAV%max_num_matches))
    lvec1=.false.
    OUTLOOP: do i=1,SAV%max_num_matches
       SAV%tol(i,:) = saved_tolerances(i,:)
       if_reduce: if(reduce)then
          tf = find_tf(comb_trans_1(i,:,:),comb_trans_2(i,:,:))
          if(abs(abs(det(comb_trans_1(i,:,:)))-1._real32).lt.1.E-6_real32) exit if_reduce
          if(verbose.ge.1)then
             write(0,*) i
             write(0,'( 3( 3(F7.3,1X), /) )') tf
          end if
          if(any( (/ (maxval(abs(reduce_vec_gcd(tf(j,:3)))),j=1,3) /)&
               .eq.0._real32))then
             exit if_reduce
          else
             tmat1(:,:) = &
                  reshape((/ 1, 0, 0, 0, 1, 0, 0, 0, 1 /),shape(tmat1(:,:)))
             tmat2(:,:) = nint(tf)
             do j=1,3
                dtmp1=1._real32
                if(any(abs(tf(j,:3)-nint(tf(j,:3))).gt.1.E-6_real32))then
                   dtmp1=get_vec_multiple(tf(j,:3),reduce_vec_gcd(tf(j,:3)))
                end if
                if(abs(dtmp1-nint(dtmp1)).gt.1.E-6_real32)then
                   dtmp1=get_frac_denom(1._real32/dtmp1)
                end if
                tmat1(j,:) = tmat1(j,:3)*nint(dtmp1)
                tmat2(j,:) = nint(tf(j,:3)*dtmp1)
             end do
             if(abs(det(tmat1)).gt.abs(det(comb_trans_1(i,:,:)))) exit if_reduce
             SAV%tf1(i,:,:) = tmat1(:,:)
             SAV%tf2(i,:,:) = tmat2(:,:)
             lvec1(i)=.true.
             cycle OUTLOOP
          end if
       end if if_reduce
       SAV%tf1(i,:,:) = nint(comb_trans_1(i,:,:))
       SAV%tf2(i,:,:) = nint(comb_trans_2(i,:,:))
    end do OUTLOOP
    SAV%tol(:,1) = SAV%tol(:,1)*100._real32
    SAV%tol(:,3) = SAV%tol(:,3)*100._real32
    write(*,*) "Total number of matches saved:",SAV%nfit


!!!-----------------------------------------------------------------------------
!!! Print the set of best matches
!!!-----------------------------------------------------------------------------
    if(verbose.gt.0)then
       do i=1,SAV%nfit
          write(*,'(/,A,I0,2X,A,I0)') &
               "Fit number: ",i,&
               "Area increase: ",&
               nint(get_area(real(SAV%tf1(i,1,:),real32),real(SAV%tf1(i,2,:),real32)))
          write(*,'("   Transmat 1:    Transmat 2:")')
          write(*,'((/,1X,3(3X,A1),3X,3(3X,A1)))') SAV%abc,SAV%abc
          write(*,'(3(/,2X,3(I3," "),3X,3(I3," ")))') &
               SAV%tf1(i,1,1:3),SAV%tf2(i,1,1:3),&
               SAV%tf1(i,2,1:3),SAV%tf2(i,2,1:3),&
               SAV%tf1(i,3,1:3),SAV%tf2(i,3,1:3)
          write(*,'(" vector mismatch (%) = ",F0.9)') SAV%tol(i,1)
          write(*,'(" angle mismatch (°)  = ",F0.9)') SAV%tol(i,2)*180/pi
          write(*,'(" area mismatch (%)   = ",F0.9)') SAV%tol(i,3)
          write(*,*) "reduced:",lvec1(i)
          write(*,*)
       end do
    end if


  end subroutine lattice_matching
!!!#############################################################################

  
!!!!#############################################################################
!!!! Apply the elastic constants to determine strain energy
!!!!#############################################################################
!!!! elastic tensor form (Voight notation):
!!!!  C1111  C1122  C1133  C1123  C1113  C1112
!!!!         C2222  C2233  C2223  C2213  C2212
!!!!                C3333  C3323  C3313  C3312
!!!!                       C2323  C2313  C2312
!!!!                              C1313  C1312
!!!!                                     C1212
!!!! (VASP swaps the final three columns to 12 23 13
!  function compensate_strains(tfmat,w_elastic_tensor,up_elastic_tensor)
!    implicit none
!    integer :: i
!    real(real32), dimension(3) :: strain_vec
!
!    integer, intent(in) :: axis
!    real(real32), dimension(3,3), intent(in) :: lat1,lat2
!    real(real32), dimension(6,6), intent(in) :: elastic_tensor
!
!    
!    ident = 0._real32
!    do i=1,3
!       ident(i,i) = 1._real32
!    end do
!
!    do i=1,3
!       strain_ratio(i) = sum(lw_elastic_tensor(i,:))/sum(up_elastic_tensor(i,:))
!    end do
!    do i=1,3
!       strain_ratio(3+i) = sum(lw_elastic_tensor(i,4:))/sum(up_elastic_tensor(i,4:))
!    end do
!
!    
!
!    
!    
!  end function compensate_strains
!!!!#############################################################################
!
!  
!!!!#############################################################################
!!!! Apply the elastic constants to determine strain energy
!!!!#############################################################################
!!!! elastic tensor form (Voight notation):
!!!!  C1111  C1122  C1133  C1123  C1113  C1112
!!!!         C2222  C2233  C2223  C2213  C2212
!!!!                C3333  C3323  C3313  C3312
!!!!                       C2323  C2313  C2312
!!!!                              C1313  C1312
!!!!                                     C1212
!!!! (VASP swaps the final three columns to 12 23 13
!  function tester(lw_lat,up_lat,lw_tfmat,up_tfmat,lw_elastic,up_elastic) result(stress_vec)
!    implicit none
!    integer :: i
!    real(real32), dimension(6) :: strain_vec, stress_vec
!
!    integer, intent(in) :: axis
!    real(real32), dimension(2,3), intent(in) :: lw_tfmat,up_tfmat
!    real(real32), dimension(3,3), intent(in) :: lw_lat,up_lat
!    real(real32), dimension(6,6), intent(in) :: lw_elastic,up_elastic
!
!
!    ! turn lw_elastic and up_elastic into 3x3x3x3 matrices
!    ! apply the lw_tfmat to lw_elastic and same for up
!    ! ... this should reduce them to a WHAT? sized matrix
!    ! ... this can then be reduced to the Voigt notation, to reduce size and number of variables
!    ! apply transformations to lw_lat and up_lat
!    ! ... then determine the transformation matrix between lw_lat and up_lat.
!    ! now compare the two sets of elements and make them equal.
!    ! apply this to one of the stress matrices. Then we have their ratio timesed by the tfmat already, so we have the amount one will expand or compress. The inverse should be the other.
!
!
!    lw_tflat = matmul(lw_lat,lw_tfmat)
!    up_tflat = matmul(up_lat,up_tfmat)
!    
!    ident = 0._real32
!    do i=1,3
!       ident(i,i) = 1._real32
!    end do
!    
!    strain_mat = matmul(lat1,inverse(lat2))-ident
!    do i=1,3
!       strain_vec(i) = strain_mat(i,i)
!    end do
!    strain_vec(4) = 2._real32*strain_mat(2,3)
!    strain_vec(5) = 2._real32*strain_mat(3,1)
!    strain_vec(6) = 2._real32*strain_mat(1,2)
!
!    stress_vec = matmul(strain_vec,elastic_tensor)
!
!    
!    
!  end function tester
!!!!#############################################################################
!
!  
!!!!#############################################################################
!!!! Apply the elastic constants to determine strain energy
!!!!#############################################################################
!!!! elastic tensor form (Voight notation):
!!!!  C1111  C1122  C1133  C1123  C1113  C1112
!!!!         C2222  C2233  C2223  C2213  C2212
!!!!                C3333  C3323  C3313  C3312
!!!!                       C2323  C2313  C2312
!!!!                              C1313  C1312
!!!!                                     C1212
!!!! (VASP swaps the final three columns to 12 23 13
!  function get_stress(lat1,lat2,axis,elastic_tensor) result(stress_vec)
!    implicit none
!    integer :: i
!    real(real32), dimension(6) :: strain_vec, stress_vec
!
!    integer, intent(in) :: axis
!    real(real32), dimension(3,3), intent(in) :: lat1,lat2
!    real(real32), dimension(6,6), intent(in) :: elastic_tensor
!
!    
!    ident = 0._real32
!    do i=1,3
!       ident(i,i) = 1._real32
!    end do
!    
!    strain_mat = matmul(lat1,inverse(lat2))-ident
!    do i=1,3
!       strain_vec(i) = strain_mat(i,i)
!    end do
!    strain_vec(4) = 2._real32*strain_mat(2,3)
!    strain_vec(5) = 2._real32*strain_mat(3,1)
!    strain_vec(6) = 2._real32*strain_mat(1,2)
!
!    stress_vec = matmul(strain_vec,elastic_tensor)
!
!    
!    
!  end function get_stress
!!!!#############################################################################

  
end module lat_compare
