!!!#############################################################################
!!! Code written by Isiah Edward Mikel Rudkin and Ned Thaddeus Taylor
!!! Code part of the ARTEMIS group (Hepplestone research group).
!!! Think Hepplestone, think HRG.
!!!#############################################################################
module plane_matching
  use constants
  use misc_linalg, only: cross,modu,get_angle,get_area,find_tf,&
       reduce_vec_gcd,gcd
  implicit none
  !! importance of vector, angle, and area
  double precision, dimension(3) :: vaa_weighting=(/1.D0,5.D0,2.5D0/)

  type :: pm_tol_type
     integer :: maxsize,maxfit,nstore
     double precision :: maxlen=20.D0
     double precision :: maxarea=400.D0
     double precision :: vec,ang,area
     double precision :: ang_weight = 10.D0
     double precision :: area_weight = 100.D0
  end type pm_tol_type


!!!updated 2021/11/11


contains
!!!#############################################################################
!!! Subroutine that sorts mainlooplist into ascending order based on ...
!!! ... total tolerance size.              
!!!#############################################################################
  subroutine datasort(list_in,tol_in)
    implicit none
    double precision, dimension(:,:,:) :: list_in
    double precision, allocatable, dimension(:,:,:) :: list_out
    double precision, dimension(:) :: tol_in
    double precision, allocatable, dimension(:) :: tol_out
    integer :: a,dummylocation,len

    len=size(list_in(:,1,1))
    allocate(list_out(len,size(list_in(1,:,1)),size(list_in(1,1,:))))
    allocate(tol_out(len))

    LOOP: do a=1,len
       dummylocation = minloc(tol_in,dim=1)
       tol_out(a) = tol_in(dummylocation)
       list_out(a,:,:) = list_in(dummylocation,:,:)
       tol_in(dummylocation) = INF
    end do LOOP
    tol_in = tol_out
    list_in = list_out

  end subroutine datasort
!!!#############################################################################


!!!#############################################################################
!!! Subroutine that sorts saved_tolerances into ascending order based on ...
!!! ... total tolerance size  
!!!#############################################################################
  subroutine datasortmain(list_in,mat1_in,mat2_in,trans1_in,trans2_in)
    implicit none
    integer :: len
    integer :: a,dummylocation
    double precision, dimension(:,:,:) :: mat1_in,mat2_in
    double precision, allocatable, dimension(:,:,:) :: mat1_out,mat2_out
    double precision, dimension(:,:,:) :: trans1_in,trans2_in
    double precision, allocatable, dimension(:,:,:) :: trans1_out,trans2_out
    double precision, dimension(:) :: list_in
    double precision, allocatable, dimension(:) :: list_out


    len = size(list_in)
    allocate(mat1_out(len,size(mat1_in(1,:,1)),size(mat1_in(1,1,:))))
    allocate(mat2_out(len,size(mat2_in(1,:,1)),size(mat2_in(1,1,:))))
    allocate(trans1_out(len,size(trans1_in(1,:,1)),size(trans1_in(1,1,:))))
    allocate(trans2_out(len,size(trans2_in(1,:,1)),size(trans2_in(1,1,:))))
    allocate(list_out(len)) 

    LOOP: do a=1,len
       dummylocation = minloc(list_in,dim=1)
       list_out(a) = list_in(dummylocation)
       mat1_out(a,:,:) = mat1_in(dummylocation,:,:)
       mat2_out(a,:,:) = mat2_in(dummylocation,:,:)
       trans1_out(a,:,:) = trans1_in(dummylocation,:,:)
       trans2_out(a,:,:) = trans2_in(dummylocation,:,:)
       list_in(dummylocation) = INF
    end do LOOP
    mat1_in = mat1_out
    mat2_in = mat2_out
    list_in = list_out
    trans1_in = trans1_out
    trans2_in = trans2_out

  end subroutine datasortmain
!!!#############################################################################


!!!#############################################################################
!!! Subroutine that sorts list into ascending order based on ...
!!! ... individual tolerance size.              
!!!#############################################################################
  subroutine datasort_tols(list_in,tol_in)
    implicit none
    integer :: i,j,len,ntol_features
    double precision, allocatable,dimension(:) :: vtmp1
    double precision, dimension(:,:,:) :: list_in
    double precision, allocatable, dimension(:,:,:) :: list_out
    double precision, dimension(:,:) :: tol_in
    double precision, allocatable, dimension(:,:) :: tol_out,tmp_store

    ntol_features = size(tol_in(1,:)) 
    len = size(list_in(:,1,1))
    allocate(list_out(len,size(list_in(1,:,1)),size(list_in(1,1,:))))
    allocate(tol_out(len,ntol_features))
    allocate(vtmp1(ntol_features))

    allocate(tmp_store(size(list_out(1,:,1)),size(list_out(1,1,:))))

    tol_out = tol_in
    list_out = list_in
    do i=1,len
       do j=i+1,len
          !if( all(tol_out(j,:) .le. tol_out(i,:) ) )then
          if( dot_product(tol_out(j,:),vaa_weighting).le.&
               dot_product(tol_out(i,:),vaa_weighting) )then
             vtmp1 = tol_out(i,:)
             tol_out(i,:) = tol_out(j,:)
             tol_out(j,:) = vtmp1
             tmp_store(:,:) = list_out(i,:,:)
             list_out(i,:,:) = list_out(j,:,:)
             list_out(j,:,:) = tmp_store(:,:)
          end if
       end do
    end do
    tol_in = tol_out
    list_in = list_out

  end subroutine datasort_tols
!!!#############################################################################


!!!#############################################################################
!!! Subroutine that sorts saved_tolerances into ascending order based on ...
!!! ... total tolerance size  
!!!#############################################################################
  subroutine datasortmain_tols(tol,mat1,mat2,trans1,trans2)
    implicit none
    integer :: i,j,len
    double precision, dimension(3) :: vtmp1
    double precision, dimension(2,2) :: dmat1
    double precision, dimension(3,3) :: dmat2
    double precision, dimension(:,:,:) :: mat1,mat2
    double precision, dimension(:,:,:) :: trans1,trans2
    double precision, dimension(:,:) :: tol


    len=size(tol,dim=1)

    do i=1,len
       do j=i+1,len,1
          !if( all(tol(j,:) .le. tol(i,:) ) )then
          if( dot_product(tol(j,:),vaa_weighting).le.&
               dot_product(tol(i,:),vaa_weighting) )then
             vtmp1=tol(i,:)
             tol(i,:)=tol(j,:)
             tol(j,:)=vtmp1
             
             dmat1 = mat1(j,:,:)
             mat1(j,:,:) = mat1(i,:,:)
             mat1(i,:,:) = dmat1
             
             dmat1 = mat2(j,:,:)
             mat2(j,:,:) = mat2(i,:,:)
             mat2(i,:,:) = dmat1
             
             dmat2 = trans1(j,:,:)
             trans1(j,:,:) = trans1(i,:,:)
             trans1(i,:,:) = dmat2
             
             dmat2 = trans2(j,:,:)
             trans2(j,:,:) = trans2(i,:,:)
             trans2(i,:,:) = dmat2
          end if
       end do
    end do


  end subroutine datasortmain_tols
!!!#############################################################################

    
!!!#############################################################################
!!! Function that checks if the matching planes we have found are a ...
!!! ...duplicate of any others already saved to the list
!!!#############################################################################
  function is_duplicate(list1,list2,lat1,lat2,sym1,sym2) result(outval)
    implicit none
    integer :: i,len
    logical :: outval
    double precision, dimension(:,:,:) :: list1,list2 ! The lists of already saved matrices
    double precision, dimension(:,:) :: lat1,lat2 ! The pair of matrices we want to check
    double precision, allocatable, dimension(:,:) :: dummy1,dummy2
    double precision, allocatable, dimension(:,:) :: tmplat1,tmplat2
    double precision, dimension(:,:,:), optional :: sym1,sym2


    len = size(list1(:,1,1))
    outval = .false.
    allocate( dummy1( size( lat1(:,1)), size(lat1(1,:)) ) )
    allocate( dummy2( size( lat2(:,1)), size(lat2(1,:)) ) )
    allocate( tmplat1( size( lat1(:,1)), size(lat1(1,:)) ) )
    allocate( tmplat2( size( lat2(:,1)), size(lat2(1,:)) ) )

    dummy1 = dble(find_tf(lat1,lat2))
    LOOP: do i=1,len
       if(all(abs(list1(i,:,:)).lt.1.D-5)) cycle LOOP
       tmplat1(:,:) = list1(i,:,:)
       tmplat2(:,:) = list2(i,:,:)
       dummy2 = dble(find_tf(tmplat1,tmplat2))

       if ( all(abs( dummy1(:,:)-dummy2(:,:) ) .lt. 1.D-5) ) then
          outval = .true.
 !         write(0,*) "error"
          exit LOOP
!       else 
!          write(0,*) i
!          write(0,'(2(F15.4,2X))') (dummy1(j,:),j=1,2)
!          write(0,*)
!          write(0,'(2(F15.4,2X))') (dummy2(j,:),j=1,2)!(dummy2(j,:),j=1,2)
!          write(0,*)
       end if
    end do LOOP
    

  end function is_duplicate
!!!#############################################################################


!!!#############################################################################
!!! Function that checks for duplicates of the miller vector.
!!! Outputs the boolean logic .true. if the vector can be reduced.
!!!#############################################################################
  function is_unique(miller,sym) result(outval)
    implicit none
    integer :: i,j
    double precision :: tol
    logical :: outval
    integer, dimension(3) :: miller
    double precision, dimension(3) :: vec_in,vec_out,vec_tmp1,vec_tmp2
    double precision, dimension(:,:,:) :: sym
 
!    if(dot_product(vec_out-vec_in,vec_out-vec_in).lt.1.D-5)
!    if(all(abs(vec_out-vec_in).lt.1.D-5))
!    any(vec_in.eq.3.D0)
!    all(vec_in.eq.3.D0)

    outval = .true.
    vec_in = dble(miller)
    vec_out = reduce_vec_gcd(vec_in)

    if (all(miller.eq.0)) then
       outval = .false.
    else if (all(abs(vec_out-vec_in).lt.1.D-5)) then
       outval = .true.
    else 
       outval = .false.
    end if
    if(.not.outval) return


    tol=1.D-5
    if(all(vec_in.le.0.D0))then
       outval=.false.
       return
    end if
    signloop1: do j=1,3
       if(abs(vec_in(j)).lt.tol) cycle signloop1
       vec_in=sign(1.D0,vec_in(j))*vec_in
       exit signloop1
    end do signloop1

    symloop1: do i=1,size(sym(:,1,1),dim=1)
       vec_out=matmul(vec_in,sym(i,:3,:3))
       if(all(abs(vec_out-vec_in).lt.tol)) cycle symloop1
       vec_tmp1(:)=abs(vec_in(:))-abs(vec_out(:))
       vec_tmp2(:)=vec_in(:)-vec_out(:)
       symloop2: do j=1,3
          if(vec_tmp1(j).gt.tol.or.&
               (abs(vec_tmp1(j)).lt.tol.and.vec_tmp2(j).lt.-tol))then
             outval=.false.
             exit symloop1
          elseif(vec_tmp1(j).lt.-tol)then
             cycle symloop1
          end if
       end do symloop2
    end do symloop1



  end function is_unique
!!!#############################################################################


!!!#############################################################################
!!! Checks whether vec1 is a unique vector after symmetry transformation
!!! This is used to check that the following match is caught if lat2's a=b
!!! lat1:
!!!   1 0
!!!   0 1
!!! lat2:
!!!   1 0     or    0 1
!!!   0 1           1 0
!!!#############################################################################
!!! WARNING!!! NOT CURRENTLY USED, CHECK USE OF IT
  function is_unique_set(vec1,vec2,sym) result(outval)
    implicit none
    integer :: i,j
    double precision :: tol
    integer, dimension(2) :: vec1,vec2
    double precision, dimension(3) :: vec_in,vec_out,vec_tmp1,vec_tmp2
    double precision, dimension(:,:,:) :: sym
    logical :: outval
    

    tol=1.D-5
    outval=.true.
    vec_in=(/ dble(vec1(1)), dble(vec1(2)), 0.D0/)
    !vec_in1=(/ dble(vec1(1)), dble(vec1(2)), 0.D0/)
    !vec_in2=(/ dble(vec2(1)), dble(vec2(2)), 0.D0/)

    symloop1: do i=1,size(sym(:,1,1),dim=1)
       ! matmul inmat with sym
       ! then compare to mat_checklist
       vec_out=matmul(vec_in,sym(i,:3,:3))
       if(all(abs(vec_out-vec_in).lt.tol)) cycle symloop1
       vec_tmp1(:)=abs(vec_in(:))-abs(vec_out(:))
       vec_tmp2(:)=vec_in(:)-vec_out(:)
       symloop2: do j=1,3
          if(vec_tmp1(j).gt.tol.or.&
               (abs(vec_tmp1(j)).lt.tol.and.vec_tmp2(j).lt.-tol))then
             outval=.false.
             exit symloop2
          elseif(vec_tmp1(j).lt.-tol)then
             cycle symloop1
          end if
       end do symloop2
    end do symloop1

    !tol=1.D-5
    !outval=.true.
    !vec_in=(/ dble(vec1(1)), dble(vec1(2)), 0.D0/)
    !
    !symloop1: do i=1,size(sym(:,1,1),dim=1)
    !   vec_out=matmul(vec_in,sym(i,:3,:3))
    !   if(all(abs(vec_out-vec_in).lt.tol)) cycle symloop1
    !   vec_tmp1(:)=abs(vec_in(:))-abs(vec_out(:))
    !   vec_tmp2(:)=vec_in(:)-vec_out(:)
    !   symloop2: do j=1,3
    !      if(vec_tmp1(j).gt.tol.or.&
    !           (abs(vec_tmp1(j)).lt.tol.and.vec_tmp2(j).lt.-tol))then
    !         outval=.false.
    !         exit symloop2
    !      elseif(vec_tmp1(j).lt.-tol)then
    !         cycle symloop1
    !      end if
    !   end do symloop2
    !end do symloop1


  end function is_unique_set
!!!#############################################################################


!!!#############################################################################
!!! This function needs to compare previous lattice matches by performing ...
!!! ... symmetry operations on them to see if they are identical matches
!!!#############################################################################
  function is_unique_match(sym1,sym2,check_set,test_list,lw_check,up_check,up_list) result(lunique)
    implicit none
    integer :: i,isym,jsym
    integer :: nlist,matched_loc
    double precision :: tol
    logical :: lunique
    double precision, dimension(2,2) :: mat1,mat2,tf
    double precision, dimension(2,4) :: inmat
    double precision, allocatable, dimension(:,:,:) :: tf_testlist,mat_testlist

    double precision, dimension(:,:,:), intent(in) :: sym1,sym2

    double precision, dimension(2,4), optional, intent(in) :: check_set
    double precision, dimension(:,:,:), intent(inout), optional :: test_list
    double precision, dimension(4), optional, intent(in) :: lw_check,up_check
    double precision, dimension(:,:), optional, intent(inout) :: up_list

    !logical :: ltest_print
    !logical, optional, intent(in) :: ltest

    !ltest_print=.false.
    !if(present(ltest)) ltest_print=ltest
    
    !! test set
    !double precision, dimension(2,2) :: test1,test2
    !test1(1,:) = [ 0, 1 ]
    !test1(2,:) = [ 3, 0 ]
    !test2(1,:) = [ 1, 0 ]
    !test2(2,:) = [ 0, -2 ]


!!!------------------------------------------------------------------------
!!! initialises tolerance and output
!!!------------------------------------------------------------------------
    tol=1.D-5
    lunique = .true.


!!!------------------------------------------------------------------------
!!! checks for whether input matrices or vectors.
!!! converts either into inmat
!!!------------------------------------------------------------------------
    if(present(check_set))then
       inmat = check_set
    else
       inmat(1,1:2) = lw_check(1:2)
       inmat(2,1:2) = lw_check(3:4)
       inmat(1,3:4) = up_check(1:2)
       inmat(2,3:4) = up_check(3:4)
    end if


!!!------------------------------------------------------------------------
!!! checks for whether input list contains lw_tfmat also.
!!! if not, uses inmat(:2,:2) for it
!!!------------------------------------------------------------------------
    if(present(test_list))then
       nlist=size(test_list(:,1,1))
       allocate(mat_testlist, source=test_list)
    else
       nlist=size(up_list(:,1))
       allocate(mat_testlist(nlist,2,4))
       do i=1,nlist
          mat_testlist(i,:2,:2) = inmat(:2,:2)
       end do

       mat_testlist(:,1,3:4) = up_list(:,1:2)
       mat_testlist(:,2,3:4) = up_list(:,3:4)
    end if


!!!------------------------------------------------------------------------
!!! finds tfmat between the list of stored matches
!!!------------------------------------------------------------------------
    allocate(tf_testlist(nlist,2,2))
    do i=1,nlist
       tf_testlist(i,:2,:2) = find_tf(&
            mat_testlist(i,:2,:2),&
            transpose(mat_testlist(i,:2,3:4)))
    end do


!!!------------------------------------------------------------------------
!!! loop to apply symmetries to determine whether input set is unique ...
!!! ... when compared against the list
!!!------------------------------------------------------------------------
    matched_loc = 0
    sym_loop1: do isym=1,size(sym1(:,1,1))
       !mat1 = matmul(inmat(:2,:2),transpose(sym1(isym,:2,:2)))
       mat1 = matmul(inmat(:2,:2),(sym1(isym,:2,:2)))
       do jsym=1,size(sym2(:,1,1))
          !mat2 = matmul(inmat(:2,3:4),transpose(sym2(jsym,:2,:2)))
          mat2 = matmul(inmat(:2,3:4),(sym2(jsym,:2,:2)))
          tf = find_tf(mat1,transpose(mat2))
          !if(ltest_print)then
          !!if(any(ISNAN(tf)))then
          !!if(all(abs(inmat(:2,:2)-test1).lt.tol))then
          !   !if(all(abs(inmat(:2,:2)-test1).lt.tol).and.&
          !   !      all(abs(inmat(:2,3:4)-test2).lt.tol))then
          !   write(0,*) isym,jsym
          !
          !   write(0,'(2(2X,F7.3))') sym1(isym,:2,:2)
          !   write(0,*)
          !   write(0,'(2(2X,F7.3))') sym2(jsym,:2,:2)
          !   write(0,*) "mat1"
          !   write(0,'(2(2X,F7.3))') mat1
          !   write(0,*) "mat2"
          !   write(0,'(2(2X,F7.3))') mat2!inmat(:2,3:4)!mat2
          !   write(0,*) "tf"
          !   write(0,'(2(2X,F7.3))') tf
          !   write(0,*)
          !
          !!   !if(isym.eq.1) stop
          !!   !if(jsym.eq.size(sym2(:,1,1))) stop
          !!   !if(isym.eq.size(sym1(:,1,1))) stop
          !!   !stop
          !end if

          do i=1,nlist
             !if(ltest_print)then
             !!if(all(abs(inmat(:2,:2)-test1).lt.tol).and.&
             !!     all(abs(inmat(:2,3:4)-test2).lt.tol))then
             !   if(i.eq.1)then
             !      write(0,*) "#############################"
             !      write(0,*) "check",i
             !      write(0,'(2(2X,F7.3))') tf_testlist(i,:2,:2)
             !      write(0,*)
             !      write(0,*) "#############################"
             !   end if
             !end if
             !!------------------------------------------------------------
             !! determine whether set is the same as this one in list
             !!------------------------------------------------------------
             if(all(abs(tf - tf_testlist(i,:,:)).lt.tol))then
                lunique = .false.
                matched_loc = i
                exit sym_loop1
                
             end if
          end do

       end do
    end do sym_loop1


!!!------------------------------------------------------------------------
!!! saves the smallest match if successful
!!!------------------------------------------------------------------------
    if(.not.lunique)then
       if(abs(get_area(inmat(:2,:2),inmat(:2,3:4))).lt.&
            abs(&
            get_area(mat_testlist(matched_loc,:2,:2),&
            mat_testlist(matched_loc,:2,3:4))))then
          mat_testlist(matched_loc,:2,:4) = inmat(:2,:4)
          if(present(test_list))then
             test_list = mat_testlist
          else
             up_list(matched_loc,1:2) = inmat(1,3:4)
             up_list(matched_loc,3:4) = inmat(2,3:4)
          end if
       end if
    end if


    !if(all(abs(inmat(:2,:2)-test1).lt.tol).and.&
    !     all(abs(inmat(:2,3:4)-test2).lt.tol)) stop

  end function is_unique_match
!!!#############################################################################





!!!#############################################################################
!!!#############################################################################
!!!                        M A I N   S E C T I O N
!!!#############################################################################
!!!#############################################################################




  
!!!#############################################################################
!!! Function to match lattices of two input planes.
!!! Matching plane (ab) from each lattice.
!!!#############################################################################
  subroutine cell_match(&
       tol,lat1,lat2,&
       transforms1,transforms2,&
       ntransforms,matched_tols,sym1,sym2)
    implicit none
    integer :: i,j,l,m,total_list_count,nvec1,nvec2, k
    real :: tol_up_ang,tol_dw_ang,tol_up_vec,tol_dw_vec
    double precision :: tiny
    double precision :: reference_mag,considered_mag
    double precision :: reference_angle,considered_angle
    type(pm_tol_type) :: tol
    double precision, dimension(3) :: lat1_veca,lat1_vecb,lat2_veca,lat2_vecb
    double precision, dimension(tol%maxfit) :: MAIN_LOOP_LIST_TOLERANCES
    !double precision, dimension(:) :: MAIN_LOOP_LIST_TOLERANCES
    integer, dimension(2,6) :: tmpmat
    double precision, dimension(2,2) :: tf,mat1,mat2
    double precision, dimension(2,3) :: considered_vectors
    double precision, dimension(3,3) :: lat1,lat2
    double precision, dimension(1000,3) :: tmp_tolerances
    double precision, allocatable, dimension(:,:) :: matched_tols
    double precision, dimension(tol%maxfit,2,4) :: MAIN_LOOP_LIST

    integer :: ntransforms
    !! The 2x2 transformation matrices output by the code.
    !! allocated when we know how many fits.
    integer, allocatable, dimension(:,:,:) :: transforms1,transforms2
    integer, allocatable, dimension(:,:) :: numstore_1,numstore_2
    integer, allocatable, dimension(:,:) :: iarrtmp1
    double precision, allocatable, dimension(:,:) :: latstore_1,latstore_2
    double precision, allocatable, dimension(:,:) :: darrtmp1
    double precision, dimension(:,:,:), optional :: sym1,sym2
  

!!! Layout of each of the 1000 cells:
!!!
!!! (int num of latvec1a, int no. of latvec1b), (int num of latvec2a, int num of latvec2b)
!!! (int num of latvec1a, int no. of latvec1b), (int num of latvec2a, int num of latvec2b)
!!!                                        =
!!! (Vector in the first plane) that matches in magnitude to (Vector in the second plane)
!!! (Vector in the first plane) that matches in magnitude to (Vector in the second plane)


  !! Number of entries in each of the lists.
  integer :: len_list_1a,len_list_1b 
  integer :: len_list_final !Length of final list of compatible vector pairs after angle check

  !! list of vec combins (of 2a and 2b) that fit vec lat1_a, mag of tol on fit
  double precision, dimension(1000,3) :: list_1a 
  !! layout: int num of lat2_a, int num of lat2_b, tol
  !! list of vec combins (of 2a and 2b) that fit vec lat1_b, mag of tol on fit
  double precision, dimension(1000,3) :: list_1b 
  !! layout: integer number of lat2_a, integer number of lat2_b, tol

  double precision, dimension(1000,5) :: list_angle_fits
  !Layout:
  ! First 2 components(1-2); integer number of lat2_a, integer number of lat2_b
  ! Next 2 components(3-4); integer number of lat2_a, integer number of lat2_b
  ! Last Component (5); Total weighted tolerance on everything


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Setting up tolerances !!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  tiny = 1.D-5
  tol_up_ang = 1.E0 + real(tol%ang)/(2.E0*pi)
  tol_dw_ang = 1.E0 - real(tol%ang)/(2.E0*pi)
  tol_up_vec = 1.E0 + real(tol%vec)!/100.D0
  tol_dw_vec = 1.E0 - real(tol%vec)!/100.D0

  if(allocated(matched_tols)) deallocate(matched_tols)
  allocate(matched_tols(tol%maxfit,3))
  MAIN_LOOP_LIST_TOLERANCES(:) = INF
  tmp_tolerances(:,:) = INF
  matched_tols(:,:) = INF


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Assign the vectors a and b for the first and second lattice,  !!!
!!! -   1a and 1b refer to the a and b lattice vectors for lattice 1
!!!       These vectors form the planes we want to match.         !!!
!!! -   2a and 2b refer to the a and b lattice vectors for lattice 2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  lat1_veca = lat1(1,:)
  lat1_vecb = lat1(2,:)
  lat2_veca = lat2(1,:)
  lat2_vecb = lat2(2,:)

 
!!!------------------------------------------------------------------------
!!! set up the vectors on lower plane
!!!------------------------------------------------------------------------
  nvec1=0
  allocate(numstore_1((2*(tol%maxsize+1))**2,2))
  allocate(latstore_1((2*(tol%maxsize+1))**2,3))
  vecmakeloop1: do l=0,tol%maxsize
     pmloop1: do i=1,-1,-2
        vecmakeloop2: do m=0,tol%maxsize
           if (l.eq.0 .and. m.eq.0) cycle vecmakeloop2
           pmloop2: do j=1,-1,-2
              nvec1=nvec1+1
              numstore_1(nvec1,:) = (/ i*l, j*m /)
              latstore_1(nvec1,:) = dble(i*l) * lat1_veca + dble(j*m) * lat1_vecb
              if(abs(modu(latstore_1(nvec1,:))).gt.tol%maxlen)then
                 nvec1=nvec1-1
                 cycle pmloop1
              end if
           end do pmloop2
        end do vecmakeloop2
     end do pmloop1
  end do vecmakeloop1
  allocate(iarrtmp1(nvec1,2))
  allocate(darrtmp1(nvec1,3))
  iarrtmp1(:nvec1,:)=numstore_1(:nvec1,:)
  call move_alloc(iarrtmp1,numstore_1)
  darrtmp1(:nvec1,:)=latstore_1(:nvec1,:)
  call move_alloc(darrtmp1,latstore_1)
  

!!!------------------------------------------------------------------------
!!! set up the vectors on upper plane
!!!------------------------------------------------------------------------
  nvec2=0
  allocate(numstore_2((2*(tol%maxsize+1))**2,2))
  allocate(latstore_2((2*(tol%maxsize+1))**2,3))
  vecmakeloop3: do l=0,tol%maxsize
     pmloop3: do i=1,-1,-2
        vecmakeloop4: do m=0,tol%maxsize
           if (l.eq.0 .and. m.eq.0) cycle vecmakeloop4
           pmloop4: do j=1,-1,-2
              nvec2=nvec2+1
              numstore_2(nvec2,:) = (/ i*l, j*m /)
              latstore_2(nvec2,:) = dble(i*l) * lat2_veca + dble(j*m) * lat2_vecb
              if(modu(latstore_2(nvec2,:)).gt.tol%maxlen)then
                 nvec2=nvec2-1
                 cycle vecmakeloop3
              end if
           end do pmloop4
        end do vecmakeloop4
     end do pmloop3
  end do vecmakeloop3
  allocate(iarrtmp1(nvec2,2))
  allocate(darrtmp1(nvec2,3))
  iarrtmp1(:nvec2,:)=numstore_2(:nvec2,:)
  call move_alloc(iarrtmp1,numstore_2)
  darrtmp1(:nvec2,:)=latstore_2(:nvec2,:)
  call move_alloc(darrtmp1,latstore_2)


!!!------------------------------------------------------------------------
!!! lower lattice vector 1 loop
!!!------------------------------------------------------------------------
  total_list_count = 0
  MAINLOOP1: do l=1,nvec1
     tmpmat(1,:2) = numstore_1(l,:2)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Finding the best fit options for the first lattice vector. !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     len_list_1a = 0
     reference_mag = modu(latstore_1(l,:)) !! mag of lattice vector to fit to
     loop102: do j=1,nvec2                      !! matcing to 1st lattice vector
        if(len_list_1a.ge.1000) exit loop102
        if(gcd( [ numstore_1(l,1),numstore_1(l,2),&
             numstore_2(j,1),numstore_2(j,2) ]).ne.1) cycle loop102
        considered_mag = modu(latstore_2(j,:))  !! Get magnitude of the vector
        !! Checking the fit (if too big or too small)
        if ( ( considered_mag .ge. (tol_dw_vec * reference_mag) ) .and. &
             ( considered_mag .le. (tol_up_vec * reference_mag) ) ) then
           len_list_1a = len_list_1a + 1
           list_1a(len_list_1a,1:2) = numstore_2(j,1:2)
           list_1a(len_list_1a,3) = abs((considered_mag - reference_mag)/reference_mag)
        end if
     end do loop102


!!!------------------------------------------------------------------------
!!! lower lattice vector 2 loop
!!!------------------------------------------------------------------------
     MAINLOOP2: do m=1,nvec1
        tmpmat(2,:2) = numstore_1(m,:2)
        if(all(latstore_1(l,:).eq.latstore_1(m,:))) cycle MAINLOOP2
        if(get_area(latstore_1(l,:),latstore_1(m,:)).gt.tol%maxarea) cycle MAINLOOP2
        if(all(cross(latstore_1(l,:),latstore_1(m,:)).eq.0.D0)) cycle MAINLOOP2
        reference_angle = get_angle(latstore_1(l,:),latstore_1(m,:))
        if (abs(reference_angle) .lt. tiny) cycle MAINLOOP2 
        
        !!! CHANGE IT TO TAKE IN A 2x2 MATRIX LATER !!!
        if(modu(latstore_1(l,:)).gt.modu(latstore_1(m,:))) cycle MAINLOOP2
        if(dot_product(latstore_1(l,:),latstore_1(m,:)).gt.&
             (0.5D0*dot_product(latstore_1(l,:),latstore_1(l,:))))&
             cycle MAINLOOP2
        !SHOULD I REMOVE THIS? WHAT DOES IT WANT TO DO?
        !if(.not.is_unique_set(numstore_1(l,:),numstore_1(m,:),sym1)) &
        !     cycle MAINLOOP2


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Finding the best fit options for the second lattice vector. !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        len_list_1b = 0
        reference_mag = modu(latstore_1(m,:))
        loop103: do j=1,nvec2
           if(len_list_1b.ge.1000) exit loop103
           if(gcd( [ numstore_1(m,1),numstore_1(m,2),&
                numstore_2(j,1),numstore_2(j,2) ]).ne.1) cycle loop103
           considered_mag = modu(latstore_2(j,:))
           if ( ( considered_mag .ge. (tol_dw_vec * reference_mag) ) .and. &
                ( considered_mag .le. (tol_up_vec * reference_mag) ) ) then
              len_list_1b = len_list_1b + 1
              list_1b(len_list_1b,1:2) = numstore_2(j,1:2)
              list_1b(len_list_1b,3) = abs((considered_mag - reference_mag)/reference_mag)
           end if
        end do loop103


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Checking the angle between all possible sets of vectors.   !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
        len_list_final = 0
        loop109: do i=1, len_list_1a
           tmpmat(1,3:4) = nint(list_1a(i,:2))
           loop110: do j=1, len_list_1b
              if(len_list_final.ge.1000) exit loop109
              considered_vectors(1,:) = list_1a(i,1)*lat2_veca + list_1a(i,2)*lat2_vecb
              considered_vectors(2,:) = list_1b(j,1)*lat2_veca + list_1b(j,2)*lat2_vecb
              considered_angle = &
                   get_angle(considered_vectors(1,:),considered_vectors(2,:))
              !if(.not.is_unique_set(nint(list_1a(i,:2)),nint(list_1b(j,:2)),sym2)) &
              !     cycle loop110
              !!--------------------------------------------------------------------------
              !! checks angle of lat2 with the reference angle from lat1
              !! if they are too disimilar, cycle
              !!--------------------------------------------------------------------------
              if (considered_angle .lt. tol_dw_ang*reference_angle) then
                 cycle loop110
              else if (considered_angle .gt. tol_up_ang*reference_angle) then
                 cycle loop110
              else
                 tmpmat(2,3:4) = nint(list_1b(j,:2))
                 !write(0,'(A,4X,"[",I3,I3,",",I3,I3,"]",6X,"[",I3,I3,",",I3,I3,"]",4X,I2,2X,I2,2X,F0.3)') &
                 !     "HERE",numstore_1(l,:2),numstore_1(m,:2),nint(list_1a(i,:2)),nint(list_1b(j,:2)),&
                 !     total_list_count,len_list_final, considered_angle
                 if(total_list_count.ne.0)then
                    if(.not.is_unique_match( sym1, sym2, &
                         check_set = dble(tmpmat),&
                         test_list = MAIN_LOOP_LIST(:total_list_count,:2,:4)))&
                         cycle loop110
                 end if
                 if(len_list_final.ne.0)then
                    if(.not.is_unique_match( sym1, sym2, &
                         check_set = dble(tmpmat),&
                         up_list = list_angle_fits(:len_list_final,:4)))&
                         cycle loop110
                 end if
                 !write(0,*) "PAST HERE", list_angle_fits(len_list_final,:)

                 len_list_final = len_list_final + 1
                 list_angle_fits(len_list_final,1:2) = list_1a(i,1:2)
                 list_angle_fits(len_list_final,3:4) = list_1b(j,1:2)
                 tmp_tolerances(len_list_final,1) = &
                      max(list_1a(i,3),list_1b(j,3))
                 tmp_tolerances(len_list_final,2) = &
                      abs(considered_angle-reference_angle)
                 tmp_tolerances(len_list_final,3) = abs(1.D0 - &
                      get_area(considered_vectors(1,:),considered_vectors(2,:))&
                      /get_area(latstore_1(l,:),latstore_1(m,:)))
                 list_angle_fits(len_list_final,5) = &
                      tol%ang_weight * abs(considered_angle-reference_angle) + &
                      list_1a(i,3) + list_1b(i,3) + &
                      tol%area_weight*get_area(latstore_1(l,:),latstore_1(m,:))
              end if
           end do loop110
        end do loop109
    
        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Searching for the best (max_matches) fits and forcing the !!!
!!! output list down to that size                             !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        loop112: do i=1, len_list_final
           mat1(1,:2)=dble(numstore_1(l,:2))
           mat1(2,:2)=dble(numstore_1(m,:2))
           mat2(1,:2)=dble(list_angle_fits(i,1:2))
           mat2(2,:2)=dble(list_angle_fits(i,3:4))
           tf=find_tf(mat1,mat2)
           do j=1,tol%maxfit
              if(all(abs(tf-find_tf(dble(MAIN_LOOP_LIST(j,:2,1:2)),&
                   dble(MAIN_LOOP_LIST(j,:2,3:4)))).lt.1.D-6))then
                 cycle loop112
              end if
           end do

           !!-------------------------------------------------------------------
           !! Filling in the first (tol%maxfit) places.
           !!-------------------------------------------------------------------
           if (total_list_count .lt. tol%maxfit) then
              total_list_count = total_list_count + 1
              MAIN_LOOP_LIST(total_list_count,1,1:2) = numstore_1(l,1:2)
              MAIN_LOOP_LIST(total_list_count,2,1:2) = numstore_1(m,1:2)
              MAIN_LOOP_LIST(total_list_count,1,3:4) = list_angle_fits(i,1:2)
              MAIN_LOOP_LIST(total_list_count,2,3:4) = list_angle_fits(i,3:4)
              MAIN_LOOP_LIST_TOLERANCES(total_list_count) = list_angle_fits(i,5)
              matched_tols(total_list_count,:) = tmp_tolerances(i,:)
              cycle loop112
           end if

           !!-------------------------------------------------------------------
           !! Sorts the data into order of tolerance
           !!-------------------------------------------------------------------
           call datasort_tols(MAIN_LOOP_LIST,matched_tols)
           if ( all(tmp_tolerances(i,:) .le. matched_tols(tol%maxfit,:)) ) then
              MAIN_LOOP_LIST(tol%maxfit,1,1:2) = numstore_1(l,1:2)
              MAIN_LOOP_LIST(tol%maxfit,2,1:2) = numstore_1(m,1:2)
              MAIN_LOOP_LIST(tol%maxfit,1,3:4) = list_angle_fits(i,1:2)
              MAIN_LOOP_LIST(tol%maxfit,2,3:4) = list_angle_fits(i,3:4)
              MAIN_LOOP_LIST_TOLERANCES(tol%maxfit) = list_angle_fits(i,5)
              matched_tols(tol%maxfit,:) = tmp_tolerances(i,:)
           end if
        end do LOOP112
        !call datasort(MAIN_LOOP_LIST,MAIN_LOOP_LIST_TOLERANCES)
        call datasort_tols(MAIN_LOOP_LIST,matched_tols)


     end do MAINLOOP2
  end do MAINLOOP1
  ntransforms = total_list_count


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Sorts the data from main loop list into transformation matrices for output !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (allocated(transforms1)) deallocate(transforms1) !! From previous iterations
  if (allocated(transforms2)) deallocate(transforms2) !! ""

  allocate(transforms1(ntransforms,2,2))
  allocate(transforms2(ntransforms,2,2))

  loop114: do i=1,ntransforms
     transforms1(i,1:2,1:2) = nint(MAIN_LOOP_LIST(i,1:2,1:2))
     transforms2(i,1:2,1:2) = nint(MAIN_LOOP_LIST(i,1:2,3:4))
  end do loop114


end subroutine cell_match
!!!#############################################################################
  
 
end module plane_matching
