module artemis__interface_identifier
  !! Module for identifying interface locations in crystal structures.
  !!
  !! Provides methods for interface identification using cumulative atomic
  !! density (CAD) and density of neighbours (DON) approaches.
  use coreutils__kind, only: real32
  use coreutils__array, only: swap, sort1D
  use artemis__misc_linalg, only: simeq,get_area,uvec
  use artemis__misc_maths, only: gauss_array,get_turn_points,overlap_indiv_points,&
       running_avg,mean
  use atomstruc, only: basis_type, species_type
  implicit none

  private

  integer, parameter :: nstep_default=1000
  !! Default number of steps for density calculations.

  type intf_info_type
     !! Derived type storing interface information.
     integer :: axis
     !! Axis perpendicular to the interface.
     real(real32), dimension(2) :: loc
     !! Locations of the two interfaces along the axis.
  end type intf_info_type
  
  type den_of_neigh_type
     !! Derived type storing density of neighbours for each atom.
     real(real32), allocatable, dimension(:,:) :: atom
     !! Density of neighbours array (atom, step).
  end type den_of_neigh_type
  type den_of_spec_type
     !! Derived type storing species-dependent density of neighbours.
     real(real32), allocatable, dimension(:,:,:) :: atom
     !! Species density array (atom, species, step).
  end type den_of_spec_type


  public :: nstep_default
  public :: intf_info_type,den_of_neigh_type
  public :: get_interface
  public :: gen_DON,gen_DONsim
  public :: get_layered_axis
  public :: gen_single_DOS,gen_single_DON


contains
!###############################################################################
  function get_interface(basis, axis) result(intf)
    !! Get the interface location using the cumulative atomic density (CAD) method.
    implicit none

    ! Arguments
    type(basis_type), intent(in) :: basis
    !! Input basis type containing atomic positions.
    integer, optional, intent(in) :: axis
    !! Axis perpendicular to the interface. Determined automatically if absent.

    ! Local variables
    integer :: nstep
    !! Number of steps for density calculation.
    real(real32) :: dist_max
    !! Maximum distance for density calculation.
    type(intf_info_type) :: intf
    !! Interface information result.
    type(den_of_spec_type), allocatable, dimension(:) :: DOS
    !! Species-dependent density of neighbours.


    dist_max = 12._real32
    DOS = gen_DOS(basis%lat,basis,dist_max)
    nstep = size(DOS(1)%atom(1,1,:))

    intf%axis = 0
    if(present(axis)) intf%axis = axis
    if(intf%axis.eq.0)then
       intf%axis = get_intf_axis_DOS(DOS, basis%lat, basis, dist_max)
    end if

    intf%loc=get_intf_CAD(basis%lat, basis, intf%axis, nstep)

    if(intf%loc(1).gt.intf%loc(2)) call swap(intf%loc(1),intf%loc(2))

  end function get_interface
!###############################################################################


!###############################################################################
  function gen_DOS(lat,bas,dist_max,scale_dist,norm) result(DOS)
    !! Generate species-dependent density of neighbours (DOS).
    implicit none

    ! Arguments
    real(real32), optional, intent(in) :: dist_max
    !! Maximum interatomic distance for DOS calculation.
    logical, optional, intent(in) :: scale_dist,norm
    !! Whether to scale by distance and whether to normalise.
    type(basis_type), intent(in) :: bas
    !! Input basis type containing atomic positions.
    real(real32), dimension(3,3), intent(in) :: lat
    !! Lattice vectors.

    ! Local variables
    integer :: i,j,k,is,ia,js,ja,count1
    !! Loop indices and counter.
    integer :: nstep,nsize
    !! Number of steps and neighbour list size.
    real(real32) :: rdist_max,rtmp1,rtmp2
    !! Local maximum distance and temporary reals.
    logical :: lscale_dist,lnorm
    !! Local flags for distance scaling and normalisation.
    real(real32) :: gauss_tol,DON_sigma,dist
    !! Gaussian tolerance, sigma, and distance value.
    integer, dimension(3) :: ncell
    !! Number of periodic cell images in each direction.
    real(real32), dimension(3) :: vrtmp1,vrtmp2
    !! Temporary 3-vectors.
    real(real32), dimension(3) :: vtmp1,vtmp2,vtmp3
    !! Temporary 3-vectors.
    real(real32), allocatable, dimension(:) :: distance
    !! Distance grid array.
    type(den_of_spec_type), allocatable, dimension(:) :: DOS
    !! Species-dependent density of neighbours result.
    real(real32), allocatable, dimension(:) :: dist_list
    !! List of interatomic distances.


    if(present(scale_dist))then
       lscale_dist=scale_dist
    else
       lscale_dist=.true.
    end if
    if(present(norm))then
       lnorm=norm
    else
       lnorm=.false.
    end if

    nstep=nstep_default
    allocate(DOS(bas%nspec))
    do is=1,bas%nspec
       allocate(DOS(is)%atom(bas%spec(is)%num,bas%nspec,nstep))
    end do

    allocate(distance(nstep))
    rdist_max=12._real32
    if(present(dist_max)) rdist_max=dist_max
    do i=1,nstep
       distance(i)=real(i,real32)*rdist_max/real(nstep,real32)
    end do

    !! should now consider lattice vector addition for obtuse cells.
    !! in obtuse cells, ncells may need to be larger than just individual
    !! distance due to similar paths
    !ncell = 2
    !i = -1
    !rtmp1_old = 1E6
    !ncell_loop1: do while ( i < 10 )
    !   i = i + 1
    !   j = -1
    !   ncell_loop2: do while ( j < 10 )
    !      j = j + 1
    !      k = -1
    !      ncell_loop3: do while ( k < 10 )
    !         k = k + 1
    !         if(i.eq.0.and.j.eq.0.and.k.eq.0) cycle ncell_loop3
    !         rtmp1 = norm2(i*lat(1,:) + j*lat(2,:) + k*lat(3,:))
    !         if(rtmp1.gt.rtmp1_old)then
    !            rtmp1_old = 1E6
    !            exit ncell_loop3
    !         end if
    !         if(rtmp1.le.rdist_max)then
    !            ncell(1) = max(ncell(1),i+1)
    !            ncell(2) = max(ncell(2),j+1)
    !            ncell(3) = max(ncell(3),k+1)
    !         end if
    !         rtmp1_old = rtmp1
    !      end do ncell_loop3
    !   end do ncell_loop2
    !end do ncell_loop1

    ncell = 0
    ncell_loop1: do i=1,3
       rtmp1 = norm2(lat(i,:))
       ncell(i) = max(ncell(i),ceiling(rdist_max/norm2(lat(i,:))))!maxval(ceiling( rdist_max/abs(lat(i,:)) ))
       do j=1,3
          if(i.eq.j) cycle
          rtmp2 = dot_product(lat(i,:),lat(j,:))
          if(sign(1._real32,rtmp1).eq.sign(1._real32,rtmp2)) cycle
          !vrtmp1 = uvec(lat(i,:)) * dot_product(uvec(lat(i,:)),lat(j,:))
          !vrtmp1 = uvec(lat(i,:)) * lat(j,:)
          vrtmp1 = merge(lat(j,:), (/0._real32, 0._real32, 0._real32/), mask = abs(lat(i,:)).gt.1.E-5_real32)
          rtmp1 = norm2(vrtmp1)
          if(abs(rtmp1).lt.1.E-5_real32) cycle
          k = 0
          vrtmp2 = lat(i,:)
          rtmp2 = norm2(vrtmp2)
          do while ( rtmp2 .le. rtmp1)
             k = k + 1
             rtmp1 = rtmp2
             vrtmp2 = lat(i,:) + real(k,real32)*vrtmp1
             rtmp2 = norm2(vrtmp2)
          end do
          if(abs(rtmp1).lt.1.E-5_real32) cycle
          ncell(i) = max(ncell(i), ceiling(rdist_max/rtmp1))
          ncell(j) = max(ncell(j), (k-1)*ceiling(rdist_max/rtmp1))
       end do
    end do ncell_loop1
    !iloop1: do i=1,3
    !   ncell(i) = ceiling( rdist_max/norm2(lat(i,:)) )
    !   jloop1: do j=i+1,3
    !      if(i.eq.j) cycle
    !      itmp1 = ceiling(rdist_max/dot_product(lat(i,:),lat(j,:)))
    !      if(ncell(i).lt.itmp1) ncell(i) = itmp1
    !      if(ncell(j).lt.itmp1) ncell(j) = itmp1
    !   end do jloop1
    !end do iloop1
    !write(0,*) "ncell",ncell
    nsize = bas%natom*(2*ncell(1)+1) * (2*ncell(2)+1) * (2*ncell(3)+1) - 1
    allocate(dist_list(nsize))

    gauss_tol=16.E0!38._real32
    DON_sigma=0.5E-1
    specloop1: do is=1,bas%nspec
       DOS(is)%atom(:,:,:)=0._real32
       atomloop1: do ia=1,bas%spec(is)%num

          specloop2: do js=1,bas%nspec
             count1=0
             dist_list = 0.0
             atomloop2: do ja=1,bas%spec(js)%num
                vtmp1(:3) = bas%spec(is)%atom(ia,:3) - bas%spec(js)%atom(ja,:3)
                do i=-ncell(1),ncell(1),1
                   vtmp2(1) = vtmp1(1) + real(i,real32)
                   do j=-ncell(2),ncell(2),1
                      vtmp2(2) = vtmp1(2) + real(j,real32)
                      kloop1: do k=-ncell(3),ncell(3),1
                         if(is.eq.js.and.ia.eq.ja)then
                            if(i.eq.0.and.j.eq.0.and.k.eq.0)then
                               cycle kloop1
                            end if
                         end if
                         count1=count1+1
                         !if(count1.gt.nsize)then
                         !   write(0,'("ERROR: Internal error in mod_intf_identifier.f90")')
                         !   write(0,'(2X,"Internal error in gen_DOS subroutine.")')
                         !   write(0,'(2X,"dist_list size allocated too small")')
                         !   stop
                         !end if
                         vtmp2(3) = vtmp1(3) + real(k,real32)
                         vtmp3 = matmul(vtmp2,lat)
                         dist_list(count1) = norm2(vtmp3)


                      end do kloop1
                   end do
                end do
             end do atomloop2
             
             DOS(is)%atom(ia,js,:) = &
                  gauss_array(distance,&
                  dist_list(:count1),DON_sigma,gauss_tol,lnorm)

             
          end do specloop2

          if(lscale_dist)then
             do i=minloc(abs(distance(:)-2._real32),dim=1),nstep
                !dist=abs(1._real32/distance(i))**2._real32
                dist=exp(-abs(distance(i)-2.E0))
                DOS(is)%atom(ia,:,i)=DOS(is)%atom(ia,:,i)*dist
             end do
          end if

       end do atomloop1
    end do specloop1


  end function gen_DOS
!###############################################################################


!###############################################################################
  function gen_DON(lat,bas,dist_max,scale_dist,norm) result(DON)
    !! Generate density of neighbours (DON) by summing DOS over all species.
    implicit none

    ! Arguments
    real(real32), optional, intent(in) :: dist_max
    !! Maximum interatomic distance.
    logical, optional, intent(in) :: scale_dist,norm
    !! Whether to scale by distance and whether to normalise.
    type(basis_type), intent(in) :: bas
    !! Input basis type containing atomic positions.
    real(real32), dimension(3,3), intent(in) :: lat
    !! Lattice vectors.

    ! Local variables
    integer :: i,is,ia,nstep
    !! Loop indices and number of steps.
    logical :: lscale_dist,lnorm
    !! Local flags for distance scaling and normalisation.
    type(den_of_spec_type), allocatable, dimension(:) :: DOS
    !! Species-dependent density of neighbours.
    type(den_of_neigh_type), allocatable, dimension(:) :: DON
    !! Density of neighbours result.


    if(present(scale_dist))then
       lscale_dist=scale_dist
    else
       lscale_dist=.true.
    end if
    if(present(norm))then
       lnorm=norm
    else
       lnorm=.false.
    end if

    if(present(dist_max))then
       DOS=gen_DOS(lat,bas,dist_max,scale_dist=lscale_dist,norm=lnorm)
    else
       DOS=gen_DOS(lat,bas,scale_dist=lscale_dist,norm=lnorm)
    end if
    nstep=size(DOS(1)%atom(1,1,:))
    allocate(DON(bas%nspec))
    do is=1,bas%nspec
       allocate(DON(is)%atom(bas%spec(is)%num,nstep))
    end do
       
    do is=1,bas%nspec
       do ia=1,bas%spec(is)%num
          do i=1,nstep
             DON(is)%atom(ia,i) = sum(DOS(is)%atom(ia,:,i))
          end do
       end do
    end do

  end function gen_DON
!###############################################################################


!###############################################################################
  function gen_DONsim(DON,dist_max,cutoff,avg_mthd) result(intf_atoms)
    !! Find unique atoms by comparing density of neighbours.
    implicit none

    ! Arguments
    type(den_of_neigh_type), dimension(:), intent(in) :: DON
    !! Input density of neighbours for each species.
    real(real32), optional, intent(in) :: dist_max,cutoff
    !! Maximum distance and similarity cutoff distance.
    integer, optional, intent(in) :: avg_mthd
    !! Averaging method for similarity cutoff (1=midpoint, 2=mean, 3=max jump).

    ! Local variables
    integer :: i,is,ia,ja,cutloc,itmp1,udef_avg_mthd
    !! Loop indices, cutoff location, temporary integer, and averaging method.
    integer :: nspec,natom,nstep
    !! Number of species, atoms, and steps.
    real(real32) :: avg,rdist_max,rcutoff,maxjump
    !! Average similarity, local max distance, cutoff, and maximum jump.
    integer, allocatable, dimension(:) :: intf_list,sumspec
    !! Interface atom list and cumulative species count.
    real(real32), allocatable, dimension(:) :: newf,simi,distance
    !! Overlap function, similarity array, and distance grid.
    integer, allocatable, dimension(:,:) :: intf_atoms
    !! Interface atoms result (atom index, species/atom number).

    type(den_of_neigh_type), allocatable, dimension(:) :: sim
    !! Similarity per atom.
    type(den_of_spec_type), allocatable, dimension(:) :: similarity
    !! Pairwise similarity between atoms.

    
!-------------------------------------------------------------------------------
! Initialise variables based on input values
!-------------------------------------------------------------------------------
    nstep=size(DON(1)%atom(1,:))
    allocate(distance(nstep))
    rdist_max=12._real32
    if(present(dist_max)) rdist_max=dist_max
    do i=1,nstep
       distance(i)=real(i,real32)*rdist_max/real(nstep,real32)
    end do
    rcutoff=4._real32
    if(present(cutoff)) rcutoff=min(rcutoff,cutoff)
    cutloc=minloc(abs(distance(:)-rcutoff),dim=1)


!-------------------------------------------------------------------------------
! Allocate arrays
!-------------------------------------------------------------------------------
    natom=0
    nspec=size(DON)
    allocate(sim(nspec))
    allocate(similarity(nspec))
    do is=1,nspec
       natom=natom+size(DON(is)%atom(:,1))       
       allocate(sim(is)%atom(&
            size(DON(is)%atom(:,1)),&
            nstep))
       allocate(similarity(is)%atom(&
            size(DON(is)%atom(:,1)),&
            size(DON(is)%atom(:,1)),&
            nstep))
    end do
    allocate(simi(natom))


!-------------------------------------------------------------------------------
! Compute the similarity between the DONs of every atom with each other atom
!-------------------------------------------------------------------------------
    natom=0
    specloop1: do is=1,nspec
       atomloop1: do ia=1,size(DON(is)%atom(:,1))
          natom=natom+1
          atomloop2: do ja=1,size(DON(is)%atom(:,1))
             newf = &
                  overlap_indiv_points(&
                  [DON(is)%atom(ia,:)],&
                  [DON(is)%atom(ja,:)])
             similarity(is)%atom(ia,ja,:)=real(newf,real32)
             deallocate(newf)
          end do atomloop2
          do i=1,nstep
             sim(is)%atom(ia,i)=sum(&
                  similarity(is)%atom(ia,:,i))/size(DON(is)%atom(:,1))
          end do
          simi(natom)=sum(sim(is)%atom(ia,:cutloc))
       end do atomloop1
    end do specloop1


!-------------------------------------------------------------------------------
! Set up intf_list values
!-------------------------------------------------------------------------------
    allocate(sumspec(nspec))
    sumspec(1)=size(DON(1)%atom(:,1))
    do is=2,nspec
       sumspec(is)=sumspec(is-1)+size(DON(is)%atom(:,1))
    end do
    allocate(intf_list(natom))
    do i=1,natom
       intf_list(i)=i
    end do


!-------------------------------------------------------------------------------
! Define similarity cutoff using user-defined averaging method
!-------------------------------------------------------------------------------
    call sort1D(simi,intf_list)
    udef_avg_mthd=1
    if(present(avg_mthd)) udef_avg_mthd=avg_mthd
    select case(udef_avg_mthd)
    case(1,2)
       select case(udef_avg_mthd)
       case(1)
          avg=(minval(simi)+maxval(simi))/2.0
       case(2)
          avg=mean(simi)
       end select
       do i=1,natom
          if(simi(i).ge.avg) exit
          itmp1=i
       end do
    case(3)
       maxjump=0._real32
       do i=2,natom
          if(simi(i)-simi(i-1).gt.maxjump)then
             maxjump = simi(i) - simi(i-1)
             avg=(simi(i)+simi(i-1))/2
             itmp1=i-1
          end if
       end do
    end select
    

!-------------------------------------------------------------------------------
! Save the species and atom numbers of the interfacial atoms
!-------------------------------------------------------------------------------
    allocate(intf_atoms(itmp1,2))
    intf_atoms=0
    do i=1,itmp1
       intf_atoms(i,1)=minloc(sumspec(:)-intf_list(i),&
            mask=sumspec(:)-intf_list(i).ge.0,dim=1)
       intf_atoms(i,2)=intf_list(i)
       if(intf_atoms(i,1).gt.1)then
          intf_atoms(i,2)=intf_atoms(i,2)-sumspec(intf_atoms(i,1)-1)
       end if
    end do


  end function gen_DONsim
!###############################################################################


!###############################################################################
! returns interface ions from two supplied bases
!###############################################################################
!  subroutine get_intf_atoms(lat1,bas1,lat2,bas2)
!    implicit none
!    real(real32), dimension(3,3) :: lat1,lat2
!    integer, allocatable, dimension(:,:) :: intf_list1,intf_list2
!
!
!    intf_list1=gen_DONsim(gen_DON(lat1,bas1),cutoff=4.0)
!    intf_list2=gen_DONsim(gen_DON(lat2,bas2),cutoff=4.0)
!
!
!  end subroutine get_intf_atoms
!###############################################################################


!###############################################################################
  function get_intf_axis_DOS(DOS,lat,bas,dist_max,cutoff,lprint) result(axis)
    !! Determine the axis perpendicular to the interface using the DOS method.
    implicit none

    ! Arguments
    real(real32), optional, intent(in) :: dist_max,cutoff
    !! Maximum distance and similarity cutoff distance.
    logical, optional :: lprint
    !! Whether to print progress information.
    type(basis_type) :: bas
    !! Input basis type containing atomic positions.
    real(real32), dimension(3,3) :: lat
    !! Lattice vectors.
    type(den_of_spec_type), allocatable, dimension(:) :: DOS
    !! Species-dependent density of neighbours.

    ! Local variables
    integer :: axis
    !! Interface axis result.
    integer :: i,is,ia,ja,l,m,n,ks,cutloc,nstep,itmp1
    !! Loop indices, cutoff location, number of steps, and temporary integer.
    real(real32) :: rdist_max,rcutoff,power,rtmp1
    !! Local max distance, cutoff, decay power, and temporary real.
    real(real32), dimension(3) :: dir_disim
    !! Directional dissimilarity along each axis.
    real(real32), dimension(3) :: vtmp1,vtmp2,vtmp3
    !! Temporary 3-vectors.
    real(real32), allocatable, dimension(:) :: sim_dist,distance
    !! Similarity-distance array and distance grid.
    type(den_of_neigh_type), allocatable, dimension(:,:) :: intf_func
    !! Interface function per axis and species.


!-------------------------------------------------------------------------------
! Initialise variables
!-------------------------------------------------------------------------------
    if(present(lprint))then
       if(lprint) write(*,'(1X,"Determining axis perpendicular to interface")')
    end if
    power=1.E0
    nstep=size(DOS(1)%atom(1,1,:))
    rdist_max=12.0
    if(present(dist_max)) rdist_max=dist_max
    allocate(distance(nstep))
    do i=1,nstep
       distance(i)=real(i,real32)*rdist_max/real(nstep,real32)
    end do
    rcutoff=4.0
    if(present(cutoff)) rcutoff=min(rcutoff,cutoff)
    cutloc=minloc(abs(distance(:)-rcutoff),dim=1)
    allocate(intf_func(3,bas%nspec))
    do i=1,3
       do is=1,bas%nspec
          allocate(intf_func(i,is)%atom(bas%spec(is)%num,2))
       end do
    end do
    allocate(sim_dist(nstep))


!-------------------------------------------------------------------------------
! Interface axis identifier loop
! Loops over each axis
!-------------------------------------------------------------------------------
    do i=1,3
       distloop2: do is=1,bas%nspec
          do ia=1,bas%spec(is)%num
             itmp1=0
             sim_dist=0._real32
             !!-----------------------------------------------------------------
             !! identifies the similarity (scaled by inverse distance) ...
             !! ... between an atom and each other atom of the same species.
             !! This shows how similar an atom is to its local environment
             !!-----------------------------------------------------------------
             do ja=1,bas%spec(is)%num
                vtmp1(:3) = bas%spec(is)%atom(ia,:3) - &
                     bas%spec(is)%atom(ja,:3)
                do l=-1,1,1
                   vtmp2(1) = vtmp1(1) + real(l,real32)
                   do m=-1,1,1
                      vtmp2(2) = vtmp1(2) + real(m,real32)
                      nloop3: do n=-1,1,1
                         vtmp2(3) = vtmp1(3) + real(n,real32)
                         vtmp3 = matmul(vtmp2,lat)
                         !rtmp1=table_func(vtmp3(i),0.8_real32)
                         !rtmp1=exp(-abs(vtmp3(i))*power)
                         rtmp1=exp(-norm2(vtmp3)*power)
                         if(rtmp1.lt.1.D-3) cycle nloop3
                         itmp1=itmp1+1

                         do ks=1,bas%nspec
                            sim_dist = sim_dist + &
                                 sqrt(overlap_indiv_points(&
                                 [DOS(is)%atom(ia,ks,:)],&
                                 [DOS(is)%atom(ja,ks,:)]))*rtmp1
                         end do
                      end do nloop3
                   end do
                end do
             end do
             !!-----------------------------------------------------------------
             !! saves similarity up to the cutoff for each atom and its location
             !!-----------------------------------------------------------------
             intf_func(i,is)%atom(ia,1)=bas%spec(is)%atom(ia,i)*norm2(lat(i,:))
             intf_func(i,is)%atom(ia,2)=sum(sim_dist(:cutloc))!/bas%spec(is)%num!/itmp1


          end do
       end do distloop2
       dir_disim(i)=0._real32
       !!-----------------------------------------------------------------------
       !! finds max difference between points within the cell along a direction
       !!-----------------------------------------------------------------------
       do is=1,bas%nspec
          do ia=1,bas%spec(is)%num
             do ja=ia+1,bas%spec(is)%num
                if( abs( &
                     intf_func(i,is)%atom(ia,1) - &
                     intf_func(i,is)%atom(ja,1)).lt.1._real32)then
                   if( abs( &
                        intf_func(i,is)%atom(ia,2) - &
                        intf_func(i,is)%atom(ja,2)).gt.dir_disim(i) )then
                      dir_disim(i) = &
                           abs(intf_func(i,is)%atom(ia,2) - &
                           intf_func(i,is)%atom(ja,2))
                   end if
                end if

             end do
          end do

       end do
    end do


!-------------------------------------------------------------------------------
! Define the interface axis as the one with the greatest difference
!-------------------------------------------------------------------------------
    axis=minloc(dir_disim,dim=1)
    if(present(lprint))then
       if(lprint) write(*,*) "Interface located along axis",axis
    end if


  end function get_intf_axis_DOS
!###############################################################################


!###############################################################################
  function get_intf_axis_CAD(lat,bas) result(axis)
    !! Determine the axis perpendicular to the interface using the CAD method.
    implicit none

    ! Arguments
    type(basis_type), intent(in) :: bas
    !! Input basis type containing atomic positions.
    real(real32), dimension(3,3), intent(in) :: lat
    !! Lattice vectors.

    ! Local variables
    integer :: i,j,is,iaxis
    !! Loop indices.
    integer :: pntl,pntr,nstep
    !! Left and right point indices and number of steps.
    real(real32) :: sigma,gauss_tol,area
    !! Gaussian sigma, tolerance, and cross-sectional area.
    integer, dimension(3) :: abc
    !! Axis permutation indices.
    real(real32), dimension(3) :: vtmp1,vtmp2,axis_vec
    !! Temporary 3-vectors and per-axis CADD maximum.
    real(real32), allocatable, dimension(:) :: rangevec
    !! Range of derivative per species.
    real(real32), allocatable, dimension(:) :: dist,multiCADD
    !! Distance grid and combined CADD.
    real(real32), allocatable, dimension(:,:) :: CAD,deriv
    !! Cumulative atomic density and its derivative.
    real(real32), allocatable, dimension(:,:,:) :: CADD
    !! Cumulative atomic density derivative coefficients.
    integer :: axis
    !! Interface axis result.



!-------------------------------------------------------------------------------
! Initialise variables
!-------------------------------------------------------------------------------
    nstep=nstep_default
    allocate(dist(nstep))
    dist=0._real32

    sigma=2._real32
    gauss_tol=16._real32
    allocate(rangevec(bas%nspec))
    allocate(deriv(bas%nspec,nstep))
    allocate(CAD(bas%nspec,nstep))
    allocate(CADD(bas%nspec,nstep,3))  !!CADD(spec,nstep,nth order deriv)
    allocate(multiCADD(nstep))
    abc = [1,2,3]

!-------------------------------------------------------------------------------
! Cycle over the 3 axes
!-------------------------------------------------------------------------------
    do iaxis=1,3
       do i=1,nstep
          dist(i)=(i-1)*norm2(lat(iaxis,:))/nstep
       end do
       abc = cshift(abc,1,1)
       area = get_area(lat(abc(1),:),lat(abc(2),:))
       CAD=0._real32
       CADD=0._real32
       !!--------------------------------------------------------------------------
       !! set up CAD and CADD
       !!--------------------------------------------------------------------------
       do is=1,bas%nspec
          !!-----------------------------------------------------------------------
          !! extend cell 1 cell above and below the interface
          !!-----------------------------------------------------------------------
          do j=-1,1,1
             CAD(is,:) = CAD(is,:) + gauss_array(&
                  dist(:),&
                  (bas%spec(is)%atom(:,iaxis)+real(j,real32))*norm2(lat(iaxis,:)),&
                  sigma,gauss_tol,.false.)
          end do
          !!-----------------------------------------------------------------------
          !! generate cumulative atomic density (CAD) from atomic density
          !!-----------------------------------------------------------------------
          do i=nstep,1,-1
             CAD(is,i)=sum(CAD(is,:i))
          end do
          !!-----------------------------------------------------------------------
          !! generate cumulative atomic density derivative (CADD)
          !!-----------------------------------------------------------------------
          do i=1,nstep
             pntl=i-1
             pntr=i+1
             do j=-1,1,1
                vtmp1(j+2)=real(i+j-1,real32)*norm2(lat(iaxis,:))/nstep
             end do
             vtmp2=0._real32
             vtmp2(2)=CAD(is,i)
             if(i.eq.1)then
                vtmp2(1)=0._real32
                vtmp2(3)=CAD(is,pntr)
                !pntl=nstep-1
             elseif(i.eq.nstep)then
                vtmp2(1)=CAD(is,pntl)
                vtmp2(3)=CAD(is,nstep)+CAD(is,1)
                !pntr=1+1
             else
                vtmp2(1)=CAD(is,pntl)
                vtmp2(3)=CAD(is,pntr)
             end if

             CADD(is,i,:)=simeq(vtmp1,vtmp2)
             deriv(is,i)=(CADD(is,i,2)+2.0*CADD(is,i,1)*dist(i))/bas%spec(is)%num
          end do
          rangevec(is)=range(deriv(is,:))
       end do


       !!--------------------------------------------------------------------------
       !! multiply the CADDs of each species into an overal CADD (multiCADD)
       !!--------------------------------------------------------------------------
       multiCADD=1._real32
       do is=1,bas%nspec
          if(rangevec(is).lt.maxval(rangevec)*5.D-2) cycle
          multiCADD(:) = multiCADD(:)*CADD(is,:,1)
       end do

       axis_vec(iaxis) = maxval(multiCADD/area)
    end do


!-------------------------------------------------------------------------------
! Determine the interfacial axis
!-------------------------------------------------------------------------------
    axis = maxloc(axis_vec,dim=1)


    

  end function get_intf_axis_CAD
!###############################################################################


!###############################################################################
  function get_intf_CAD(lat,bas,axis,num_step,lprint) result(intf_loc)
    !! Find interface location using the cumulative atomic density (CAD) method.
    implicit none

    ! Arguments
    integer :: axis
    !! Axis perpendicular to the interface.
    integer, optional, intent(in) :: num_step
    !! Number of steps for density calculation.
    type(basis_type) :: bas
    !! Input basis type containing atomic positions.
    real(real32), dimension(3,3) :: lat
    !! Lattice vectors.
    logical, optional :: lprint
    !! Whether to print debug output.

    ! Local variables
    integer :: i,j,is
    !! Loop indices.
    integer :: pntl,pntr,nstep
    !! Left and right point indices and number of steps.
    real(real32) :: sigma, gauss_tol
    !! Gaussian sigma and tolerance.
    real(real32), dimension(2) :: intf_loc
    !! Interface location result.
    real(real32), dimension(3) :: vtmp1,vtmp2
    !! Temporary 3-vectors.
    integer, allocatable, dimension(:) :: ivec1
    !! Turning point indices.
    real(real32), allocatable, dimension(:) :: rangevec
    !! Range of derivative per species.
    real(real32), allocatable, dimension(:) :: dist,multiCADD
    !! Distance grid and combined CADD.
    real(real32), allocatable, dimension(:,:) :: CAD,deriv
    !! Cumulative atomic density and its derivative.
    real(real32), allocatable, dimension(:,:,:) :: CADD
    !! Cumulative atomic density derivative coefficients.
    real(real32) :: diff
    !! Distance difference between interface candidates.


!-------------------------------------------------------------------------------
! Initialise variables
!-------------------------------------------------------------------------------
    nstep=nstep_default
    if(present(num_step)) nstep=num_step
    allocate(dist(nstep))
    dist=0._real32
    do i=1,nstep
       dist(i)=(i-1)*norm2(lat(axis,:))/nstep
    end do

    sigma=2._real32
    gauss_tol=16._real32
    allocate(rangevec(bas%nspec))
    allocate(deriv(bas%nspec,nstep))
    allocate(CAD(bas%nspec,nstep))
    allocate(CADD(bas%nspec,nstep,3))  !!CADD(spec,nstep,nth order deriv)
    CAD=0._real32
    CADD=0._real32

   
!-------------------------------------------------------------------------------
! Set up CAD and CADD
!-------------------------------------------------------------------------------
    do is=1,bas%nspec
       !!-----------------------------------------------------------------------
       !! extend cell 1 cell above and below the interface
       !!-----------------------------------------------------------------------
       do j=-1,1,1
          CAD(is,:) = CAD(is,:) + gauss_array(&
               dist(:),&
               (bas%spec(is)%atom(:,axis)+real(j,real32))*norm2(lat(axis,:)),&
               sigma,gauss_tol,.false.)
       end do
       !!-----------------------------------------------------------------------
       !! generate cumulative atomic density (CAD) from atomic density
       !!-----------------------------------------------------------------------
       do i=nstep,1,-1
          CAD(is,i)=sum(CAD(is,:i))
       end do
       !!-----------------------------------------------------------------------
       !! generate cumulative atomic density derivative (CADD)
       !!-----------------------------------------------------------------------
       do i=1,nstep
          pntl=i-1
          pntr=i+1
          do j=-1,1,1
             vtmp1(j+2)=real(i+j-1,real32)*norm2(lat(axis,:))/nstep
          end do
          vtmp2=0._real32
          vtmp2(2)=CAD(is,i)
          if(i.eq.1)then
             vtmp2(1)=0._real32
             vtmp2(3)=CAD(is,pntr)
             !pntl=nstep-1
          elseif(i.eq.nstep)then
             vtmp2(1)=CAD(is,pntl)
             vtmp2(3)=CAD(is,nstep)+CAD(is,1)
             !pntr=1+1
          else
             vtmp2(1)=CAD(is,pntl)
             vtmp2(3)=CAD(is,pntr)
          end if

          CADD(is,i,:)=simeq(vtmp1,vtmp2)
          deriv(is,i)=(CADD(is,i,2)+2.0*CADD(is,i,1)*dist(i))/bas%spec(is)%num
       end do
       rangevec(is)=range(deriv(is,:))
    end do
    

!-------------------------------------------------------------------------------
! Multiply the CADDs of each species into an overall CADD (multiCADD)
!-------------------------------------------------------------------------------
    allocate(multiCADD(nstep))
    multiCADD=1._real32
    do is=1,bas%nspec
       if(rangevec(is).lt.maxval(rangevec)*5.D-2) cycle
       multiCADD(:)=multiCADD(:)*CADD(is,:,1)
    end do


!-------------------------------------------------------------------------------
! Identify whether system is likely a planar defect
!-------------------------------------------------------------------------------
    if(count(abs(multiCADD).lt.1.D-8).gt.0.9*nstep)then
       write(*,'(1X,"System has same species-split density across system")')
       write(*,'(1X,"Likely a planar defect")')
       write(*,'(1X,"Use another interface identifier method...")')
    end if


!-------------------------------------------------------------------------------
! Smooth the multiCADD by applying a running average
!-------------------------------------------------------------------------------
    multiCADD=running_avg(multiCADD,window=9,lperiodic=.true.)
    multiCADD=abs(multiCADD)
    if(present(lprint))then
       if(lprint) then
          do i=1,nstep
             write(52,*) dist(i),multiCADD(i)
          end do
       end if
    end if


!-------------------------------------------------------------------------------
! Find the turning points of the multiCADD and attribute them to
! the two interfaces
!-------------------------------------------------------------------------------
    ivec1 = get_turn_points([multiCADD(:)],window=8,lperiodic=.true.)

    intf_loc(1)=dist(ivec1(size(ivec1)))
    do i = size(ivec1) - 1, 1, -1
       diff = abs(intf_loc(1)-dist(ivec1(i)))
       ! map back into the original space if greater than the size of the cell
       if(abs(diff).gt.0.5*norm2(lat(axis,:)))then
          diff = diff - sign(1._real32,diff) * norm2(lat(axis,:))
       end if
       if(abs(diff).gt.2._real32)then
          intf_loc(2)=dist(ivec1(i))
          exit
       end if
    end do


  end function get_intf_CAD
!###############################################################################


!###############################################################################
  function get_layered_axis(lat,bas,lprint) result(axis)
    !! Determine whether a structure is layered and return the layered axis.
    implicit none

    ! Arguments
    type(basis_type), intent(in) :: bas
    !! Input basis type containing atomic positions.
    real(real32), dimension(3,3), intent(in) :: lat
    !! Lattice vectors.
    logical, optional, intent(in) :: lprint
    !! Whether to print progress information.

    ! Local variables
    integer :: i,is,j,nstep,diffcount,axis
    !! Loop indices, number of steps, diff count, and axis result.
    !integer, dimension(3) :: nturns
    real(real32) :: sigma, gauss_tol
    !! Gaussian sigma and tolerance.
    logical :: udef_lprint
    !! Local print flag.
    real(real32), dimension(3) :: diff
    !! Ratio of max to min atomic density per axis.
    real(real32), dimension(3,2) :: minmax
    !! Min and max atomic density per axis.
    real(real32), allocatable, dimension(:) :: AD,dist
    !! Atomic density and distance grid.
    !integer, allocatable, dimension(:) :: ivec1


!-------------------------------------------------------------------------------
! Initialise variables
!-------------------------------------------------------------------------------
    sigma=0.5_real32
    gauss_tol=16._real32
    if(present(lprint))then
       udef_lprint=lprint
    else
       udef_lprint=.false.
    end if


!-------------------------------------------------------------------------------
! Cycle over axes to find atomic density (AD) mapping along each axis
!-------------------------------------------------------------------------------
    axis_loop1: do i=1,3
       if(allocated(dist)) deallocate(dist)
       nstep=nint(norm2(lat(i,:))/0.001_real32)
       allocate(dist(nstep))
       dist=0._real32
       do j=1,nstep
          dist(j)=(j-1)*norm2(lat(i,:))/nstep
       end do
       
       if(allocated(AD)) deallocate(AD)       
       allocate(AD(nstep))       
       AD=0._real32
       do is=1,bas%nspec
          do j=-1,1,1
             AD(:) = AD(:) + gauss_array(&
                  dist(:),&
                  (bas%spec(is)%atom(:,i)+real(j,real32))*norm2(lat(i,:)),&
                  sigma,gauss_tol,.false.)
          end do
       end do
       !ivec1=get_turn_points(dble(AD(:)),window=8,lperiodic=.true.)
       !nturns(i)=size(ivec1,dim=1)
       minmax(i,1)=minval(AD)
       minmax(i,2)=maxval(AD)
       diff(i)=minmax(i,2)/minmax(i,1)
    end do axis_loop1


!-------------------------------------------------------------------------------
! Check each axis
!-------------------------------------------------------------------------------
    axis=0
    select case(count(diff.gt.huge(0._real32)))
    case(1)
       axis=maxloc(diff(:),dim=1)
       if(udef_lprint) write(0,'("Found a 2D system along ",I0)') axis
    case(2)
       axis=minloc(diff(:),dim=1)
       if(udef_lprint) write(0,'("Found a 1D system along ",I0)') axis
       axis=-1
    case(3)
       if(udef_lprint) write(0,*) "Found a 0D system"
       axis=-2
    case default
       axis_loop2: do i=1,3
! ADD A TOLERANCE FOR 'COULD BE LAYERED'
          diffcount=count(diff(i).gt.5._real32*diff(:))
          if(diffcount.eq.2)then
             axis=i
             exit axis_loop2
          end if
       end do axis_loop2
       if(udef_lprint)then
          if(axis.eq.0)then
             write(0,'("System is not layered")')
          else
             write(0,'("Found a 2D system along ",I0)') axis
          end if
       end if
    end select


  end function get_layered_axis
!###############################################################################


!###############################################################################
! Locate two interfaces (currently unused)
!###############################################################################
!  function locate_two_intfs(func,ivec1,lmax) result(intf_loc)
!    implicit none
!    integer :: loc,i
!    real(real32) :: rtmp1
!    logical :: luse_max
!    integer, dimension (:) :: ivec1
!    integer, dimension(2) :: intf_loc
!    real(real32), dimension(:) :: func
!    logical, optional :: lmax
!
!
!    intf_loc=0
!    luse_max=.false.
!    if(present(lmax)) luse_max=lmax
!
!    if(luse_max)then
!       rtmp1=tiny(0.0)
!       loc=maxloc((/ ( func(ivec1(i)),i=1,size(ivec1) ) /),dim=1)
!    else
!       rtmp1=huge(0.0)
!       loc=minloc((/ ( func(ivec1(i)),i=1,size(ivec1) ) /),dim=1)
!    end if
!    intf_loc(1)=ivec1(loc)
!
!
!    do i=1,size(ivec1)
!       if(i.eq.loc) cycle
!       if(luse_max.and.func(ivec1(i)).lt.rtmp1)then
!          rtmp1=func(ivec1(i))
!          intf_loc(2)=ivec1(i)
!       elseif(.not.luse_max.and.func(ivec1(i)).gt.rtmp1)then
!          rtmp1=func(ivec1(i))
!          intf_loc(2)=ivec1(i)
!       end if
!    end do
!
!  end function locate_two_intfs
!###############################################################################


!###############################################################################
  function gen_single_DOS(lat,bas,ispec,iatom,dist_max,weight_dist) result(DOS)
    !! Generate species-dependent density of neighbours for a single atom.
    implicit none

    ! Arguments
    integer, intent(in) :: ispec,iatom
    !! Species and atom indices.
    real(real32), dimension(3,3), intent(in) :: lat
    !! Lattice vectors.
    real(real32), optional, intent(in) :: dist_max
    !! Maximum interatomic distance.
    logical, optional, intent(in) :: weight_dist
    !! Whether to apply distance weighting.

    ! Local variables
    integer :: i,j,k,js,ja,count1
    !! Loop indices and counter.
    integer :: nstep
    !! Number of steps.
    real(real32) :: rdist_max
    !! Local maximum distance.
    real(real32) :: gauss_tol,DON_sigma,dist,dist_cutoff,rtmp1
    !! Gaussian tolerance, sigma, distance, cutoff, and temporary real.
    type(basis_type) :: bas
    !! Input basis type containing atomic positions.
    logical :: lweight
    !! Local distance weighting flag.
    real(real32), dimension(3) :: vtmp1,vtmp2,vtmp3
    !! Temporary 3-vectors.
    real(real32), allocatable, dimension(:) :: distance
    !! Distance grid array.
    real(real32), allocatable, dimension(:,:) :: DOS
    !! Species-dependent density of neighbours result.
    real(real32), allocatable, dimension(:) :: dist_list
    !! List of interatomic distances.


    nstep=nstep_default
    allocate(dist_list(bas%natom*27-1))
    allocate(DOS(bas%nspec,nstep))

    allocate(distance(nstep))
    rdist_max=12._real32
    if(present(dist_max)) rdist_max=dist_max
    do i=1,nstep
       distance(i)=real(i,real32)*rdist_max/real(nstep,real32)
    end do

    gauss_tol=16.E0!38._real32
    DON_sigma=0.5E-1
    dist_cutoff=dist_max+sqrt(2*gauss_tol*DON_sigma**2)

    DOS(:,:)=0._real32
    specloop1: do js=1,bas%nspec
       count1=0
       dist_list = 0.0
       atomloop1: do ja=1,bas%spec(js)%num
          vtmp1(:3) = bas%spec(ispec)%atom(iatom,:3) - bas%spec(js)%atom(ja,:3)
          do i=-1,1,1
             vtmp2(1) = vtmp1(1) + real(i,real32)
             do j=-1,1,1
                vtmp2(2) = vtmp1(2) + real(j,real32)
                kloop1: do k=-1,1,1
                   if(ispec.eq.js.and.iatom.eq.ja)then
                      if(i.eq.0.and.j.eq.0.and.k.eq.0)then
                         cycle kloop1
                      end if
                   end if
                   vtmp2(3) = vtmp1(3) + real(k,real32)
                   vtmp3 = matmul(vtmp2,lat)
                   rtmp1=norm2(vtmp3)
                   if(rtmp1.gt.dist_cutoff) cycle kloop1
                   count1=count1+1
                   dist_list(count1) = norm2(vtmp3)

                end do kloop1
             end do
          end do
       end do atomloop1

       DOS(js,:) = &
            gauss_array(distance,&
            dist_list(:count1),DON_sigma,gauss_tol,.false.)

    end do specloop1


    if(present(weight_dist))then
       lweight=weight_dist
    else
       lweight=.true.
    end if

    if(lweight)then
       do i=minloc(abs(distance(:)-2.E0),dim=1),nstep
          dist=exp(-abs(distance(i)-2.E0))
          DOS(:,i)=DOS(:,i)*dist
       end do
    end if


  end function gen_single_DOS
!###############################################################################


!###############################################################################
  function gen_single_DON(lat,bas,ispec,iatom,dist_max) result(DON)
    !! Generate density of neighbours for a single atom.
    implicit none

    ! Arguments
    integer, intent(in) :: ispec,iatom
    !! Species and atom indices.
    real(real32), optional, intent(in) :: dist_max
    !! Maximum interatomic distance.

    ! Local variables
    integer :: i,nstep
    !! Loop index and number of steps.
    type(basis_type) :: bas
    !! Input basis type containing atomic positions.
    real(real32), dimension(3,3) :: lat
    !! Lattice vectors.
    real(real32), allocatable, dimension(:) :: DON
    !! Density of neighbours result.
    real(real32), allocatable, dimension(:,:) :: DOS
    !! Species-dependent density of neighbours.


    if(present(dist_max))then
       DOS=gen_single_DOS(lat,bas,ispec,iatom,dist_max)
    else
       DOS=gen_single_DOS(lat,bas,ispec,iatom)
    end if

    nstep=size(DOS(1,:))
    allocate(DON(nstep))
       
    do i=1,nstep
       DON(i) = sum(DOS(:,i))
    end do

  end function gen_single_DON
!###############################################################################


end module artemis__interface_identifier
