!!!#############################################################################
!!! Code written by Ned Thaddeus Taylor and Francis Huw Davies
!!! Code part of the ARTEMIS group (Hepplestone research group).
!!! Think Hepplestone, think HRG.
!!!#############################################################################
!!!module contains lattice- and basis-related functions and subroutines.
!!!module includes the following functions and subroutines:
!!! MATNORM          (normalises a 3x3 matrix)
!!! min_dist         (min distance between a point in a cell and nearest atom)
!!! get_surface_normal (return the vector normal to the surface of the plane ...
!!!                   ... constructed by the other two vectors)
!!! get_atom_height  (get the value of the atom along that axis)
!!! get_min_bulk_bond
!!! shifter          (shifts the basis along the cell by an amount)
!!! shift_region     (shifts the basis within a region along the cell by amount)
!!! vacuumer         (adds a vacuum gap to the location specified)
!!! set_vacuum
!!! ortho_axis       (makes specified axis perpendicular to plane of other two)
!!! transformer      (applies a transformation matrix to a lattice and basis)
!!! change_basis     (convert basis into direct coords wrt another lattice)
!!! region_rot       (rotates a region specified along an axis about that axis)
!!! normalise_basis  (convert basis coordinates to be within val-> val-1)
!!! centre_of_geom   (prints centre of geom of a molecule
!!! centre_of_mass   (prints centre of mass of a molecule)
!!! primitive_lat    (reorientates the lattice to the primitive lattice)
!!! reducer
!!! mkNiggli_lat
!!! reduced_check
!!! planecutter      (generates transformation mat to obtain miller plane)
!!! bas_merge        (merges two supplied bases)
!!! bas_lat_merge    (merges two supplied bases and lattices)
!!! split_bas
!!! get_bulk
!!! get_centre_atom
!!! get_wyckoff      (returns an array of the similar atoms)
!!! get_shortest_bond
!!!#############################################################################
module artemis__geom_utils
  use artemis__constants, only: real32, pi
  use artemis__geom_rw, only: basis_type,geom_write
  use artemis__sym, only: confine_type, gldfnd, tol_sym_default
  use artemis__misc, only: swap, sort2D
  use misc_linalg, only: cross,outer_product,cross_matrix,uvec,modu,&
       get_vol,det,inverse,inverse_3x3,LUinv,reduce_vec_gcd,get_vec_multiple,&
       proj,GramSchmidt,LLL_reduce
  implicit none

  type wyck_atom_type
     integer, allocatable, dimension(:) :: atom
  end type wyck_atom_type
  type wyck_spec_type
     type(wyck_atom_type), allocatable, dimension(:) :: spec
  end type wyck_spec_type
  type bond_type
     real(real32) :: length
     integer, dimension(2,2) :: atoms
  end type bond_type

  
  interface get_closest_atom
     procedure get_closest_atom_1D,get_closest_atom_3D
  end interface get_closest_atom



contains

!###############################################################################
  function compare_stoichiometry(basis1, basis2) result(output)
    !! Check if two basis structures have the same stoichiometry ratio
    !!
    !! This function compares the stoichiometry ratios of two basis structures
    !! It returns true if the relative proportions of all atomic species are 
    !! identical and all species names match between both structures
    implicit none
    type(basis_type), intent(in) :: basis1, basis2
    logical :: output
   
    integer :: is, js, total_atoms1, total_atoms2
    real(real32) :: ratio1, ratio2, tol
    logical :: found_match
   
    ! Set tolerance for floating-point comparisons
    tol = 1.E-5_real32
   
    ! Initialize output to true, will set to false if any condition fails
    output = .true.
   
    ! Check if both basis have the same number of species
    if (basis1%nspec /= basis2%nspec) then
       output = .false.
       return
    end if
   
    ! Get total number of atoms in each basis
    total_atoms1 = sum(basis1%spec(:)%num)
    total_atoms2 = sum(basis2%spec(:)%num)
   
    ! Compare each species in basis1 with corresponding species in basis2
    do is = 1, basis1%nspec
       found_match = .false.
      
       ! Find matching species in basis2
       do js = 1, basis2%nspec
          ! Check if species names match
          if (basis1%spec(is)%name == basis2%spec(js)%name) then
             found_match = .true.
            
             ! Calculate and compare stoichiometry ratios
             ratio1 = real(basis1%spec(is)%num, real32) / real(total_atoms1, real32)
             ratio2 = real(basis2%spec(js)%num, real32) / real(total_atoms2, real32)
            
             ! Check if ratios are equal within tolerance
             if (abs(ratio1 - ratio2) .gt. tol) then
                output = .false.
                return
             end if
             
             exit  ! Found matching species, continue to next species in basis1
          end if
       end do
      
       ! If no matching species found in basis2, stoichiometry can't be the same
       if (.not. found_match) then
          output = .false.
          return
       end if
    end do

  end function compare_stoichiometry
!###############################################################################


!!!#############################################################################
!!! Normalises a 3x3 matrix to the form:
!!! a 0 0
!!! b c 0
!!! d e f
!!! NEW NAME: lat_low
!!!#############################################################################
  function MATNORM(lat) result(nlat)
    implicit none
    real(real32), dimension(3,3) :: lat, nlat
    nlat(1,1)=sqrt(lat(1,1)**2+lat(1,2)**2+lat(1,3)**2)
    nlat(1,2)=0.0
    nlat(1,3)=0.0

    nlat(2,1)=(lat(1,1)*lat(2,1)+lat(1,2)*lat(2,2)+lat(1,3)*lat(2,3))/nlat(1,1)
    nlat(2,2)=&
         sqrt((lat(1,2)*lat(2,3)-lat(1,3)*lat(2,2))**2+&
         (lat(1,3)*lat(2,1)-lat(1,1)*lat(2,3))**2+&
         (lat(1,1)*lat(2,2)-lat(1,2)*lat(2,1))**2)/nlat(1,1)
    nlat(2,3)=0.0

    nlat(3,1)=(lat(1,1)*lat(3,1)+lat(1,2)*lat(3,2)+lat(1,3)*lat(3,3))/nlat(1,1)
    nlat(3,2)=(&
         lat(2,1)*lat(3,1)+&
         lat(2,2)*lat(3,2)+&
         lat(2,3)*lat(3,3)-&
         nlat(2,1)*nlat(3,1))/nlat(2,2)
    nlat(3,3)=sqrt(&
         lat(3,1)**2+lat(3,2)**2+&
         lat(3,3)**2-nlat(3,1)**2-nlat(3,2)**2)
  end function MATNORM
!!!#############################################################################


!!!#############################################################################
!!! Finds distance between a location in a cell and ...
!!! ... the nearest atom to that point either above ...
!!! ... or below
!!!#############################################################################
  function min_dist(bas,axis,loc,above)
    implicit none
    integer :: is,axis
    real(real32) :: min_dist,pos
    real(real32), intent(in) :: loc
    type(basis_type) :: bas
    logical :: labove
    logical,optional :: above


    pos=loc
    labove=.false.
    if(present(above)) labove=above
    aboveloop: if(labove)then
       min_dist=huge(0._real32)
       if(all( (/ (bas%spec(is)%atom(:,axis),is=1,bas%nspec) /).lt.pos))&
            pos=pos-1._real32
    else
       min_dist=-huge(0._real32)
       if(all( (/ (bas%spec(is)%atom(:,axis),is=1,bas%nspec) /).gt.pos))&
            pos=pos-1._real32
    end if aboveloop


    do is=1,bas%nspec
       if(.not.labove.and.maxval(bas%spec(is)%atom(:,axis)-pos,&
            mask=(bas%spec(is)%atom(:,axis)-pos.le.0._real32)).gt.min_dist) then
          min_dist=maxval(bas%spec(is)%atom(:,axis)-pos,&
               mask=(bas%spec(is)%atom(:,axis)-pos.le.0._real32))
       elseif(labove.and.minval(bas%spec(is)%atom(:,axis)-pos,&
            mask=(bas%spec(is)%atom(:,axis)-pos.ge.0._real32)).lt.min_dist) then
          min_dist=minval(bas%spec(is)%atom(:,axis)-pos,&
               mask=(bas%spec(is)%atom(:,axis)-pos.ge.0._real32))
       end if
    end do

  end function min_dist
!!!#############################################################################


!!!#############################################################################
!!! Get the value of the atom along that axis
!!!#############################################################################
  function get_atom_height(bas,atom,axis) result(val)
    implicit none
    integer :: i,axis,atom,sum_atom
    real(real32) :: val
    type(basis_type) :: bas

    val=0._real32
    sum_atom=0
    do i=1,bas%nspec
       if(atom.le.sum_atom+bas%spec(i)%num)then
          val=bas%spec(i)%atom(atom-sum_atom,axis)
          return
       end if
       sum_atom=sum_atom+bas%spec(i)%num
    end do


  end function get_atom_height
!!!#############################################################################


!!!#############################################################################
!!! returns minimum bond within bulk
!!!#############################################################################
  function get_min_bulk_bond(basis) result(min_bond)
    implicit none
    type(basis_type), intent(in) :: basis

    integer :: is,ia,js,ja
    real(real32) :: dtmp1,min_bond
    real(real32), dimension(3) :: vdtmp1


    min_bond=huge(0._real32)
    if(basis%natom.eq.1)then
       min_bond = min( &
            modu(basis%lat(1,:3)), &
            modu(basis%lat(2,:3)), &
            modu(basis%lat(3,:3)) &
       )
       return
    end if

    do is = 1, basis%nspec
       do ia = 1, basis%spec(is)%num

          do js=1,basis%nspec
             atmloop: do ja=1,basis%spec(js)%num
                if(is.eq.js.and.ia.eq.ja) cycle atmloop
                vdtmp1 = basis%spec(js)%atom(ja,:3) - basis%spec(is)%atom(ia,:3)
                vdtmp1 = &
                     vdtmp1(1)*basis%lat(1,:3) + &
                     vdtmp1(2)*basis%lat(2,:3) + &
                     vdtmp1(3)*basis%lat(3,:3)
                dtmp1 = modu(vdtmp1)
                if(dtmp1.lt.min_bond) min_bond = dtmp1
             end do atmloop
          end do

       end do
    end do

  end function get_min_bulk_bond
!!!#############################################################################


!!!#############################################################################
!!! returns minimum bond for a specified atom
!!!#############################################################################
  function get_min_bond(basis,is,ia,axis,labove,tol) result(vsave)
    implicit none
    integer :: js,ja
    integer :: iaxis
    real(real32) :: dtmp1,min_bond,dtol
    logical :: ludef_above
    real(real32), dimension(3) :: vdtmp1, vsave

    integer, intent(in) :: is,ia
    type(basis_type), intent(in) :: basis

    integer, intent(in), optional :: axis
    real(real32), intent(in), optional :: tol
    logical, intent(in), optional :: labove

    if(present(tol))then
       dtol = tol
    else
       dtol = 1.E-5_real32
    end if

    if(present(labove))then
       ludef_above=labove
    else
       ludef_above=.false.
    end if

    if(present(axis))then
       iaxis=axis
    else
       iaxis=0
    end if

    min_bond=huge(0._real32)
    
    do js=1,basis%nspec
       atmloop: do ja=1,basis%spec(js)%num
          if(is.eq.js.and.ia.eq.ja) cycle atmloop
          vdtmp1 = basis%spec(js)%atom(ja,:3) - basis%spec(is)%atom(ia,:3)
          if(iaxis.gt.0)then
             if(abs(vdtmp1(iaxis)).lt.dtol) cycle atmloop
             if(ludef_above)then
                vdtmp1(iaxis) = 1._real32 + vdtmp1(iaxis)
             else
                vdtmp1(iaxis) = vdtmp1(iaxis) - 1._real32
             end if
          end if
          vdtmp1 = &
               vdtmp1(1)*basis%lat(1,:3) + &
               vdtmp1(2)*basis%lat(2,:3) + &
               vdtmp1(3)*basis%lat(3,:3)
          dtmp1 = modu(vdtmp1)
          if(dtmp1.lt.min_bond)then
             min_bond = dtmp1
             vsave = vdtmp1
          end if
       end do atmloop
    end do


  end function get_min_bond
!!!#############################################################################


!!!#############################################################################
!!! returns minimum bond for a specified atom
!!!#############################################################################
  function get_min_dist(lat,bas,loc,lignore_close,axis,labove,lreal,tol) &
       result(vsave)
    implicit none
    integer :: js,ja
    integer :: iaxis
    real(real32) :: dtmp1,min_bond,dtol
    logical :: ludef_above,ludef_real
    real(real32), dimension(3) :: vdtmp1,vdtmp2,vsave

    logical, intent(in) :: lignore_close
    type(basis_type), intent(in) :: bas
    real(real32), dimension(3), intent(in) :: loc
    real(real32), dimension(3,3), intent(in) :: lat

    integer, intent(in), optional :: axis
    real(real32), intent(in), optional :: tol
    logical, intent(in), optional :: labove, lreal

    !! CORRECT tol TO ACCOUNT FOR LATTICE SIZE
    if(present(tol))then
       dtol = tol
    else
       dtol = 1.E-5_real32
    end if

    if(present(labove))then
       ludef_above=labove
    else
       ludef_above=.false.
    end if

    if(present(lreal))then
       ludef_real=lreal
    else
       ludef_real=.true.
    end if

    if(present(axis))then
       iaxis=axis
    else
       iaxis=0
    end if

    min_bond=huge(0._real32)
    vsave = 0._real32
    do js=1,bas%nspec
       atmloop: do ja=1,bas%spec(js)%num
          vdtmp1 = bas%spec(js)%atom(ja,:3) - loc
          if(lignore_close.and.modu(vdtmp1).lt.dtol) cycle atmloop
          if(iaxis.gt.0)then
             if(abs(vdtmp1(iaxis)).lt.dtol) cycle atmloop
             if(ludef_above)then
                vdtmp1(iaxis) = 1._real32 + vdtmp1(iaxis)
             else
                vdtmp1(iaxis) = vdtmp1(iaxis) - 1._real32
             end if
          end if
          vdtmp2 = &
               vdtmp1(1)*lat(1,:3) + &
               vdtmp1(2)*lat(2,:3) + &
               vdtmp1(3)*lat(3,:3)
          dtmp1 = modu(vdtmp2)
          if(dtmp1.lt.min_bond)then
             min_bond = dtmp1
             if(ludef_real)then
                vsave = vdtmp1
             else
                vsave = vdtmp2
             end if
          end if
       end do atmloop
    end do


  end function get_min_dist
!!!#############################################################################


!!!#############################################################################
!!! Shifts the basis along a, b or c by amount 'shift'
!!!#############################################################################
  subroutine shifter(basis,axis,shift,renormalise)
    implicit none
    type(basis_type), intent(inout) :: basis
    integer, intent(in) :: axis
    real(real32), intent(in) :: shift
    logical, optional, intent(in) ::renormalise
    integer :: i,j
    logical :: renormalise_

    renormalise_=.false.
    if(present(renormalise)) renormalise_ = renormalise

    do i=1,basis%nspec
       do j=1,basis%spec(i)%num
         basis%spec(i)%atom(j,axis) = basis%spec(i)%atom(j,axis) + shift
          if(renormalise_) basis%spec(i)%atom(j,axis) = basis%spec(i)%atom(j,axis) - &
               floor(basis%spec(i)%atom(j,axis))
       end do
    end do

  end subroutine shifter
!!!#############################################################################


!!!#############################################################################
!!! Shifts basis in a region by an amount
!!!#############################################################################
  subroutine shift_region(bas,region_axis,region_lw,region_up,shift_axis,shift,renorm)
    implicit none
    integer :: is,ia,shift_axis,region_axis
    real(real32) :: shift,region_lw,region_up
    type(basis_type) :: bas
    logical, optional ::renorm
    logical :: lrenorm

    lrenorm=.false.
    if(present(renorm)) lrenorm=renorm

    do is=1,bas%nspec
       do ia=1,bas%spec(is)%num
          if(bas%spec(is)%atom(ia,region_axis).ge.region_lw.and.&
               bas%spec(is)%atom(ia,region_axis).le.region_up)then
             bas%spec(is)%atom(ia,shift_axis) = &
                  bas%spec(is)%atom(ia,shift_axis) + shift
             if(lrenorm) bas%spec(is)%atom(ia,shift_axis) = &
                  bas%spec(is)%atom(ia,shift_axis) - &
                  floor(bas%spec(is)%atom(ia,shift_axis))
          end if
       end do
    end do

  end subroutine shift_region
!!!#############################################################################


!!!#############################################################################
!!! Return the surface normal vector
!!!#############################################################################
  function get_surface_normal(lat,axis) result(normal)
    implicit none
    real(real32) :: component
    integer, dimension(3) :: order=(/1,2,3/)
    real(real32), dimension(3) :: normal

    integer, intent(in) :: axis
    real(real32), dimension(3,3), intent(in) :: lat

    order = cshift(order,3-axis)
    normal = cross([lat(order(1),:)],[lat(order(2),:)])
    component = dot_product(lat(3,:),normal) / modu(normal)**2._real32
    normal = normal * component

    return
  end function get_surface_normal
!!!#############################################################################


!!!#############################################################################
!!! Adjusts the amount of vacuum at a location ...
!!! ... within a cell and adjusts the basis accordingly
!!!#############################################################################
  subroutine vacuumer(lat,bas,axis,loc,add,tol)
    implicit none
    integer :: is,ia
    real(real32) :: rtol,rloc,ortho_scale
    real(real32) :: cur_vac,inc,diff,mag_old,mag_new
    real(real32),dimension(3) :: normal

    integer, intent(in) :: axis
    real(real32), intent(in) :: add,loc
    type(basis_type), intent(inout) :: bas
    real(real32),dimension(3,3), intent(inout) :: lat

    real(real32), optional, intent(in) :: tol


    !! get surface normal vector
    normal = get_surface_normal(lat,axis)
    ortho_scale = modu(lat(axis,:))/modu(normal)


    rtol = 1.E-5_real32
    inc = add
    if(present(tol)) rtol = tol
    cur_vac = min_dist(bas,axis,loc,.true.) - min_dist(bas,axis,loc,.false.)
    cur_vac = cur_vac * modu(lat(axis,:))
    diff = cur_vac + inc
    if(diff.lt.0._real32)then
       write(0,*) "WARNING! Removing vacuum entirely"
    end if

    mag_old = modu(lat(axis,:))
    mag_new = ( mag_old + inc ) / mag_old
    lat(axis,:) = lat(axis,:) * mag_new
    inc = inc / modu(lat(axis,:)) * ortho_scale
    rtol = rtol / mag_old
    rloc = loc / mag_new + rtol


    do is=1,bas%nspec
       do ia=1,bas%spec(is)%num
          bas%spec(is)%atom(ia,axis) = bas%spec(is)%atom(ia,axis) / mag_new
          if(bas%spec(is)%atom(ia,axis).gt.rloc) then
             bas%spec(is)%atom(ia,axis) = bas%spec(is)%atom(ia,axis) + inc
          end if
       end do
    end do


  end subroutine vacuumer
!!!#############################################################################


!!!#############################################################################
!!! Adjusts the amount of vacuum at a location ...
!!! ... within a cell and adjusts the basis accordingly
!!!#############################################################################
  subroutine set_vacuum(basis,axis,loc,vac,tol)
    implicit none
    integer :: is,ia
    real(real32) :: rtol,rloc,ortho_scale
    real(real32) :: cur_vac,diff,mag_old,mag_new
    real(real32),dimension(3) :: normal

    integer, intent(in) :: axis
    real(real32), intent(in) :: vac,loc
    type(basis_type), intent(inout) :: basis

    real(real32), optional, intent(in) :: tol


    !! get surface normal vector
    normal = get_surface_normal(basis%lat,axis)
    ortho_scale = modu(basis%lat(axis,:))/modu(normal)


    rtol = 0._real32
    if(present(tol)) rtol = tol
    if(vac.lt.0._real32)then
       write(0,*) "WARNING! Removing vacuum entirely"
    end if
    cur_vac = min_dist(basis,axis,loc,.true.) - min_dist(basis,axis,loc,.false.)
    cur_vac = cur_vac * modu(normal)
    diff = ( vac - cur_vac ) * ortho_scale

    mag_old = modu(basis%lat(axis,:))
    mag_new = ( mag_old + diff ) / mag_old
    basis%lat(axis,:) = basis%lat(axis,:) * mag_new
    diff = diff / modu(basis%lat(axis,:))
    rtol = rtol / mag_old
    rloc = loc / mag_new + rtol



    do is=1,basis%nspec
       do ia=1,basis%spec(is)%num
          basis%spec(is)%atom(ia,axis) = basis%spec(is)%atom(ia,axis) / mag_new
          if(basis%spec(is)%atom(ia,axis).gt.rloc) then
             basis%spec(is)%atom(ia,axis) = basis%spec(is)%atom(ia,axis) + diff
          end if
       end do
    end do


  end subroutine set_vacuum
!!!#############################################################################


!!!#############################################################################
!!! Takes a lattice and makes the defined axis orthogonal to the other two
!!! WARNING! THIS IS FOR SLAB STRUCTURES! IT REMOVES PERIODICITY ALONG THAT AXIS
!!!#############################################################################
  subroutine ortho_axis(basis,axis)
    implicit none
    type(basis_type), intent(inout) :: basis
    integer :: axis
    real(real32) :: ortho_comp
    integer, dimension(3) :: order
    real(real32), dimension(3) :: ortho_vec
    real(real32), dimension(3,3) :: lat


    order = [ 1, 2, 3 ]
    order = cshift( order, 3 - axis )
    lat = basis%lat

    ortho_vec=cross( [ lat(order(1),:) ] , [ lat(order(2),:) ] )
    ortho_comp=dot_product([ lat(3,:) ],ortho_vec)/modu(ortho_vec)**2._real32
    ortho_vec=ortho_vec*ortho_comp

    lat(3,:)=ortho_vec
    call basis%change_lattice(lat)

    return
  end subroutine ortho_axis
!!!#############################################################################


!!!#############################################################################
!!! Applies a transformation matrix to a lattice ...
!!! ... and extends the basis where needed
!!!#############################################################################
  subroutine transformer(basis, tfmat, map)
    implicit none
    integer :: i,j,k,l,m,n,is,ia
    integer :: satom,dim
    real(real32) :: tol,vol_inc
    logical :: lmap
    type(basis_type), intent(inout) :: basis
    type(basis_type) :: sbas
    integer, dimension(3) :: latmin,latmax
    real(real32), dimension(3):: translvec,tolvec
    integer, allocatable, dimension(:) :: tmp_map_atom
    integer, allocatable, dimension(:,:,:) :: new_map
    real(real32), allocatable, dimension(:,:) :: tmpbas
    real(real32), dimension(3,3) :: tfmat,invmat

    integer, allocatable, dimension(:,:,:), optional, intent(inout) :: map

    vol_inc = abs(det(basis%lat))
    if(vol_inc.lt.0.5_real32)then
       write(0,'(1X,"ERROR: Internal error in transformer function")')
       write(0,'(2X,"transformer in mod_geom_utils.f90 been supplied a&
            & lattice with almost zero determinant")')
       write(0,'(2X,"determinant = ",F0.9)') vol_inc
       write(0,'(3(1X,F7.2))') basis%lat
       stop
    end if
    call basis%normalise(ceil_val = 1._real32, floor_coords = .true., round_coords = .false.)
    vol_inc=abs(det(tfmat))
    sbas%lat=matmul(tfmat,basis%lat)
    invmat=inverse_3x3(tfmat)
    translvec=0._real32
    dim=size(basis%spec(1)%atom(1,:))
    
    
    !!--------------------------------------------------------------------------
    !! If map is present, sets up new map
    !!--------------------------------------------------------------------------
    lmap = .false.
    if_map: if(present(map))then
       if(all(map.eq.-1))then
          exit if_map
       end if
       lmap = .true.
       allocate(new_map(&
            basis%nspec,&
            ceiling(vol_inc)*maxval(basis%spec(:)%num,dim=1),2))
       new_map=0
       if(all(map.eq.0))then
          do is=1,basis%nspec
             map(is,:basis%spec(is)%num,1) = is
             do ia=1,basis%spec(is)%num
                map(is,ia,2) = ia
             end do
          end do
       end if
    end if if_map
    
    
    !!--------------------------------------------------------------------------
    !! Convert tolerance from Å to a fraction of each direction
    !!--------------------------------------------------------------------------
    tol = 1.E-3_real32 !! in Å
    do i=1,3
       tolvec(i)=tol/modu(sbas%lat(i,:))
    end do
    if(vol_inc.lt.minval(tolvec))then
       write(0,'(1X,"ERROR: Internal error in transformer function")')
       write(0,'(2X,"transformer in mod_geom_utils.f90 been supplied a&
            & transformation matrix with almost zero determinant")')
       write(0,'(2X,"determinant = ",F0.9)') vol_inc
       write(0,'(3(1X,F7.2))') tfmat
       stop
    end if
    
    
    !!--------------------------------------------------------------------------
    !! extends basis from min to sum(:) for each column
    !!--------------------------------------------------------------------------
    !! latmin is how far down you have to move to get it back to 0.
    !! hence, latmin=-nint(max*vol)
    !! latmax is how far you still have left to go.
    !! hence, latmax=vol-nint(max*vol)
    !! 
    !! This is why the sum works.
    !! Distance above origin using a is the sum of all positive a's in tfmat
    !! Distance below origin using a is the sum of all negative a's in tfmat
    !!
    !!     ____b
    !!    /   /
    !!  a/   /
    !!  /___/
    !! o 
    !!   
    !!    /\
    !!  a/  \b
    !!  /   /
    !! o\  /
    !!   \/
    !!----------------------------------
    !latmin(i)=(minval(invmat(i,:))-ceiling(minval(invmat(i,:))))*vol
    !latmax(i)=(maxval(invmat(i,:))-floor(minval(invmat(i,:))))*vol
    !latmin(i)=(min(minval(invmat(i,:)),0._real32)-ceiling(minval(invmat(i,:))))*vol
    !latmax(i)=(ceiling(maxval(invmat(i,:3)))-maxval(invmat(i,:3)) )*vol
    do i=1,3
       latmin(i)=floor(sum(tfmat(:3,i),mask=tfmat(:3,i).lt.0._real32))-1
       latmax(i)=ceiling(sum(tfmat(:3,i),mask=tfmat(:3,i).gt.0._real32))+1
    end do
    
    
    !!--------------------------------------------------------------------------
    !! transform the basis
    !!--------------------------------------------------------------------------
    do i=1,basis%nspec
       do j=1,basis%spec(i)%num
          basis%spec(i)%atom(j,:3)=matmul(basis%spec(i)%atom(j,:3),invmat)
       end do
    end do
    
    
    !!--------------------------------------------------------------------------
    !! generates atoms to fill the supercell
    !!--------------------------------------------------------------------------
    allocate(sbas%spec(basis%nspec))
    sbas%sysname = basis%sysname
    sbas%nspec = 0
    sbas%natom = 0
    spec_loop1: do is = 1, basis%nspec
       if(allocated(tmpbas)) deallocate(tmpbas)
       allocate(tmpbas(basis%spec(is)%num*(&
            (abs(latmax(3))+abs(latmin(3))+1)*&
            (abs(latmax(2))+abs(latmin(2))+1)*&
            (abs(latmax(1))+abs(latmin(1))+1)),3))
       satom=0
       if(lmap)then
          allocate(tmp_map_atom(ceiling(vol_inc)*basis%spec(is)%num))
       end if
       do ia = 1, basis%spec(is)%num
          do n=latmin(3),latmax(3)!,1
             translvec(3)=real(n, real32)
             do m=latmin(2),latmax(2)!,1
                translvec(2)=real(m, real32)
                inloop: do l=latmin(1),latmax(1)!,1
                   translvec(1)=real(l, real32)
                   tmpbas(satom+1,:3) = &
                        basis%spec(is)%atom(ia,:3) + matmul(translvec,invmat)
                   !!tmpbas(satom+1,:3)=&
                   !!     matmul((basis%spec(is)%atom(ia,:3)+translvec),invmat)
                   !where(abs(tmpbas(satom+1,:3)-nint(tmpbas(satom+1,k))).lt.tol)
                   !   tmpbas(satom+1,:3)=nint(tmpbas(satom+1,:3))
                   !end where
                   !if(any(tmpbas(satom+1,:).ge.1._real32).or.&
                   !     any(tmpbas(satom+1,:).lt.0._real32)) cycle
                   !if(any(tmpbas(satom+1,:).ge.1._real32+tol).or.&
                   !     any(tmpbas(satom+1,:).lt.0._real32-tol)) cycle
                   if(any(tmpbas(satom+1,:).ge.1._real32-tol).or.&
                        any(tmpbas(satom+1,:).lt.0._real32-tol)) cycle inloop !??? cycle inloop or spec_loop1?
                   tmpbas(satom+1,:3) = tmpbas(satom+1,:3) - &
                        real(floor(tmpbas(satom+1,:3)),real32)
                   do k=1,satom
                      if(all(mod(abs(tmpbas(satom+1,:3)-tmpbas(k,:3)),1._real32).le.&
                           tol)) cycle inloop
                   end do
                   if(lmap) tmp_map_atom(satom+1)=map(is,ia,2)
                   satom=satom+1
                end do inloop
             end do
          end do
       end do
       if(satom.eq.0)then
          if(lmap) deallocate(tmp_map_atom)
          cycle spec_loop1
       end if
       sbas%nspec=sbas%nspec+1
       sbas%spec(sbas%nspec)%num=satom
       sbas%natom=sbas%natom+satom
       sbas%spec(sbas%nspec)%name=basis%spec(is)%name
       allocate(sbas%spec(sbas%nspec)%atom(satom,dim))
       sbas%spec(sbas%nspec)%atom(1:satom,:3)=tmpbas(1:satom,:3)
       if(dim.eq.4) sbas%spec(sbas%nspec)%atom(1:satom,4)=1._real32
       deallocate(tmpbas)
       deallocate(basis%spec(is)%atom)
       if(lmap)then
          new_map(sbas%nspec,:satom,1) = is
          new_map(sbas%nspec,:satom,2) = tmp_map_atom(:satom)
          deallocate(tmp_map_atom)
       end if
    end do spec_loop1
    
    
    !!--------------------------------------------------------------------------
    !! check to see if successfully generated correct number of atoms
    !!--------------------------------------------------------------------------
    if(all(abs(tfmat-nint(tfmat)).lt.tol))then
       if(nint(basis%natom*vol_inc).ne.sbas%natom)then
          write(0,'(1X,"ERROR: Internal error in transformer function")')
          write(0,'(2X,"Transformer in mod_geom_utils.f90 has failed to &
               &generate enough atoms when extending the cell")')
          write(0,'(2X,"Generated ",I0," atoms, whilst expecting ",I0," atoms")') &
               sbas%natom,nint(basis%natom*vol_inc)
          write(0,*) basis%natom,nint(vol_inc)
          write(0,'(3(1X,F7.2))') tfmat
          open(60,file="broken_cell.vasp")
          call geom_write(60,sbas)
          close(60)
          stop
       end if
    end if
    
    
    !!--------------------------------------------------------------------------
    !! saves new lattice and basis to original set
    !!--------------------------------------------------------------------------
    basis%lat = sbas%lat
    deallocate(basis%spec)
    allocate(basis%spec(sbas%nspec))
    basis%sysname=sbas%sysname
    basis%nspec=sbas%nspec
    basis%natom=sbas%natom
    do i=1,sbas%nspec
       allocate(basis%spec(i)%atom(sbas%spec(i)%num,dim))
       basis%spec(i)=sbas%spec(i)
    end do
    
    
    !!--------------------------------------------------------------------------
    !! sets up the new map, if map supplied
    !!--------------------------------------------------------------------------
    if(lmap)then
       deallocate(map)
       call move_alloc(new_map,map)
    end if



  end subroutine transformer
!!!#############################################################################


!!!#############################################################################
!!! Convert basis from direct coords in one lattice ...
!!! ... into direct coords wrt another lattice
!!!#############################################################################
  function change_basis(vec,old_lat,new_lat)
    implicit none
    real(real32), dimension(3) :: change_basis,vec
    real(real32), dimension(3,3), intent(in) :: old_lat,new_lat
    real(real32), dimension(3,3) :: inew_lat
    inew_lat=inverse_3x3(new_lat)
    change_basis=matmul(transpose(inew_lat),matmul(old_lat,vec))
  end function change_basis
!!!#############################################################################


!!!#############################################################################
!!! rotates a region along an axis about that axis
!!!#############################################################################
  subroutine region_rot(bas,lat,angle,axis,bound1,bound2,tvec)
    implicit none
    integer :: axis,i,j
    real(real32) :: angle,bound1,bound2
    real(real32), dimension(3) :: u,centre
    real(real32), dimension(3,3) :: rotmat,ident,lat,invlat
    type(basis_type) :: bas
    real(real32), optional, dimension(3) :: tvec

    centre=(/0.5,0.5,0.0/)
    if(present(tvec)) centre=tvec
    ident=0._real32
    do i=1,3
       ident(i,i)=1._real32
    end do

!!! DEFINE ROTMAT BEFORE THIS
    u=0._real32
    u(axis)=-1._real32
    rotmat=&
         (cos(angle)*ident)+&
         (sin(angle))*cross_matrix(u)+&
         (1-cos(angle))*outer_product(u,u)


!!! Transform the rotation matrix into direct space
    invlat=LUinv(lat)
    rotmat=matmul(lat,rotmat)
    rotmat=matmul(rotmat,invlat)


!!! Rotate the basis within the bounds
    do i=1,bas%nspec
       do j=1,bas%spec(i)%num
          if(bas%spec(i)%atom(j,axis).lt.bound1.or.&
               bas%spec(i)%atom(j,axis).gt.bound2) cycle
          bas%spec(i)%atom(j,:3)=&
               matmul(rotmat,bas%spec(i)%atom(j,:3)-centre)+centre
       end do
    end do


    return
  end subroutine region_rot
!!!#############################################################################


!!!#############################################################################
!!! finds the centre of geometry of the supplied basis
!!!#############################################################################
  function centre_of_geom(bas) result(centre)
    implicit none
    integer :: is,ia,j
    real(real32), dimension(3) :: centre
    type(basis_type) :: bas

    centre=0._real32
    do is=1,bas%nspec
       do ia=1,bas%spec(is)%num
          do j=1,3
             centre(j)=centre(j)+bas%spec(is)%atom(ia,j)
          end do
       end do
    end do

    centre=centre/bas%natom

    return
  end function centre_of_geom
!!!#############################################################################


!!!#############################################################################
!!! finds the centre of mass of the supplied basis
!!!#############################################################################
  function centre_of_mass(bas) result(centre)
    implicit none
    integer :: is,ia,j
    real(real32) :: tot_mass
    real(real32), dimension(3) :: centre
    type(basis_type) :: bas

    centre=0._real32
    tot_mass=0._real32
    do is=1,bas%nspec
       tot_mass=tot_mass+bas%spec(is)%mass*bas%spec(is)%num
       do ia=1,bas%spec(is)%num
          do j=1,3
             centre(j)=centre(j)+bas%spec(is)%atom(ia,j)*bas%spec(is)%mass
          end do
       end do
    end do

    centre=centre/tot_mass

    return
  end function centre_of_mass
!!!#############################################################################


!###############################################################################
  subroutine primitive_lat(basis)
    !! Reorientate lattice to the primitive lattice of its type
    !!
    !! NEED TO SET UP TO WORK FOR THE EXTRA SWAPPINGS OF A, B AND C
    implicit none

    ! Arguments
    type(basis_type), intent(inout) :: basis
    !! Structure data

    ! Local variables
    integer :: i,j
    !! Loop indices
    real(real32) :: rtmp1
    !! Temporary variable
    real(real32), dimension(3) :: scal
    !! Scaling factors
    real(real32), dimension(3,3) :: lat, plat
    !! Lattice matrices
    real(real32), dimension(3,3) :: tmat1, tmat2
    !! Temporary matrices
    real(real32), dimension(3,3,4) :: special
    !! Special lattice matrices


    !!---------------------------------------------------------------
    !! makes all lattice vectors unity
    !!---------------------------------------------------------------
    call reducer(basis)
    lat  = basis%lat
    plat = lat
    do i = 1, 3
       scal(i) = modu(lat(i,:))
       lat(i,:) = lat(i,:) / scal(i)
    end do


    !!---------------------------------------------------------------
    !! sets up the special set of primitive lattices
    !!---------------------------------------------------------------
    special(:,:,1) = transpose( reshape( (/&
         1._real32, 0._real32, 0._real32,&
         0._real32, 1._real32, 0._real32,&
         0._real32, 0._real32, 1._real32/), shape(lat) ) )
    special(:,:,2) = transpose( reshape( (/&
         1._real32, 0._real32, 0._real32,&
         -0.5_real32, sqrt(3._real32)/2._real32, 0._real32,&
         0._real32, 0._real32, 1.0_real32/), shape(lat) ) )
    special(:,:,3) = transpose( reshape( (/&
         0.0_real32, 1._real32, 1._real32,&
         1._real32, 0._real32, 1._real32,&
         1._real32, 1._real32, 0.0_real32/), shape(lat) ) )
    special(:,:,3) = special(:,:,3) / sqrt(2._real32)
    special(:,:,4) = transpose( reshape( (/&
         -1._real32,  1._real32,  1._real32,&
         1._real32, -1._real32,  1._real32,&
         1._real32,  1._real32, -1._real32/), shape(lat) ) )
    special(:,:,4) = special(:,:,4) / sqrt(3._real32)


    !!---------------------------------------------------------------
    !! cycles special set to find primitive lattice of supplied lat
    !!---------------------------------------------------------------
    tmat1 = matmul(lat,transpose(lat))
    checkloop: do i = 1, 4
       !tfmat=matmul(lat,inverse_3x3(special(:,:,i)))
       !tfmat=matmul(tfmat,transpose(tfmat))
       tmat2 = matmul(special(:,:,i),transpose(special(:,:,i)))
       rtmp1 = tmat2(1,1) / tmat1(1,1)
       !if(all(abs(tfmat-nint(tfmat)).lt.1.E-8_real32))then
       if(all(abs(tmat1*rtmp1-tmat2).lt.1.E-6_real32))then
          do j = 1, 3
             plat(j,:) = scal(j) * special(j,:,i)
          end do
          exit checkloop
       end if
    end do checkloop

    basis%lat = plat

  end subroutine primitive_lat
!###############################################################################


!###############################################################################
  subroutine reducer(basis, tmptype, verbose)
    !! Reduce the cell using Buerger's algorithm
    implicit none

    ! Arguments
    type(basis_type), intent(inout) :: basis
    !! Structure data
    integer, intent(in), optional :: tmptype
    !! Cell type
    integer, intent(in), optional :: verbose
    !! Verbosity level

    ! Local variables
    integer :: cell_type
    !! Cell type
    integer :: i,j,k,count,limit
    !! Loop indices
    real(real32), dimension(3,3) :: newlat,transmat,S,tmp_mat
    !! Lattice matrices
    real(real32) :: tiny,pi2
    !! Constants
    integer :: verbose_
    !! Verbosity level
    logical :: lreduced
    !! Boolean whether cell is reduced


    !---------------------------------------------------------------------------
    ! set up inital variable values
    !---------------------------------------------------------------------------
    verbose_ = 0
    if(present(verbose)) verbose_ = verbose
    cell_type=2
    if(present(tmptype)) cell_type=tmptype
    S=0._real32
    count=0
    limit=100
    lreduced=.false.
    tiny = 1.E-5_real32 * (get_vol(basis%lat))**(1._real32/3._real32)
    pi2 = 2._real32*atan(1._real32)
    transmat = 0._real32
    do i = 1, 3
       transmat(i,i) = 1._real32
    end do
    newlat = basis%lat


    !---------------------------------------------------------------------------
    ! perform checks on the other main conditions defined by Niggli
    !---------------------------------------------------------------------------
    find_reduced: do while(.not.lreduced)
       count = count + 1
       call mkNiggli_lat(basis%lat,newlat,transmat,S)
       lreduced = reduced_check(newlat, cell_type, S, verbose_)
       if(lreduced) exit
       if(verbose_.gt.1) then
          write(*,*)
          write(*,*) count
          write(*,*) "###############"
          write(*,*) (transmat(i,:),i=1,3)
          write(*,*)
          write(*,*) (newlat(i,:),i=1,3)
       end if
       if(count.gt.limit) then
          write(0,'("FAILED to find the reduced cell within ",I0," steps")') count
          exit
       end if


       !! A1 & A2 
       do i=1,2
          j=i+1
          if(S(i,i)-S(j,j).gt.tiny) then
             call swap(transmat(i,:),transmat(j,:))
             transmat=-transmat
             if(i.eq.2) cycle find_reduced
             call mkNiggli_lat(basis%lat,newlat,transmat,S)
          end if
       end do


       !! A3
       i=1;j=1;k=1
       if(S(2,3).lt.0) i=-1
       if(S(1,3).lt.0) j=-1
       if(S(1,2).lt.0) k=-1
       if(i*j*k.gt.0) then
          tmp_mat=reshape((/i,0,0,  0,j,0,  0,0,k/),shape(tmp_mat))
          transmat=matmul(transpose(tmp_mat),transmat)
          call mkNiggli_lat(basis%lat,newlat,transmat,S)
       end if


       !! A4
       i=1;j=1;k=1
       if(S(1,2).lt.0) i=-1
       if(S(1,3).lt.0) j=-1
       if(S(1,2).lt.0) k=-1
       if(i*j*k.gt.0) then
          tmp_mat=reshape((/i,0,0,  0,j,0,  0,0,k/),shape(tmp_mat))
          transmat=matmul(transpose(tmp_mat),transmat)
          call mkNiggli_lat(basis%lat,newlat,transmat,S)
       end if


       tmp_mat=reshape((/1,0,0,  0,1,0,  0,0,1/),shape(tmp_mat))
       !! A5
       if(abs(2*S(2,3)).gt.S(2,2)+tiny.or.&
            (abs(2*S(2,3)-S(2,2)).le.tiny.and.2*S(1,3).lt.S(1,2)).or.&
            (abs(2*S(2,3)+S(2,2)).le.tiny.and.S(1,2).lt.0._real32))then
          tmp_mat(2,3)=((-1)**(cell_type+1))*floor((2*S(2,3)+S(2,2))/(2*S(2,2)))
          transmat=matmul(transpose(tmp_mat),transmat)
          cycle find_reduced
          !       elseif(cell_type.eq.1.and.S(2,3).lt.0._real32)then
          !          tmp_mat(2,3)=1._real32
          !          transmat=matmul(transpose(tmp_mat),transmat)
          !          cycle find_reduced
       end if


       !! A6
       if(abs(2*S(1,3)).gt.S(1,1)+tiny.or.&
            (abs(2*S(1,3)-S(1,1)).le.tiny.and.2*S(2,3).lt.S(1,2)).or.&
            (abs(2*S(1,3)+S(1,1)).le.tiny.and.S(1,2).lt.0._real32))then
          tmp_mat(1,3)=((-1)**(cell_type+1))*floor((2*S(1,3)+S(1,1))/(2*S(1,1)))
          transmat=matmul(transpose(tmp_mat),transmat)
          cycle find_reduced
          !       elseif(cell_type.eq.1.and.S(1,3).lt.0._real32)then
          !          tmp_mat(1,3)=1._real32
          !          transmat=matmul(transpose(tmp_mat),transmat)
          !          cycle find_reduced
       end if


       !! A7
       if(abs(2*S(1,2)).gt.S(1,1)+tiny.or.&
            (abs(2*S(1,2)-S(1,1)).le.tiny.and.2*S(2,3).lt.S(1,3)).or.&
            (abs(2*S(1,2)+S(1,1)).le.tiny.and.S(1,3).lt.0._real32))then
          tmp_mat(1,2)=((-1)**(cell_type+1))*floor((2*S(1,2)+S(1,1))/(2*S(1,1)))
          transmat=matmul(transpose(tmp_mat),transmat)
          cycle find_reduced
          !       elseif(cell_type.eq.1.and.S(1,2).lt.0._real32)then
          !          tmp_mat(1,2)=1._real32
          !          transmat=matmul(transpose(tmp_mat),transmat)
          !          cycle find_reduced
       end if


       !! A8
       if(cell_type.eq.2.and.&
            2*(S(2,3)+S(1,3)+S(1,2))+S(1,1)+S(2,2).lt.tiny.or.&
            (abs(2*(S(2,3)+S(1,3)+S(1,2))+S(1,1)+S(2,2)).le.tiny.and.&
            2*(S(1,1)+2*S(1,3))+2*S(1,2).gt.tiny))then
          tmp_mat(1,3)=((-1)**(cell_type+1))*floor( &
               ( 2*(S(2,3) + S(1,3) + S(1,2)) + S(1,1) + S(2,2)  )/&
               (2*(2*S(1,2) + S(1,1) + S(2,2) ))  )
          tmp_mat(2,3)=tmp_mat(1,3)
          transmat=matmul(transpose(tmp_mat),transmat)
          cycle find_reduced     
       end if

       lreduced=.true.
    end do find_reduced


    if(abs(det(transmat)+1._real32).le.tiny)then
       tmp_mat=reshape((/-1,0,0,  0,-1,0,  0,0,-1/),shape(tmp_mat))
       transmat=matmul(transpose(tmp_mat),transmat)
    end if
    call mkNiggli_lat(basis%lat,newlat,transmat,S)
    lreduced = reduced_check(newlat, cell_type, S, verbose_)
    if(verbose_.gt.1) then
       write(*,*) lreduced
       write(*,*) (transmat(i,:),i=1,3)
    end if


    !---------------------------------------------------------------------------
    ! Renormalise the lattice and basis into the new lattice
    !---------------------------------------------------------------------------
    basis%lat = newlat
    do i = 1, basis%nspec
       do j = 1, basis%spec(i)%num
          basis%spec(i)%atom(j,:3) = &
               matmul( basis%spec(i)%atom(j,:3), inverse_3x3(transmat) )
          basis%spec(i)%atom(j,:3) = &
               basis%spec(i)%atom(j,:3) - floor( basis%spec(i)%atom(j,:3) )
       end do
    end do

  end subroutine reducer
!###############################################################################


!!!#############################################################################
!!! Subroutine to set up the required dot products of the lattice
!!!#############################################################################
!!! a = lat(:,1),  b=lat(:,2),  c=lat(:,3)
!!! S(1,1) = a.a,   S(2,2) = b.b,   S(3,3) = c.c
!!! S(1,2) = a.b,   S(1,3) = a.c,   S(2,3) = b.c
  subroutine mkNiggli_lat(lat,newlat,transmat,S)
    implicit none
    real(real32), dimension(3,3) :: lat,newlat,transmat,S
    real(real32), dimension(3) :: a,b,c


    newlat=matmul(transmat,lat)
    a=newlat(1,:);b=newlat(2,:);c=newlat(3,:)

    S(1,1)=dot_product(a,a)
    S(2,2)=dot_product(b,b)
    S(3,3)=dot_product(c,c)
    S(2,3)=dot_product(b,c)
    S(1,3)=dot_product(a,c)
    S(1,2)=dot_product(a,b)

    return
  end subroutine mkNiggli_lat
!!!#############################################################################


!!!#############################################################################
!!! Function to check whether cell satisfies all the main ...
  ! ... Niggli conditions (1928)
!!!#############################################################################
!!! tiny = tolerance to satisfy conditions
!!! lat = lattice being checked
!!! a = lat(:,1),  b=lat(:,2),  c=lat(:,3)
!!! S(1,1) = a.a,   S(2,2) = b.b,   S(3,3) = c.c
!!! S(1,2) = a.b,   S(1,3) = a.c,   S(2,3) = b.c
!!! Type I  = Sij (i!=j) are all positive (angles <90)
!!! Type II = Sij (i!=j) are all negative or any zero (angles >=90)
!!! Cell is reduced if, and only if, all conditions are ...
!!! ... satisfied (Niggli 1928)
  function reduced_check(lat, cell_type, S, verbose) result(check)
    implicit none
    real(real32), dimension(3,3), intent(in) :: lat
    real(real32), dimension(3,3), intent(out) :: S
    integer, intent(in) :: cell_type
    integer :: verbose

    real(real32) :: tiny,alpha,beta,gamma,pi2
    real(real32), dimension(3) :: a,b,c
    logical :: check



    pi2 = 2._real32*atan(1._real32)
    check = .false.
    tiny = 1.E-3_real32

    a = lat(1,:); b = lat(2,:); c = lat(3,:)
    S(1,1) = dot_product(a,a)
    S(2,2) = dot_product(b,b)
    S(3,3) = dot_product(c,c)
    S(2,3) = dot_product(b,c)
    S(1,3) = dot_product(a,c)
    S(1,2) = dot_product(a,b)

    alpha=acos(S(2,3)/sqrt(S(2,2)*S(3,3)))
    beta=acos(S(1,3)/sqrt(S(1,1)*S(3,3)))
    gamma=acos(S(1,2)/sqrt(S(1,1)*S(2,2)))

    if(S(2,2)-S(3,3).lt.tiny.and.S(1,1)-S(2,2).lt.tiny) then
       check=.true.
    else
       check=.false.
       return
    end if
    if(cell_type.eq.1.and.&
         alpha.le.pi2.and.beta.le.pi2.and.gamma.le.pi2.and.&
         S(1,2)-0.5_real32*S(1,1).lt.tiny.and.&
         S(1,3)-0.5_real32*S(1,1).lt.tiny.and.&
         S(2,3)-0.5_real32*S(2,2).lt.tiny) then !Type I
       check=.true.
       if(verbose.gt.0) write(0,*) "Found Type I reduced Niggli cell"
    elseif(cell_type.eq.2.and.&
         alpha.ge.pi2-tiny.and.beta.ge.pi2-tiny.and.gamma.ge.pi2-tiny.and.&
         abs(S(1,2))-0.5_real32*S(1,1).lt.tiny.and.&
         abs(S(1,3))-0.5_real32*S(1,1).lt.tiny.and.&
         abs(S(2,3))-0.5_real32*S(2,2).lt.tiny.and.&
         (abs(S(2,3))+abs(S(1,3))+abs(S(1,2)))-0.5_real32*(S(1,1)+S(2,2)).lt.tiny) then !Type II
       if(abs(S(1,2))-0.5_real32*S(1,1).le.tiny.and.S(1,3).gt.tiny) return
       if(abs(S(1,3))-0.5_real32*S(1,1).le.tiny.and.S(1,2).gt.tiny) return
       if(abs(S(2,3))-0.5_real32*S(2,2).le.tiny.and.S(1,2).gt.tiny) return
       if((abs(S(2,3))+abs(S(1,3))+abs(S(1,2)))-0.5_real32*(S(1,1)+S(2,2)).gt.tiny.and.&
            S(1,1)-(2._real32*abs(S(1,3))+abs(S(1,2))).gt.tiny) return
       check=.true.
       if(verbose.gt.1) write(0,*) "Found Type II reduced Niggli cell"
    else
       check=.false.
    end if

  end function reduced_check
!!!#############################################################################


!!!#############################################################################
!!! planecutter
!!!#############################################################################
  function planecutter(lat, plane) result(tfmat)
    implicit none
    real(real32), dimension(3,3), intent(in) :: lat
    real(real32), dimension(3), intent(in) :: plane

    integer :: i,j,itmp1
    real(real32) :: tol
    integer, dimension(3) :: order
    real(real32), dimension(3) :: plane_,tvec1
    real(real32), dimension(3,3) :: lat_,b,tfmat,invlat,reclat



!!!-----------------------------------------------------------------------------
!!! Initialise variables and matrices
!!!-----------------------------------------------------------------------------
    tol    = 1.E-4_real32
    plane_ = plane
    lat_   = lat
    invlat = inverse(lat_)
    reclat = transpose(invlat)
    if(all(plane_.le.0._real32)) plane_ = -plane_
    plane_ = reduce_vec_gcd(plane_)
    order  = [ 1, 2, 3 ]


!!!-----------------------------------------------------------------------------
!!! Align the normal vector such that all non-zero values are left of all zeros
!!!-----------------------------------------------------------------------------
    do i=1,2
       if(plane_(i).eq.0)then
          if(all(plane_(i:).eq.0._real32)) exit
          itmp1=maxloc(plane_(i+1:),mask=plane_(i+1:).ne.0,dim=1)+i
          call swap(order(i),order(itmp1))
          call swap(plane_(i),plane_(itmp1))
          call swap(lat_(:,i),lat_(:,itmp1))
          call swap(lat_(i,:),lat_(itmp1,:))
          call swap(reclat(:,i),reclat(:,itmp1))
          call swap(reclat(i,:),reclat(itmp1,:))
       end if
    end do
    !plane_=matmul(plane_,reclat)


!!!-----------------------------------------------------------------------------
!!! Perform Lenstra-Lenstra-Lovász reduction
!!!-----------------------------------------------------------------------------
    b(1,:) = [ -plane_(2),plane_(1),0._real32 ]
    b(2,:) = [ -plane_(3),0._real32,plane_(1) ]
    b(3,:) = plane_
    tfmat = b
    b(:2,:) = LLL_reduce(b(:2,:))


!!!-----------------------------------------------------------------------------
!!! Checking whether b1 and b2 are still perpendicular to b3 and have size ...
!!! ... greater than zero
!!!-----------------------------------------------------------------------------
    if(dot_product(b(1,:),b(3,:)).gt.tol)then
       write(0,'("ERROR: Internatl error in planecutter")')
       write(0,'(2X,"Error in planecutter subroutine in mod_geom_utils.f90")')
       write(0,'(2X,"b1 not perpendicular to b3")')
       write(0,'(2X,"b1 = ",3(1X,F0.3))') b(1,:)
       write(0,'(2X,"b3 = ",3(1X,F0.3))') b(3,:)
       write(0,'(2X,"b1·b3 = ",F0.3)') dot_product(b(1,:),b(3,:))
       write(0,'("Inform developers of this issue")')
       write(0,'("Stopping...")')
       stop
    elseif(dot_product(b(2,:),b(3,:)).gt.tol)then
       write(0,'("ERROR: Internatl error in planecutter")')
       write(0,'(2X,"Error in planecutter subroutine in mod_geom_utils.f90")')
       write(0,'(2X,"b2 not perpendicular to b3")')
       write(0,'(2X,"b2 = ",3(1X,F6.2))') b(2,:)
       write(0,'(2X,"b3 = ",3(1X,F6.2))') b(3,:)
       write(0,'(2X,"b1·b3 = ",F6.3)') dot_product(b(2,:),b(3,:))
       write(0,'("Inform developers of this issue")')
       write(0,'("Stopping...")')
       stop
    elseif(dot_product(b(1,:),b(1,:)).lt.tol)then
       write(0,'("ERROR: Internatl error in planecutter")')
       write(0,'(2X,"Error in planecutter subroutine in mod_geom_utils.f90")')
       write(0,'(2X,"b1 has zero size")')
       write(0,'(2X,"b1 = ",3(1X,F6.2))') b(1,:)
       write(0,'("Inform developers of this issue")')
       write(0,'("Stopping...")')
       stop
    elseif(dot_product(b(2,:),b(2,:)).lt.tol)then
       write(0,'("ERROR: Internatl error in planecutter")')
       write(0,'(2X,"Error in planecutter subroutine in mod_geom_utils.f90")')
       write(0,'(2X,"b2 has zero size")')
       write(0,'(2X,"b2 = ",3(1X,F6.2))') b(2,:)
       write(0,'("Inform developers of this issue")')
       write(0,'("Stopping...")')
       stop
    end if

    !b = matmul(b,lat_)
    

!!!-----------------------------------------------------------------------------
!!! Fix normal vector and lattice
!!!-----------------------------------------------------------------------------
    do i=1,3
       if(i.eq.order(i)) cycle
       call swap(lat_(i,:),lat_(order(i),:))
       call swap(lat_(:,i),lat_(:,order(i)))
       call swap(b(:,i),b(:,order(i)))
       call swap(order(order(i)),order(i))
    end do


!!!-----------------------------------------------------------------------------
!!! Convert the new lattice to direct coordinates
!!! Make it such that it is a fully integerised transformation matrix
!!!-----------------------------------------------------------------------------
    !b=matmul(b,invlat)
    where(abs(b(:,:)).lt.tol)
       b(:,:)=0._real32
    end where
    !write(0,'(3(2X,F9.3))') (b(j,:),j=1,3)
    !write(0,*) 
    reduce_loop: do i=1,3
       b(i,:)=reduce_vec_gcd(b(i,:))
       if(any(abs(b(i,:)-nint(b(i,:))).gt.tol))then
          write(0,'("Issue with plane ",3(1X,I0))') nint(plane)
          write(0,*) plane_
          write(0,'("row ",I0," of the following matrix")') i
          write(0,'(3(2X,F9.3))') (b(j,:),j=1,3)
          write(0,'(1X,"ERROR: Internal error in planecutter function")')
          write(0,'(2X,"Planecutter in mod_geom_utils.f90 is unable to find a&
               & perpendicular plane")')
          b=0._real32
          exit
       end if
    end do reduce_loop
    if(det(b).lt.0._real32)then
       tvec1  = b(2,:)
       b(2,:) = b(1,:)
       b(1,:) = tvec1
    end if
    if(abs(det(b)).lt.tol)then
       write(0,'(1X,"ERROR: Internal error in planecutter function")')
       write(0,'(2X,"Planecutter in mod_geom_utils.f90 has generated a 0&
            & determinant matrix")')
       write(0,'(3(2X,F9.3))') (b(j,:),j=1,3)
       b=0._real32
       !stop
    end if
    tfmat=b


    return
  end function planecutter
!!!#############################################################################


!###############################################################################
  function basis_merge(basis1,basis2,length,map1,map2) result(output)
    !! Merge two supplied bases
    !!
    !! Merge two bases assuming that the lattice is the same
    implicit none

    ! Arguments
    type(basis_type) :: output
    !! Output merged basis.
    class(basis_type), intent(in) :: basis1, basis2
    !! Input bases to merge.
    integer, intent(in), optional :: length
    !! Number of dimensions for atomic positions (default 3).
    integer, allocatable, dimension(:,:,:), optional, intent(inout) :: map1,map2
    !! Maps for atoms in the two bases.

    ! Local variables
    integer :: i, j, k, itmp, dim
    !! Loop counters.
    logical :: lmap
    !! Boolean for map presence.
    integer, allocatable, dimension(:) :: match
    !! Array to match species.
    integer, allocatable, dimension(:,:,:) :: new_map
    !! New map for merged basis.



    !---------------------------------------------------------------------------
    ! set up number of species
    !---------------------------------------------------------------------------
    dim=3
    if(present(length)) dim=length

    allocate(match(basis2%nspec))
    match=0
    output%nspec=basis1%nspec
    do i = 1, basis2%nspec
       if(.not.any(basis2%spec(i)%name.eq.basis1%spec(:)%name))then
          output%nspec=output%nspec+1
       end if
    end do
    allocate(output%spec(output%nspec))
    output%spec(:basis1%nspec)%num=basis1%spec(:)%num
    output%spec(:basis1%nspec)%name=basis1%spec(:)%name


    write(output%sysname,'(A,"+",A)') &
         trim(basis1%sysname),trim(basis2%sysname)
    k=basis1%nspec
    spec1check: do i = 1, basis2%nspec
       do j = 1, basis1%nspec
          if(basis2%spec(i)%name.eq.basis1%spec(j)%name)then
             output%spec(j)%num=output%spec(j)%num+basis2%spec(i)%num
             match(i)=j
             cycle spec1check
          end if
       end do
       k=k+1
       match(i)=k
       output%spec(k)%num=basis2%spec(i)%num
       output%spec(k)%name=basis2%spec(i)%name
    end do spec1check


    !---------------------------------------------------------------------------
    ! if map is present, sets up new map
    !---------------------------------------------------------------------------
    lmap = .false.
    if_map: if(present(map1).and.present(map2))then
       if(all(map1.eq.-1)) exit if_map
       lmap = .true.
       allocate(new_map(&
            output%nspec,&
            maxval(output%spec(:)%num,dim=1),2))
       new_map = 0
    end if if_map


    !---------------------------------------------------------------------------
    ! set up atoms in merged basis
    !---------------------------------------------------------------------------
    do i = 1, basis1%nspec
       allocate(output%spec(i)%atom(output%spec(i)%num,dim))
       output%spec(i)%atom(:,:)=0._real32
       output%spec(i)%atom(1:basis1%spec(i)%num,:3)=basis1%spec(i)%atom(:,:3)
       if(lmap) new_map(i,:basis1%spec(i)%num,:)=map1(i,:basis1%spec(i)%num,:)
    end do
    do i = 1, basis2%nspec
       if(match(i).gt.basis1%nspec)then
          allocate(output%spec(match(i))%atom(output%spec(match(i))%num,dim))
          output%spec(match(i))%atom(:,:)=0._real32
          output%spec(match(i))%atom(:,:3)=basis2%spec(i)%atom(:,:3)
          if(lmap) new_map(match(i),:basis2%spec(i)%num,:) = &
               map2(i,:basis2%spec(i)%num,:)
       else
          itmp=basis1%spec(match(i))%num
          output%spec(match(i))%atom(itmp+1:basis2%spec(i)%num+itmp,:3) = &
               basis2%spec(i)%atom(:,:3)   
          if(lmap) new_map(match(i),itmp+1:basis2%spec(i)%num+itmp,:) = &
               map2(i,:basis2%spec(i)%num,:)      
       end if
    end do
    output%natom=sum(output%spec(:)%num)
    output%lat = basis1%lat


    if(lmap) call move_alloc(new_map,map1)

    return
  end function basis_merge
!###############################################################################


!###############################################################################
  function basis_stack(basis1,basis2,axis,offset,length,map1,map2) result(output)
    !! Merge two supplied bases
    !!
    !! Merge two bases assuming that the lattice is the same
    implicit none

    ! Arguments
    type(basis_type) :: output
    !! Output merged basis.
    class(basis_type), intent(in) :: basis1, basis2
    !! Input bases to merge.
    integer, intent(in), optional :: length
    !! Number of dimensions for atomic positions (default 3).
    integer, intent(in) :: axis
    !! Axis for the offset.
    real(real32), dimension(3), intent(in) :: offset
    !! Offset for the merged basis.
    integer, allocatable, dimension(:,:,:), optional, intent(inout) :: map1,map2
    !! Maps for atoms in the two bases.

    ! Local variables
    integer :: i, j, k, length_
    !! Loop counters.
    real(real32) :: loc, c1_ratio, c2_ratio, zgap, add
    !! Lattice parameters.
    type(basis_type) :: basis1_, basis2_
    integer, dimension(3) :: order
    !! Order of axes.
    real(real32), dimension(3) :: unit_vec
    !! Unit vector for the axis.
    real(real32), dimension(3) :: offset_
    !! Offset for the merged basis.
    integer, allocatable, dimension(:) :: match
    !! Array to match species.
    real(real32), dimension(3,3) :: output_lat
    !! Output lattice.


    !---------------------------------------------------------------------------
    ! copy basis1 and basis2
    !---------------------------------------------------------------------------
    call basis1_%copy(basis1)
    call basis2_%copy(basis2)


    !---------------------------------------------------------------------------
    ! set up number of species
    !---------------------------------------------------------------------------
    length_ = 3
    if(present(length)) length_ = length

    allocate(match(basis2_%nspec))
    match=0
    output%nspec=basis1_%nspec
    do i = 1, basis2_%nspec
       if(.not.any(basis2_%spec(i)%name.eq.basis1_%spec(:)%name))then
          output%nspec=output%nspec+1
       end if
    end do
    allocate(output%spec(output%nspec))
    output%spec(:basis1_%nspec)%num=basis1_%spec(:)%num
    output%spec(:basis1_%nspec)%name=basis1_%spec(:)%name


    write(output%sysname,'(A,"+",A)') &
         trim(basis1_%sysname),trim(basis2_%sysname)
    k=basis1_%nspec
    spec1check: do i = 1, basis2_%nspec
       do j = 1, basis1_%nspec
          if(basis2_%spec(i)%name.eq.basis1_%spec(j)%name)then
             output%spec(j)%num=output%spec(j)%num+basis2_%spec(i)%num
             match(i)=j
             cycle spec1check
          end if
       end do
       k=k+1
       match(i)=k
       output%spec(k)%num=basis2_%spec(i)%num
       output%spec(k)%name=basis2_%spec(i)%name
    end do spec1check


    !-----------------------------------------------------------------------------
    ! Shifts cells to 
    !-----------------------------------------------------------------------------
    loc=0.D0
    basis1_%lat=MATNORM(basis1_%lat)
    add = -min_dist(basis1_,axis,loc,.true.)
    call shifter(basis1_,axis,add,.true.)

    basis2_%lat=MATNORM(basis2_%lat)
    add = -min_dist(basis2_,axis,loc,.true.)
    call shifter(basis2_,axis,add,.true.)



    !---------------------------------------------------------------------------
    ! handle offset
    !---------------------------------------------------------------------------
    loc = 1._real32
    call set_vacuum(basis1_,axis,loc,offset(axis))
    call set_vacuum(basis2_,axis,loc,offset(axis))
  
    order = [ 1, 2, 3 ]
    order = cshift(order,3-axis)
    do k = 1, 2
       offset_(order(k)) = offset(order(k)) / modu(basis1_%lat(order(k),:))
    end do
    unit_vec = uvec(basis1_%lat(order(3),:))
    zgap = offset_(order(3)) / unit_vec(order(3))



    !---------------------------------------------------------------------------
    ! makes supercell
    !---------------------------------------------------------------------------
    output_lat(order(1),:) = basis1_%lat(order(1),:)
    output_lat(order(2),:) = basis1_%lat(order(2),:)
    unit_vec = uvec(basis1_%lat(axis,:))
    output_lat(axis,:) = basis1_%lat(axis,:) + modu(basis2_%lat(axis,:)) * unit_vec
    c1_ratio = modu(basis1_%lat(axis,:)) / modu(output_lat(axis,:))
    c2_ratio = modu(basis2_%lat(axis,:)) / modu(output_lat(axis,:))


!!!-----------------------------------------------------------------------------
!!! merge list of atomic types and respective numbers for both structures
!!!-----------------------------------------------------------------------------
    do i=1,basis1_%nspec
       basis1_%spec(i)%atom(:,axis) = basis1_%spec(i)%atom(:,axis) * c1_ratio
    end do
    do i=1,basis2_%nspec
       basis2_%spec(i)%atom(:,axis) = basis2_%spec(i)%atom(:,axis)*c2_ratio + c1_ratio
       do k=1,2
          basis2_%spec(i)%atom(:,order(k)) = basis2_%spec(i)%atom(:,order(k)) + offset_(order(k))
       end do
    end do

    if(present(map1).and.present(map2))then
       output = basis_merge(basis1_,basis2_,map1=map1,map2=map2)
    else
       output = basis_merge(basis1_,basis2_)
    end if
    output%lat = output_lat
    call output%normalise(ceil_val = 1._real32, floor_coords = .true.)

  end function basis_stack
!###############################################################################


!!!#############################################################################
!!! splits basis into an array of bases
!!!#############################################################################
  function split_bas(inbas,loc_vec,axis,lall_same_nspec,map1,map2) result(bas_arr)
    implicit none
    integer :: i,is,ia,itmp1,nregions,axis,nspec
    logical :: lsame
    logical :: lmap,lmove
    type(basis_type) :: tbas
    real(real32), allocatable, dimension(:,:) :: dloc_vec
    logical, optional :: lall_same_nspec

    type(basis_type),intent(in) :: inbas
    real(real32), dimension(:,:), intent(in) :: loc_vec
    type(basis_type), allocatable, dimension(:) :: bas_arr

    type map_type   
       integer, allocatable, dimension(:,:,:) :: spec       
    end type map_type
    type(map_type), dimension(2) :: map
    integer, allocatable, dimension(:,:,:), optional, intent(inout) :: map1,map2


    !!--------------------------------------------------------------------------
    !! If map1 is present, sets up map2
    !!--------------------------------------------------------------------------
    if(present(map1).and.present(map2))then
       lmap=.true.
    else
       lmap=.false.
    end if


    !!--------------------------------------------------------------------------
    !! Sets up region locations
    !!--------------------------------------------------------------------------
    lsame=.true.
    if(present(lall_same_nspec)) lsame=lall_same_nspec
    nregions=size(loc_vec(:,1),dim=1)
    allocate(dloc_vec(nregions,2))
    dloc_vec(:,:)=loc_vec(:,:)-floor(loc_vec(:,:))
    where(dloc_vec(:,2).lt.dloc_vec(:,1))
       dloc_vec(:,2)=dloc_vec(:,2)+1._real32
    end where
    allocate(bas_arr(nregions))

    
    !!--------------------------------------------------------------------------
    !! Loops over regions and assigns atoms to them
    !!--------------------------------------------------------------------------
    regionloop1: do i=1,nregions
       bas_arr(i)%natom = 0
       bas_arr(i)%nspec = inbas%nspec
       write(bas_arr(i)%sysname,'(A,"_region_",I0)') trim(inbas%sysname),i
       allocate(bas_arr(i)%spec(inbas%nspec))
       if(lmap) allocate(map(i)%spec(bas_arr(i)%nspec,maxval(inbas%spec(:)%num),2))

       specloop1: do is=1,inbas%nspec
          bas_arr(i)%spec(is)%name=inbas%spec(is)%name
          bas_arr(i)%spec(is)%num = &
               count(inbas%spec(is)%atom(:,axis) - &
               floor(inbas%spec(is)%atom(:,axis)-dloc_vec(i,1))&
               .lt.dloc_vec(i,2))
          bas_arr(i)%natom = bas_arr(i)%natom + bas_arr(i)%spec(is)%num
          allocate(bas_arr(i)%spec(is)%atom(bas_arr(i)%spec(is)%num,3))
          itmp1=0
          atomloop1: do ia=1,inbas%spec(is)%num
             if(inbas%spec(is)%atom(ia,axis) - &
                  floor(inbas%spec(is)%atom(ia,axis)-dloc_vec(i,1))&
                  .lt.dloc_vec(i,2))then
                itmp1=itmp1+1
                bas_arr(i)%spec(is)%atom(itmp1,:3) = inbas%spec(is)%atom(ia,:3)
                if(lmap) map(i)%spec(is,itmp1,:) = map1(is,ia,:)
             end if

          end do atomloop1
       end do specloop1
    end do regionloop1


    !!--------------------------------------------------------------------------
    !! Revmoes null species from regions, if specified
    !!--------------------------------------------------------------------------
    if(.not.lsame)then
       do i=1,nregions
          nspec=0
          lmove=.false.
          tbas%nspec=count(bas_arr(i)%spec(:)%num.gt.0)
          tbas%natom=bas_arr(i)%natom
          tbas%sysname=bas_arr(i)%sysname
          tbas%lat = inbas%lat
          allocate(tbas%spec(tbas%nspec))

          if(lmap.and.i.eq.1)then
             if(allocated(map1)) deallocate(map1)
             allocate(map1(tbas%nspec,size(map(1)%spec(1,:,1),dim=1),2))
             map1=0
          elseif(lmap.and.i.eq.2)then
             if(allocated(map2)) deallocate(map2)
             allocate(map2(tbas%nspec,size(map(2)%spec(1,:,1),dim=1),2))
             map2=0
          end if

          specloop2: do is=1,bas_arr(i)%nspec
             if(bas_arr(i)%spec(is)%num.eq.0)then
                lmove=.true.
                cycle specloop2
             else
                tbas%spec(nspec+1)%num=bas_arr(i)%spec(is)%num
                tbas%spec(nspec+1)%name=bas_arr(i)%spec(is)%name
                if(lmap.and.i.eq.1) map1(nspec+1,:,:) = map(i)%spec(is,:,:)
                if(lmap.and.i.eq.2) map2(nspec+1,:,:) = map(i)%spec(is,:,:)
                call move_alloc(bas_arr(i)%spec(is)%atom,tbas%spec(nspec+1)%atom)
                lmove=.false.                
             end if
             nspec=nspec+1
          end do specloop2
          call bas_arr(i)%copy(tbas)
          deallocate(tbas%spec)
       end do
    end if

  end function split_bas
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
    real(real32) :: scale,projection,dtmp1
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
          projection=1.E2_real32
          trans_loop: do j=1,ntrans+3
             dtmp1 = dot_product(trans(j,:),trans(ntrans+i,:))
             if(dtmp1.lt.tol_sym) cycle trans_loop

             do k=1,i-1,1
                if(modu(abs(cross( [ trans(j,:) ], [ dmat1(k,:) ]))).lt.1.E-8_real32) cycle trans_loop
             end do

             dtmp1 = modu( [ trans(j,:) ] )
             if(dtmp1.lt.projection)then
                projection=dtmp1
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
                     [ tol_sym,tol_sym,tol_sym ])) cycle atcheck
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
    !! next line necessary as FCC and BCC do not conform to Niggli reduced ...
    !! ... cell definitions.
    call primitive_lat(basis)


    
  end subroutine get_primitive_cell
!!!#############################################################################

!!!#############################################################################
!!! returns the bulk basis and lattice of 
!!!#############################################################################
  subroutine get_bulk(lat,bas,axis,bulk_lat,bulk_bas)
    implicit none
    integer :: is,ia,ja,len,itmp1
    integer :: minspecloc,minatomloc,nxtatomloc
    real(real32), dimension(3) :: transvec
    real(real32), dimension(2,2) :: regions
    logical, allocatable, dimension(:) :: atom_mask
    type(basis_type), allocatable, dimension(:) :: splitbas

    integer, intent(in) :: axis
    type(basis_type), intent(in) :: bas
    real(real32), dimension(3,3), intent(in):: lat
    type(basis_type), intent(out) :: bulk_bas
    real(real32), dimension(3,3), intent(out) :: bulk_lat


    minspecloc = minloc(bas%spec(:)%num,mask=bas%spec(:)%num.ne.0,dim=1)
    if(bas%spec(minspecloc)%num.eq.1)then
       write(0,'("ERROR: Internal error in get_bulk")')
       write(0,'(2X,"get_bulk subroutine in mod_geom_utils.f90 unable cannot &
            &find enough atoms to reproduce a bulk from")')
       stop
    end if
    minatomloc = minloc(bas%spec(minspecloc)%atom(:,axis),dim=1)



    allocate(atom_mask(bas%spec(minspecloc)%num))
    atom_mask = .true.



    itmp1 = minatomloc
    where(bas%spec(minspecloc)%atom(:,axis).le.&
         bas%spec(minspecloc)%atom(itmp1,axis))
       atom_mask(:) = .false.
    end where
    region_loop1: do

       atom_mask(itmp1) = .false.
       where(bas%spec(minspecloc)%atom(:,axis).lt.&
            bas%spec(minspecloc)%atom(itmp1,axis))
          atom_mask(:) = .false.
       end where
       nxtatomloc = minloc(&
            abs(&
            bas%spec(minspecloc)%atom(:,axis)-&
            bas%spec(minspecloc)%atom(itmp1,axis)),&
            mask = atom_mask,&
            dim=1)
       write(0,*) minatomloc,nxtatomloc
       
       transvec = bas%spec(minspecloc)%atom(nxtatomloc,:) - &
            bas%spec(minspecloc)%atom(minatomloc,:)

       regions(1,1:2) = [ &
            bas%spec(minspecloc)%atom(minatomloc,axis), &
            bas%spec(minspecloc)%atom(nxtatomloc,axis) ]
       regions(2,1:2) = [ &
            bas%spec(minspecloc)%atom(nxtatomloc,axis), &
            bas%spec(minspecloc)%atom(minatomloc,axis) ]
       splitbas = split_bas(bas,regions,axis)


       spec_loop1: do is=1,bas%nspec
          atom_loop1: do ia=1,splitbas(1)%spec(is)%num

             atom_loop2: do ja=1,splitbas(2)%spec(is)%num

                if( all( abs( ( splitbas(1)%spec(is)%atom(ia,:3) + transvec ) - &
                     splitbas(2)%spec(is)%atom(ja,:3) ).lt.1.E-5_real32 ) )then
                   write(0,*) ia,ja
                   cycle atom_loop1


                end if

             end do atom_loop2
             itmp1 = nxtatomloc
             cycle region_loop1

          end do atom_loop1
       end do spec_loop1
       exit region_loop1

    end do region_loop1


    len=size(bas%spec(1)%atom(1,:))
    allocate(bulk_bas%spec(bas%nspec))
    bulk_bas%nspec=bas%nspec
    bulk_bas%sysname=bas%sysname


    bulk_lat = lat
    bulk_lat(axis,:) = matmul(transvec,lat)


    call bulk_bas%copy(splitbas(1))
    call bulk_bas%change_lattice(bulk_lat)

  end subroutine get_bulk
!!!#############################################################################


!!!#############################################################################
!!! returns the atom closest to the centre of the region for a species
!!!#############################################################################
  function get_centre_atom(bas,spec,axis,lw,up) result(iatom)
    implicit none
    integer :: ia
    integer :: iatom
    real(real32) :: dtmp1,dtmp2,centre
    real(real32) :: dlw,dup
    integer, intent(in) :: spec,axis
    real(real32), intent(in) :: lw,up
    type(basis_type), intent(in) :: bas


    iatom=0
    dtmp1 = 1._real32
    if(lw.gt.up)then
       dlw = lw
       dup = 1._real32 + up
    else
       dlw = lw
       dup = up
    end if
    centre = (dlw + dup)/2._real32
    do ia=1,bas%spec(spec)%num
       dtmp2=bas%spec(spec)%atom(ia,axis)&
            -ceiling(bas%spec(spec)%atom(ia,axis)-dup)
       if(abs(dtmp2-centre).lt.dtmp1)then
          dtmp1=abs(dtmp2-centre)
          iatom=ia
       end if
    end do


  end function get_centre_atom
!!!#############################################################################


!!!#############################################################################
!!! returns the atom closest to the location
!!!#############################################################################
  function get_closest_atom_1D(bas,axis,loc,species,above,below) result(atom)
    implicit none
    integer :: is,ia
    integer :: is_start,is_end
    real(real32) :: dtmp1,dtmp2
    logical :: labove,lbelow
    integer, intent(in) :: axis
    real(real32), intent(in) :: loc
    integer, dimension(2) :: atom
    type(basis_type), intent(in) :: bas

    integer, optional, intent(in) :: species
    logical, optional, intent(in) :: above,below

    
    atom=[0,0]
    dtmp1 = 1._real32
    if(present(species))then
       is_start=species
       is_end=species
    else
       is_start=1
       is_end=bas%nspec
    end if
    labove=.false.
    lbelow=.false.
    if(present(above)) labove=above
    if(present(below)) lbelow=below

    do is=is_start,is_end
       atom_loop1: do ia=1,bas%spec(is)%num
          if(labove.and.bas%spec(is)%atom(ia,axis).lt.loc)then
             cycle atom_loop1
          elseif(lbelow.and.bas%spec(is)%atom(ia,axis).gt.loc)then
             cycle atom_loop1
          end if
          dtmp2=bas%spec(is)%atom(ia,axis)&
               -ceiling(bas%spec(is)%atom(ia,axis)-(loc+0.5_real32))
          if(abs(dtmp2-loc).lt.dtmp1)then
             dtmp1=abs(dtmp2-loc)
             atom=[is,ia]
          end if
       end do atom_loop1
    end do
       



  end function get_closest_atom_1D
!!!-----------------------------------------------------
!!!-----------------------------------------------------
  function get_closest_atom_3D(lat,bas,loc,species) result(atom)
    implicit none
    integer :: is,ia
    integer :: is_start,is_end
    real(real32) :: dtmp1,dtmp2
    real(real32), dimension(3) :: vtmp1
    real(real32), dimension(3), intent(in) :: loc
    real(real32), dimension(3,3), intent(in) :: lat
    integer, dimension(2) :: atom
    type(basis_type), intent(in) :: bas

    integer, optional, intent(in) :: species

    
    atom=[0,0]
    dtmp1 = 1._real32
    if(present(species))then
       is_start=species
       is_end=species
    else
       is_start=1
       is_end=bas%nspec
    end if

    spec_loop1: do is=is_start,is_end
       atom_loop1: do ia=1,bas%spec(is)%num
          vtmp1 = bas%spec(is)%atom(ia,:) - loc
          vtmp1 = vtmp1 - ceiling(vtmp1 - 0.5_real32)
          vtmp1 = matmul(vtmp1,lat)
          dtmp2 = modu(vtmp1)
          if(dtmp2.lt.dtmp1)then
             dtmp1=dtmp2
             atom=[is,ia]
          end if
       end do atom_loop1
    end do spec_loop1
    

  end function get_closest_atom_3D
!!!#############################################################################


!!!#############################################################################
!!! returns the wyckoff atom for each
!!!#############################################################################
  function get_wyckoff(bas,axis) result(wyckoff)
    implicit none
    integer :: is,ia,ja,itmp1,itmp2!ref_atom
    integer :: minspecloc,minatomloc,nxtatomloc
    real(real32) :: up_loc,lw_loc,up_loc2,lw_loc2
    real(real32), dimension(3) :: transvec,tmp_vec1,tmp_vec2,tmp_vec3,tvec
    logical, allocatable, dimension(:) :: atom_mask
    type(wyck_spec_type) :: wyckoff
    integer, intent(in) :: axis
    type(basis_type), intent(in) :: bas

    type l_bulk_type
       logical, allocatable, dimension(:) :: atom
    end type l_bulk_type
    type(l_bulk_type), allocatable, dimension(:) :: l_bulk_atoms

    
    
!!!-----------------------------------------------------------------------------
!!! Finds the species with the minimum number of atoms
!!! Finds upper and lower locations for "slab" and finds atom nearest to the ...
!!! ... centre of that region
!!!-----------------------------------------------------------------------------
    minspecloc = minloc(bas%spec(:)%num,mask=bas%spec(:)%num.ne.0,dim=1)
    minatomloc = minloc(bas%spec(minspecloc)%atom(:,axis),dim=1)
    nxtatomloc = maxloc(bas%spec(minspecloc)%atom(:,axis),dim=1)
    lw_loc = bas%spec(minspecloc)%atom(minatomloc,axis)
    up_loc = bas%spec(minspecloc)%atom(nxtatomloc,axis)
    minatomloc = &
         maxloc(bas%spec(minspecloc)%atom(:,axis)-(lw_loc+up_loc)/2._real32,dim=1,&
         mask=bas%spec(minspecloc)%atom(:,axis)-(lw_loc+up_loc)/2._real32.le.0._real32)
    allocate(atom_mask(bas%spec(minspecloc)%num))
    atom_mask = .true.



!!! INSTEAD OF STARTING FROM BOTTOM, START FROM CLOSEST BELOW MIDDLE AND CLOSEST ABOVE MIDDLE
!!! THEN WORK YOUR WAY OUT FROM THAT GOING 1 BELOW, THEN 1 ABOVE, etc.


!!!-----------------------------------------------------------------------------
!!! Set up lower atom location
!!!-----------------------------------------------------------------------------
    itmp1 = minatomloc
    lw_loc = bas%spec(minspecloc)%atom(minatomloc,axis)
    up_loc = 1._real32
    allocate(l_bulk_atoms(bas%nspec))
    do is=1,bas%nspec
       allocate(l_bulk_atoms(is)%atom(bas%spec(is)%num))
       l_bulk_atoms(is)%atom(:) = .false.
    end do

    
!!!-----------------------------------------------------------------------------
!!! Loops over atoms in cell until it finds a reproducible set to define ...
!!! ... as the bulk
!!!-----------------------------------------------------------------------------
    region_loop1: do
       !!-----------------------------------------------------------------------
       !! Mask of whether an atom has been checked for bulk limits or not
       !!-----------------------------------------------------------------------
       atom_mask(itmp1) = .false.
       where(bas%spec(minspecloc)%atom(:,axis).lt.&
            bas%spec(minspecloc)%atom(itmp1,axis))
          atom_mask(:) = .false.
       end where
       if(all(.not.atom_mask))then
          write(0,'("ERROR: Internal error in get_wyckoff")')
          write(0,'(2X,"Error in subroutine get_wyckoff in mod_geom_utils.f90")')
          write(0,'(2X,"No bulk found")')
          write(0,'(2X,"Exiting subroutine...")')
          return
       end if


       !!-----------------------------------------------------------------------
       !! Defines next atom to use as uper limit for bulk cell
       !!-----------------------------------------------------------------------
       nxtatomloc = minloc(&
            abs(&
            bas%spec(minspecloc)%atom(:,axis)-&
            bas%spec(minspecloc)%atom(itmp1,axis)),&
            mask = atom_mask,&
            dim=1)
       transvec = bas%spec(minspecloc)%atom(nxtatomloc,:) - &
            bas%spec(minspecloc)%atom(minatomloc,:)


       !!-----------------------------------------------------------------------
       !! Checks atoms within a region to see if they reproduce layer above
       !!-----------------------------------------------------------------------
       up_loc = lw_loc + transvec(axis)
       !if(lw_loc.eq.up_loc) up_loc = up_loc + 1.E-8_real32  !! IS THIS NEEDED?
       if(lw_loc.gt.up_loc)then
          write(0,'("ERROR: Internal error in get_wyckoff")')
          write(0,'(2X,"Error in subroutine get_wyckoff in mod_geom_utils.f90")')
          write(0,'(2X,"Region size is negative")')
          write(0,'(2X,"Stopping...")')
          stop
       end if
       spec_loop1: do is=1,bas%nspec
          l_bulk_atoms(is)%atom(:)=.false.
          if(bas%spec(is)%num.eq.0) cycle spec_loop1
          atom_loop1: do ia=1,bas%spec(is)%num
             if(bas%spec(is)%atom(ia,axis).lt.lw_loc.or.&
                  bas%spec(is)%atom(ia,axis).ge.up_loc)then
                cycle atom_loop1
             else
                l_bulk_atoms(is)%atom(ia)=.true.
             end if
             atom_loop2: do ja=1,bas%spec(is)%num
                tmp_vec2 = ( bas%spec(is)%atom(ia,:3) + transvec ) - &
                     bas%spec(is)%atom(ja,:3)
                !! SAME ISSUE HERE AS BELOW
                !! NEED TO TAKE INTO ACCOUNT THAT THEY WORK IN UNISON
                tmp_vec2 = tmp_vec2 - ceiling( tmp_vec2 - 0.5_real32 )


                if( all( abs(tmp_vec2).lt.1.E-5_real32 ) )then
                   cycle atom_loop1
                end if

             end do atom_loop2
             itmp1 = nxtatomloc
             cycle region_loop1

          end do atom_loop1
          if(all(.not.l_bulk_atoms(is)%atom(:)))then
             itmp1 = nxtatomloc
             cycle region_loop1
          end if
       end do spec_loop1


       !!-----------------------------------------------------------------------
       !! Checks the layer above
       !!-----------------------------------------------------------------------
       lw_loc2 = lw_loc + transvec(axis)
       up_loc2 = lw_loc2 + transvec(axis)
       spec_loop2: do is=1,bas%nspec
          if(bas%spec(is)%num.eq.0) cycle spec_loop2
          atom_loop3: do ia=1,bas%spec(is)%num
             if(bas%spec(is)%atom(ia,axis).lt.lw_loc2.or.&
                  bas%spec(is)%atom(ia,axis).ge.up_loc2) cycle atom_loop3
             tmp_vec1 = bas%spec(is)%atom(ia,:3) + transvec
             if( all(bas%spec(is)%atom(:,axis).lt.tmp_vec1(axis)-1.E-5_real32) ) cycle atom_loop3
             atom_loop4: do ja=1,bas%spec(is)%num
                tmp_vec2 = tmp_vec1 - bas%spec(is)%atom(ja,:3)
                tmp_vec2 = tmp_vec2 - ceiling( tmp_vec2 - 0.5_real32 )
                if( all( abs(tmp_vec2).lt.1.E-5_real32 ) )then
                   cycle atom_loop3
                end if
             end do atom_loop4
             itmp1 = nxtatomloc
             cycle region_loop1

          end do atom_loop3
          if(all(.not.l_bulk_atoms(is)%atom(:)))then
             itmp1 = nxtatomloc
             cycle region_loop1
          end if
       end do spec_loop2

       !!-----------------------------------------------------------------------
       !! If it gets to this point, then it has found a bulk cell and exits
       !!-----------------------------------------------------------------------
       exit region_loop1


    end do region_loop1



!!!-----------------------------------------------------------------------------
!!! Using the bulk definition, loop runs through checking which atom maps ...
!!! ... onto which through the bulk translation.
!!! Defines each atom's cell centre wyckoff atom
!!!-----------------------------------------------------------------------------
    allocate(wyckoff%spec(bas%nspec))
    do is=1,bas%nspec
       allocate(wyckoff%spec(is)%atom(bas%spec(is)%num))
       wyckoff%spec(is)%atom(:)=0
       atom_loop5: do ia=1,bas%spec(is)%num
          if(l_bulk_atoms(is)%atom(ia))then
             wyckoff%spec(is)%atom(ia) = ia
          end if
          !write(0,*) is,ia,l_bulk_atoms(is)%atom(ia)
          tmp_vec2 = bas%spec(is)%atom(ia,:3)
          atom_loop6: do ja=1,bas%spec(is)%num
             if(.not.l_bulk_atoms(is)%atom(ja)) cycle atom_loop6
             tmp_vec3 = tmp_vec2 - bas%spec(is)%atom(ja,:3)
             itmp1 = nint(tmp_vec3(axis)/transvec(axis))
             tvec = itmp1*transvec
             tvec = tvec - ceiling(tvec-1._real32)
             !tmp_vec3 = tmp_vec3/transvec
             !tmp_vec3 = reduce_vec_gcd(tmp_vec3)
             itmp2 = nint(get_vec_multiple(tvec,tmp_vec3))

             if(itmp1.eq.0) cycle atom_loop6
             tmp_vec3 = tmp_vec3 - tvec!itmp1*tvec
             tmp_vec3 = tmp_vec3 - ceiling(tmp_vec3 - 0.5_real32)
             !THIS IS WHERE WE NEED TO MAKE IT RIGHT
             !! FIND THE GCD AND DIVIDE
             if(all(abs(tmp_vec3).lt.1.E-5_real32))then
                if(wyckoff%spec(is)%atom(ja).ne.0)then
                   wyckoff%spec(is)%atom(ia) = wyckoff%spec(is)%atom(ja)
                else
                   wyckoff%spec(is)%atom(ia) = ja
                end if
                cycle atom_loop5
             end if
          end do atom_loop6
       end do atom_loop5


       if(any(wyckoff%spec(is)%atom(:).eq.0))then
          write(0,'("ERROR: Internal error in get_wyckoff")')
          write(0,'(2X,"Error in subroutine get_wyckoff in mod_geom_utils.f90")')
          write(0,'(2X,"Not all wyckoff atoms found")')
          do ia=1,bas%spec(is)%num
             write(0,*) is,ia,wyckoff%spec(is)%atom(ia)
          end do
          write(0,'(2X,"Stopping...")')
          stop
       end if


    end do



  end function get_wyckoff
!!!#############################################################################


!!!#############################################################################
!!! identify the shortest bond in the crystal, takes in crystal basis
!!!#############################################################################
  function get_shortest_bond(basis) result(bond)
    implicit none
    type(basis_type), intent(in) :: basis

    integer :: is,js,ia,ja,ja_start
    real(real32) :: dist,min_bond
    type(bond_type) :: bond
    real(real32), dimension(3) :: vec
    integer, dimension(2,2) :: atoms
    
    min_bond = huge(0._real32)
    atoms = 0
    do is = 1, basis%nspec
       do js = is, basis%nspec
          do ia = 1, basis%spec(is)%num
             if(is.eq.js)then
                ja_start = ia + 1
             else
                ja_start = 1
             end if
             do ja=ja_start,basis%spec(js)%num
                vec = basis%spec(is)%atom(ia,:3) - basis%spec(js)%atom(ja,:3)
                vec = vec - ceiling(vec - 0.5_real32)
                vec = matmul(vec,basis%lat)
                dist = modu(vec)
                if(dist.lt.min_bond)then
                   min_bond = dist
                   atoms(1,:) = [ is, ia ]
                   atoms(2,:) = [ js, ja ]
                end if
             end do
          end do
       end do
    end do
    bond%length = min_bond
    bond%atoms = atoms

  end function get_shortest_bond
!!!#############################################################################

  
!!!#############################################################################
!!! shares strain between two lattices
!!!#############################################################################
  subroutine share_strain(lat1,lat2,bulk_mod1,bulk_mod2,axis,lcompensate)
    implicit none
    integer :: i
    integer :: iaxis
    real(real32) :: area1,area2,delta1,delta2
    integer, dimension(3) :: abc=(/1,2,3/)
    real(real32), dimension(3) :: strain

    real(real32), intent(in) :: bulk_mod1,bulk_mod2
    real(real32), dimension(3,3), intent(inout) :: lat1,lat2

    integer, optional, intent(in) :: axis
    logical, optional, intent(in) :: lcompensate

    iaxis=3
    if(present(axis)) iaxis=axis
 
    abc=cshift(abc,3-iaxis)
    area1 = modu(cross(lat1(abc(1),:),lat1(abc(2),:)))
    area2 = modu(cross(lat2(abc(1),:),lat2(abc(2),:)))
    delta1 = - (1._real32 - area2/area1)/(1._real32 + (area2/area1)*(bulk_mod1/bulk_mod2))
    delta2 = - (1._real32 - area1/area2)/(1._real32 + (area1/area2)*(bulk_mod2/bulk_mod1))
    write(0,*) "areas", area1,area2
    write(0,*) "deltas", delta1,delta2
    write(0,*) "modulus", bulk_mod1,bulk_mod2
    do i=1,3
       if(i.eq.iaxis) cycle
       strain(:) = lat1(i,:)-lat2(i,:)
       lat1(i,:) = lat1(i,:) * (1._real32 + delta1)
       lat2(i,:) = lat1(i,:)
    end do
    
    if(present(lcompensate))then
       if(lcompensate)then
          lat1(abc(3),:) =  lat1(abc(3),:) * (1._real32 - delta1/(1._real32 + delta1))  
          lat2(abc(3),:) =  lat2(abc(3),:) * (1._real32 - delta2/(1._real32 + delta2))
       end if
    end if

  end subroutine share_strain
!!!#############################################################################

end module artemis__geom_utils
