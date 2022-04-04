!!!#############################################################################
!!! INTERFACES CARD SUBROUTINES
!!! Code written by Ned Thaddeus Taylor and Isiah Edward Mikel Rudkin
!!! Code part of the ARTEMIS group (Hepplestone research group).
!!! Think Hepplestone, think HRG.
!!!#############################################################################
module interface_subroutines
  use io
  use misc_linalg,          only: uvec,modu,get_area,inverse
  use inputs
  use interface_identifier, only: intf_info_type,&
       get_interface,get_layered_axis,gen_DON
  use edit_geom,            only: planecutter,primitive_lat,ortho_axis,&
       shift_region,set_vacuum,transformer,shifter,&
       get_min_bulk_bond,clone_bas,bas_lat_merge,get_shortest_bond,bond_type
  use mod_sym,              only: term_arr_type,confine_type,gldfnd,&
       get_terminations,print_terminations,setup_ladder,get_primitive_cell
  use swapping,              only: rand_swapper
  use shifting !!! CHANGE TO SHIFTER?
  implicit none
  integer, private :: intf=0

  
  type term_list_type
     integer :: term
     double precision :: loc
  end type term_list_type
  private :: term_list_type

  type(bulk_DON_type), dimension(2) :: bulk_DON


!!!updated  2022/01/17


contains
!!!#############################################################################
!!! Generates and prints terminations parallel to the supplied miller plane
!!!#############################################################################
  subroutine gen_terminations(lat,bas,miller_plane,axis,directory,&
       thickness,udef_layer_sep)
    implicit none
    character(len=200) :: dirname
    type(term_arr_type) :: term
    double precision, dimension(3,3) :: tfmat

    integer, intent(in) :: axis
    type(bas_type), intent(in) :: bas
    integer, dimension(3), intent(in) :: miller_plane
    double precision, dimension(3,3), intent(in) :: lat

    integer, optional, intent(in) :: thickness
    double precision, optional, intent(in) :: udef_layer_sep
    character(len=*), optional, intent(in) :: directory


    write(6,'(1X,"Using supplied plane...")')
    tfmat=planecutter(lat,dble(miller_plane))
    call transformer(lat,bas,tfmat)
    !call err_abort_print_struc(lat,bas,"check.vasp","stop")


    if(present(udef_layer_sep)) then
       term = get_terminations(lat,bas,axis,lprint=.true.,layer_sep=udef_layer_sep)
    else
       term = get_terminations(lat,bas,axis,lprint=.true.,layer_sep=layer_sep)
    end if
    term%arr(:)%hmin = term%arr(:)%hmin - 1.D-8
    call setup_ladder(lat,bas,axis,term)


    if(present(directory))then
       dirname = directory
    else
       dirname = "DTERMINATIONS"
    end if
    if(present(thickness))then
       call print_terminations(term,lat,bas,trim(dirname),thickness,vacuum,&
            lortho = lortho)
    else
       call print_terminations(term,lat,bas,trim(dirname),lortho = lortho)
    end if


    return
  end subroutine gen_terminations
!!!#############################################################################


!!!#############################################################################
!!! generate interfaces
!!!#############################################################################
  subroutine gen_interfaces_restart(lat,bas)
    implicit none
    integer :: is,ia,js,ja
    double precision :: dtmp1,min_bond,min_bond1,min_bond2
    type(bas_type) :: bas
    type(intf_info_type) :: intf
    double precision, dimension(3) :: vtmp1
    double precision, dimension(3,3) :: lat


    call system('mkdir -p '//trim(adjustl(dirname)))
    call chdir(dirname)

    min_bond1=huge(0.D0)
    min_bond2=huge(0.D0)
    if(any(udef_intf_loc.lt.0.D0))then
       if(ludef_axis)then
          intf=get_interface(lat,bas,axis)
       else
          intf=get_interface(lat,bas)
       end if
       intf%loc=intf%loc/modu(lat(intf%axis,:))
       write(6,*) "interface axis:",intf%axis
       write(6,*) "interface loc:",intf%loc
       !! write interface location to a file for user to refer back to
       open(unit=10,file="interface_location.dat")
       write(10,'(1X,"AXIS = ",I0)') intf%axis
       write(10,'(1X,"INTF_LOC = ",2(2X,F9.6))') intf%loc
       close(10)
    else
       intf%axis = axis
       intf%loc = udef_intf_loc
    end if
    specloop1: do is=1,bas%nspec
       atomloop1: do ia=1,bas%spec(is)%num

          specloop2: do js=1,bas%nspec
             atomloop2: do ja=1,bas%spec(js)%num
                if(is.eq.js.and.ia.eq.ja) cycle atomloop2
                if( &
                     ( bas%spec(is)%atom(ia,intf%axis).gt.intf%loc(1).and.&
                     bas%spec(is)%atom(ia,intf%axis).lt.intf%loc(2) ).and.&
                     ( bas%spec(js)%atom(ja,intf%axis).gt.intf%loc(1).and.&
                     bas%spec(js)%atom(ja,intf%axis).lt.intf%loc(2) ) )then
                   vtmp1 = (bas%spec(is)%atom(ia,:3)-bas%spec(js)%atom(ja,:3))
                   vtmp1 = matmul(vtmp1,lat)
                   dtmp1 = modu(vtmp1)
                   if(dtmp1.lt.min_bond1) min_bond1 = dtmp1
                elseif( &
                     ( bas%spec(is)%atom(ia,intf%axis).lt.intf%loc(1).or.&
                     bas%spec(is)%atom(ia,intf%axis).gt.intf%loc(2) ).and.&
                     ( bas%spec(js)%atom(ja,intf%axis).lt.intf%loc(1).or.&
                     bas%spec(js)%atom(ja,intf%axis).gt.intf%loc(2) ) )then
                   vtmp1 = (bas%spec(is)%atom(ia,:3)-bas%spec(js)%atom(ja,:3))
                   vtmp1 = matmul(vtmp1,lat)
                   dtmp1 = modu(vtmp1)
                   if(dtmp1.lt.min_bond2) min_bond2 = dtmp1
                end if

             end do atomloop2
          end do specloop2
    
       end do atomloop1
    end do specloop1

    min_bond = ( min_bond1 + min_bond2 )/2.D0
    write(6,'(1X,"Avg min bulk bond: ",F0.3," Å")') min_bond
    write(6,'(1X,"Trans-interfacial scaling factor:",F0.3)') c_scale
    call gen_shifts_and_swaps(lat,bas,intf%axis,intf%loc,min_bond,&
         ishift,nshift,&
         iswap,swap_den,nswap)


  end subroutine gen_interfaces_restart
!!!#############################################################################


!!!#############################################################################
!!! generate interfaces
!!!#############################################################################
  subroutine gen_interfaces(tolerance,inlw_lat,inup_lat,inlw_bas,inup_bas)
    implicit none
    integer :: i,j,iterm,jterm,ntrans,ifit,iunique,old_natom,itmp1,old_intf
    integer :: lw_ncells,up_ncells,istep
    integer :: lw_layered_axis,up_layered_axis
    integer :: intf_start,intf_end
    integer :: lw_term_start,lw_term_end,up_term_start,up_term_end
    integer :: natom_check
    double precision :: avg_min_bond,tmp_vac,dtmp1
    double precision :: lw_height,up_height
    character(3) :: abc
    character(1024) :: pwd,intf_dir,dirpath,msg
    logical :: ludef_lw_surf,ludef_up_surf
    type(bas_type) :: sbas
    type(bas_type) :: inlw_bas,inup_bas
    type(bas_type) :: lw_bas,up_bas,tlw_bas,tup_bas
    type(tol_type) :: tolerance
    type(confine_type) :: confine
    type(latmatch_type) :: SAV
    type(term_arr_type) :: lw_term,up_term
    integer, dimension(3) :: ivtmp1
    double precision, dimension(2) :: intf_loc
    double precision, dimension(3) :: init_offset=[0.D0,0.D0,2.D0]
    double precision, dimension(3,3) :: mtmp1,DONup_lat
    double precision, dimension(3,3) :: tfmat,slat,inlw_lat,inup_lat
    double precision, dimension(3,3) :: lw_lat,up_lat,tlw_lat,tup_lat
    double precision, allocatable, dimension(:) :: vtmp1
    integer, allocatable, dimension(:,:,:) :: lw_map,up_map,t1lw_map,t1up_map,t2lw_map,t2up_map
    double precision, allocatable, dimension(:,:) :: trans
    type(term_list_type), allocatable, dimension(:) :: lw_list,up_list

    
!!!-----------------------------------------------------------------------------
!!! determines the primitive and niggli reduced cell for each bulk
!!!-----------------------------------------------------------------------------
    if(lw_use_pricel)then
       call get_primitive_cell(inlw_lat,inlw_bas)
    end if
    if(up_use_pricel)then
       call get_primitive_cell(inup_lat,inup_bas)
    end if

    
!!!-----------------------------------------------------------------------------
!!! investigates individual bulks and their bondlengths
!!!-----------------------------------------------------------------------------
    avg_min_bond = &
         ( get_min_bulk_bond(inlw_lat,inlw_bas) + &
         get_min_bulk_bond(inup_lat,inup_bas) )/2.D0
    write(6,'(1X,"Avg min bulk bond: ",F0.3," Å")') avg_min_bond
    write(6,'(1X,"Trans-interfacial scaling factor:",F0.3)') c_scale
    if(ishift.eq.-1) nshift=1
    

!!!-----------------------------------------------------------------------------
!!! gets bulk DONs, if ISHIFT = 4
!!!-----------------------------------------------------------------------------
    allocate(lw_map(inlw_bas%nspec,maxval(inlw_bas%spec(:)%num,dim=1),2))
    allocate(up_map(inup_bas%nspec,maxval(inup_bas%spec(:)%num,dim=1),2))    
    if(ishift.eq.4.or.ishift.eq.0)then
       lw_map=0
       bulk_DON(1)%spec=gen_DON(inlw_lat,inlw_bas,&
            dist_max=max_bondlength,&
            scale_dist=.false.,&
            norm=.true.)
       if(all(abs(bulk_DON(1)%spec(1)%atom(:,:)).lt.1.D0))then
          call err_abort("ISSUE WITH THE LOWER BULK DON!!!")
       end if
       open(unit=13,file="lw_DON.dat")
       do j=1,1000
          write(13,*) &
               (j-1)*max_bondlength/1000,&
               bulk_DON(1)%spec(1)%atom(1,j)
       end do
       close(13)
       !call exit()
       up_map=0
       bulk_DON(2)%spec=gen_DON(inup_lat,inup_bas,&
            dist_max=max_bondlength,&
            scale_dist=.false.,&
            norm=.true.)
       if(all(abs(bulk_DON(2)%spec(1)%atom(:,:)).lt.1.D0))then
          call err_abort("ISSUE WITH THE UPPER BULK DON!!!")
       end if
    else
       lw_map=-1
       up_map=-1       
    end if


!!!-----------------------------------------------------------------------------
!!! checks whether system appears layered
!!!-----------------------------------------------------------------------------
    lw_layered_axis=get_layered_axis(inlw_lat,inlw_bas)
    if(.not.lw_layered.and.lw_layered_axis.gt.0)then
       ivtmp1=0
       ivtmp1(lw_layered_axis)=1
       if(ludef_lw_layered)then
          write(msg,'("Lower crystal appears layered along axis ",I0,"\n&
               &Partial layer terminations will be generated\n&
               &We suggest using LW_MILLER =",3(1X,I1))') lw_layered_axis,ivtmp1
          call print_warning(trim(msg))
       else
          write(msg,'("Lower crystal has been identified as layered\nalong",3(1X,I1),"\n&
               &Confining crystal to this plane and\nstoichiometric terminations.\n&
               &If you don''t want this, set\nLW_LAYERED = .FALSE.")') &
               ivtmp1
          call print_warning(trim(msg))
          lw_mplane=ivtmp1
          lw_layered=.true.
       end if
    elseif(lw_layered.and.lw_layered_axis.gt.0.and.all(lw_mplane.eq.0))then
       lw_mplane(lw_layered_axis)=1
    end if

    up_layered_axis=get_layered_axis(inup_lat,inup_bas)
    if(.not.up_layered.and.up_layered_axis.gt.0)then
       ivtmp1=0
       ivtmp1(up_layered_axis)=1
       if(ludef_up_layered)then
          write(msg,'("Upper crystal appears layered along axis ",I0,"\n&
               &Partial layer terminations will be generated\n&
               &We suggest using UP_MILLER =",3(1X,I1))') up_layered_axis,ivtmp1
          call print_warning(trim(msg))
       else
          write(msg,'("Upper crystal has been identified as layered\nalong",3(1X,I1),"\n&
               &Confining crystal to this plane and\nstoichiometric terminations.\n&
               &If you don''t want this, set\nUP_LAYERED = .FALSE.")') &
               ivtmp1
          call print_warning(trim(msg))
          up_mplane=ivtmp1
          up_layered=.true.
       end if
    elseif(up_layered.and.up_layered_axis.gt.0.and.all(up_mplane.eq.0))then
       up_mplane(up_layered_axis)=1
    end if


!!!-----------------------------------------------------------------------------
!!! Finds and stores the best matches between the materials
!!!-----------------------------------------------------------------------------
    tmp_vac=14.D0
    call getcwd(pwd)
    old_intf = -1
    intf=0
    abc="abc"
    inlw_lat=primitive_lat(inlw_lat)
    inup_lat=primitive_lat(inup_lat)
    if(any(lw_mplane.ne.0))then
       if(imatch.ne.0)then
          abc="ab"
          tfmat=planecutter(inlw_lat,dble(lw_mplane))
          call transformer(inlw_lat,inlw_bas,tfmat,lw_map)
          SAV=get_best_match(&
               tolerance,&
               inlw_lat,inup_lat,&
               inlw_bas,inup_bas,&
               trim(abc),"abc",lprint_matches,ierror,imatch=imatch)
       elseif(any(up_mplane.ne.0))then
          SAV=get_best_match(&
               tolerance,&
               inlw_lat,inup_lat,&
               inlw_bas,inup_bas,&
               trim(abc),"abc",lprint_matches,ierror,imatch=imatch,&
               plane1=lw_mplane,plane2=up_mplane,nmiller=nmiller)
       else
          SAV=get_best_match(&
               tolerance,&
               inlw_lat,inup_lat,&
               inlw_bas,inup_bas,&
               trim(abc),"abc",lprint_matches,ierror,imatch=imatch,&
               plane1=lw_mplane,nmiller=nmiller)
       end if
    elseif(any(up_mplane.ne.0))then
       SAV=get_best_match(&
            tolerance,&
            inlw_lat,inup_lat,&
            inlw_bas,inup_bas,&
            trim(abc),"abc",lprint_matches,ierror,imatch=imatch,&
            plane2=up_mplane,nmiller=nmiller)
    else
       SAV=get_best_match(&
            tolerance,&
            inlw_lat,inup_lat,&
            inlw_bas,inup_bas,&
            trim(abc),"abc",lprint_matches,ierror,imatch=imatch,&
            nmiller=nmiller)
    end if
    if(min(tolerance%nstore,SAV%nfit).eq.0)then
       write(0,'("No matches found.")')
       write(0,'("Exiting...")')
       call exit()
    else
       write(0,'(1X,"Number of matches found: ",I0)')&
            min(tolerance%nstore,SAV%nfit)
    end if
    write(6,'(1X,"Maximum number of generated interfaces will be: ",I0)')&
         nterm*nshift*tolerance%nstore
    if(.not.lgen_interfaces)then
       write(0,'(1X,"Told not to generate interfaces, just find matches.")')
       write(0,'("Exiting...")')
       call exit()
    end if

       
!!!-----------------------------------------------------------------------------
!!! Saves current directory and moves to new directory
!!!-----------------------------------------------------------------------------
    call system('mkdir -p '//trim(adjustl(dirname)))
    call chdir(dirname)
    call getcwd(intf_dir)

    if(iintf.gt.0)then
       intf_start=iintf
       intf_end=iintf
       write(6,'(1X,"Generating only interfaces for match ",I0)') iintf
    else
       intf_start=1
       intf_end=min(tolerance%nstore,SAV%nfit)
    end if
    iunique=0
!!!-----------------------------------------------------------------------------
!!! Applies the best match transformations
!!!-----------------------------------------------------------------------------
    intf_loop: do ifit=intf_start,intf_end
       write(6,'("Fit number: ",I0)') ifit
       call clone_bas(inlw_bas,lw_bas,inlw_lat,lw_lat)
       call clone_bas(inup_bas,up_bas,inup_lat,up_lat)
       if(allocated(t1lw_map)) deallocate(t1lw_map)
       if(allocated(t1up_map)) deallocate(t1up_map)
       allocate(t1lw_map,source=lw_map)
       allocate(t1up_map,source=up_map)
       

       !!-----------------------------------------------------------------------
       !! Applies the best match transformations
       !!-----------------------------------------------------------------------
       call transformer(lw_lat,lw_bas,dble(SAV%tf1(ifit,:,:)),t1lw_map)
       call transformer(up_lat,up_bas,dble(SAV%tf2(ifit,:,:)),t1up_map)


       !!-----------------------------------------------------------------------
       !! Determines the cell change for the upper lattice to get the new DON
       !!-----------------------------------------------------------------------
       if(ishift.eq.4)then
          !! Issue with using this method when large deformations result in large
          !! angle changes. REMOVING IT FOR NOW AND RETURNING TO CALCULATING DONS
          !! FOR THE SUPERCELL.
          t1up_map=0 !TEMPORARY TO USE SUPERCELL DONS.
          !do i=1,2
          !   mtmp1(i,:) = &
          !        ( modu(lw_lat(i,:)) )*uvec(up_lat(i,:))
          !end do
          !mtmp1(3,:) = up_lat(3,:)
          !DONup_lat = matmul(mtmp1,inverse(dble(SAV%tf2(ifit,:,:))))
          !if(ierror.eq.1)then
          !   write(0,*) "#####################################"
          !   write(0,*) "ifit", ifit
          !   write(0,*) "undeformed lattice"
          !   write(0,'(3(2X,F6.2))') (mtmp1(i,:),i=1,3)
          !   write(0,*)
          !   write(0,*) "deformed lattice"
          !   write(0,'(3(2X,F8.4))') (DONup_lat(i,:),i=1,3)
          !   write(0,*)
          !end if
          deallocate(bulk_DON(2)%spec)
          bulk_DON(2)%spec=gen_DON(up_lat,up_bas,&
               dist_max=max_bondlength,&
               scale_dist=.false.,&
               norm=.true.)
          !call err_abort_print_struc(DONup_lat,inup_bas,"bulk_up_term.vasp",&
          !     "",.false.)
       end if


       !!-----------------------------------------------------------------------
       !! Finds smallest thickness of the lower slab and increases to ...
       !!user-defined thickness
       !! SHOULD MAKE IT LATER MAKE DIFFERENT SETS OF THICKNESSES
       !!-----------------------------------------------------------------------
       confine%l=.false.
       confine%axis=axis
       confine%laxis=.false.
       confine%laxis(axis)=.true.
       old_natom=lw_bas%natom
       if(allocated(trans)) deallocate(trans)
       allocate(trans(minval(lw_bas%spec(:)%num+2),3))
       call gldfnd(confine,lw_bas,lw_bas,trans,ntrans)
       tfmat(:,:)=0.D0
       tfmat(1,1)=1.D0
       tfmat(2,2)=1.D0
       if(ntrans.eq.0)then
          tfmat(3,3)=1.D0
       else
          itmp1=minloc(abs(trans(:ntrans,axis)),dim=1,&
               mask=abs(trans(:ntrans,axis)).gt.1.D-3/modu(lw_lat(axis,:)))
          tfmat(3,:)=trans(itmp1,:)
       end if
       if(all(abs(tfmat(3,:)).lt.1.D-5)) tfmat(3,3) = 1.D0
       call transformer(lw_lat,lw_bas,tfmat,t1lw_map)

       
       !!-----------------------------------------------------------------------
       !! Finds all terminations parallel to the surface plane
       !!-----------------------------------------------------------------------
       if(allocated(lw_term%arr)) deallocate(lw_term%arr)
       lw_term=get_terminations(lw_lat,lw_bas,axis,&
            lprint=lprint_terms,layer_sep=lw_layer_sep)
       lw_term%arr(:)%hmin = lw_term%arr(:)%hmin - 1.D-8


       !!-----------------------------------------------------------------------
       !! Sort out ladder rungs (checks whether the material is centrosymmetric)
       !!-----------------------------------------------------------------------
       call setup_ladder(lw_lat,lw_bas,axis,lw_term)


       !!-----------------------------------------------------------------------
       !! Defines height of lower slab from user-defined values
       !!-----------------------------------------------------------------------
       ludef_lw_surf = .false.
       lw_term_start = 1
       lw_term_end = min(lw_term%nterm,nterm)
       if(all(lw_surf.ne.0))then
          ludef_lw_surf = .true.
          lw_list=get_term_list(lw_term)
          !do iterm=1,size(lw_list)
          !   write(0,*) lw_list(iterm)
          !end do
          lw_term_start = lw_surf(1)
          lw_term_end = lw_surf(1)

          if(allocated(vtmp1)) deallocate(vtmp1)
          allocate(vtmp1(size(lw_list)))
          lw_height = lw_term%arr(lw_term_start)%hmin
          do i=lw_thickness,2,-1
             vtmp1 = lw_list(:)%loc - lw_height
             vtmp1 = vtmp1 - ceiling( vtmp1 - 1.D0 )
             itmp1 = minloc( vtmp1(:), dim=1,&
                  mask=&
                  vtmp1(:).gt.0.and.&
                  lw_list(:)%term.eq.lw_surf(1))
             lw_height = lw_height + vtmp1(itmp1)
          end do
          vtmp1 = lw_list(:)%loc - lw_height
          vtmp1 = vtmp1 - ceiling( vtmp1 - 1.D0 )
          itmp1 = minloc( vtmp1(:), dim=1,&
               mask=&
               vtmp1(:).gt.0.and.&
               lw_list(:)%term.eq.lw_surf(2))
          lw_height = lw_height + vtmp1(itmp1) - lw_term%arr(lw_term_start)%hmin

          lw_ncells = ceiling(lw_height)
          lw_height = lw_height/dble(lw_ncells)

          if(.not.lw_term%lmirror)then
             dtmp1 = lw_term%arr(lw_surf(2))%hmax - lw_term%arr(lw_surf(2))%hmin
             if(dtmp1.lt.0.D0) dtmp1 = dtmp1 + 1.D0
             lw_height = lw_height + dtmp1
          end if
       else
          lw_ncells = int((lw_thickness-1)/lw_term%nstep) + 1
       end if


       !!-----------------------------------------------------------------------
       !! Extends lower slab to user-defined thickness
       !!-----------------------------------------------------------------------
       tfmat(:,:)=0.D0
       tfmat(1,1)=1.D0
       tfmat(2,2)=1.D0
       tfmat(3,3)=lw_ncells
       call transformer(lw_lat,lw_bas,tfmat,t1lw_map)
       if(mod(real(old_natom*lw_ncells)/real(lw_bas%natom),1.0).gt.1.D-5)then
          write(0,'(1X,"ERROR: Internal error in interfaces subroutine")')
          write(0,'(2X,"gldfnd subroutine did not reproduce a sensible &
               &primitive cell for lower crystal")')
          write(0,'(2X,"Generated ",I0," atoms, from the original ",&
               &I0," atoms")') &
               lw_bas%natom/itmp1,old_natom
          if(ierror.eq.1)then
             call chdir(intf_dir)
             call err_abort_print_struc(lw_lat,lw_bas,&
                  "broken_primitive.vasp",&
                  "As IPRINT = 1, code is now exiting...")
          end if
          write(0,'(2X,"Skipping this lattice match...")')
          cycle intf_loop
       end if


       !!-----------------------------------------------------------------------
       !! Finds smallest thickness of the upper slab and increases to ...
       !! user-defined thickness
       !! SHOULD MAKE IT LATER MAKE DIFFERENT SETS OF THICKNESSES
       !!-----------------------------------------------------------------------
       old_natom=up_bas%natom
       deallocate(trans)
       allocate(trans(minval(up_bas%spec(:)%num+2),3))
       call gldfnd(confine,up_bas,up_bas,trans,ntrans)
       tfmat(:,:)=0.D0
       tfmat(1,1)=1.D0
       tfmat(2,2)=1.D0
       if(ntrans.eq.0)then
          tfmat(3,3)=1.D0
       else
          itmp1=minloc(abs(trans(:ntrans,axis)),dim=1,&
               mask=abs(trans(:ntrans,axis)).gt.1.D-3/modu(lw_lat(axis,:)))
          tfmat(3,:)=trans(itmp1,:)
       end if
       if(all(abs(tfmat(3,:)).lt.1.D-5)) tfmat(3,3) = 1.D0
       call transformer(up_lat,up_bas,tfmat,t1up_map)

       
       !!-----------------------------------------------------------------------
       !! Finds all up_lat unique terminations parallel to the surface plane
       !!-----------------------------------------------------------------------
       if(allocated(up_term%arr)) deallocate(up_term%arr)
       up_term=get_terminations(up_lat,up_bas,axis,&
            lprint=lprint_terms,layer_sep=up_layer_sep)
       up_term%arr(:)%hmin = up_term%arr(:)%hmin - 1.D-8


       !!-----------------------------------------------------------------------
       !! Sort out ladder rungs (checks whether the material is centrosymmetric)
       !!-----------------------------------------------------------------------
       call setup_ladder(up_lat,up_bas,axis,up_term)


       !!-----------------------------------------------------------------------
       !! Defines height of upper slab from user-defined values
       !!-----------------------------------------------------------------------
       ludef_up_surf = .false.
       up_term_start = 1
       up_term_end = min(up_term%nterm,nterm)
       if(all(up_surf.ne.0))then
          ludef_up_surf = .true.
          up_list=get_term_list(up_term)
          !do iterm=1,size(up_list)
          !   write(0,*) up_list(iterm)
          !end do
          up_term_start = up_surf(1)
          up_term_end = up_surf(1)

          if(allocated(vtmp1)) deallocate(vtmp1)
          allocate(vtmp1(size(up_list)))
          up_height = up_term%arr(up_term_start)%hmin
          do i=up_thickness,2,-1
             vtmp1 = up_list(:)%loc - up_height
             vtmp1 = vtmp1 - ceiling( vtmp1 - 1.D0 )
             itmp1 = minloc( vtmp1(:), dim=1,&
                  mask=&
                  vtmp1(:).gt.0.and.&
                  up_list(:)%term.eq.up_surf(1))
             up_height = up_height + vtmp1(itmp1)
          end do
          vtmp1 = up_list(:)%loc - up_height
          vtmp1 = vtmp1 - ceiling( vtmp1 - 1.D0 )
          itmp1 = minloc( vtmp1(:), dim=1,&
               mask=&
               vtmp1(:).gt.0.and.&
               up_list(:)%term.eq.up_surf(2))
          up_height = up_height + vtmp1(itmp1) - up_term%arr(up_term_start)%hmin

          if(.not.up_term%lmirror)then
             dtmp1 = up_term%arr(up_surf(2))%hmax - up_term%arr(up_surf(2))%hmin
             if(dtmp1.lt.0.D0) dtmp1 = dtmp1 + 1.D0
             up_height = up_height + dtmp1
          end if

          up_ncells = ceiling(up_height)
          up_height = up_height/dble(up_ncells)
       else
          up_ncells = int((up_thickness-1)/up_term%nstep) + 1
       end if


       !!-----------------------------------------------------------------------
       !! Extends upper slab to user-defined thickness
       !!-----------------------------------------------------------------------
       tfmat(:,:)=0.D0
       tfmat(1,1)=1.D0
       tfmat(2,2)=1.D0
       up_ncells = int((up_thickness-1)/up_term%nstep)+1
       tfmat(3,3)=up_ncells
       call transformer(up_lat,up_bas,tfmat,t1up_map)
       if(mod(real(old_natom*up_ncells)/real(up_bas%natom),1.0).gt.1.D-5)then
          write(0,'(1X,"ERROR: Internal error in interfaces subroutine")')
          write(0,'(2X,"gldfnd subroutine did not reproduce a sensible &
               &primitive cell for upper crystal")')
          write(0,'(2X,"Generated ",I0," atoms, from the original ",&
               &I0," atoms")') &
               up_bas%natom/up_thickness,old_natom
          write(0,'(2X,"Skipping this lattice match...")')
          if(ierror.eq.1)then
             call chdir(intf_dir)
             call err_abort_print_struc(up_lat,up_bas,&
                  "broken_primitive.vasp",&
                  "As IPRINT = 1, code is now exiting...")
          end if
          cycle intf_loop
       end if


       !!-----------------------------------------------------------------------
       !! Readjust termination plane locations and print
       !!-----------------------------------------------------------------------
       lw_term%arr(:)%hmin = lw_term%arr(:)%hmin/dble(lw_ncells)
       lw_term%arr(:)%hmax = lw_term%arr(:)%hmax/dble(lw_ncells)
       !lw_term%arr(:)%add = lw_term%arr(:)%add/dble(lw_ncells)
       lw_term%tol = lw_term%tol/dble(lw_ncells)
       up_term%arr(:)%hmin = up_term%arr(:)%hmin/dble(up_ncells)
       up_term%arr(:)%hmax = up_term%arr(:)%hmax/dble(up_ncells)
       !up_term%arr(:)%add = up_term%arr(:)%add/dble(up_ncells)
       up_term%tol = up_term%tol/dble(up_ncells)
       write(6,'(1X,"Number of unique terminations: ",I0,2X,I0)') &
            lw_term%nterm,up_term%nterm

       !!-----------------------------------------------------------------------
       !! Cycle over terminations of both materials and generates interfaces ...
       !! ... composed of all of the possible combinations of the two
       !!-----------------------------------------------------------------------
       lw_term_loop: do iterm=lw_term_start,lw_term_end
          call clone_bas(lw_bas,tlw_bas,lw_lat,tlw_lat)
          if(allocated(t2lw_map)) deallocate(t2lw_map)
          allocate(t2lw_map,source=t1lw_map)
          !!--------------------------------------------------------------------
          !! Shifts lower material to specified termination
          !!--------------------------------------------------------------------
          call shifter(tlw_bas,lw_term%axis,-lw_term%arr(iterm)%hmin,.true.)
          tfmat=0.D0
          istep = 0
          do j=1,3
             tfmat(j,j)=1.D0
             if(j.eq.lw_term%axis)then
                if(ludef_lw_surf)then
                   tfmat(j,j) = lw_height + lw_term%tol/8.D0
                elseif(lw_term%lmirror)then
                   istep = lw_thickness - (lw_ncells-1)*lw_term%nstep
                   if(istep.ne.0)then
                      dtmp1 = (lw_ncells-1) + lw_term%arr(iterm)%ladder(istep)
                      dtmp1 = dtmp1/(lw_ncells)
                      tfmat(j,j) = dtmp1 + lw_term%tol/8.D0
                      tfmat(j,j) = tfmat(j,j) + &
                           (lw_term%arr(iterm)%hmax - lw_term%arr(iterm)%hmin)
                   end if
                   !dtmp1 = dble(lw_thickness-1)+lw_term%arr(iterm)%add
                   !if(dtmp1.eq.0.D0) dtmp1=1.D0
                   !tfmat(j,j) = dtmp1 + lw_term%tol/8.D0 !tfmat(j,j)+(&
                   !!tfmat(j,j) = tfmat(j,j) + (&
                   !!     !lw_term%arr(iterm)%hmax-&
                   !!     !lw_term%arr(iterm)%hmin+&
                   !!     lw_term%arr(iterm)%add)+lw_term%tol/8.D0
                else
                   tfmat(j,j)=tfmat(j,j)+(&
                        lw_term%arr(iterm)%hmax-&
                        lw_term%arr(iterm)%hmin)+lw_term%tol/8.D0
                end if
             end if
          end do
          natom_check = lw_bas%natom
          if(iterm.lt.lw_term_end)then
             do j=1,max(1,up_term%nstep-istep)
                natom_check = natom_check - sum(lw_term%arr(:iterm)%natom)
             end do
          end if
          call transformer(tlw_lat,tlw_bas,tfmat,t2lw_map)
          if(tlw_bas%natom.ne.natom_check)then
             write(msg, '("NUMBER OF ATOMS IN LOWER SLAB! &
                  &Expected ",I0," but generated ",I0," instead")') &
                  natom_check,tlw_bas%natom
             call err_abort(trim(msg),fmtd=.true.)
                call err_abort_print_struc(tlw_lat,tlw_bas,"lw_term.vasp",&
                     trim(msg),.true.)
          end if
          !!--------------------------------------------------------------------
          !! Applied slab_cuber to orthogonalise lower material
          !!-------------------------------------------------------------------
          ortho_check1: do j=1,2 !! MAKE THIS GLOBAL, NOT JUST FOR AXIS 3!!!
             if(abs(dot_product(tlw_lat(j,:),tlw_lat(3,:))).gt.1.D-5)then
                call ortho_axis(tlw_lat,tlw_bas,axis)
                exit ortho_check1
             end if
          end do ortho_check1
          call set_vacuum(tlw_lat,tlw_bas,lw_term%axis,1.D0,tmp_vac)
          !call err_abort_print_struc(tlw_lat,tlw_bas,"check.vasp","stop")


          
          !!--------------------------------------------------------------------
          !! Cycles over terminations of upper material
          !!--------------------------------------------------------------------
          up_term_loop: do jterm=up_term_start,up_term_end
             call clone_bas(up_bas,tup_bas,up_lat,tup_lat)
             if(allocated(t2up_map)) deallocate(t2up_map)
             allocate(t2up_map,source=t1up_map)
             tfmat=0.D0
             istep=0
             !!-----------------------------------------------------------------
             !! Shifts upper material to specified termination
             !!-----------------------------------------------------------------
             call shifter(tup_bas,up_term%axis,-up_term%arr(jterm)%hmin,.true.)
             !! NEED TO SORT OUT LAYER MIN AND MAX, ISSUES OF ROUNDING !!
             !! COULD USE THE CHECK NEAREST ATOM CODE AND USE THAT VALUE !!
             do j=1,3
                tfmat(j,j)=1.D0
                if(j.eq.up_term%axis)then
                   if(ludef_up_surf)then
                      tfmat(j,j) = up_height + up_term%tol/8.D0
                   elseif(up_term%lmirror)then
                      istep = up_thickness - (up_ncells-1)*up_term%nstep
                      if(istep.ne.0)then
                         dtmp1 = (up_ncells-1) + up_term%arr(jterm)%ladder(istep)
                         dtmp1 = dtmp1/(up_ncells)
                         tfmat(j,j) = dtmp1 + up_term%tol/8.D0
                         tfmat(j,j) = tfmat(j,j) + &
                              (up_term%arr(jterm)%hmax - up_term%arr(jterm)%hmin)
                      end if
                      !dtmp1 = dble(up_thickness-1)+up_term%arr(jterm)%add
                      !if(dtmp1.eq.0.D0) dtmp1=1.D0
                      !tfmat(j,j) = dtmp1 + up_term%tol/8.D0 !tfmat(j,j)-(&
                      !!tfmat(j,j) = tfmat(j,j) + (&
                      !!     !up_term%arr(jterm)%hmax-&
                      !!     !up_term%arr(jterm)%hmin+&
                      !!     up_term%arr(jterm)%add)+up_term%tol/8.D0
                      !!     !up_mirror%loc)!+up_term%tol/4.D0
                   else
                      tfmat(j,j)=tfmat(j,j)+(&
                           up_term%arr(jterm)%hmax-&
                           up_term%arr(jterm)%hmin)+up_term%tol/8.D0
                   end if
                end if
             end do
             natom_check = up_bas%natom
             if(jterm.lt.up_term_end)then
                do j=1,max(1,up_term%nstep-istep)
                   natom_check = natom_check - sum(up_term%arr(:jterm)%natom)
                end do
             end if
             call transformer(tup_lat,tup_bas,tfmat,t2up_map)
             if(tup_bas%natom.ne.natom_check)then
                write(msg, '("NUMBER OF ATOMS IN UPPER SLAB! &
                     &Expected ",I0," but generated ",I0," instead")') &
                     natom_check,tup_bas%natom
                !call err_abort(trim(msg),fmtd=.true.)
                call err_abort_print_struc(tup_lat,tup_bas,"up_term.vasp",&
                     trim(msg),.true.)
             end if
             !!-----------------------------------------------------------------
             !! Applied slab_cuber to orthogonalise upper material
             !!-----------------------------------------------------------------
             ortho_check2: do j=1,2 !! MAKE THIS GLOBAL, NOT JUST FOR AXIS 3!!!
                if(abs(dot_product(tup_lat(j,:),tup_lat(3,:))).gt.1.D-5)then
                   call ortho_axis(tup_lat,tup_bas,axis)
                   exit ortho_check2
                end if
             end do ortho_check2
             call ortho_axis(tup_lat,tup_bas,axis)
             call set_vacuum(tup_lat,tup_bas,up_term%axis,1.D0,tmp_vac)
 
             
             !!-----------------------------------------------------------------
             !! Checks stoichiometry
             !!-----------------------------------------------------------------
             if(tlw_bas%nspec.ne.inlw_bas%nspec.or.any(&
                  (inlw_bas%spec(1)%num*tlw_bas%spec(:)%num)&
                  /tlw_bas%spec(1)%num.ne.inlw_bas%spec(:)%num))then
                write(6,'("WARNING: This lower surface termination is not &
                     &stoichiometric")')
                if(lw_layered)then
                   write(6,'(2X,"As lower structure is layered, stoichiometric &
                        &surfaces are required.")')
                   write(6,'(2X,"Skipping this termination...")')
                   cycle lw_term_loop
                end if
             end if
             if(tup_bas%nspec.ne.inup_bas%nspec.or.any(&
                  (inup_bas%spec(1)%num*tup_bas%spec(:)%num)&
                  /tup_bas%spec(1)%num.ne.inup_bas%spec(:)%num))then
                write(6,'("WARNING: This upper surface termination is not &
                     &stoichiometric")')
                if(up_layered)then
                   write(6,'(2X,"As upper structure is layered, stoichiometric &
                        &surfaces are required.")')
                   write(6,'(2X,"Skipping this termination...")')
                   cycle up_term_loop
                end if
             end if


             !!-----------------------------------------------------------------
             !! Merge the two bases and lattices and define the interface loc
             !!-----------------------------------------------------------------
             call bas_lat_merge(&
                  slat,sbas,&
                  tlw_lat,tup_lat,&
                  tlw_bas,tup_bas,axis,init_offset(:),&
                  t2lw_map,t2up_map)
             intf_loc(1) = ( modu(tlw_lat(axis,:)) + 0.5D0*init_offset(axis) - &
                  tmp_vac)/modu(slat(axis,:))
             intf_loc(2) = ( modu(tlw_lat(axis,:)) + modu(tup_lat(axis,:)) + &
                  1.5D0*init_offset(axis) - 2.D0*tmp_vac )/modu(slat(axis,:))
             if(ierror.ge.1)then
                write(0,*) "interface:",intf_loc
                if(ierror.eq.1.and.iunique.eq.icheck_intf-1)then
                   call chdir(intf_dir)
                   call err_abort_print_struc(tlw_lat,tlw_bas,"lw_term.vasp",&
                        "",.false.)
                   call err_abort_print_struc(tup_lat,tup_bas,"up_term.vasp",&
                        "As IPRINT = 1 and ICHECK has been set, &
                        &code is now exiting...")
                elseif(ierror.eq.2.and.iunique.eq.icheck_intf-1)then
                   call chdir(intf_dir)
                   call err_abort_print_struc(slat,sbas,"test_intf.vasp",&
                        "As IPRINT = 2 and ICHECK has been set, &
                        &code is now exiting...")
                end if
             end if


             !!-----------------------------------------------------------------
             !! Saves current directory and moves to new directory
             !!-----------------------------------------------------------------
             if(intf.gt.old_intf)then
                iunique=iunique+1
                if(ishift.gt.0.and.nshift.gt.1) &
                     write(6,'(1X,"Generating shifts for unique interface ",&
                     &I0,":")') iunique
                write(dirpath,'(A,I0.2)') trim(adjustl(subdir_prefix)),iunique
                call system('mkdir -p '//trim(adjustl(dirpath)))
             else
                write(dirpath,'(A,I0.2)') trim(adjustl(subdir_prefix)),iunique
             end if
             call chdir(dirpath)
             old_intf = intf

             
             !!-----------------------------------------------------------------
             !! Writes information of current match to file in save directory
             !!-----------------------------------------------------------------
             call  output_intf_data(SAV, ifit, lw_term, iterm, up_term, jterm)


             !!-----------------------------------------------------------------
             !! Generates shifts and swaps and prints the subsequent structures
             !!-----------------------------------------------------------------
             call gen_shifts_and_swaps(slat,sbas,axis,intf_loc,avg_min_bond,&
                  ishift,nshift,&
                  iswap,swap_den,nswap,t2lw_map)

             if(intf.ge.nintf) exit intf_loop
             !call chdir(dirname)
             call chdir(intf_dir)

          end do up_term_loop
       end do lw_term_loop
       !!-----------------------------------------------------------------------
       !! Returns to working directory
       !!-----------------------------------------------------------------------
       call chdir(intf_dir)

    end do intf_loop

    call chdir(pwd)


    return
  end subroutine gen_interfaces
!!!#############################################################################


!!!#############################################################################
!!! Takes input interface structure and generates a set of shifts and swaps.
!!! Prints these new structures to POSCARs.
!!!#############################################################################
!!! ISWAP METHOD NOT YET SET UP
  subroutine gen_shifts_and_swaps(lat,bas,axis,intf_loc,bond,&
       ishift,nshift,&
       iswap,swap_den,nswap,&
       map)
    implicit none
    integer :: shift_unit=10
    integer :: ounit,iaxis,k,l
    integer :: ngen_swaps,nswaps_per_cell
    double precision :: dtmp1
    type(bas_type) :: tbas
    type(bond_type) :: min_bond
    character(1024) :: filename,dirpath,pwd1,pwd2,msg
    integer, dimension(3) :: abc
    double precision, dimension(2) :: intf_loc
    double precision, dimension(3) :: toffset
    double precision, dimension(3,3) :: tlat
    type(bas_type), allocatable, dimension(:) :: bas_arr
    double precision, allocatable, dimension(:,:) :: output_shifts

    integer, intent(in) :: axis
    integer, intent(in) :: nshift,nswap
    integer, intent(in) :: ishift,iswap
    double precision, intent(in) :: bond,swap_den
    type(bas_type), intent(in) :: bas
    double precision, dimension(3,3), intent(in) :: lat

    integer, dimension(:,:,:), optional, intent(in) :: map


!!!-----------------------------------------------------------------------------
!!! Sets up shift axis
!!!-----------------------------------------------------------------------------
    abc = [ 1, 2, 3 ]
    abc = cshift(abc,axis)


!!!-----------------------------------------------------------------------------
!!! Sets up and moves to appropriate directories
!!!-----------------------------------------------------------------------------
    call getcwd(pwd1)
    if(ishift.gt.0.or.nshift.gt.1)then
       call system('mkdir -p '//trim(adjustl(shiftdir)))
       call chdir(shiftdir)
    end if
    call getcwd(pwd2)
    open(unit=shift_unit,file="shift_vals.txt")
    write(shift_unit,&
         '("# interface_num    shift (a,b,c) units=(direct,direct,Å)")')


!!!-----------------------------------------------------------------------------
!!! Generates sets of shifts based on shift version
!!!-----------------------------------------------------------------------------
    if(ishift.eq.0) allocate(output_shifts(nshift,3))
    select case(ishift)
    case(1)
       output_shifts(1,:3)=0.D0
       do k=2,nshift
          do iaxis=1,2
             call random_number(output_shifts(k,iaxis))
          end do
       end do
    case(2)
       output_shifts = get_fit_shifts(&
            lat=lat,bas=bas,&
            bond=bond,&
            axis=axis,&
            intf_loc=intf_loc,&
            depth=intf_depth,&
            nstore=nshift)
    case(3)
       output_shifts = get_descriptive_shifts(&
            lat=lat,bas=bas,&
            bond=bond,&
            axis=axis,&
            intf_loc=intf_loc,&
            depth=intf_depth,c_scale=c_scale,&
            nstore=nshift,lprint=lprint_shifts)
    case(4)
       if(present(map))then
          output_shifts = get_shifts_DON(&
               lat=lat,bas=bas,&
               axis=axis,&
               intf_loc=intf_loc,&
               nstore=nshift,c_scale=c_scale,offset=offset(1,:3),&
               lprint=lprint_shifts,bulk_DON=bulk_DON,bulk_map=map,&
               max_bondlength=max_bondlength)
       else
          output_shifts = get_shifts_DON(&
               lat=lat,bas=bas,&
               axis=axis,&
               intf_loc=intf_loc,&
               nstore=nshift,c_scale=c_scale,offset=offset(1,:3),&
               lprint=lprint_shifts,&
               max_bondlength=max_bondlength)
       end if
       if(size(output_shifts(:,1)).eq.0)then
          write(0,'(2X,"No shifts were identified with ISHIFT = 4 for this lattice match")')
          write(0,'(2X,"We suggest increasing MBOND_MAXLEN to find shifts")')
          write(0,'("Skipping interface...")')
          return
       end if
    case default
      ! nshift=1 !!! SORT THIS OUT !!! RESET NSHIFT DUE TO ISHIFT
       if(.not.allocated(output_shifts)) allocate(output_shifts(1,3))
       output_shifts(:,:) = offset
       do iaxis=1,2
          output_shifts(1,iaxis) = output_shifts(1,iaxis)!/modu(lat(iaxis,:))
       end do
    end select
    if(ishift.gt.0)then
       output_shifts(:,axis) = output_shifts(:,axis)*modu(lat(axis,:))
    end if


!!!-----------------------------------------------------------------------------
!!! Prints number of shifts to terminal
!!!-----------------------------------------------------------------------------
    write(6,'(3X,"Number of unique shifts structures: ",I0)') nshift


!!!-----------------------------------------------------------------------------
!!! Determines number of swaps across the interface
!!!-----------------------------------------------------------------------------
    nswaps_per_cell=nint(swap_den*get_area(lat(abc(1),:),lat(abc(2),:)))
    if(iswap.ne.0)then
       write(6,&
            '(" Generating ",I0," swaps per structure ")') nswaps_per_cell
    end if


!!!-----------------------------------------------------------------------------
!!! Prints each unique shift structure
!!!-----------------------------------------------------------------------------
    shift_loop: do k=1,nshift
       call clone_bas(bas,tbas,lat,tlat)
       toffset=output_shifts(k,:3)
       do iaxis=1,2
          call shift_region(tbas,axis,&
               intf_loc(1),intf_loc(2),&
               shift_axis=iaxis,shift=toffset(iaxis),renorm=.true.)
       end do
       dtmp1=modu(tlat(axis,:))
       call set_vacuum(&
            lat=tlat,bas=tbas,&
            axis=axis,loc=maxval(intf_loc(:)),&
            vac=toffset(axis))
       dtmp1=minval(intf_loc(:))*dtmp1/modu(tlat(axis,:))
       call set_vacuum(&
            lat=tlat,bas=tbas,&
            axis=axis,loc=dtmp1,&
            vac=toffset(axis))
       min_bond = get_shortest_bond(tlat,tbas)
       if(min_bond%length.le.1.5D0)then
          write(msg,'("Smallest bond in the interface structure is\nless than 1.5 Å.")')
          call print_warning(trim(msg))
          write(6,'(2X,"bond length: ",F9.6)') min_bond%length
          write(6,'(2X,"atom 1:",I4,2X,I4)') min_bond%atoms(1,:)
          write(6,'(2X,"atom 2:",I4,2X,I4)') min_bond%atoms(2,:)
       end if


       !!-----------------------------------------------------------------------
       !! prints shift vector to shift_vals.txt
       !!-----------------------------------------------------------------------
       write(shift_unit,'(2X,I0.2,15X,"(",2(" ",F9.6,", ")," ",F9.6," )")') &
            k,toffset(:)


       !!-----------------------------------------------------------------------
       !! Merges lower and upper materials
       !! Writes interfaces to output directories
       !!-----------------------------------------------------------------------
       intf=intf+1
       ounit=100+intf
       if(ishift.gt.0.or.nshift.gt.1)then
          write(dirpath,'(A,I0.2)') trim(adjustl(subdir_prefix)),k
          call system('mkdir -p '//trim(adjustl(dirpath)))
          write(filename,'(A,"/",A)') trim(adjustl(dirpath)),trim(out_filename)
       else
          filename = trim(out_filename)
       end if
       write(6,'(2X,"Writing interface ",I0,"...")') intf
       open(unit=ounit,file=trim(adjustl(filename)))
       call geom_write(ounit,tlat,tbas)
       close(ounit)
       if(intf.ge.nintf) return


       !!-----------------------------------------------------------------------
       !! Performs swaps within the shifted structures if requested
       !!-----------------------------------------------------------------------
       if_swap: if(iswap.ne.0)then
          bas_arr = rand_swapper(tlat,tbas,axis,swap_depth,&
               nswaps_per_cell,nswap,intf_loc,iswap,seed,sigma=swap_sigma,&
               require_mirror=lswap_mirror)
          ngen_swaps = nswap
          LOOPswaps: do l=1,nswap
             if (bas_arr(l)%nspec.eq.0) then
                ngen_swaps = l - 1
                exit LOOPswaps
             end if
          end do LOOPswaps
          if(ngen_swaps.eq.0)then
             exit if_swap
          end if
          call chdir(dirpath)
          call system('mkdir -p '//trim(adjustl(swapdir)))
          call chdir(swapdir)
          write(6,'(3X,"Number of unique swap structures: ",I0)') ngen_swaps
          do l=1,ngen_swaps
             write(dirpath,'(A,I0.2)') trim(adjustl(subdir_prefix)),l
             call system('mkdir -p '//trim(adjustl(dirpath)))
             write(filename,'(A,"/",A)') &
                  trim(adjustl(dirpath)),trim(out_filename)
             ounit=100+l
             write(6,'(3X,"Writing swap ",I0,"...")') l
             open(unit=ounit,file=trim(adjustl(filename)))
             call geom_write(ounit,tlat,bas_arr(l))
             close(ounit)
          end do
          deallocate(bas_arr)
          call chdir(pwd2)
       end if if_swap


    end do shift_loop
    call chdir(pwd1)
    close(unit=shift_unit)


  end subroutine gen_shifts_and_swaps
!!!#############################################################################


!!!#############################################################################
!!! changes terminations to one long list of the top surfaces of the crystal
!!!#############################################################################
  function get_term_list(term) result(list)
    implicit none
    integer :: i,j,itmp1,nlist,loc
    type(term_arr_type), intent(in) :: term
    
    type(term_list_type) :: tmp_element
    type(term_list_type), allocatable, dimension(:) :: list


    if(.not.allocated(term%arr(1)%ladder))then
       nlist=term%nterm
       allocate(list(nlist))
       list(:)%loc = term%arr(:)%hmin
       do i=1,term%nterm
          list(i)%term = i
       end do
    else
       nlist = term%nstep*term%nterm
       allocate(list(nlist))
       itmp1=0
       do i=1,term%nterm
          do j=1,term%nstep
             itmp1=itmp1+1
             list(itmp1)%loc = term%arr(i)%hmin+term%arr(i)%ladder(j)
             list(itmp1)%loc = list(itmp1)%loc - &
                  ceiling( list(itmp1)%loc - 1.D0 )
             list(itmp1)%term=i
          end do
       end do
    end if

    !! sort the list now
    do i=1,nlist
       loc=minloc(list(i:nlist)%loc,dim=1)+i-1
       tmp_element=list(i)
       list(i)=list(loc)
       list(loc)=tmp_element
    end do



  end function get_term_list
!!!#############################################################################

  
!!!#############################################################################
!!! write structure data in each structure directory
!!!#############################################################################
  subroutine output_intf_data(SAV, ifit, lw_term, ilw_term, up_term, iup_term)
    implicit none
    integer :: unit

    integer, intent(in) :: ifit, ilw_term, iup_term
    type(term_arr_type), intent(in) :: lw_term, up_term
    type(latmatch_type), intent(in) :: SAV

    
    unit=99
    open(unit=unit, file="struc_dat.txt")
    write(unit,'("Lattice match:")')
    write(unit,'((1X,3(3X,A1),3X,3(3X,A1)),3(/,2X,3(I3," "),3X,3(I3," ")))') &
         SAV%abc,SAV%abc,&
         SAV%tf1(ifit,1,1:3),SAV%tf2(ifit,1,1:3),&
         SAV%tf1(ifit,2,1:3),SAV%tf2(ifit,2,1:3),&
         SAV%tf1(ifit,3,1:3),SAV%tf2(ifit,3,1:3)
    write(unit,'(" vector mismatch (%) = ",F0.9)') SAV%tol(ifit,1)
    write(unit,'(" angle mismatch (°)  = ",F0.9)') SAV%tol(ifit,2)
    write(unit,'(" area mismatch (%)   = ",F0.9)') SAV%tol(ifit,3)
    write(unit,*)
    write(unit,'(" Lower crystal Miller plane: ",3(I3," "))') SAV%tf1(ifit,3,1:3)
    write(unit,'(" Lower termination")')
    write(unit,'(1X,"Term.",3X,"Min layer loc",3X,"Max layer loc",3X,"no. atoms")')
    write(unit,'(1X,I3,8X,F7.5,9X,F7.5,8X,I3)') &
            ilw_term,lw_term%arr(ilw_term)%hmin,lw_term%arr(ilw_term)%hmax,lw_term%arr(ilw_term)%natom
    write(unit,*)
    write(unit,'(" Upper crystal Miller plane: ",3(I3," "))') SAV%tf2(ifit,3,1:3)
    write(unit,'(" Upper termination")')
    write(unit,'(1X,"Term.",3X,"Min layer loc",3X,"Max layer loc",3X,"no. atoms")')
    write(unit,'(1X,I3,8X,F7.5,9X,F7.5,8X,I3)') &
            iup_term,up_term%arr(iup_term)%hmin,up_term%arr(iup_term)%hmax,up_term%arr(iup_term)%natom
    write(unit,*)
    close(unit)
    
    return
  end subroutine output_intf_data
!!!#############################################################################


end module interface_subroutines
