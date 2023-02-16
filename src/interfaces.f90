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
       shift_region,set_vacuum,transformer,shifter,reducer,&
       get_min_bulk_bond,clone_bas,bas_lat_merge,get_shortest_bond,bond_type,&
       share_strain, normalise_basis
  use mod_sym,              only: term_arr_type,confine_type,gldfnd,&
       get_terminations,get_primitive_cell
  use swapping,              only: rand_swapper
  use shifting !!! CHANGE TO SHIFTER?
  implicit none
  integer, private :: intf=0
  double precision, private, parameter :: tmp_vac = 14.D0

  
  type term_list_type
     integer :: term
     double precision :: loc
  end type term_list_type
  private :: term_list_type

  type(bulk_DON_type), dimension(2) :: bulk_DON


!!!updated  2023/02/16


contains
!!!#############################################################################
!!! Generates and prints terminations parallel to the supplied miller plane
!!!#############################################################################
  subroutine gen_terminations(lat,bas,miller_plane,axis,directory,&
       thickness,udef_layer_sep)
    implicit none
    integer :: unit
    integer :: itmp1,j,iterm,term_start,term_end,iterm_step
    integer :: old_natom,ncells,thickness_val,ntrans
    double precision :: height
    character(len=1024) :: dirname,filename,pwd
    logical :: ludef_surf,lignore
    type(bas_type) :: tmp_bas1,tmp_bas2
    type(confine_type) :: confine
    type(term_arr_type) :: term
    double precision, dimension(3,3) :: tfmat,tmp_lat1,tmp_lat2

    integer, allocatable, dimension(:,:,:) :: bas_map,t1bas_map
    double precision, allocatable, dimension(:,:) :: trans

    integer, intent(in) :: axis
    type(bas_type), intent(in) :: bas
    integer, dimension(3), intent(in) :: miller_plane
    double precision, dimension(3,3), intent(in) :: lat

    integer, optional, intent(in) :: thickness
    double precision, optional, intent(in) :: udef_layer_sep
    character(len=*), optional, intent(in) :: directory


    !! copy lattice and basis for manipulating
    call clone_bas(bas,tmp_bas1,lat,tmp_lat1)
    allocate(bas_map(tmp_bas1%nspec,maxval(tmp_bas1%spec(:)%num,dim=1),2))
    bas_map=-1


    write(6,'(1X,"Using supplied plane...")')
    tfmat=planecutter(tmp_lat1,dble(miller_plane))
    call transformer(tmp_lat1,tmp_bas1,tfmat,bas_map)
    !call err_abort_print_struc(lat,bas,"check.vasp","stop")

    !!-----------------------------------------------------------------------
    !! Finds smallest thickness of the slab and increases to ...
    !! ... user-defined thickness
    !!-----------------------------------------------------------------------
    confine%l=.false.
    confine%axis=axis
    confine%laxis=.false.
    confine%laxis(axis)=.true.
    old_natom=tmp_bas1%natom
    if(allocated(trans)) deallocate(trans)
    allocate(trans(minval(tmp_bas1%spec(:)%num+2),3))
    call gldfnd(confine,tmp_bas1,tmp_bas1,trans,ntrans)
    tfmat(:,:)=0.D0
    tfmat(1,1)=1.D0
    tfmat(2,2)=1.D0
    if(ntrans.eq.0)then
       tfmat(3,3)=1.D0
    else
       itmp1=minloc(abs(trans(:ntrans,axis)),dim=1,&
            mask=abs(trans(:ntrans,axis)).gt.1.D-3/modu(tmp_lat1(axis,:)))
       tfmat(3,:)=trans(itmp1,:)
    end if
    if(all(abs(tfmat(3,:)).lt.1.D-5)) tfmat(3,3) = 1.D0
    call transformer(tmp_lat1,tmp_bas1,tfmat,bas_map)

    !! get the terminations
    if(present(udef_layer_sep)) then
       term = get_terminations(tmp_lat1,tmp_bas1,axis,&
            lprint=.true.,layer_sep=udef_layer_sep)
    else
       term = get_terminations(tmp_lat1,tmp_bas1,axis,&
            lprint=.true.,layer_sep=layer_sep)
    end if

    !! set thickness if provided by user
    if(present(thickness))then
       thickness_val = thickness
    else
       thickness_val = 1
    end if

    !! make directory and change to that directory
    if(present(directory))then
       dirname = directory
    else
       dirname = "DTERMINATIONS"
    end if
    call system('mkdir -p '//trim(adjustl(dirname)))
    call getcwd(pwd)
    call chdir(dirname)

    !! determine tolerance for layer separations (termination tolerance)
    !! ... this is different from layer_sep
    call set_layer_tol(term)

    !! determine required extension and perform that
    call set_slab_height(tmp_lat1,tmp_bas1,bas_map,term,lw_surf,old_natom,&
         height,thickness_val,ncells,&
         term_start,term_end,iterm_step,ludef_surf,&
         dirname,"lw",lignore)

    
    !! loop over terminations and write them
    do iterm=term_start,term_end,iterm_step
       call clone_bas(tmp_bas1,tmp_bas2,tmp_lat1,tmp_lat2)
       if(allocated(t1bas_map)) deallocate(t1bas_map)
       allocate(t1bas_map,source=bas_map)
       call prepare_slab(tmp_lat2,tmp_bas2,bas_map,term,iterm,&
            thickness_val,ncells,height,ludef_surf,lw_surf(2),&
            "lw",lignore,lortho,vacuum)

       !!-----------------------------------------------------------------------
       !! Prints structure
       !!-----------------------------------------------------------------------
       unit=100+iterm
       write(filename,'("POSCAR_term",I0)') iterm
       open(unit,file=trim(filename))
       call geom_write(unit,tmp_lat2,tmp_bas2)
       close(unit)
    end do

    !! return to parent directory
    call chdir(pwd)


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
    integer :: iterm_step,jterm_step
    integer :: lw_ncells,up_ncells,istep
    integer :: lw_layered_axis,up_layered_axis
    integer :: intf_start,intf_end
    integer :: lw_term_start,lw_term_end,up_term_start,up_term_end
    double precision :: avg_min_bond
    double precision :: lw_height,up_height
    character(3) :: abc
    character(1024) :: pwd,intf_dir,dirpath,msg
    logical :: ludef_lw_surf,ludef_up_surf,lcycle
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
    !double precision, dimension(3,3) :: mtmp1,DONup_lat
    double precision, dimension(3,3) :: tfmat,slat,inlw_lat,inup_lat
    double precision, dimension(3,3) :: lw_lat,up_lat,tlw_lat,tup_lat
    integer, allocatable, dimension(:,:,:) :: lw_map,t1lw_map,t2lw_map
    integer, allocatable, dimension(:,:,:) :: up_map,t1up_map,t2up_map
    double precision, allocatable, dimension(:,:) :: trans


!!!-----------------------------------------------------------------------------
!!! determines the primitive and niggli reduced cell for each bulk
!!!-----------------------------------------------------------------------------
    write(6,*)
    if(lw_use_pricel)then
       write(6,'(1X,"Using primitive cell for lower material")')
       call get_primitive_cell(inlw_lat,inlw_bas)
    else
       write(6,'(1X,"Using supplied cell for lower material")')
       call reducer(inlw_lat,inlw_bas)
       inlw_lat=primitive_lat(inlw_lat)
    end if
    if(up_use_pricel)then
       write(6,'(1X,"Using primitive cell for upper material")')
       call get_primitive_cell(inup_lat,inup_bas)
    else
       write(6,'(1X,"Using supplied cell for upper material")')
       call reducer(inup_lat,inup_bas)
       inup_lat=primitive_lat(inup_lat)
    end if
    write(6,*)
    

    
!!!-----------------------------------------------------------------------------
!!! investigates individual bulks and their bondlengths
!!!-----------------------------------------------------------------------------
    avg_min_bond = &
         ( get_min_bulk_bond(inlw_lat,inlw_bas) + &
         get_min_bulk_bond(inup_lat,inup_bas) )/2.D0
    write(6,'(1X,"Avg min bulk bond: ",F0.3," Å")') avg_min_bond
    write(6,'(1X,"Trans-interfacial scaling factor: ",F0.3)') c_scale
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
    call getcwd(pwd)
    old_intf = -1
    intf=0
    abc="abc"
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


       !!-----------------------------------------------------------------------
       !! Sort out ladder rungs (checks whether the material is centrosymmetric)
       !!-----------------------------------------------------------------------
       !call setup_ladder(lw_lat,lw_bas,axis,lw_term)
       if(sum(lw_term%arr(:)%natom)*lw_term%nstep.ne.lw_bas%natom)then
          write(msg, '("ERROR: Number of atoms in lower layers not correct: "&
               &I0,2X,I0)') sum(lw_term%arr(:)%natom)*lw_term%nstep,lw_bas%natom
          call err_abort(trim(msg),fmtd=.true.)
       end if
       call set_layer_tol(lw_term)


       !!-----------------------------------------------------------------------
       !! Defines height of lower slab from user-defined values
       !!-----------------------------------------------------------------------
       call set_slab_height(lw_lat,lw_bas,t1lw_map,lw_term,lw_surf, old_natom,&
            lw_height,lw_thickness,lw_ncells,&
            lw_term_start,lw_term_end,iterm_step,ludef_lw_surf,&
            intf_dir,"lw",lcycle)
       if(lcycle) cycle intf_loop


       !!-----------------------------------------------------------------------
       !! Finds smallest thickness of the upper slab and increases to ...
       !! ... user-defined thickness
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


       !!-----------------------------------------------------------------------
       !! Sort out ladder rungs (checks whether the material is centrosymmetric)
       !!-----------------------------------------------------------------------
       !call setup_ladder(up_lat,up_bas,axis,up_term)
       if(sum(up_term%arr(:)%natom)*up_term%nstep.ne.up_bas%natom)then
          write(msg, '("ERROR: Number of atoms in upper layers not correct: "&
               &I0,2X,I0)') sum(up_term%arr(:)%natom)*up_term%nstep,up_bas%natom
          call err_abort(trim(msg),fmtd=.true.)
       end if
       call set_layer_tol(up_term)


       !!-----------------------------------------------------------------------
       !! Defines height of upper slab from user-defined values
       !!-----------------------------------------------------------------------
       call set_slab_height(up_lat,up_bas,t1up_map,up_term,up_surf,old_natom,&
            up_height,up_thickness,up_ncells,&
            up_term_start,up_term_end,jterm_step,ludef_up_surf,&
            intf_dir,"up",lcycle)
       if(lcycle) cycle intf_loop


       !!-----------------------------------------------------------------------
       !! Print termination plane locations
       !!-----------------------------------------------------------------------
       write(6,'(1X,"Number of unique terminations: ",I0,2X,I0)') &
            lw_term%nterm,up_term%nterm

       !!-----------------------------------------------------------------------
       !! Cycle over terminations of both materials and generates interfaces ...
       !! ... composed of all of the possible combinations of the two
       !!-----------------------------------------------------------------------
       lw_term_loop: do iterm=lw_term_start,lw_term_end,iterm_step
          call clone_bas(lw_bas,tlw_bas,lw_lat,tlw_lat)
          if(allocated(t2lw_map)) deallocate(t2lw_map)
          allocate(t2lw_map,source=t1lw_map)
          !!--------------------------------------------------------------------
          !! Shifts lower material to specified termination
          !!--------------------------------------------------------------------
          call prepare_slab(tlw_lat,tlw_bas,t2lw_map,lw_term,iterm,&
               lw_thickness,lw_ncells,lw_height,ludef_lw_surf,lw_surf(2),&
               "lw",lcycle)
          if(lcycle) cycle lw_term_loop

          
          !!--------------------------------------------------------------------
          !! Cycles over terminations of upper material
          !!--------------------------------------------------------------------
          up_term_loop: do jterm=up_term_start,up_term_end,jterm_step
             call clone_bas(up_bas,tup_bas,up_lat,tup_lat)
             if(allocated(t2up_map)) deallocate(t2up_map)
             allocate(t2up_map,source=t1up_map)
             call prepare_slab(tup_lat,tup_bas,t2up_map,up_term,jterm,&
                  up_thickness,up_ncells,up_height,ludef_up_surf,up_surf(2),&
                  "up",lcycle)
             if(lcycle) cycle up_term_loop

             
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
             !! Use the bulk moduli to determine the strain sharing
             !!-----------------------------------------------------------------
             if(lw_bulk_modulus.ne.0.E0.and.up_bulk_modulus.ne.0.E0)then
                call share_strain(tlw_lat,tup_lat,&
                     lw_bulk_modulus,up_bulk_modulus,lcompensate=.not.lc_fix)
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
             call  output_intf_data(SAV, ifit, lw_term, iterm, up_term, jterm,&
                  lw_use_pricel,up_use_pricel)


             !!-----------------------------------------------------------------
             !! Generates shifts and swaps and prints the subsequent structures
             !!-----------------------------------------------------------------
             call gen_shifts_and_swaps(slat,sbas,axis,intf_loc,avg_min_bond,&
                  ishift,nshift,&
                  iswap,swap_den,nswap,t2lw_map)

             if(intf.ge.nintf) exit intf_loop
             !call chdir(dirname)
             call chdir(intf_dir)

             if(ludef_up_surf) exit up_term_loop
          end do up_term_loop
          if(ludef_lw_surf) exit lw_term_loop
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
!!! sets the maximum height of the slab
!!!#############################################################################
  subroutine set_slab_height(lat, bas, map, term, surf, old_natom,&
       height, thickness, ncells,&
       term_start, term_end, term_step, ludef_surf,&
       intf_dir, lwup_in, lcycle)
    implicit none
    integer :: i,itmp1
    double precision :: dtmp1
    character(2) :: lwup
    character(5) :: lowerupper
    character(1024) :: msg
    double precision, dimension(3,3) :: tfmat
    double precision, allocatable, dimension(:) :: vtmp1
    type(term_list_type), allocatable, dimension(:) :: list

    integer, intent(in) :: thickness, old_natom
    integer, intent(inout) :: term_start, term_end, ncells
    integer, intent(out) :: term_step
    double precision, intent(inout) :: height
    character(2), intent(in) :: lwup_in
    character(1024), intent(in) :: intf_dir
    logical, intent(inout) :: ludef_surf
    logical, intent(out) :: lcycle
    type(bas_type), intent(inout) :: bas
    type(term_arr_type), intent(inout) :: term
    integer, dimension(2), intent(in) :: surf
    double precision, dimension(3,3), intent(inout) :: lat

    integer, allocatable, dimension(:,:,:), intent(inout) :: map
    

    !!--------------------------------------------------------------------
    !! Initialise variables
    !!--------------------------------------------------------------------
    lwup=to_lower(lwup_in)
    if(lwup.eq."lw") lowerupper="LOWER"
    if(lwup.eq."up") lowerupper="UPPER"

    lcycle = .false.
    height = 0.D0


    !!-----------------------------------------------------------------------
    !! Defines height of slab from user-defined values
    !!-----------------------------------------------------------------------
    ludef_surf = .false.
    term_start = 1
    term_end = min(term%nterm,nterm)
    if(all(surf.ne.0))then
       if(any(surf.gt.term%nterm))then
          write(msg, '(A2,"_SURFACE VALUES INVALID!\nOne or more value &
               &exceeds the maximum number of terminations in the &
               structure.\n&
               &  Supplied values: ",I0,1X,I0,"\n&
               &  Maximum allowed: ",I0)') lwup, surf, term%nterm
          call err_abort(trim(msg),fmtd=.true.)
       end if
       ludef_surf = .true.
       list = get_term_list(term)
       !! set term_start to first surface value
       term_start = surf(1)
       !! set term_end to first surface value as a user-defined surface ...
       !! ... should not be cycled over.
       !! it is just one, potentially assymetric, slab to be explored.
       term_end = surf(1)

       !! determines the maximum number of cells required
       allocate(vtmp1(size(list)))
       height = term%arr(term_start)%hmin
       do i=thickness,2,-1
          vtmp1 = list(:)%loc - height
          vtmp1 = vtmp1 - ceiling( vtmp1 - 1.D0 )
          itmp1 = minloc( vtmp1(:), dim=1,&
               mask=&
               vtmp1(:).gt.0.and.&
               list(:)%term.eq.surf(1))
          height = height + vtmp1(itmp1)
       end do
       vtmp1 = list(:)%loc - height
       !vtmp1 = vtmp1 - ceiling( vtmp1 - 1.D0 )
       where(vtmp1.lt.-1.D-5)
          vtmp1 = vtmp1 - ceiling( vtmp1 + 1.D-5 - 1.D0 )
       end where
       itmp1 = minloc( vtmp1(:), dim=1,&
            mask=&
            vtmp1(:).ge.-1.D-5.and.&
            list(:)%term.eq.surf(2))
       !!write(0,*) "temp",itmp1
       !!write(0,*) "temp",list(:)%loc
       !!write(0,*) "SURFACES",surf
       !write(0,*) "look",term%arr(term_start)%hmin, term_start
       !write(0,*) vtmp1(itmp1),itmp1
       !write(0,*) list(:)%loc
       !write(0,*) list(:)%loc-height
       !write(0,*) vtmp1
       !write(0,*) list(:)%term
       !write(0,*) "height check1", height
       height = height + vtmp1(itmp1) - term%arr(term_start)%hmin
       !write(0,*) "height check2", height

       !write(0,*) "mirror?",term%lmirror
       !! if there is no mirror, we need to remove extra layers in the cell
       !if(.not.term%lmirror)then
          ! get thickness of top/surface layer
          dtmp1 = term%arr(surf(2))%hmax - term%arr(surf(2))%hmin
          if(dtmp1.lt.-1.D-5) dtmp1 = dtmp1 + 1.D0
          height = height + dtmp1 !(1.D0 - dtmp1)
       !end if

       !write(0,*) "HEIGHT", height
       ncells = ceiling(height)
       height = height/dble(ncells)
    end if
    !write(0,*) "ncells",ncells
    !write(0,*) "height",height

    
    !!-----------------------------------------------------------------------
    !! Define termination iteration counter
    !!-----------------------------------------------------------------------
    if(term_end.lt.term_start)then
       term_step = -1
    else
       term_step = 1
    end if

    
    !!-----------------------------------------------------------------------
    !! Extend slab to user-defined thickness
    !!-----------------------------------------------------------------------
    !write(0,*) "HERE",term%nstep,thickness
    !write(0,*) thickness-1, (thickness-1)/term%nstep,int((thickness-1)/term%nstep)+1
    if(.not.ludef_surf) ncells = int((thickness-1)/term%nstep)+1
    !write(0,*) ncells
    tfmat(:,:)=0.D0
    tfmat(1,1)=1.D0
    tfmat(2,2)=1.D0
    tfmat(3,3)=ncells
    !write(0,*) "test0",ncells
    call transformer(lat,bas,tfmat,map)
    !write(0,*) "test1"
    if(mod(real(old_natom*ncells)/real(bas%natom),1.0).gt.1.D-5)then
       write(0,'(1X,"ERROR: Internal error in interfaces subroutine")')
       write(0,'(2X,"gldfnd subroutine did not reproduce a sensible &
            &primitive cell for ",A5," crystal")') lowerupper
       write(0,'(2X,"Generated ",I0," atoms, from the original ",&
            &I0," atoms")') &
            bas%natom/itmp1,old_natom
       if(ierror.eq.1)then
          call chdir(intf_dir)
          call err_abort_print_struc(lat,bas,&
               "broken_primitive.vasp",&
               "As IPRINT = 1, code is now exiting...")
       end if
       write(0,'(2X,"Skipping this lattice match...")')
       lcycle=.true.
    end if

    
    !!-----------------------------------------------------------------------
    !! Readjust termination plane locations
    !! ... i.e. divide all termination values by the number of cells
    !!-----------------------------------------------------------------------
    term%arr(:)%hmin = term%arr(:)%hmin/dble(ncells)
    term%arr(:)%hmax = term%arr(:)%hmax/dble(ncells)
    term%tol = term%tol/dble(ncells)


  end subroutine set_slab_height
!!!#############################################################################


!!!#############################################################################
!!! Set the tolerance for layer definitions
!!!#############################################################################
  subroutine set_layer_tol(term)
    implicit none
    integer :: i
    double precision :: dtmp1

    type(term_arr_type), intent(inout) :: term
    

    do i=1,term%nterm
       if(i.eq.1)then
          dtmp1 = abs(term%arr(i)%hmin - &
               (term%arr(term%nterm)%hmax+term%arr(i)%ladder(term%nstep)-1.D0)&
               )/4.D0
       else
          dtmp1 = abs(term%arr(i)%hmin-term%arr(i-1)%hmax)/4.D0
       end if
       if(dtmp1.lt.term%tol)then
          term%tol = dtmp1
       end if
    end do

    !! add the tolerances to the edges of the layers
    !! this ensures that all atoms in the layers are captured
    term%arr(:)%hmin = term%arr(:)%hmin - term%tol
    term%arr(:)%hmax = term%arr(:)%hmax + term%tol


  end subroutine set_layer_tol
!!!#############################################################################


!!!#############################################################################
!!! Prepares lattice and basis to specified termination
!!!#############################################################################
!!! Supply a supercell that can be cut down to the size of the slab ...
!!! ... i.e. the input structure must be larger or equal to the desired output
  subroutine prepare_slab(lat, bas, map, term, iterm, thickness, ncells, &
       height, ludef_surf, udef_top_iterm, lwup_in, lcycle, &
       ludef_ortho, udef_vacuum)
    implicit none
    integer :: j, j_start, istep, natom_check
    double precision :: vacuum, dtmp1
    character(2) :: lwup
    character(5) :: lowerupper
    character(1024) :: msg
    logical :: lortho
    integer, dimension(3) :: abc=(/1,2,3/)
    double precision, dimension(3,3) :: tfmat
    integer, allocatable, dimension(:) :: iterm_list

    integer, intent(in) :: iterm, udef_top_iterm, thickness, ncells
    double precision, intent(in) :: height
    character(2), intent(in) :: lwup_in
    logical, intent(in) :: ludef_surf
    logical, intent(out) :: lcycle
    type(bas_type), intent(inout) :: bas
    type(term_arr_type), intent(in) :: term
    double precision, dimension(3,3), intent(inout) :: lat

    integer, allocatable, dimension(:,:,:), intent(inout) :: map
    logical, optional, intent(in) :: ludef_ortho
    double precision, optional, intent(in) :: udef_vacuum

    !!--------------------------------------------------------------------
    !! Initialise variables
    !!--------------------------------------------------------------------
    lwup=to_lower(lwup_in)
    if(lwup.eq."lw") lowerupper="LOWER"
    if(lwup.eq."up") lowerupper="UPPER"
    lcycle = .false.
    dtmp1=0.D0
    tfmat=0.D0
    istep = thickness - (ncells-1)*term%nstep
    natom_check = bas%natom

    if(present(ludef_ortho))then
       lortho = ludef_ortho
    else
       lortho = .true.
    end if

    if(present(udef_vacuum))then
       vacuum = udef_vacuum
    else
       vacuum = tmp_vac
    end if


    !!--------------------------------------------------------------------
    !! Set up list for checking expected number of atoms
    !!--------------------------------------------------------------------
    allocate(iterm_list(term%nterm))
    do j=1,term%nterm
       iterm_list(j) = j
    end do
    iterm_list=cshift(iterm_list,term%nterm-iterm+1)
    if(ludef_surf)then
       j_start = udef_top_iterm - iterm + 1
       if(j_start.le.0) j_start = j_start + term%nterm
       j_start = j_start + 1 !+ (istep-1)*term%nterm/term%nstep
    else
       !! handle ladder steps that are equivalent
       j_start = 2 !+ (istep-1)*term%nterm/term%nstep
    end if


    !!--------------------------------------------------------------------
    !! Shift lower material to specified termination
    !!--------------------------------------------------------------------
    call shifter(bas,term%axis,-term%arr(iterm)%hmin,.true.)
    !open(20,file="test.vasp")
    !call geom_write(20,lat,bas)
    !close(20)


    !!--------------------------------------------------------------------
    !! Determine cell reduction to specified termination
    !!--------------------------------------------------------------------
    !write(0,*) "LUDEF_SURF?", ludef_surf
    do j=1,3
       tfmat(j,j)=1.D0
       if(j.eq.term%axis)then
          if(ludef_surf)then
             tfmat(j,j) = height !+ term%tol*2.D0
          else!if(term%lmirror)then
             if(istep.ne.0)then
                dtmp1 = (ncells-1) + term%arr(iterm)%ladder(istep)
                dtmp1 = dtmp1/(ncells)
                tfmat(j,j) = dtmp1 !+ term%tol*2.D0
                tfmat(j,j) = tfmat(j,j) + &
                     (term%arr(iterm)%hmax - term%arr(iterm)%hmin)
             end if
          !else
          !   tfmat(j,j) = tfmat(j,j) + (&
          !        term%arr(iterm)%hmax - &
          !        term%arr(iterm)%hmin) + term%tol*2.D0
          end if
       end if
    end do


    !!--------------------------------------------------------------------
    !! Apply transformation and shift cell back to bottom of layer
    !! ... i.e. account for the tolerance that has been added to layer ...
    !! ... hmin and hmax
    !!--------------------------------------------------------------------
    call transformer(lat,bas,tfmat,map)
    call shifter(bas,term%axis,-term%tol/tfmat(term%axis,term%axis),.true.)


    !!--------------------------------------------------------------------
    !! Check number of atoms is expected
    !!--------------------------------------------------------------------
    if(term%nterm.gt.1.or.term%nstep.gt.1)then
       do j=1,max(0,term%nstep-istep),1
          natom_check = natom_check - sum(term%arr(:)%natom)
       end do
       do j=j_start,term%nterm,1
          natom_check = natom_check - term%arr(iterm_list(j))%natom
       end do
    end if
    if(bas%natom.ne.natom_check)then
       write(msg, '("NUMBER OF ATOMS IN '//to_upper(lowerupper)//' SLAB! &
            &Expected ",I0," but generated ",I0," instead")') &
            natom_check,bas%natom
       if(tfmat(term%axis,term%axis).gt.1.D0)then
          write(0,'("THE TRANSFORMATION IS GREATER THAN ONE ",F0.9)') &
               tfmat(term%axis,term%axis)
       end if
       !call err_abort(trim(msg),fmtd=.true.)
       call err_abort_print_struc(lat,bas,lwup//"_term.vasp",&
            trim(msg),.true.)
       lcycle = .true.
    end if


    !!--------------------------------------------------------------------
    !! Apply slab_cuber to orthogonalise lower material
    !!--------------------------------------------------------------------
    call set_vacuum(lat,bas,term%axis,1.D0-term%tol/tfmat(term%axis,term%axis),vacuum)
    !call err_abort_print_struc(lat,bas,"check.vasp","stop")
    abc=cshift(abc,3-term%axis)
    if(lortho)then
       ortho_check: do j=1,2
          if(abs(dot_product(lat(abc(j),:),lat(axis,:))).gt.1.D-5)then
             call ortho_axis(lat,bas,term%axis)
             exit ortho_check
          end if
       end do ortho_check
    end if
    call normalise_basis(bas,dtmp=0.9999D0,lfloor=.true.)


  end subroutine prepare_slab
!!!#############################################################################


!!!#############################################################################
!!! write structure data in each structure directory
!!!#############################################################################
  subroutine output_intf_data(SAV, ifit, lw_term, ilw_term, up_term, iup_term, lw_pricel,up_pricel)
    implicit none
    integer :: unit

    integer, intent(in) :: ifit, ilw_term, iup_term
    logical, intent(in) :: lw_pricel,up_pricel
    type(term_arr_type), intent(in) :: lw_term, up_term
    type(latmatch_type), intent(in) :: SAV


    
    unit=99
    open(unit=unit, file="struc_dat.txt")
    write(unit,'("Lower material primitive cell used: ",L1)') lw_pricel
    write(unit,'("Upper material primitive cell used: ",L1)') lw_pricel
    write(unit,*)
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
