!!!#############################################################################
!!! ARTEMIS
!!! Code written by Ned Thaddeus Taylor and Francis Huw Davies
!!! Code part of the ARTEMIS group (Hepplestone research group).
!!! Think Hepplestone, think HRG.
!!!#############################################################################
program artemis_executable
  use artemis
  use inputs
  implicit none


  integer :: i, j, unit
  integer, dimension(:), allocatable :: match_and_term_idx_list, idx_list(:)
  character(len=256) :: filepath, filename
  type(artemis_generator_type) :: generator
  type(basis_type), allocatable, dimension(:) :: structures



!!!-----------------------------------------------------------------------------
!!! set up global variables
!!!-----------------------------------------------------------------------------
  call set_global_vars()
  generator%tol_sym = tol_sym


!!!-----------------------------------------------------------------------------
!!! checks what task has been called and starts the appropriate codes
!!!-----------------------------------------------------------------------------
!!!  SEARCH  = Substitutions, Extension, Additions & Rotations for Creating Heterostructures
!!!  ASPECT  = Additions, Substitutions & Positional Editing of Crystals Tool
!!!  ARTEMIS = Ab initio Restructuring Tool  Enabling Modelling of  Interface Structures
!!!  ARTIE   = Alloying & Rotating Tool for Intermixed structure Editing ??? 
  select case(task)
  case(0) ! cell_edit/ASPECT
     write(*,'(1X,"task ",I0," set",/,1X,"Performing Cell Edits")') task
     if(lsurf_gen)then
        write(0,'(1X,"Finding terminations for lower material.")')

        call generator%set_tolerance( &
             tolerance = tolerance &
        )
        call generator%set_materials( &
             structure_lw = struc1_bas, &
             use_pricel_lw = lw_use_pricel &
        )
        call generator%set_surface_properties( &
             miller_lw = lw_mplane, &
             is_layered_lw = lw_layered, &
             require_stoichiometry_lw = lw_require_stoich, &
             vacuum_gap = vacuum_gap, &
             layer_separation_cutoff = layer_sep &
        )

        structures = generator%get_terminations(1, &
             surface = lw_surf, &
             num_layers = lw_num_layers, &
             thickness = lw_thickness, &
             orthogonalise = lortho, &
             print_termination_info = lprint_terms, &
             verbose = verbose &
        )
        filepath = "DTERMINATIONS"
        call system("mkdir -p " // trim(filepath))
        do i = 1, size(structures)
           write(filename, '(A,"/POSCAR_term",I0)') &
                trim(adjustl(filepath)), i
           open(newunit=unit, status='replace', file=trim(filename))
           call geom_write(unit, structures(i))
           close(unit)
        end do
        write(0,'(1X,"Terminations printed.",/,1X,"Exiting...")')
        stop
     end if
     call edit_structure(&
          lat=struc1_lat,bas=struc1_bas,&
          ofile=out_filename,edits=edits,&
          lnorm=lnorm_lat)

  case(1) ! interfaces/ARTEMIS/SEARCH
     write(*,'(1X,"task ",I0," set",/,1X,"Performing Interface Generation")') task
     generator%max_num_structures = max_num_structures
     generator%axis = axis
     call generator%set_tolerance( &
          tolerance = tolerance &
     )
     call generator%set_match_method( &
          method = match_method, &
          max_num_matches = max_num_matches, & ! this is maxfit/nstore
          max_num_terms = max_num_terms, &
          max_num_planes = max_num_planes, &
          compensate_normal = compensate_normal &
     )
     call generator%set_shift_method( &
          method = shift_method, &
          num_shifts = num_shifts, &
          shifts = shifts, &
          interface_depth = interface_depth, &
          separation_scale = separation_scale, &
          depth_method = depth_method, &
          bondlength_cutoff = bondlength_cutoff &
     )
     call generator%set_swap_method( &
          method = swap_method, &
          num_swaps = num_swaps, &
          swap_density = swap_density, &
          swap_depth = swap_depth, &
          swap_sigma = swap_sigma, &
          require_mirror_swaps = require_mirror_swaps &
     )
     call generator%set_materials( &
          structure_lw = struc1_bas, structure_up = struc2_bas, &
          use_pricel_lw = lw_use_pricel, use_pricel_up = up_use_pricel, &
          elastic_constants_lw = [ lw_bulk_modulus ], &
          elastic_constants_up = [ up_bulk_modulus ] &
     )
     call generator%set_surface_properties( &
          miller_lw = lw_mplane, miller_up = up_mplane, &
          is_layered_lw = lw_layered, is_layered_up = up_layered, &
          require_stoichiometry_lw = lw_require_stoich, &
          require_stoichiometry_up = up_require_stoich, &
          layer_separation_cutoff = [ lw_layer_sep, up_layer_sep ], &
          vacuum_gap = vacuum_gap &
     )
     if(.not.ludef_lw_layered) call generator%reset_is_layered_lw()
     if(.not.ludef_up_layered) call generator%reset_is_layered_up()

     !!-------------------------------------------------------------------------
     !! surface generator
     !!-------------------------------------------------------------------------
     if(lsurf_gen)then
        if(all(lw_mplane.eq.0))then
           write(*,'("No Miller plane defined for lower material.")')
           write(*,'("Skipping...")')
        else
           write(*,'(1X,"Finding terminations for lower material.")')
           structures = generator%get_terminations(1, &
                surface = lw_surf, &
                num_layers = lw_num_layers, &
                thickness = lw_thickness, &
                orthogonalise = lortho, &
                print_termination_info = lprint_terms, &
                verbose = verbose &
           )
           filepath = "DTERMINATIONS/DLW_TERMS"
           call system("mkdir -p " // trim(filepath))
           do i = 1, size(structures)
              write(filename, '(A,"/POSCAR_term",I0)') &
                   trim(adjustl(filepath)), i
              open(newunit=unit, status='replace', file=trim(filename))
              call geom_write(unit, structures(i))
              close(unit)
           end do
        end if
        if(all(up_mplane.eq.0))then
           write(*,'("No Miller plane defined for upper material.")')
           write(*,'("Skipping...")')
        else
           write(*,'(1X,"Finding terminations for upper material.")')
           structures = generator%get_terminations(2, &
                surface = up_surf, &
                num_layers = up_num_layers, &
                thickness = up_thickness, &
                orthogonalise = lortho, &
                print_termination_info = lprint_terms, &
                verbose = verbose &
           )
           filepath = "DTERMINATIONS/DUP_TERMS"
           call system("mkdir -p " // trim(filepath))
           do i = 1, size(structures)
              write(filename, '(A,"/POSCAR_term",I0)') &
                   trim(adjustl(filepath)), i
              open(newunit=unit, status='replace', file=trim(filename))
              call geom_write(unit, structures(i))
              close(unit)
           end do
        end if
        write(*,'(1X,"Terminations printed.",/,1X,"Exiting...")')
        stop
     end if
     

     !!-------------------------------------------------------------------------
     !! interface generator
     !!-------------------------------------------------------------------------
     if(irestart.eq.0)then
        call generator%generate( &
             surface_lw = lw_surf, surface_up = up_surf, &
             thickness_lw = lw_thickness, thickness_up = up_thickness, &
             num_layers_lw = lw_num_layers, num_layers_up = up_num_layers, &
             reduce_matches = reduce_matches, &
             print_lattice_match_info = lprint_matches, &
             print_termination_info = lprint_terms, &
             print_shift_info = lprint_shifts, &
             break_on_fail = break_on_fail, &
             icheck_term_pair = icheck_term_pair, &
             interface_idx = interface_idx, &
             seed = clock, &
             verbose = verbose &
        )
     else
        call generator%restart(struc1_bas)
     end if
     allocate(match_and_term_idx_list(0))
     do i = 1, generator%num_structures
        write(filepath, '(A,"/",A,I0.2)') &
             trim(adjustl(dirname)), &
             trim(adjustl(subdir_prefix)), &
             generator%structure_data(i)%match_and_term_idx
        if(all(generator%structure_data(1:i-1:1)%match_and_term_idx.ne.generator%structure_data(i)%match_and_term_idx))then
           call system("mkdir -p " // trim(filepath))
           call generator%write_match_and_term_data(i, &
                directory = trim(filepath), &
                filename = "struc_data.txt" &
           )
           match_and_term_idx_list = [ match_and_term_idx_list, generator%structure_data(i)%match_and_term_idx ]
        end if
        if(generator%structure_data(i)%shift_idx.gt.0)then
           write(filepath, '(A,"/",A,"/",A,I0.2)') &
                trim(adjustl(filepath)), trim(adjustl(shiftdir)), &
                trim(adjustl(subdir_prefix)), &
                generator%structure_data(i)%shift_idx
        end if
        if(generator%structure_data(i)%swap_idx.gt.0)then
           write(filepath, '(A,"/",A,"/",A,I0.2)') &
                trim(adjustl(filepath)), &
                trim(adjustl(swapdir)), &
                trim(adjustl(subdir_prefix)), &
                generator%structure_data(i)%swap_idx
        end if
        call system("mkdir -p " // trim(filepath))
        write(filename, '(A,"/",A)') trim(filepath), "POSCAR"
        open(newunit=unit, status='replace', file=trim(filename))
        call geom_write(unit, generator%structures(i))
        close(unit)
     end do
     ! get all indices with the same match_idx
     ! write the shift data associated with all of them
     do i = 1, size(match_and_term_idx_list)
        idx_list = pack([(j, j=1, generator%num_structures)], &
                          generator%structure_data(:)%match_and_term_idx .eq. match_and_term_idx_list(i) )
        if(size(idx_list).eq.0) cycle
        write(filepath, '(A,"/",A,I0.2,"/",A)') &
             trim(dirname), &
             trim(subdir_prefix), &
             generator%structure_data(idx_list(1))%match_and_term_idx, &
             trim(shiftdir)
        call generator%write_shift_data(idx_list, &
             directory = trim(filepath), &
             filename = "shift_data.txt" &
        )
     end do


  case(2) ! defects/ARTIE
     write(*,'(1X,"task ",I0," set",/,1X,"Performing Defect Generation")') task
     

  case default
     write(*,'(1X,"No task selected.")')
     write(*,'(1X,"Exiting code...")')
     call exit()
  end select

 

end program artemis_executable

