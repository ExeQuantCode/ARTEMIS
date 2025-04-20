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


  integer :: i, unit
  character(len=256) :: filename
  type(artemis_generator_type) :: generator
  type(basis_type), allocatable, dimension(:) :: structures



!!!-----------------------------------------------------------------------------
!!! set up global variables
!!!-----------------------------------------------------------------------------
  call set_global_vars()


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
               is_layered_lw = lw_layered &
        )

        structures = generator%get_terminations(1, &
             num_layers = lw_num_layers, &
             thickness = lw_thickness &
        )
        do i = 1, size(structures)
           write(filename, '(A,I0,A)') "term_", i, ".vasp"
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
     call generator%set_tolerance( &
          tolerance = tolerance &
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
          layer_separation_cutoff = [ lw_layer_sep, up_layer_sep ] &
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
                num_layers = lw_num_layers, &
                thickness = lw_thickness &
           )
           do i = 1, size(structures)
              write(filename, '(A,I0,A)') "lw_term_", i, ".vasp"
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
                num_layers = up_num_layers, &
                thickness = up_thickness &
           )
           do i = 1, size(structures)
              write(filename, '(A,I0,A)') "up_term_", i, ".vasp"
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
             print_lattice_match_info = lprint_matches, &
             print_termination_info = lprint_terms, &
             print_shift_info = lprint_shifts &
        )
        call generator%write_structures(directory = "DINTERFACES", prefix= "")
     else
        call generator%restart(struc1_bas)
     end if


  case(2) ! defects/ARTIE
     write(*,'(1X,"task ",I0," set",/,1X,"Performing Defect Generation")') task
     

  case default
     write(*,'(1X,"No task selected.")')
     write(*,'(1X,"Exiting code...")')
     call exit()
  end select

 

end program artemis_executable

