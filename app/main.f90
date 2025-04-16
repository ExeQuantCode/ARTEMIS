!!!#############################################################################
!!! ARTEMIS
!!! Code written by Ned Thaddeus Taylor and Francis Huw Davies
!!! Code part of the ARTEMIS group (Hepplestone research group).
!!! Think Hepplestone, think HRG.
!!!#############################################################################
program artemis_executable
  use artemis
  implicit none


  type(artemis_termination_generator_type) :: term_gen
  type(artemis_interface_generator_type) :: intf_gen



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
     write(6,'(1X,"task ",I0," set",/,1X,"Performing Cell Edits")') task
     if(lsurf_gen)then
        write(0,'(1X,"Finding terminations for lower material.")')
        term_gen%layer_separation_cutoff = layer_sep
        call term_gen%generate(struc1_bas,lw_mplane,axis,&
             num_layers = lw_num_layers, &
             thickness = lw_thickness &
        )
        call term_gen%write_structures(directory = "DTERMINATIONS", prefix= "term_")
        write(0,'(1X,"Terminations printed.",/,1X,"Exiting...")')
        stop
     end if
     call edit_structure(&
          lat=struc1_lat,bas=struc1_bas,&
          ofile=out_filename,edits=edits,&
          lnorm=lnorm_lat)

  case(1) ! interfaces/ARTEMIS/SEARCH
     write(6,'(1X,"task ",I0," set",/,1X,"Performing Interface Generation")') task

     !!-------------------------------------------------------------------------
     !! surface generator
     !!-------------------------------------------------------------------------
     if(lsurf_gen)then
        
        call system('mkdir -p DTERMINATIONS')
        call chdir("DTERMINATIONS")
        
        if(all(lw_mplane.eq.0))then
           write(6,'("No Miller plane defined for lower material.")')
           write(6,'("Skipping...")')
        else
           write(6,'(1X,"Finding terminations for lower material.")')
           term_gen%layer_separation_cutoff = lw_layer_sep
           call term_gen%generate(struc1_bas,lw_mplane,axis,&
                num_layers = lw_num_layers, &
                thickness = lw_thickness &
           )
           call term_gen%write_structures(directory = "DTERMINATIONS", prefix= "lw_")
        end if
        if(all(up_mplane.eq.0))then
           write(6,'("No Miller plane defined for upper material.")')
           write(6,'("Skipping...")')
        else
           write(6,'(1X,"Finding terminations for upper material.")')
           term_gen%layer_separation_cutoff = up_layer_sep
           call term_gen%generate(struc2_bas,up_mplane,axis,&
                num_layers = up_num_layers, &
                thickness = up_thickness &
           )
           call term_gen%write_structures(directory = "DTERMINATIONS", prefix= "up_")
        end if
        write(6,'(1X,"Terminations printed.",/,1X,"Exiting...")')
        stop
     end if
     

     !!-------------------------------------------------------------------------
     !! interface generator
     !!-------------------------------------------------------------------------
     if(irestart.eq.0)then
        call gen_interfaces(tolerance,&
             struc1_lat,struc2_lat,&
             struc1_bas,struc2_bas)
     else
        call gen_interfaces_restart(struc1_lat,struc1_bas)
     end if


  case(2) ! defects/ARTIE
     write(6,'(1X,"task ",I0," set",/,1X,"Performing Defect Generation")') task
     

  case default
     write(6,'(1X,"No task selected.")')
     write(6,'(1X,"Exiting code...")')
     call exit()
  end select

 

end program artemis_executable

