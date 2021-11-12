!!!#############################################################################
!!! ARTEMIS
!!! Code written by Ned Thaddeus Taylor and Francis Huw Davies
!!! Code part of the ARTEMIS group (Hepplestone research group).
!!! Think Hepplestone, think HRG.
!!!#############################################################################
program artemis
  use inputs
  use interface_subroutines
  implicit none


!!!updated 2021/11/12
 

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
        call gen_terminations(struc1_lat,struc1_bas,lw_mplane,axis,&
             thickness=lw_thickness)
        write(0,'(1X,"Terminations printed.",/,1X,"Exiting...")')
        stop
     end if
     call edit_structure(&
          lat=struc1_lat,bas=struc1_bas,&
          ofile=out_filename,edits=edits)

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
           call gen_terminations(struc1_lat,struc1_bas,lw_mplane,axis,&
                directory="DLW_TERMS",thickness=lw_thickness,udef_layer_sep=lw_layer_sep)
        end if
        if(all(up_mplane.eq.0))then
           write(6,'("No Miller plane defined for upper material.")')
           write(6,'("Skipping...")')
        else
           write(6,'(1X,"Finding terminations for upper material.")')
           call gen_terminations(struc2_lat,struc2_bas,up_mplane,axis,&
                directory="DUP_TERMS",thickness=up_thickness,udef_layer_sep=up_layer_sep)
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

 

end program artemis

