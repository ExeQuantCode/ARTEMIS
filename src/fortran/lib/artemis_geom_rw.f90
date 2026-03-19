module artemis__geom_rw
  !! Thin re-export wrapper for backwards compatibility.
  !! All types and IO routines are now provided by the atomstruc library.
  use atomstruc, only: &
       igeom_input, igeom_output, &
       basis_type, species_type, &
       geom_read, geom_write, &
       get_element_properties
  implicit none

  private

  public :: igeom_input, igeom_output
  public :: basis_type, species_type
  public :: geom_read, geom_write
  public :: get_element_properties

end module artemis__geom_rw

