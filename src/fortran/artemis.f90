module artemis
  !! Public API module for the ARTEMIS library.
  !!
  !! Re-exports the main types and procedures used by consumers of the
  !! ARTEMIS interface-generation library.
  use atomstruc, only: basis_type, &
       geom_write, geom_read
  use artemis__structure_cache, only: &
       store_last_generated_structures, &
       retrieve_last_generated_structures
  use artemis__interface_identifier, only: intf_info_type
  use artemis__interface_translations, only: get_interface_translations
  use artemis__generator, only: artemis_generator_type
  use artemis__io_utils, only: artemis__suppress_warnings
  implicit none


  ! allow the identify_interface procedure to be called externally

end module artemis
