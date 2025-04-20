module artemis
  use artemis__geom_rw, only: basis_type, &
       geom_write, geom_read
  use artemis__structure_cache, only: &
       store_last_generated_structures, &
       retrieve_last_generated_structures
  use artemis__interface_identifier, only: intf_info_type
  use artemis__generator, only: artemis_generator_type
  implicit none


  ! allow the identify_interface procedure to be called externally

end module artemis