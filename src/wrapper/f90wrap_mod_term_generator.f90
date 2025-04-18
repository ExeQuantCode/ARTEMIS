! Module artemis__termination_generator defined in file ../src/fortran/lib/mod_term_generator.f90

subroutine f90wrap_artemis_termination_generator_type__get__layer_sepace78(this, f90wrap_layer_separation_cutoff)
    use artemis__termination_generator, only: artemis_termination_generator_type
    implicit none
    type artemis_termination_generator_type_ptr_type
        type(artemis_termination_generator_type), pointer :: p => NULL()
    end type artemis_termination_generator_type_ptr_type
    integer, intent(in)   :: this(2)
    type(artemis_termination_generator_type_ptr_type) :: this_ptr
    real(4), intent(out) :: f90wrap_layer_separation_cutoff
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_layer_separation_cutoff = this_ptr%p%layer_separation_cutoff
end subroutine f90wrap_artemis_termination_generator_type__get__layer_sepace78

subroutine f90wrap_artemis_termination_generator_type__set__layer_sepae7ef(this, f90wrap_layer_separation_cutoff)
    use artemis__termination_generator, only: artemis_termination_generator_type
    implicit none
    type artemis_termination_generator_type_ptr_type
        type(artemis_termination_generator_type), pointer :: p => NULL()
    end type artemis_termination_generator_type_ptr_type
    integer, intent(in)   :: this(2)
    type(artemis_termination_generator_type_ptr_type) :: this_ptr
    real(4), intent(in) :: f90wrap_layer_separation_cutoff
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%layer_separation_cutoff = f90wrap_layer_separation_cutoff
end subroutine f90wrap_artemis_termination_generator_type__set__layer_sepae7ef

subroutine f90wrap_term_gen__artemis_termination293d(this)
    use artemis__termination_generator, only: artemis_termination_generator_type
    implicit none
    
    type artemis_termination_generator_type_ptr_type
        type(artemis_termination_generator_type), pointer :: p => NULL()
    end type artemis_termination_generator_type_ptr_type
    type(artemis_termination_generator_type_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_term_gen__artemis_termination293d

subroutine f90wrap_term_gen__artemis_terminationdf16(this)
    use artemis__termination_generator, only: artemis_termination_generator_type
    implicit none
    
    type artemis_termination_generator_type_ptr_type
        type(artemis_termination_generator_type), pointer :: p => NULL()
    end type artemis_termination_generator_type_ptr_type
    type(artemis_termination_generator_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_term_gen__artemis_terminationdf16

subroutine f90wrap_term_gen__generate__binding__2af7(this, basis, miller_plane, axis, surface, &
    num_layers, thickness, orthogonalise, normalise, break_on_fail, n0)
    use artemis__termination_generator, only: artemis_termination_generator_type
    use artemis__geom_rw, only: basis_type
    implicit none
    
    type artemis_termination_generator_type_ptr_type
        type(artemis_termination_generator_type), pointer :: p => NULL()
    end type artemis_termination_generator_type_ptr_type
    type basis_type_ptr_type
        type(basis_type), pointer :: p => NULL()
    end type basis_type_ptr_type
    type(artemis_termination_generator_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    type(basis_type_ptr_type) :: basis_ptr
    integer, intent(in), dimension(2) :: basis
    integer, dimension(3), intent(in) :: miller_plane
    integer, intent(in) :: axis
    integer, intent(in), optional, dimension(n0) :: surface
    integer, intent(in), optional :: num_layers
    real(4), intent(in), optional :: thickness
    logical, intent(in), optional :: orthogonalise
    logical, intent(in), optional :: normalise
    logical, intent(in), optional :: break_on_fail
    integer :: n0
    !f2py intent(hide), depend(surface) :: n0 = shape(surface,0)
    this_ptr = transfer(this, this_ptr)
    basis_ptr = transfer(basis, basis_ptr)
    call this_ptr%p%generate(basis=basis_ptr%p, miller_plane=miller_plane, axis=axis, surface=surface, &
        num_layers=num_layers, thickness=thickness, orthogonalise=orthogonalise, normalise=normalise, &
        break_on_fail=break_on_fail)
end subroutine f90wrap_term_gen__generate__binding__2af7

! End of module artemis__termination_generator defined in file ../src/fortran/lib/mod_term_generator.f90

