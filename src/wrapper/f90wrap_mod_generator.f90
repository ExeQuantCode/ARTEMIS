! Module artemis__generator defined in file ../src/fortran/lib/mod_intf_generator.f90

subroutine f90wrap_artemis_gen_type__get__num_structures(this, f90wrap_num_structures)
    use artemis__generator, only: artemis_generator_type
    implicit none
    type artemis_generator_type_ptr_type
        type(artemis_generator_type), pointer :: p => NULL()
    end type artemis_generator_type_ptr_type
    integer, intent(in)   :: this(2)
    type(artemis_generator_type_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_num_structures
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_num_structures = this_ptr%p%num_structures
end subroutine f90wrap_artemis_gen_type__get__num_structures

subroutine f90wrap_artemis_gen_type__set__num_structures(this, f90wrap_num_structures)
    use artemis__generator, only: artemis_generator_type
    implicit none
    type artemis_generator_type_ptr_type
        type(artemis_generator_type), pointer :: p => NULL()
    end type artemis_generator_type_ptr_type
    integer, intent(in)   :: this(2)
    type(artemis_generator_type_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_num_structures
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%num_structures = f90wrap_num_structures
end subroutine f90wrap_artemis_gen_type__set__num_structures

subroutine f90wrap_artemis_gen_type__get__max_num_structures(this, f90wrap_max_num_structures)
    use artemis__generator, only: artemis_generator_type
    implicit none
    type artemis_generator_type_ptr_type
        type(artemis_generator_type), pointer :: p => NULL()
    end type artemis_generator_type_ptr_type
    integer, intent(in)   :: this(2)
    type(artemis_generator_type_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_max_num_structures
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_max_num_structures = this_ptr%p%max_num_structures
end subroutine f90wrap_artemis_gen_type__get__max_num_structures

subroutine f90wrap_artemis_gen_type__set__max_num_structures(this, f90wrap_max_num_structures)
    use artemis__generator, only: artemis_generator_type
    implicit none
    type artemis_generator_type_ptr_type
        type(artemis_generator_type), pointer :: p => NULL()
    end type artemis_generator_type_ptr_type
    integer, intent(in)   :: this(2)
    type(artemis_generator_type_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_max_num_structures
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%max_num_structures = f90wrap_max_num_structures
end subroutine f90wrap_artemis_gen_type__set__max_num_structures


subroutine f90wrap_artemis_gen_type__get__structure_lw(this, f90wrap_structure_lw)
    use artemis__generator, only: artemis_generator_type
    use artemis__geom_rw, only: basis_type
    implicit none
    type artemis_generator_type_ptr_type
        type(artemis_generator_type), pointer :: p => NULL()
    end type artemis_generator_type_ptr_type
    type basis_type_ptr_type
        type(basis_type), pointer :: p => NULL()
    end type basis_type_ptr_type
    integer, intent(in)   :: this(2)
    type(artemis_generator_type_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_structure_lw(2)
    type(basis_type_ptr_type) :: structure_lw_ptr
    
    this_ptr = transfer(this, this_ptr)
    structure_lw_ptr%p => this_ptr%p%structure_lw
    f90wrap_structure_lw = transfer(structure_lw_ptr,f90wrap_structure_lw)
end subroutine f90wrap_artemis_gen_type__get__structure_lw

subroutine f90wrap_artemis_gen_type__set__structure_lw(this, f90wrap_structure_lw)
    use artemis__generator, only: artemis_generator_type
    use artemis__geom_rw, only: basis_type
    implicit none
    type artemis_generator_type_ptr_type
        type(artemis_generator_type), pointer :: p => NULL()
    end type artemis_generator_type_ptr_type
    type basis_type_ptr_type
        type(basis_type), pointer :: p => NULL()
    end type basis_type_ptr_type
    integer, intent(in)   :: this(2)
    type(artemis_generator_type_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_structure_lw(2)
    type(basis_type_ptr_type) :: structure_lw_ptr
    
    this_ptr = transfer(this, this_ptr)
    structure_lw_ptr = transfer(f90wrap_structure_lw,structure_lw_ptr)
    this_ptr%p%structure_lw = structure_lw_ptr%p
end subroutine f90wrap_artemis_gen_type__set__structure_lw

subroutine f90wrap_artemis_gen_type__get__structure_up(this, f90wrap_structure_up)
    use artemis__generator, only: artemis_generator_type
    use artemis__geom_rw, only: basis_type
    implicit none
    type artemis_generator_type_ptr_type
        type(artemis_generator_type), pointer :: p => NULL()
    end type artemis_generator_type_ptr_type
    type basis_type_ptr_type
        type(basis_type), pointer :: p => NULL()
    end type basis_type_ptr_type
    integer, intent(in)   :: this(2)
    type(artemis_generator_type_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_structure_up(2)
    type(basis_type_ptr_type) :: structure_up_ptr
    
    this_ptr = transfer(this, this_ptr)
    structure_up_ptr%p => this_ptr%p%structure_up
    f90wrap_structure_up = transfer(structure_up_ptr,f90wrap_structure_up)
end subroutine f90wrap_artemis_gen_type__get__structure_up

subroutine f90wrap_artemis_gen_type__set__structure_up(this, f90wrap_structure_up)
    use artemis__generator, only: artemis_generator_type
    use artemis__geom_rw, only: basis_type
    implicit none
    type artemis_generator_type_ptr_type
        type(artemis_generator_type), pointer :: p => NULL()
    end type artemis_generator_type_ptr_type
    type basis_type_ptr_type
        type(basis_type), pointer :: p => NULL()
    end type basis_type_ptr_type
    integer, intent(in)   :: this(2)
    type(artemis_generator_type_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_structure_up(2)
    type(basis_type_ptr_type) :: structure_up_ptr
    
    this_ptr = transfer(this, this_ptr)
    structure_up_ptr = transfer(f90wrap_structure_up,structure_up_ptr)
    this_ptr%p%structure_up = structure_up_ptr%p
end subroutine f90wrap_artemis_gen_type__set__structure_up

subroutine f90wrap_artemis_gen_type__array__elastic_co4c3f(this, nd, dtype, dshape, dloc)
    use artemis__generator, only: artemis_generator_type
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type artemis_generator_type_ptr_type
        type(artemis_generator_type), pointer :: p => NULL()
    end type artemis_generator_type_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(artemis_generator_type_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 11
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%elastic_constants_lw)) then
        dshape(1:1) = shape(this_ptr%p%elastic_constants_lw)
        dloc = loc(this_ptr%p%elastic_constants_lw)
    else
        dloc = 0
    end if
end subroutine f90wrap_artemis_gen_type__array__elastic_co4c3f

subroutine f90wrap_artemis_gen_type__array__elastic_coedb6(this, nd, dtype, dshape, dloc)
    use artemis__generator, only: artemis_generator_type
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type artemis_generator_type_ptr_type
        type(artemis_generator_type), pointer :: p => NULL()
    end type artemis_generator_type_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(artemis_generator_type_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 11
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%elastic_constants_up)) then
        dshape(1:1) = shape(this_ptr%p%elastic_constants_up)
        dloc = loc(this_ptr%p%elastic_constants_up)
    else
        dloc = 0
    end if
end subroutine f90wrap_artemis_gen_type__array__elastic_coedb6

subroutine f90wrap_artemis_gen_type__get__use_pricel_lw(this, f90wrap_use_pricel_lw)
    use artemis__generator, only: artemis_generator_type
    implicit none
    type artemis_generator_type_ptr_type
        type(artemis_generator_type), pointer :: p => NULL()
    end type artemis_generator_type_ptr_type
    integer, intent(in)   :: this(2)
    type(artemis_generator_type_ptr_type) :: this_ptr
    logical, intent(out) :: f90wrap_use_pricel_lw
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_use_pricel_lw = this_ptr%p%use_pricel_lw
end subroutine f90wrap_artemis_gen_type__get__use_pricel_lw

subroutine f90wrap_artemis_gen_type__set__use_pricel_lw(this, f90wrap_use_pricel_lw)
    use artemis__generator, only: artemis_generator_type
    implicit none
    type artemis_generator_type_ptr_type
        type(artemis_generator_type), pointer :: p => NULL()
    end type artemis_generator_type_ptr_type
    integer, intent(in)   :: this(2)
    type(artemis_generator_type_ptr_type) :: this_ptr
    logical, intent(in) :: f90wrap_use_pricel_lw
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%use_pricel_lw = f90wrap_use_pricel_lw
end subroutine f90wrap_artemis_gen_type__set__use_pricel_lw

subroutine f90wrap_artemis_gen_type__get__use_pricel_up(this, f90wrap_use_pricel_up)
    use artemis__generator, only: artemis_generator_type
    implicit none
    type artemis_generator_type_ptr_type
        type(artemis_generator_type), pointer :: p => NULL()
    end type artemis_generator_type_ptr_type
    integer, intent(in)   :: this(2)
    type(artemis_generator_type_ptr_type) :: this_ptr
    logical, intent(out) :: f90wrap_use_pricel_up
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_use_pricel_up = this_ptr%p%use_pricel_up
end subroutine f90wrap_artemis_gen_type__get__use_pricel_up

subroutine f90wrap_artemis_gen_type__set__use_pricel_up(this, f90wrap_use_pricel_up)
    use artemis__generator, only: artemis_generator_type
    implicit none
    type artemis_generator_type_ptr_type
        type(artemis_generator_type), pointer :: p => NULL()
    end type artemis_generator_type_ptr_type
    integer, intent(in)   :: this(2)
    type(artemis_generator_type_ptr_type) :: this_ptr
    logical, intent(in) :: f90wrap_use_pricel_up
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%use_pricel_up = f90wrap_use_pricel_up
end subroutine f90wrap_artemis_gen_type__set__use_pricel_up

subroutine f90wrap_artemis_gen_type__array__miller_lw(this, nd, dtype, dshape, dloc)
    use artemis__generator, only: artemis_generator_type
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type artemis_generator_type_ptr_type
        type(artemis_generator_type), pointer :: p => NULL()
    end type artemis_generator_type_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(artemis_generator_type_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    dshape(1:1) = shape(this_ptr%p%miller_lw)
    dloc = loc(this_ptr%p%miller_lw)
end subroutine f90wrap_artemis_gen_type__array__miller_lw

subroutine f90wrap_artemis_gen_type__get__is_layered_lw(this, f90wrap_is_layered_lw)
    use artemis__generator, only: artemis_generator_type
    implicit none
    type artemis_generator_type_ptr_type
        type(artemis_generator_type), pointer :: p => NULL()
    end type artemis_generator_type_ptr_type
    integer, intent(in)   :: this(2)
    type(artemis_generator_type_ptr_type) :: this_ptr
    logical, intent(out) :: f90wrap_is_layered_lw
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_is_layered_lw = this_ptr%p%is_layered_lw
end subroutine f90wrap_artemis_gen_type__get__is_layered_lw

subroutine f90wrap_artemis_gen_type__set__is_layered_lw(this, f90wrap_is_layered_lw)
    use artemis__generator, only: artemis_generator_type
    implicit none
    type artemis_generator_type_ptr_type
        type(artemis_generator_type), pointer :: p => NULL()
    end type artemis_generator_type_ptr_type
    integer, intent(in)   :: this(2)
    type(artemis_generator_type_ptr_type) :: this_ptr
    logical, intent(in) :: f90wrap_is_layered_lw
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%is_layered_lw = f90wrap_is_layered_lw
end subroutine f90wrap_artemis_gen_type__set__is_layered_lw

subroutine f90wrap_artemis_gen_type__get__is_layered_up(this, f90wrap_is_layered_up)
    use artemis__generator, only: artemis_generator_type
    implicit none
    type artemis_generator_type_ptr_type
        type(artemis_generator_type), pointer :: p => NULL()
    end type artemis_generator_type_ptr_type
    integer, intent(in)   :: this(2)
    type(artemis_generator_type_ptr_type) :: this_ptr
    logical, intent(out) :: f90wrap_is_layered_up
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_is_layered_up = this_ptr%p%is_layered_up
end subroutine f90wrap_artemis_gen_type__get__is_layered_up

subroutine f90wrap_artemis_gen_type__set__is_layered_up(this, f90wrap_is_layered_up)
    use artemis__generator, only: artemis_generator_type
    implicit none
    type artemis_generator_type_ptr_type
        type(artemis_generator_type), pointer :: p => NULL()
    end type artemis_generator_type_ptr_type
    integer, intent(in)   :: this(2)
    type(artemis_generator_type_ptr_type) :: this_ptr
    logical, intent(in) :: f90wrap_is_layered_up
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%is_layered_up = f90wrap_is_layered_up
end subroutine f90wrap_artemis_gen_type__set__is_layered_up

subroutine f90wrap_artemis_gen_type__get__ludef_is_lay4aa6(this, f90wrap_ludef_is_layered_lw)
    use artemis__generator, only: artemis_generator_type
    implicit none
    type artemis_generator_type_ptr_type
        type(artemis_generator_type), pointer :: p => NULL()
    end type artemis_generator_type_ptr_type
    integer, intent(in)   :: this(2)
    type(artemis_generator_type_ptr_type) :: this_ptr
    logical, intent(out) :: f90wrap_ludef_is_layered_lw
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_ludef_is_layered_lw = this_ptr%p%ludef_is_layered_lw
end subroutine f90wrap_artemis_gen_type__get__ludef_is_lay4aa6

subroutine f90wrap_artemis_gen_type__set__ludef_is_lay87a5(this, f90wrap_ludef_is_layered_lw)
    use artemis__generator, only: artemis_generator_type
    implicit none
    type artemis_generator_type_ptr_type
        type(artemis_generator_type), pointer :: p => NULL()
    end type artemis_generator_type_ptr_type
    integer, intent(in)   :: this(2)
    type(artemis_generator_type_ptr_type) :: this_ptr
    logical, intent(in) :: f90wrap_ludef_is_layered_lw
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%ludef_is_layered_lw = f90wrap_ludef_is_layered_lw
end subroutine f90wrap_artemis_gen_type__set__ludef_is_lay87a5

subroutine f90wrap_artemis_gen_type__get__ludef_is_lay60fd(this, f90wrap_ludef_is_layered_up)
    use artemis__generator, only: artemis_generator_type
    implicit none
    type artemis_generator_type_ptr_type
        type(artemis_generator_type), pointer :: p => NULL()
    end type artemis_generator_type_ptr_type
    integer, intent(in)   :: this(2)
    type(artemis_generator_type_ptr_type) :: this_ptr
    logical, intent(out) :: f90wrap_ludef_is_layered_up
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_ludef_is_layered_up = this_ptr%p%ludef_is_layered_up
end subroutine f90wrap_artemis_gen_type__get__ludef_is_lay60fd

subroutine f90wrap_artemis_gen_type__set__ludef_is_laye6e4(this, f90wrap_ludef_is_layered_up)
    use artemis__generator, only: artemis_generator_type
    implicit none
    type artemis_generator_type_ptr_type
        type(artemis_generator_type), pointer :: p => NULL()
    end type artemis_generator_type_ptr_type
    integer, intent(in)   :: this(2)
    type(artemis_generator_type_ptr_type) :: this_ptr
    logical, intent(in) :: f90wrap_ludef_is_layered_up
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%ludef_is_layered_up = f90wrap_ludef_is_layered_up
end subroutine f90wrap_artemis_gen_type__set__ludef_is_laye6e4






subroutine f90wrap_artemis_gen_type__get__shift_method(this, f90wrap_shift_method)
    use artemis__generator, only: artemis_generator_type
    implicit none
    type artemis_generator_type_ptr_type
        type(artemis_generator_type), pointer :: p => NULL()
    end type artemis_generator_type_ptr_type
    integer, intent(in)   :: this(2)
    type(artemis_generator_type_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_shift_method
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_shift_method = this_ptr%p%shift_method
end subroutine f90wrap_artemis_gen_type__get__shift_method

subroutine f90wrap_artemis_gen_type__set__shift_method(this, f90wrap_shift_method)
    use artemis__generator, only: artemis_generator_type
    implicit none
    type artemis_generator_type_ptr_type
        type(artemis_generator_type), pointer :: p => NULL()
    end type artemis_generator_type_ptr_type
    integer, intent(in)   :: this(2)
    type(artemis_generator_type_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_shift_method
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%shift_method = f90wrap_shift_method
end subroutine f90wrap_artemis_gen_type__set__shift_method

subroutine f90wrap_artemis_gen_type__get__num_shifts(this, f90wrap_num_shifts)
    use artemis__generator, only: artemis_generator_type
    implicit none
    type artemis_generator_type_ptr_type
        type(artemis_generator_type), pointer :: p => NULL()
    end type artemis_generator_type_ptr_type
    integer, intent(in)   :: this(2)
    type(artemis_generator_type_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_num_shifts
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_num_shifts = this_ptr%p%num_shifts
end subroutine f90wrap_artemis_gen_type__get__num_shifts

subroutine f90wrap_artemis_gen_type__set__num_shifts(this, f90wrap_num_shifts)
    use artemis__generator, only: artemis_generator_type
    implicit none
    type artemis_generator_type_ptr_type
        type(artemis_generator_type), pointer :: p => NULL()
    end type artemis_generator_type_ptr_type
    integer, intent(in)   :: this(2)
    type(artemis_generator_type_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_num_shifts
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%num_shifts = f90wrap_num_shifts
end subroutine f90wrap_artemis_gen_type__set__num_shifts

subroutine f90wrap_artemis_gen_type__array__shifts(this, nd, dtype, dshape, dloc)
    use artemis__generator, only: artemis_generator_type
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type artemis_generator_type_ptr_type
        type(artemis_generator_type), pointer :: p => NULL()
    end type artemis_generator_type_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(artemis_generator_type_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 11
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%shifts)) then
        dshape(1:2) = shape(this_ptr%p%shifts)
        dloc = loc(this_ptr%p%shifts)
    else
        dloc = 0
    end if
end subroutine f90wrap_artemis_gen_type__array__shifts

subroutine f90wrap_artemis_gen_type__get__interface_depth(this, f90wrap_interface_depth)
    use artemis__generator, only: artemis_generator_type
    implicit none
    type artemis_generator_type_ptr_type
        type(artemis_generator_type), pointer :: p => NULL()
    end type artemis_generator_type_ptr_type
    integer, intent(in)   :: this(2)
    type(artemis_generator_type_ptr_type) :: this_ptr
    real(4), intent(out) :: f90wrap_interface_depth
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_interface_depth = this_ptr%p%interface_depth
end subroutine f90wrap_artemis_gen_type__get__interface_depth

subroutine f90wrap_artemis_gen_type__set__interface_depth(this, f90wrap_interface_depth)
    use artemis__generator, only: artemis_generator_type
    implicit none
    type artemis_generator_type_ptr_type
        type(artemis_generator_type), pointer :: p => NULL()
    end type artemis_generator_type_ptr_type
    integer, intent(in)   :: this(2)
    type(artemis_generator_type_ptr_type) :: this_ptr
    real(4), intent(in) :: f90wrap_interface_depth
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%interface_depth = f90wrap_interface_depth
end subroutine f90wrap_artemis_gen_type__set__interface_depth

subroutine f90wrap_artemis_gen_type__get__separation_scale(this, f90wrap_separation_scale)
    use artemis__generator, only: artemis_generator_type
    implicit none
    type artemis_generator_type_ptr_type
        type(artemis_generator_type), pointer :: p => NULL()
    end type artemis_generator_type_ptr_type
    integer, intent(in)   :: this(2)
    type(artemis_generator_type_ptr_type) :: this_ptr
    real(4), intent(out) :: f90wrap_separation_scale
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_separation_scale = this_ptr%p%separation_scale
end subroutine f90wrap_artemis_gen_type__get__separation_scale

subroutine f90wrap_artemis_gen_type__set__separation_scale(this, f90wrap_separation_scale)
    use artemis__generator, only: artemis_generator_type
    implicit none
    type artemis_generator_type_ptr_type
        type(artemis_generator_type), pointer :: p => NULL()
    end type artemis_generator_type_ptr_type
    integer, intent(in)   :: this(2)
    type(artemis_generator_type_ptr_type) :: this_ptr
    real(4), intent(in) :: f90wrap_separation_scale
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%separation_scale = f90wrap_separation_scale
end subroutine f90wrap_artemis_gen_type__set__separation_scale

subroutine f90wrap_artemis_gen_type__get__depth_method(this, f90wrap_depth_method)
    use artemis__generator, only: artemis_generator_type
    implicit none
    type artemis_generator_type_ptr_type
        type(artemis_generator_type), pointer :: p => NULL()
    end type artemis_generator_type_ptr_type
    integer, intent(in)   :: this(2)
    type(artemis_generator_type_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_depth_method
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_depth_method = this_ptr%p%depth_method
end subroutine f90wrap_artemis_gen_type__get__depth_method

subroutine f90wrap_artemis_gen_type__set__depth_method(this, f90wrap_depth_method)
    use artemis__generator, only: artemis_generator_type
    implicit none
    type artemis_generator_type_ptr_type
        type(artemis_generator_type), pointer :: p => NULL()
    end type artemis_generator_type_ptr_type
    integer, intent(in)   :: this(2)
    type(artemis_generator_type_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_depth_method
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%depth_method = f90wrap_depth_method
end subroutine f90wrap_artemis_gen_type__set__depth_method

subroutine f90wrap_artemis_gen_type__array__shift_data(this, nd, dtype, dshape, dloc)
    use artemis__generator, only: artemis_generator_type
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type artemis_generator_type_ptr_type
        type(artemis_generator_type), pointer :: p => NULL()
    end type artemis_generator_type_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(artemis_generator_type_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 11
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%shift_data)) then
        dshape(1:2) = shape(this_ptr%p%shift_data)
        dloc = loc(this_ptr%p%shift_data)
    else
        dloc = 0
    end if
end subroutine f90wrap_artemis_gen_type__array__shift_data

subroutine f90wrap_artemis_gen_type__get__swap_method(this, f90wrap_swap_method)
    use artemis__generator, only: artemis_generator_type
    implicit none
    type artemis_generator_type_ptr_type
        type(artemis_generator_type), pointer :: p => NULL()
    end type artemis_generator_type_ptr_type
    integer, intent(in)   :: this(2)
    type(artemis_generator_type_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_swap_method
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_swap_method = this_ptr%p%swap_method
end subroutine f90wrap_artemis_gen_type__get__swap_method

subroutine f90wrap_artemis_gen_type__set__swap_method(this, f90wrap_swap_method)
    use artemis__generator, only: artemis_generator_type
    implicit none
    type artemis_generator_type_ptr_type
        type(artemis_generator_type), pointer :: p => NULL()
    end type artemis_generator_type_ptr_type
    integer, intent(in)   :: this(2)
    type(artemis_generator_type_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_swap_method
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%swap_method = f90wrap_swap_method
end subroutine f90wrap_artemis_gen_type__set__swap_method

subroutine f90wrap_artemis_gen_type__get__num_swaps(this, f90wrap_num_swaps)
    use artemis__generator, only: artemis_generator_type
    implicit none
    type artemis_generator_type_ptr_type
        type(artemis_generator_type), pointer :: p => NULL()
    end type artemis_generator_type_ptr_type
    integer, intent(in)   :: this(2)
    type(artemis_generator_type_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_num_swaps
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_num_swaps = this_ptr%p%num_swaps
end subroutine f90wrap_artemis_gen_type__get__num_swaps

subroutine f90wrap_artemis_gen_type__set__num_swaps(this, f90wrap_num_swaps)
    use artemis__generator, only: artemis_generator_type
    implicit none
    type artemis_generator_type_ptr_type
        type(artemis_generator_type), pointer :: p => NULL()
    end type artemis_generator_type_ptr_type
    integer, intent(in)   :: this(2)
    type(artemis_generator_type_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_num_swaps
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%num_swaps = f90wrap_num_swaps
end subroutine f90wrap_artemis_gen_type__set__num_swaps

subroutine f90wrap_artemis_gen_type__get__swap_density(this, f90wrap_swap_density)
    use artemis__generator, only: artemis_generator_type
    implicit none
    type artemis_generator_type_ptr_type
        type(artemis_generator_type), pointer :: p => NULL()
    end type artemis_generator_type_ptr_type
    integer, intent(in)   :: this(2)
    type(artemis_generator_type_ptr_type) :: this_ptr
    real(4), intent(out) :: f90wrap_swap_density
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_swap_density = this_ptr%p%swap_density
end subroutine f90wrap_artemis_gen_type__get__swap_density

subroutine f90wrap_artemis_gen_type__set__swap_density(this, f90wrap_swap_density)
    use artemis__generator, only: artemis_generator_type
    implicit none
    type artemis_generator_type_ptr_type
        type(artemis_generator_type), pointer :: p => NULL()
    end type artemis_generator_type_ptr_type
    integer, intent(in)   :: this(2)
    type(artemis_generator_type_ptr_type) :: this_ptr
    real(4), intent(in) :: f90wrap_swap_density
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%swap_density = f90wrap_swap_density
end subroutine f90wrap_artemis_gen_type__set__swap_density

subroutine f90wrap_artemis_gen_type__get__swap_depth(this, f90wrap_swap_depth)
    use artemis__generator, only: artemis_generator_type
    implicit none
    type artemis_generator_type_ptr_type
        type(artemis_generator_type), pointer :: p => NULL()
    end type artemis_generator_type_ptr_type
    integer, intent(in)   :: this(2)
    type(artemis_generator_type_ptr_type) :: this_ptr
    real(4), intent(out) :: f90wrap_swap_depth
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_swap_depth = this_ptr%p%swap_depth
end subroutine f90wrap_artemis_gen_type__get__swap_depth

subroutine f90wrap_artemis_gen_type__set__swap_depth(this, f90wrap_swap_depth)
    use artemis__generator, only: artemis_generator_type
    implicit none
    type artemis_generator_type_ptr_type
        type(artemis_generator_type), pointer :: p => NULL()
    end type artemis_generator_type_ptr_type
    integer, intent(in)   :: this(2)
    type(artemis_generator_type_ptr_type) :: this_ptr
    real(4), intent(in) :: f90wrap_swap_depth
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%swap_depth = f90wrap_swap_depth
end subroutine f90wrap_artemis_gen_type__set__swap_depth

subroutine f90wrap_artemis_gen_type__get__swap_sigma(this, f90wrap_swap_sigma)
    use artemis__generator, only: artemis_generator_type
    implicit none
    type artemis_generator_type_ptr_type
        type(artemis_generator_type), pointer :: p => NULL()
    end type artemis_generator_type_ptr_type
    integer, intent(in)   :: this(2)
    type(artemis_generator_type_ptr_type) :: this_ptr
    real(4), intent(out) :: f90wrap_swap_sigma
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_swap_sigma = this_ptr%p%swap_sigma
end subroutine f90wrap_artemis_gen_type__get__swap_sigma

subroutine f90wrap_artemis_gen_type__set__swap_sigma(this, f90wrap_swap_sigma)
    use artemis__generator, only: artemis_generator_type
    implicit none
    type artemis_generator_type_ptr_type
        type(artemis_generator_type), pointer :: p => NULL()
    end type artemis_generator_type_ptr_type
    integer, intent(in)   :: this(2)
    type(artemis_generator_type_ptr_type) :: this_ptr
    real(4), intent(in) :: f90wrap_swap_sigma
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%swap_sigma = f90wrap_swap_sigma
end subroutine f90wrap_artemis_gen_type__set__swap_sigma

subroutine f90wrap_artemis_gen_type__get__require_mirr41cf( &
    this, f90wrap_require_mirror_swaps)
    use artemis__generator, only: artemis_generator_type
    implicit none
    type artemis_generator_type_ptr_type
        type(artemis_generator_type), pointer :: p => NULL()
    end type artemis_generator_type_ptr_type
    integer, intent(in)   :: this(2)
    type(artemis_generator_type_ptr_type) :: this_ptr
    logical, intent(out) :: f90wrap_require_mirror_swaps
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_require_mirror_swaps = this_ptr%p%require_mirror_swaps
end subroutine f90wrap_artemis_gen_type__get__require_mirr41cf

subroutine f90wrap_artemis_gen_type__set__require_mirr3bfa( &
    this, f90wrap_require_mirror_swaps)
    use artemis__generator, only: artemis_generator_type
    implicit none
    type artemis_generator_type_ptr_type
        type(artemis_generator_type), pointer :: p => NULL()
    end type artemis_generator_type_ptr_type
    integer, intent(in)   :: this(2)
    type(artemis_generator_type_ptr_type) :: this_ptr
    logical, intent(in) :: f90wrap_require_mirror_swaps
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%require_mirror_swaps = f90wrap_require_mirror_swaps
end subroutine f90wrap_artemis_gen_type__set__require_mirr3bfa

subroutine f90wrap_artemis_gen_type__get__match_method(this, f90wrap_match_method)
    use artemis__generator, only: artemis_generator_type
    implicit none
    type artemis_generator_type_ptr_type
        type(artemis_generator_type), pointer :: p => NULL()
    end type artemis_generator_type_ptr_type
    integer, intent(in)   :: this(2)
    type(artemis_generator_type_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_match_method
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_match_method = this_ptr%p%match_method
end subroutine f90wrap_artemis_gen_type__get__match_method

subroutine f90wrap_artemis_gen_type__set__match_method(this, f90wrap_match_method)
    use artemis__generator, only: artemis_generator_type
    implicit none
    type artemis_generator_type_ptr_type
        type(artemis_generator_type), pointer :: p => NULL()
    end type artemis_generator_type_ptr_type
    integer, intent(in)   :: this(2)
    type(artemis_generator_type_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_match_method
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%match_method = f90wrap_match_method
end subroutine f90wrap_artemis_gen_type__set__match_method

subroutine f90wrap_artemis_gen_type__get__max_num_matches(this, f90wrap_max_num_matches)
    use artemis__generator, only: artemis_generator_type
    implicit none
    type artemis_generator_type_ptr_type
        type(artemis_generator_type), pointer :: p => NULL()
    end type artemis_generator_type_ptr_type
    integer, intent(in)   :: this(2)
    type(artemis_generator_type_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_max_num_matches
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_max_num_matches = this_ptr%p%max_num_matches
end subroutine f90wrap_artemis_gen_type__get__max_num_matches

subroutine f90wrap_artemis_gen_type__set__max_num_matches(this, f90wrap_max_num_matches)
    use artemis__generator, only: artemis_generator_type
    implicit none
    type artemis_generator_type_ptr_type
        type(artemis_generator_type), pointer :: p => NULL()
    end type artemis_generator_type_ptr_type
    integer, intent(in)   :: this(2)
    type(artemis_generator_type_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_max_num_matches
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%max_num_matches = f90wrap_max_num_matches
end subroutine f90wrap_artemis_gen_type__set__max_num_matches

subroutine f90wrap_artemis_gen_type__get__max_num_terms(this, f90wrap_max_num_terms)
    use artemis__generator, only: artemis_generator_type
    implicit none
    type artemis_generator_type_ptr_type
        type(artemis_generator_type), pointer :: p => NULL()
    end type artemis_generator_type_ptr_type
    integer, intent(in)   :: this(2)
    type(artemis_generator_type_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_max_num_terms
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_max_num_terms = this_ptr%p%max_num_terms
end subroutine f90wrap_artemis_gen_type__get__max_num_terms

subroutine f90wrap_artemis_gen_type__set__max_num_terms(this, f90wrap_max_num_terms)
    use artemis__generator, only: artemis_generator_type
    implicit none
    type artemis_generator_type_ptr_type
        type(artemis_generator_type), pointer :: p => NULL()
    end type artemis_generator_type_ptr_type
    integer, intent(in)   :: this(2)
    type(artemis_generator_type_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_max_num_terms
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%max_num_terms = f90wrap_max_num_terms
end subroutine f90wrap_artemis_gen_type__set__max_num_terms

subroutine f90wrap_artemis_gen_type__get__max_num_planes(this, f90wrap_max_num_planes)
    use artemis__generator, only: artemis_generator_type
    implicit none
    type artemis_generator_type_ptr_type
        type(artemis_generator_type), pointer :: p => NULL()
    end type artemis_generator_type_ptr_type
    integer, intent(in)   :: this(2)
    type(artemis_generator_type_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_max_num_planes
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_max_num_planes = this_ptr%p%max_num_planes
end subroutine f90wrap_artemis_gen_type__get__max_num_planes

subroutine f90wrap_artemis_gen_type__set__max_num_planes(this, f90wrap_max_num_planes)
    use artemis__generator, only: artemis_generator_type
    implicit none
    type artemis_generator_type_ptr_type
        type(artemis_generator_type), pointer :: p => NULL()
    end type artemis_generator_type_ptr_type
    integer, intent(in)   :: this(2)
    type(artemis_generator_type_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_max_num_planes
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%max_num_planes = f90wrap_max_num_planes
end subroutine f90wrap_artemis_gen_type__set__max_num_planes

subroutine f90wrap_artemis_gen_type__get__fix_normal(this, f90wrap_fix_normal)
    use artemis__generator, only: artemis_generator_type
    implicit none
    type artemis_generator_type_ptr_type
        type(artemis_generator_type), pointer :: p => NULL()
    end type artemis_generator_type_ptr_type
    integer, intent(in)   :: this(2)
    type(artemis_generator_type_ptr_type) :: this_ptr
    logical, intent(out) :: f90wrap_fix_normal
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_fix_normal = this_ptr%p%fix_normal
end subroutine f90wrap_artemis_gen_type__get__fix_normal

subroutine f90wrap_artemis_gen_type__set__fix_normal(this, f90wrap_fix_normal)
    use artemis__generator, only: artemis_generator_type
    implicit none
    type artemis_generator_type_ptr_type
        type(artemis_generator_type), pointer :: p => NULL()
    end type artemis_generator_type_ptr_type
    integer, intent(in)   :: this(2)
    type(artemis_generator_type_ptr_type) :: this_ptr
    logical, intent(in) :: f90wrap_fix_normal
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%fix_normal = f90wrap_fix_normal
end subroutine f90wrap_artemis_gen_type__set__fix_normal

subroutine f90wrap_artemis_gen_type__get__bondlength_c21a8(this, f90wrap_bondlength_cutoff)
    use artemis__generator, only: artemis_generator_type
    implicit none
    type artemis_generator_type_ptr_type
        type(artemis_generator_type), pointer :: p => NULL()
    end type artemis_generator_type_ptr_type
    integer, intent(in)   :: this(2)
    type(artemis_generator_type_ptr_type) :: this_ptr
    real(4), intent(out) :: f90wrap_bondlength_cutoff
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_bondlength_cutoff = this_ptr%p%bondlength_cutoff
end subroutine f90wrap_artemis_gen_type__get__bondlength_c21a8

subroutine f90wrap_artemis_gen_type__set__bondlength_cbd11(this, f90wrap_bondlength_cutoff)
    use artemis__generator, only: artemis_generator_type
    implicit none
    type artemis_generator_type_ptr_type
        type(artemis_generator_type), pointer :: p => NULL()
    end type artemis_generator_type_ptr_type
    integer, intent(in)   :: this(2)
    type(artemis_generator_type_ptr_type) :: this_ptr
    real(4), intent(in) :: f90wrap_bondlength_cutoff
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%bondlength_cutoff = f90wrap_bondlength_cutoff
end subroutine f90wrap_artemis_gen_type__set__bondlength_cbd11

subroutine f90wrap_artemis_gen_type__array__layer_sepa90a5(this, nd, dtype, dshape, dloc)
    use artemis__generator, only: artemis_generator_type
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type artemis_generator_type_ptr_type
        type(artemis_generator_type), pointer :: p => NULL()
    end type artemis_generator_type_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(artemis_generator_type_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 11
    this_ptr = transfer(this, this_ptr)
    dshape(1:1) = shape(this_ptr%p%layer_separation_cutoff)
    dloc = loc(this_ptr%p%layer_separation_cutoff)
end subroutine f90wrap_artemis_gen_type__array__layer_sepa90a5

subroutine f90wrap_intf_gen__artemis_gen_type_initialise(this)
    use artemis__generator, only: artemis_generator_type
    implicit none
    
    type artemis_generator_type_ptr_type
        type(artemis_generator_type), pointer :: p => NULL()
    end type artemis_generator_type_ptr_type
    type(artemis_generator_type_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_intf_gen__artemis_gen_type_initialise

subroutine f90wrap_intf_gen__artemis_gen_type_finalise(this)
    use artemis__generator, only: artemis_generator_type
    implicit none
    
    type artemis_generator_type_ptr_type
        type(artemis_generator_type), pointer :: p => NULL()
    end type artemis_generator_type_ptr_type
    type(artemis_generator_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_intf_gen__artemis_gen_type_finalise

subroutine f90wrap_intf_gen__set_tolerance__bindind_agt(this, vector_mismatch, angle_mismatch, &
    area_mismatch, max_length, max_area, max_fit, max_extension, angle_weight, area_weight)
    use artemis__generator, only: artemis_generator_type
    implicit none
    
    type artemis_generator_type_ptr_type
        type(artemis_generator_type), pointer :: p => NULL()
    end type artemis_generator_type_ptr_type
    type(artemis_generator_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    real(4), intent(in), optional :: vector_mismatch
    real(4), intent(in), optional :: angle_mismatch
    real(4), intent(in), optional :: area_mismatch
    real(4), intent(in), optional :: max_length
    real(4), intent(in), optional :: max_area
    integer, intent(in), optional :: max_fit
    integer, intent(in), optional :: max_extension
    real(4), intent(in), optional :: angle_weight
    real(4), intent(in), optional :: area_weight
    this_ptr = transfer(this, this_ptr)
    call this_ptr%p%set_tolerance( &
        vector_mismatch=vector_mismatch, &
        angle_mismatch=angle_mismatch, &
        area_mismatch=area_mismatch, max_length=max_length, &
        max_area=max_area, max_fit=max_fit, max_extension=max_extension, &
        angle_weight=angle_weight, area_weight=area_weight)
end subroutine f90wrap_intf_gen__set_tolerance__bindind_agt

subroutine f90wrap_intf_gen__set_shift_method__binding_agt(this, method, num_shifts, shifts, &
    interface_depth, separation_scale, depth_method, n0)
    use artemis__generator, only: artemis_generator_type
    implicit none
    
    type artemis_generator_type_ptr_type
        type(artemis_generator_type), pointer :: p => NULL()
    end type artemis_generator_type_ptr_type
    type(artemis_generator_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    integer, intent(in), optional :: method
    integer, intent(in), optional :: num_shifts
    real(4), dimension(n0), intent(in), optional :: shifts
    real(4), intent(in), optional :: interface_depth
    real(4), intent(in), optional :: separation_scale
    integer, intent(in), optional :: depth_method
    integer :: n0
    !f2py intent(hide), depend(shifts) :: n0 = shape(shifts,0)
    this_ptr = transfer(this, this_ptr)
    call this_ptr%p%set_shift_method( &
        method=method, num_shifts=num_shifts, &
        shifts=shifts, interface_depth=interface_depth, &
        separation_scale=separation_scale, depth_method=depth_method)
end subroutine f90wrap_intf_gen__set_shift_method__binding_agt

subroutine f90wrap_intf_gen__set_materials__binding_agt(this, structure_lw, structure_up, &
    elastic_constants_lw, elastic_constants_up, use_pricel_lw, use_pricel_up, n0, n1)
    use artemis__generator, only: artemis_generator_type
    use artemis__geom_rw, only: basis_type
    implicit none
    
    type basis_type_ptr_type
        type(basis_type), pointer :: p => NULL()
    end type basis_type_ptr_type
    type artemis_generator_type_ptr_type
        type(artemis_generator_type), pointer :: p => NULL()
    end type artemis_generator_type_ptr_type
    type(artemis_generator_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    type(basis_type_ptr_type) :: structure_lw_ptr
    integer, intent(in), dimension(2) :: structure_lw
    type(basis_type_ptr_type) :: structure_up_ptr
    integer, intent(in), dimension(2) :: structure_up
    real(4), intent(in), optional, dimension(n0) :: elastic_constants_lw
    real(4), intent(in), optional, dimension(n1) :: elastic_constants_up
    logical, intent(in), optional :: use_pricel_lw
    logical, intent(in), optional :: use_pricel_up
    integer :: n0
    !f2py intent(hide), depend(elastic_constants_lw) :: n0 = shape(elastic_constants_lw,0)
    integer :: n1
    !f2py intent(hide), depend(elastic_constants_up) :: n1 = shape(elastic_constants_up,0)
    this_ptr = transfer(this, this_ptr)
    structure_lw_ptr = transfer(structure_lw, structure_lw_ptr)
    structure_up_ptr = transfer(structure_up, structure_up_ptr)
    call this_ptr%p%set_materials(structure_lw=structure_lw_ptr%p, structure_up=structure_up_ptr%p, &
        elastic_constants_lw=elastic_constants_lw, elastic_constants_up=elastic_constants_up, use_pricel_lw=use_pricel_lw, &
        use_pricel_up=use_pricel_up)
end subroutine f90wrap_intf_gen__set_materials__binding_agt

subroutine f90wrap_intf_gen__set_surface_properties__binding_agt(this, miller_lw, miller_up, is_layered_lw, &
    is_layered_up)
    use artemis__generator, only: artemis_generator_type
    implicit none
    
    type artemis_generator_type_ptr_type
        type(artemis_generator_type), pointer :: p => NULL()
    end type artemis_generator_type_ptr_type
    type(artemis_generator_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    integer, dimension(3), intent(in), optional :: miller_lw
    integer, dimension(3), intent(in), optional :: miller_up
    logical, intent(in), optional :: is_layered_lw
    logical, intent(in), optional :: is_layered_up
    this_ptr = transfer(this, this_ptr)
    call this_ptr%p%set_surface_properties(miller_lw=miller_lw, miller_up=miller_up, is_layered_lw=is_layered_lw, &
        is_layered_up=is_layered_up)
end subroutine f90wrap_intf_gen__set_surface_properties__binding_agt

subroutine f90wrap_intf_gen__reset_is_layered_lw__binding_agt(this)
    use artemis__generator, only: artemis_generator_type
    implicit none
    
    type artemis_generator_type_ptr_type
        type(artemis_generator_type), pointer :: p => NULL()
    end type artemis_generator_type_ptr_type
    type(artemis_generator_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    call this_ptr%p%reset_is_layered_lw()
end subroutine f90wrap_intf_gen__reset_is_layered_lw__binding_agt

subroutine f90wrap_intf_gen__reset_is_layered_up__binding_agt(this)
    use artemis__generator, only: artemis_generator_type
    implicit none
    
    type artemis_generator_type_ptr_type
        type(artemis_generator_type), pointer :: p => NULL()
    end type artemis_generator_type_ptr_type
    type(artemis_generator_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    call this_ptr%p%reset_is_layered_up()
end subroutine f90wrap_intf_gen__reset_is_layered_up__binding_agt











subroutine f90wrap_intf_gen__get_terminations__binding__agt( &
    this, identifier, miller, surface, num_layers, thickness, &
    orthogonalise, normalise, break_on_fail, &
    verbose, exit_code, &
    n_ret_structures, n0)
    use artemis__geom_rw, only: basis_type
    use artemis__generator, only: artemis_generator_type
    use artemis__structure_cache, only: store_last_generated_structures
    implicit none

    type artemis_generator_type_ptr_type
        type(artemis_generator_type), pointer :: p => NULL()
    end type artemis_generator_type_ptr_type

    type basis_type_xnum_array
        type(basis_type), dimension(:), allocatable :: items
    end type basis_type_xnum_array

    type basis_type_xnum_array_ptr_type
        type(basis_type_xnum_array), pointer :: p => NULL()
    end type basis_type_xnum_array_ptr_type
    type(artemis_generator_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    integer, intent(in):: identifier
    integer, dimension(3), intent(in), optional :: miller
    integer, intent(in), optional, dimension(n0) :: surface
    real(4), intent(in), optional :: thickness
    integer, intent(in), optional :: num_layers
    logical, intent(in), optional :: orthogonalise
    logical, intent(in), optional :: normalise
    logical, intent(in), optional :: break_on_fail
    integer, intent(in), optional :: verbose
    integer, optional, intent(out) :: exit_code
    integer :: n0
    !f2py intent(hide), depend(surface_lw) :: n0 = shape(surface_lw,0)
    type(basis_type), allocatable, dimension(:) :: local_structures
    integer, intent(out) :: n_ret_structures

    this_ptr = transfer(this, this_ptr)
    local_structures = this_ptr%p%get_terminations( &
        identifier=identifier, miller=miller, surface=surface, &
        num_layers=num_layers, thickness=thickness, &
        orthogonalise=orthogonalise, normalise=normalise, &
        break_on_fail=break_on_fail, verbose=verbose, exit_code=exit_code)

    n_ret_structures = size(local_structures, dim=1)

    ! Store local_structures in a module-level array so Python can retrieve it
    call store_last_generated_structures(local_structures)
end subroutine f90wrap_intf_gen__get_terminations__binding__agt

subroutine f90wrap_retrieve_last_generated_structures(num_structures, structures)
    use artemis__geom_rw, only: basis_type
    use artemis__structure_cache, only: retrieve_last_generated_structures
    implicit none

    type basis_type_xnum_array
        type(basis_type), dimension(:), allocatable :: items
    end type basis_type_xnum_array

    type basis_type_xnum_array_ptr_type
        type(basis_type_xnum_array), pointer :: p => NULL()
    end type basis_type_xnum_array_ptr_type
    integer, intent(in) :: num_structures
    integer, intent(inout), dimension(2) :: structures
    type(basis_type_xnum_array_ptr_type) :: structures_ptr

    structures_ptr = transfer(structures, structures_ptr)
    structures_ptr%p%items = retrieve_last_generated_structures()
    structures = transfer(structures_ptr, structures)
end subroutine f90wrap_retrieve_last_generated_structures





subroutine f90wrap_intf_gen__generate__binding__agt( &
    this, surface_lw, surface_up, &
    thickness_lw, thickness_up, &
    num_layers_lw, num_layers_up, &
    print_lattice_match_info, print_termination_info, print_shift_info, &
    break_on_fail, icheck_match, interface_idx, &
    generate_structures, &
    seed, verbose, exit_code, &
    n0, n1)
    use artemis__generator, only: artemis_generator_type
    implicit none
    
    type artemis_generator_type_ptr_type
        type(artemis_generator_type), pointer :: p => NULL()
    end type artemis_generator_type_ptr_type
    type(artemis_generator_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    integer, intent(in), optional, dimension(n0) :: surface_lw
    integer, intent(in), optional, dimension(n1) :: surface_up
    real(4), intent(in), optional :: thickness_lw
    real(4), intent(in), optional :: thickness_up
    integer, intent(in), optional :: num_layers_lw
    integer, intent(in), optional :: num_layers_up
    logical, intent(in), optional :: print_lattice_match_info
    logical, intent(in), optional :: print_termination_info
    logical, intent(in), optional :: print_shift_info
    logical, intent(in), optional :: break_on_fail
    integer, intent(in), optional :: icheck_match
    integer, intent(in), optional :: interface_idx
    logical, intent(in), optional :: generate_structures
    integer, intent(in), optional :: seed
    integer, intent(in), optional :: verbose
    integer, optional, intent(inout) :: exit_code
    integer :: n0
    !f2py intent(hide), depend(surface_lw) :: n0 = shape(surface_lw,0)
    integer :: n1
    !f2py intent(hide), depend(surface_up) :: n1 = shape(surface_up,0)
    this_ptr = transfer(this, this_ptr)
    call this_ptr%p%generate(surface_lw=surface_lw, surface_up=surface_up, thickness_lw=thickness_lw, &
        thickness_up=thickness_up, num_layers_lw=num_layers_lw, num_layers_up=num_layers_up, &
        print_lattice_match_info=print_lattice_match_info, print_termination_info=print_termination_info, &
        print_shift_info=print_shift_info, break_on_fail=break_on_fail, icheck_match=icheck_match, &
        interface_idx=interface_idx, generate_structures=generate_structures, seed=seed, verbose=verbose, &
        exit_code=exit_code)
end subroutine f90wrap_intf_gen__generate__binding__agt

subroutine f90wrap_intf_gen__restart__binding__agt(this, basis, interface_location, &
    print_shift_info, seed, verbose, exit_code)
    use artemis__generator, only: artemis_generator_type
    use artemis__geom_rw, only: basis_type
    implicit none
    
    type basis_type_ptr_type
        type(basis_type), pointer :: p => NULL()
    end type basis_type_ptr_type
    type artemis_generator_type_ptr_type
        type(artemis_generator_type), pointer :: p => NULL()
    end type artemis_generator_type_ptr_type
    type(artemis_generator_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    type(basis_type_ptr_type) :: basis_ptr
    integer, intent(in), dimension(2) :: basis
    real(4), dimension(2), intent(in), optional :: interface_location
    logical, intent(in), optional :: print_shift_info
    integer, intent(in), optional :: seed
    integer, intent(in), optional :: verbose
    integer, optional, intent(inout) :: exit_code
    this_ptr = transfer(this, this_ptr)
    basis_ptr = transfer(basis, basis_ptr)
    call this_ptr%p%restart(basis=basis_ptr%p, interface_location=interface_location, print_shift_info=print_shift_info, &
        seed=seed, verbose=verbose, exit_code=exit_code)
end subroutine f90wrap_intf_gen__restart__binding__agt

subroutine f90wrap_intf_gen__get_structures__binding__agt(this, ret_structures)
    use artemis__generator, only: artemis_generator_type
    use artemis__geom_rw, only: basis_type
    implicit none

    type artemis_generator_type_ptr_type
        type(artemis_generator_type), pointer :: p => NULL()
    end type artemis_generator_type_ptr_type

    type basis_type_xnum_array
        type(basis_type), dimension(:), allocatable :: items
    end type basis_type_xnum_array

    type basis_type_xnum_array_ptr_type
        type(basis_type_xnum_array), pointer :: p => NULL()
    end type basis_type_xnum_array_ptr_type
    type(artemis_generator_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    integer, intent(out), dimension(2) :: ret_structures
    type(basis_type_xnum_array_ptr_type) :: ret_structures_ptr

    this_ptr = transfer(this, this_ptr)
    ret_structures_ptr%p%items = this_ptr%p%get_structures()
    ret_structures = transfer(ret_structures_ptr,ret_structures)
end subroutine f90wrap_intf_gen__get_structures__binding__agt


!###############################################################################
! generated structures handling
!###############################################################################
subroutine f90wrap_artemis_gen_type__array_getitem__structures( &
     f90wrap_this, f90wrap_i, structuresitem &
)
    
    use artemis__generator, only: artemis_generator_type
    use artemis__geom_rw, only: basis_type
    implicit none
    
    type artemis_generator_type_ptr_type
        type(artemis_generator_type), pointer :: p => NULL()
    end type artemis_generator_type_ptr_type
    type basis_type_ptr_type
        type(basis_type), pointer :: p => NULL()
    end type basis_type_ptr_type
    integer, intent(in) :: f90wrap_this(2)
    type(artemis_generator_type_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_i
    integer, intent(out) :: structuresitem(2)
    type(basis_type_ptr_type) :: structures_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%structures)) then
        if (f90wrap_i < 1 .or. f90wrap_i > size(this_ptr%p%structures)) then
            call f90wrap_abort("array index out of range")
        else
            structures_ptr%p => this_ptr%p%structures(f90wrap_i)
            structuresitem = transfer(structures_ptr,structuresitem)
        endif
    else
        call f90wrap_abort("derived type array not allocated")
    end if
end subroutine f90wrap_artemis_gen_type__array_getitem__structures

subroutine f90wrap_artemis_gen_type__array_setitem__structures( &
     f90wrap_this, f90wrap_i, structuresitem &
)
    
    use artemis__generator, only: artemis_generator_type
    use artemis__geom_rw, only: basis_type
    implicit none
    
    type artemis_generator_type_ptr_type
        type(artemis_generator_type), pointer :: p => NULL()
    end type artemis_generator_type_ptr_type
    type basis_type_ptr_type
        type(basis_type), pointer :: p => NULL()
    end type basis_type_ptr_type
    integer, intent(in) :: f90wrap_this(2)
    type(artemis_generator_type_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_i
    integer, intent(in) :: structuresitem(2)
    type(basis_type_ptr_type) :: structures_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%structures)) then
        if (f90wrap_i < 1 .or. f90wrap_i > size(this_ptr%p%structures)) then
            call f90wrap_abort("array index out of range")
        else
            structures_ptr = transfer(structuresitem,structures_ptr)
            this_ptr%p%structures(f90wrap_i) = structures_ptr%p
        endif
    else
        call f90wrap_abort("derived type array not allocated")
    end if
end subroutine f90wrap_artemis_gen_type__array_setitem__structures

subroutine f90wrap_artemis_gen_type__array_len__structures( &
     f90wrap_this, f90wrap_n &
)
    
    use artemis__generator, only: artemis_generator_type
    use artemis__geom_rw, only: basis_type
    implicit none
    
    type artemis_generator_type_ptr_type
        type(artemis_generator_type), pointer :: p => NULL()
    end type artemis_generator_type_ptr_type
    integer, intent(out) :: f90wrap_n
    integer, intent(in) :: f90wrap_this(2)
    type(artemis_generator_type_ptr_type) :: this_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%structures)) then
        f90wrap_n = size(this_ptr%p%structures)
    else
        f90wrap_n = 0
    end if
end subroutine f90wrap_artemis_gen_type__array_len__structures
!###############################################################################

! End of module artemis__generator defined in file ../src/fortran/lib/mod_intf_generator.f90

