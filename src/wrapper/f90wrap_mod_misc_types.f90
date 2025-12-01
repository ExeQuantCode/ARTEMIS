! Module artemis__misc_types defined in file ../fortran/lib/mod_misc_types.f90

subroutine f90wrap_struc_data_type__get__match_idx(this, f90wrap_match_idx)
    use artemis__misc_types, only: struc_data_type
    implicit none
    type struc_data_type_ptr_type
        type(struc_data_type), pointer :: p => NULL()
    end type struc_data_type_ptr_type
    integer, intent(in)   :: this(2)
    type(struc_data_type_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_match_idx
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_match_idx = this_ptr%p%match_idx
end subroutine f90wrap_struc_data_type__get__match_idx

subroutine f90wrap_struc_data_type__set__match_idx(this, f90wrap_match_idx)
    use artemis__misc_types, only: struc_data_type
    implicit none
    type struc_data_type_ptr_type
        type(struc_data_type), pointer :: p => NULL()
    end type struc_data_type_ptr_type
    integer, intent(in)   :: this(2)
    type(struc_data_type_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_match_idx
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%match_idx = f90wrap_match_idx
end subroutine f90wrap_struc_data_type__set__match_idx

subroutine f90wrap_struc_data_type__get__shift_idx(this, f90wrap_shift_idx)
    use artemis__misc_types, only: struc_data_type
    implicit none
    type struc_data_type_ptr_type
        type(struc_data_type), pointer :: p => NULL()
    end type struc_data_type_ptr_type
    integer, intent(in)   :: this(2)
    type(struc_data_type_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_shift_idx
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_shift_idx = this_ptr%p%shift_idx
end subroutine f90wrap_struc_data_type__get__shift_idx

subroutine f90wrap_struc_data_type__set__shift_idx(this, f90wrap_shift_idx)
    use artemis__misc_types, only: struc_data_type
    implicit none
    type struc_data_type_ptr_type
        type(struc_data_type), pointer :: p => NULL()
    end type struc_data_type_ptr_type
    integer, intent(in)   :: this(2)
    type(struc_data_type_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_shift_idx
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%shift_idx = f90wrap_shift_idx
end subroutine f90wrap_struc_data_type__set__shift_idx

subroutine f90wrap_struc_data_type__get__swap_idx(this, f90wrap_swap_idx)
    use artemis__misc_types, only: struc_data_type
    implicit none
    type struc_data_type_ptr_type
        type(struc_data_type), pointer :: p => NULL()
    end type struc_data_type_ptr_type
    integer, intent(in)   :: this(2)
    type(struc_data_type_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_swap_idx
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_swap_idx = this_ptr%p%swap_idx
end subroutine f90wrap_struc_data_type__get__swap_idx

subroutine f90wrap_struc_data_type__set__swap_idx(this, f90wrap_swap_idx)
    use artemis__misc_types, only: struc_data_type
    implicit none
    type struc_data_type_ptr_type
        type(struc_data_type), pointer :: p => NULL()
    end type struc_data_type_ptr_type
    integer, intent(in)   :: this(2)
    type(struc_data_type_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_swap_idx
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%swap_idx = f90wrap_swap_idx
end subroutine f90wrap_struc_data_type__set__swap_idx

subroutine f90wrap_struc_data_type__get__from_pricel_lw(this, f90wrap_from_pricel_lw)
    use artemis__misc_types, only: struc_data_type
    implicit none
    type struc_data_type_ptr_type
        type(struc_data_type), pointer :: p => NULL()
    end type struc_data_type_ptr_type
    integer, intent(in)   :: this(2)
    type(struc_data_type_ptr_type) :: this_ptr
    logical, intent(out) :: f90wrap_from_pricel_lw
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_from_pricel_lw = this_ptr%p%from_pricel_lw
end subroutine f90wrap_struc_data_type__get__from_pricel_lw

subroutine f90wrap_struc_data_type__set__from_pricel_lw(this, f90wrap_from_pricel_lw)
    use artemis__misc_types, only: struc_data_type
    implicit none
    type struc_data_type_ptr_type
        type(struc_data_type), pointer :: p => NULL()
    end type struc_data_type_ptr_type
    integer, intent(in)   :: this(2)
    type(struc_data_type_ptr_type) :: this_ptr
    logical, intent(in) :: f90wrap_from_pricel_lw
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%from_pricel_lw = f90wrap_from_pricel_lw
end subroutine f90wrap_struc_data_type__set__from_pricel_lw

subroutine f90wrap_struc_data_type__get__from_pricel_up(this, f90wrap_from_pricel_up)
    use artemis__misc_types, only: struc_data_type
    implicit none
    type struc_data_type_ptr_type
        type(struc_data_type), pointer :: p => NULL()
    end type struc_data_type_ptr_type
    integer, intent(in)   :: this(2)
    type(struc_data_type_ptr_type) :: this_ptr
    logical, intent(out) :: f90wrap_from_pricel_up
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_from_pricel_up = this_ptr%p%from_pricel_up
end subroutine f90wrap_struc_data_type__get__from_pricel_up

subroutine f90wrap_struc_data_type__set__from_pricel_up(this, f90wrap_from_pricel_up)
    use artemis__misc_types, only: struc_data_type
    implicit none
    type struc_data_type_ptr_type
        type(struc_data_type), pointer :: p => NULL()
    end type struc_data_type_ptr_type
    integer, intent(in)   :: this(2)
    type(struc_data_type_ptr_type) :: this_ptr
    logical, intent(in) :: f90wrap_from_pricel_up
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%from_pricel_up = f90wrap_from_pricel_up
end subroutine f90wrap_struc_data_type__set__from_pricel_up

subroutine f90wrap_struc_data_type__array__term_lw_idx(this, nd, dtype, dshape, dloc)
    use artemis__misc_types, only: struc_data_type
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type struc_data_type_ptr_type
        type(struc_data_type), pointer :: p => NULL()
    end type struc_data_type_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(struc_data_type_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    dshape(1:1) = shape(this_ptr%p%term_lw_idx)
    dloc = loc(this_ptr%p%term_lw_idx)
end subroutine f90wrap_struc_data_type__array__term_lw_idx

subroutine f90wrap_struc_data_type__array__term_up_idx(this, nd, dtype, dshape, dloc)
    use artemis__misc_types, only: struc_data_type
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type struc_data_type_ptr_type
        type(struc_data_type), pointer :: p => NULL()
    end type struc_data_type_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(struc_data_type_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    dshape(1:1) = shape(this_ptr%p%term_up_idx)
    dloc = loc(this_ptr%p%term_up_idx)
end subroutine f90wrap_struc_data_type__array__term_up_idx

subroutine f90wrap_struc_data_type__array__transform_lw(this, nd, dtype, dshape, dloc)
    use artemis__misc_types, only: struc_data_type
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type struc_data_type_ptr_type
        type(struc_data_type), pointer :: p => NULL()
    end type struc_data_type_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(struc_data_type_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    dshape(1:2) = shape(this_ptr%p%transform_lw)
    dloc = loc(this_ptr%p%transform_lw)
end subroutine f90wrap_struc_data_type__array__transform_lw

subroutine f90wrap_struc_data_type__array__transform_up(this, nd, dtype, dshape, dloc)
    use artemis__misc_types, only: struc_data_type
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type struc_data_type_ptr_type
        type(struc_data_type), pointer :: p => NULL()
    end type struc_data_type_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(struc_data_type_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    dshape(1:2) = shape(this_ptr%p%transform_up)
    dloc = loc(this_ptr%p%transform_up)
end subroutine f90wrap_struc_data_type__array__transform_up

subroutine f90wrap_struc_data_type__get__approx_thickness_lw(this, f90wrap_approx_thickness_lw)
    use artemis__misc_types, only: struc_data_type
    implicit none
    type struc_data_type_ptr_type
        type(struc_data_type), pointer :: p => NULL()
    end type struc_data_type_ptr_type
    integer, intent(in)   :: this(2)
    type(struc_data_type_ptr_type) :: this_ptr
    real(4), intent(out) :: f90wrap_approx_thickness_lw
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_approx_thickness_lw = this_ptr%p%approx_thickness_lw
end subroutine f90wrap_struc_data_type__get__approx_thickness_lw

subroutine f90wrap_struc_data_type__set__approx_thickness_lw(this, f90wrap_approx_thickness_lw)
    use artemis__misc_types, only: struc_data_type
    implicit none
    type struc_data_type_ptr_type
        type(struc_data_type), pointer :: p => NULL()
    end type struc_data_type_ptr_type
    integer, intent(in)   :: this(2)
    type(struc_data_type_ptr_type) :: this_ptr
    real(4), intent(in) :: f90wrap_approx_thickness_lw
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%approx_thickness_lw = f90wrap_approx_thickness_lw
end subroutine f90wrap_struc_data_type__set__approx_thickness_lw

subroutine f90wrap_struc_data_type__get__approx_thickness_up(this, f90wrap_approx_thickness_up)
    use artemis__misc_types, only: struc_data_type
    implicit none
    type struc_data_type_ptr_type
        type(struc_data_type), pointer :: p => NULL()
    end type struc_data_type_ptr_type
    integer, intent(in)   :: this(2)
    type(struc_data_type_ptr_type) :: this_ptr
    real(4), intent(out) :: f90wrap_approx_thickness_up
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_approx_thickness_up = this_ptr%p%approx_thickness_up
end subroutine f90wrap_struc_data_type__get__approx_thickness_up

subroutine f90wrap_struc_data_type__set__approx_thickness_up(this, f90wrap_approx_thickness_up)
    use artemis__misc_types, only: struc_data_type
    implicit none
    type struc_data_type_ptr_type
        type(struc_data_type), pointer :: p => NULL()
    end type struc_data_type_ptr_type
    integer, intent(in)   :: this(2)
    type(struc_data_type_ptr_type) :: this_ptr
    real(4), intent(in) :: f90wrap_approx_thickness_up
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%approx_thickness_up = f90wrap_approx_thickness_up
end subroutine f90wrap_struc_data_type__set__approx_thickness_up

subroutine f90wrap_struc_data_type__array__mismatch(this, nd, dtype, dshape, dloc)
    use artemis__misc_types, only: struc_data_type
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type struc_data_type_ptr_type
        type(struc_data_type), pointer :: p => NULL()
    end type struc_data_type_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(struc_data_type_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 11
    this_ptr = transfer(this, this_ptr)
    dshape(1:1) = shape(this_ptr%p%mismatch)
    dloc = loc(this_ptr%p%mismatch)
end subroutine f90wrap_struc_data_type__array__mismatch

subroutine f90wrap_struc_data_type__array__shift(this, nd, dtype, dshape, dloc)
    use artemis__misc_types, only: struc_data_type
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type struc_data_type_ptr_type
        type(struc_data_type), pointer :: p => NULL()
    end type struc_data_type_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(struc_data_type_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 11
    this_ptr = transfer(this, this_ptr)
    dshape(1:1) = shape(this_ptr%p%shift)
    dloc = loc(this_ptr%p%shift)
end subroutine f90wrap_struc_data_type__array__shift

subroutine f90wrap_struc_data_type__get__swap_density(this, f90wrap_swap_density)
    use artemis__misc_types, only: struc_data_type
    implicit none
    type struc_data_type_ptr_type
        type(struc_data_type), pointer :: p => NULL()
    end type struc_data_type_ptr_type
    integer, intent(in)   :: this(2)
    type(struc_data_type_ptr_type) :: this_ptr
    real(4), intent(out) :: f90wrap_swap_density
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_swap_density = this_ptr%p%swap_density
end subroutine f90wrap_struc_data_type__get__swap_density

subroutine f90wrap_struc_data_type__set__swap_density(this, f90wrap_swap_density)
    use artemis__misc_types, only: struc_data_type
    implicit none
    type struc_data_type_ptr_type
        type(struc_data_type), pointer :: p => NULL()
    end type struc_data_type_ptr_type
    integer, intent(in)   :: this(2)
    type(struc_data_type_ptr_type) :: this_ptr
    real(4), intent(in) :: f90wrap_swap_density
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%swap_density = f90wrap_swap_density
end subroutine f90wrap_struc_data_type__set__swap_density

subroutine f90wrap_struc_data_type__array__approx_eff_swap_conc(this, nd, dtype, dshape, dloc)
    use artemis__misc_types, only: struc_data_type
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type struc_data_type_ptr_type
        type(struc_data_type), pointer :: p => NULL()
    end type struc_data_type_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(struc_data_type_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 11
    this_ptr = transfer(this, this_ptr)
    dshape(1:1) = shape(this_ptr%p%approx_eff_swap_conc)
    dloc = loc(this_ptr%p%approx_eff_swap_conc)
end subroutine f90wrap_struc_data_type__array__approx_eff_swap_conc

subroutine f90wrap_misc_types__struc_data_type_initialise(this)
    use artemis__misc_types, only: struc_data_type
    implicit none
    
    type struc_data_type_ptr_type
        type(struc_data_type), pointer :: p => NULL()
    end type struc_data_type_ptr_type
    type(struc_data_type_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_misc_types__struc_data_type_initialise

subroutine f90wrap_misc_types__struc_data_type_finalise(this)
    use artemis__misc_types, only: struc_data_type
    implicit none
    
    type struc_data_type_ptr_type
        type(struc_data_type), pointer :: p => NULL()
    end type struc_data_type_ptr_type
    type(struc_data_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_misc_types__struc_data_type_finalise

! End of module artemis__misc_types defined in file ../fortran/lib/mod_misc_types.f90

