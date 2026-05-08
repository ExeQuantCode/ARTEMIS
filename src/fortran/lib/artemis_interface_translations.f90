module artemis__interface_translations
  !! Module for detecting primitive relative shift vectors in interface cells.
  use coreutils__kind, only: real32
  use coreutils__linalg, only: inverse_3x3, cross
  use atomstruc, only: basis_type
  use artemis__geom_utils, only: split_bas
  use artemis__interface_identifier, only: intf_info_type, get_interface
  use artemis__misc_linalg, only: get_frac_denom, lcm
  use artemis__sym, only: confine_type, gldfnd
  implicit none

  private
  public :: get_interface_translations

  real(real32), parameter :: interface_translation_cart_tol = 7.5e-2_real32

contains
!###############################################################################
  subroutine get_interface_translations(basis, t1, t2)
    !! Return the two smallest independent relative in-plane shifts for an interface.
    implicit none

    type(basis_type), intent(in) :: basis
    real(real32), dimension(3), intent(out) :: t1, t2

    integer :: analytical_capacity
    integer :: lower_candidate_capacity, upper_candidate_capacity
    integer :: lower_valid_capacity, upper_valid_capacity
    integer :: lower_group_capacity, upper_group_capacity
    integer :: relative_capacity, valid_capacity
    integer :: n_lower_basis, n_upper_basis
    integer :: n_lower_candidates, n_upper_candidates
    integer :: n_lower_valid, n_upper_valid
    integer :: n_lower_group, n_upper_group
    integer :: n_relative, n_valid
    integer, dimension(2) :: plane_axes
    integer, dimension(3) :: work_axes
    real(real32) :: frac_tol
    real(real32), allocatable, dimension(:,:) :: regions
    real(real32), allocatable, dimension(:,:,:) :: frac_lower, frac_upper
    real(real32), allocatable, dimension(:,:) :: lower_basis, upper_basis
    real(real32), allocatable, dimension(:,:) :: lower_candidates, upper_candidates
    real(real32), allocatable, dimension(:,:) :: lower_valid, upper_valid
    real(real32), allocatable, dimension(:,:) :: lower_group, upper_group
    real(real32), allocatable, dimension(:,:) :: relative_candidates, valid_vectors
    real(real32), dimension(3,3) :: lat_work
    real(real32), dimension(3) :: candidate, work_t1, work_t2
    integer, allocatable, dimension(:) :: lower_spec_num, upper_spec_num
    logical :: lower_exact, upper_exact
    type(basis_type) :: basis_frac
    type(basis_type), allocatable, dimension(:) :: split_basis
    type(intf_info_type) :: interface_info


    call copy_basis_fractional(basis, basis_frac)
    interface_info = get_interface(basis)
    interface_info%loc = interface_info%loc / norm2(basis%lat(interface_info%axis, :))
    call get_plane_axes(interface_info%axis, plane_axes)
    work_axes = [ plane_axes(1), plane_axes(2), interface_info%axis ]
    lat_work(1, :) = basis%lat(work_axes(1), :)
    lat_work(2, :) = basis%lat(work_axes(2), :)
    lat_work(3, :) = basis%lat(work_axes(3), :)

    allocate(regions(2, 2))
    regions(1, 1:2) = interface_info%loc(1:2)
    regions(2, 1) = interface_info%loc(2)
    regions(2, 2) = interface_info%loc(1)
    split_basis = split_bas(basis_frac, regions, interface_info%axis)

    frac_tol = get_fractional_tolerance(lat_work, interface_translation_cart_tol)

    analytical_capacity = max(1, split_basis(2)%natom)
    call infer_inplane_translations_gldfnd(split_basis(2), interface_info%axis, &
         work_axes, frac_tol, analytical_capacity, lower_basis, n_lower_basis, &
         lower_exact)
    analytical_capacity = max(1, split_basis(1)%natom)
    call infer_inplane_translations_gldfnd(split_basis(1), interface_info%axis, &
         work_axes, frac_tol, analytical_capacity, upper_basis, n_upper_basis, &
         upper_exact)

    if (lower_exact .and. upper_exact) then
       call build_relative_candidates_from_bases(&
            lower_basis, n_lower_basis, upper_basis, n_upper_basis, &
            frac_tol, relative_candidates, n_relative &
       )
       if (n_relative.ge.2) then
          call select_translation_pair_from_group(lat_work, relative_candidates, &
               n_relative, frac_tol, work_t1, work_t2)
          call unpack_translation(work_t1, work_axes, t1)
          call unpack_translation(work_t2, work_axes, t2)
          return
       end if
    end if

    call cache_fractional_positions(split_basis(2), work_axes, frac_lower, &
         lower_spec_num)
    call cache_fractional_positions(split_basis(1), work_axes, frac_upper, &
         upper_spec_num)

    lower_candidate_capacity = get_candidate_capacity(lower_spec_num)
    upper_candidate_capacity = get_candidate_capacity(upper_spec_num)
    lower_valid_capacity = max(8, 2 * lower_candidate_capacity + 8)
    upper_valid_capacity = max(8, 2 * upper_candidate_capacity + 8)

    allocate(lower_candidates(lower_candidate_capacity, 3))
    allocate(upper_candidates(upper_candidate_capacity, 3))
    allocate(lower_valid(lower_valid_capacity, 3))
    allocate(upper_valid(upper_valid_capacity, 3))
    lower_candidates = 0._real32
    upper_candidates = 0._real32
    lower_valid = 0._real32
    upper_valid = 0._real32
    n_lower_candidates = 0
    n_upper_candidates = 0
    n_lower_valid = 0
    n_upper_valid = 0

    call collect_pair_candidates(frac_lower, lower_spec_num, lower_candidates, &
         n_lower_candidates, frac_tol)
    call collect_pair_candidates(frac_upper, upper_spec_num, upper_candidates, &
         n_upper_candidates, frac_tol)

    candidate = [ 0.5_real32, 0.0_real32, 0.0_real32 ]
    call append_unique_translation(lower_candidates, n_lower_candidates, candidate, &
         frac_tol)
    call append_unique_translation(upper_candidates, n_upper_candidates, candidate, &
         frac_tol)
    candidate = [ 0.0_real32, 0.5_real32, 0.0_real32 ]
    call append_unique_translation(lower_candidates, n_lower_candidates, candidate, &
         frac_tol)
    call append_unique_translation(upper_candidates, n_upper_candidates, candidate, &
         frac_tol)
    candidate = [ 1.0_real32, 0.0_real32, 0.0_real32 ]
    call append_unique_translation(lower_candidates, n_lower_candidates, candidate, &
         frac_tol)
    call append_unique_translation(upper_candidates, n_upper_candidates, candidate, &
         frac_tol)
    candidate = [ 0.0_real32, 1.0_real32, 0.0_real32 ]
    call append_unique_translation(lower_candidates, n_lower_candidates, candidate, &
         frac_tol)
    call append_unique_translation(upper_candidates, n_upper_candidates, candidate, &
         frac_tol)

    call collect_valid_translations(&
         frac_lower, lower_spec_num, lat_work, &
         lower_candidates, n_lower_candidates, &
         lower_valid, n_lower_valid, &
         frac_tol, interface_translation_cart_tol &
    )
    call expand_valid_translations(&
         frac_lower, lower_spec_num, lat_work, &
         lower_valid, n_lower_valid, &
         frac_tol, interface_translation_cart_tol &
    )

    call collect_valid_translations(&
         frac_upper, upper_spec_num, lat_work, &
         upper_candidates, n_upper_candidates, &
         upper_valid, n_upper_valid, &
         frac_tol, interface_translation_cart_tol &
    )
    call expand_valid_translations(&
         frac_upper, upper_spec_num, lat_work, &
         upper_valid, n_upper_valid, &
         frac_tol, interface_translation_cart_tol &
    )

    if (allocated(relative_candidates)) deallocate(relative_candidates)
    call build_relative_candidates_from_bases(&
         lower_valid, n_lower_valid, upper_valid, n_upper_valid, frac_tol, &
         relative_candidates, n_relative)
    if (n_relative.ge.2) then
       call select_translation_pair_from_group(lat_work, relative_candidates, &
            n_relative, frac_tol, work_t1, work_t2)
       call unpack_translation(work_t1, work_axes, t1)
       call unpack_translation(work_t2, work_axes, t2)
       return
    end if

    lower_group_capacity = max(64, lower_candidate_capacity)
    upper_group_capacity = max(64, upper_candidate_capacity)
    allocate(lower_group(lower_group_capacity, 3))
    allocate(upper_group(upper_group_capacity, 3))
    lower_group = 0._real32
    upper_group = 0._real32
    call build_translation_group(lower_valid, n_lower_valid, lower_group, &
         n_lower_group, frac_tol)
    call build_translation_group(upper_valid, n_upper_valid, upper_group, &
         n_upper_group, frac_tol)

    if (allocated(relative_candidates)) deallocate(relative_candidates)
    if (allocated(relative_candidates)) deallocate(relative_candidates)
    call build_relative_candidates_from_bases(&
         lower_group, n_lower_group, upper_group, n_upper_group, frac_tol, &
         relative_candidates, n_relative)
    if (n_relative.ge.2) then
       call select_translation_pair_from_group(lat_work, relative_candidates, &
            n_relative, frac_tol, work_t1, work_t2)
       call unpack_translation(work_t1, work_axes, t1)
       call unpack_translation(work_t2, work_axes, t2)
       return
    end if

    if (allocated(relative_candidates)) deallocate(relative_candidates)
    relative_capacity = max(8, n_lower_group * n_upper_group)
    valid_capacity = max(8, 8 * relative_capacity + 8)
    allocate(relative_candidates(relative_capacity, 3))
    allocate(valid_vectors(valid_capacity, 3))
    relative_candidates = 0._real32
    valid_vectors = 0._real32
    n_relative = 0
    n_valid = 0

    call collect_relative_candidates(&
         lower_group, n_lower_group, upper_group, n_upper_group, &
         relative_candidates, n_relative, frac_tol &
    )
    if (n_relative.gt.0) then
       valid_vectors(1:n_relative, :) = relative_candidates(1:n_relative, :)
       n_valid = n_relative
    end if
    call select_translation_pair(lat_work, valid_vectors, n_valid, frac_tol, &
         work_t1, work_t2)
    call unpack_translation(work_t1, work_axes, t1)
    call unpack_translation(work_t2, work_axes, t2)

  end subroutine get_interface_translations
!###############################################################################


!###############################################################################
  subroutine cache_fractional_positions(basis, work_axes, frac_atoms, spec_num)
    !! Cache basis coordinates in fractional coordinates.
    implicit none

    type(basis_type), intent(in) :: basis
    integer, dimension(3), intent(in) :: work_axes
    real(real32), allocatable, dimension(:,:,:), intent(out) :: frac_atoms
    integer, allocatable, dimension(:), intent(out) :: spec_num

    integer :: is, ia, max_atoms
    real(real32), dimension(3,3) :: invlat
    real(real32), dimension(3) :: atom_frac


    allocate(spec_num(basis%nspec))
    spec_num(:) = basis%spec(:)%num
    max_atoms = max(1, maxval(spec_num))
    allocate(frac_atoms(basis%nspec, max_atoms, 3))
    frac_atoms = 0._real32

    invlat = 0._real32
    if (basis%lcart) invlat = inverse_3x3(basis%lat)

    do is = 1, basis%nspec
       do ia = 1, spec_num(is)
          atom_frac(:) = basis%spec(is)%atom(ia, 1:3)
          if (basis%lcart) atom_frac(:) = matmul(atom_frac, invlat)
          call wrap_zero_one(atom_frac)
          frac_atoms(is, ia, 1) = atom_frac(work_axes(1))
          frac_atoms(is, ia, 2) = atom_frac(work_axes(2))
          frac_atoms(is, ia, 3) = atom_frac(work_axes(3))
       end do
    end do

  end subroutine cache_fractional_positions
!###############################################################################


!###############################################################################
  subroutine copy_basis_fractional(input_basis, output_basis)
    !! Copy a basis and convert coordinates to fractional if needed.
    implicit none

    type(basis_type), intent(in) :: input_basis
    type(basis_type), intent(out) :: output_basis

    integer :: is, ia
    real(real32), dimension(3,3) :: invlat


    call output_basis%copy(input_basis)
    if (.not. output_basis%lcart) return

    invlat = inverse_3x3(output_basis%lat)
    do is = 1, output_basis%nspec
       do ia = 1, output_basis%spec(is)%num
          output_basis%spec(is)%atom(ia, 1:3) = &
               matmul(output_basis%spec(is)%atom(ia, 1:3), invlat)
          call wrap_zero_one(output_basis%spec(is)%atom(ia, 1:3))
       end do
    end do
    output_basis%lcart = .false.

  end subroutine copy_basis_fractional
!###############################################################################


!###############################################################################
  subroutine get_plane_axes(axis, plane_axes)
    !! Return the two in-plane axes for a given interface normal.
    implicit none

    integer, intent(in) :: axis
    integer, dimension(2), intent(out) :: plane_axes


    select case(axis)
    case(1)
       plane_axes = [ 2, 3 ]
    case(2)
       plane_axes = [ 1, 3 ]
    case default
       plane_axes = [ 1, 2 ]
    end select

  end subroutine get_plane_axes
!###############################################################################


!###############################################################################
  integer function get_candidate_capacity(spec_num) result(capacity)
    !! Return a conservative candidate buffer size for one slab.
    implicit none

    integer, dimension(:), intent(in) :: spec_num
    integer :: natom


    natom = sum(spec_num)
    capacity = max(8, natom * max(1, natom) + 8)

  end function get_candidate_capacity
!###############################################################################


!###############################################################################
  real(real32) function get_fractional_tolerance(lat, cart_tol) result(frac_tol)
  !! Convert a Cartesian tolerance into a fractional in-plane tolerance.
  implicit none

  real(real32), dimension(3,3), intent(in) :: lat
  real(real32), intent(in) :: cart_tol
  real(real32) :: min_len


  min_len = min(norm2(lat(1, :)), norm2(lat(2, :)))
  frac_tol = max(cart_tol / max(min_len, 1.0e-12_real32), 1.0e-6_real32)

  end function get_fractional_tolerance
!###############################################################################


!###############################################################################
  subroutine infer_inplane_translations_gldfnd(basis, axis, work_axes, frac_tol, &
       capacity, translations, n_translations, success)
    !! Infer in-plane slab translations directly from the slab symmetry.
    implicit none

    type(basis_type), intent(in) :: basis
    integer, intent(in) :: axis, capacity
    integer, dimension(3), intent(in) :: work_axes
    real(real32), intent(in) :: frac_tol
    real(real32), allocatable, dimension(:,:), intent(out) :: translations
    integer, intent(out) :: n_translations
    logical, intent(out) :: success

    integer :: i, n_raw
    real(real32), allocatable, dimension(:,:) :: raw_translations
    real(real32), dimension(3) :: candidate
    type(confine_type) :: confine


    allocate(raw_translations(max(1, capacity), 3))
    raw_translations = 0._real32
    allocate(translations(max(1, capacity), 3))
    translations = 0._real32

    confine%l = .true.
    confine%axis = axis
    confine%laxis = .false.
    confine%laxis(axis) = .true.

    n_translations = 0
    call gldfnd(confine, basis, basis, raw_translations, n_raw, frac_tol)

    do i = 1, n_raw
       candidate = 0._real32
       candidate(1) = raw_translations(i, work_axes(1))
       candidate(2) = raw_translations(i, work_axes(2))
       call append_unique_translation(translations, n_translations, candidate, &
            frac_tol)
    end do

    success = n_translations.ge.2

  end subroutine infer_inplane_translations_gldfnd
!###############################################################################


!###############################################################################
  subroutine build_relative_candidates_from_bases(&
       lower_basis, n_lower_basis, upper_basis, n_upper_basis, frac_tol, &
       candidates, n_candidates)
    !! Build the exact relative subgroup directly from the slab generators.
    implicit none

    real(real32), dimension(:,:), intent(in) :: lower_basis, upper_basis
    integer, intent(in) :: n_lower_basis, n_upper_basis
    real(real32), intent(in) :: frac_tol
    real(real32), allocatable, dimension(:,:), intent(out) :: candidates
    integer, intent(out) :: n_candidates

    integer :: i, j, queue_head, queue_tail
    integer :: modulus, n_generators, next_x, next_y
    integer, allocatable, dimension(:,:) :: generators, queue
    logical :: success
    logical, allocatable, dimension(:,:) :: seen
    real(real32), dimension(3) :: candidate


    call integerise_translation_generators(&
         lower_basis, n_lower_basis, upper_basis, n_upper_basis, frac_tol, &
         modulus, generators, n_generators, success)
    allocate(candidates(max(1, modulus * modulus - 1), 3))
    candidates = 0._real32
    n_candidates = 0
    if (.not.success) return

    allocate(seen(0:modulus - 1, 0:modulus - 1))
    allocate(queue(max(1, modulus * modulus), 2))
    seen = .false.
    queue = 0
    queue_head = 1
    queue_tail = 1
    queue(1, :) = 0
    seen(0, 0) = .true.

    do while (queue_head.le.queue_tail)
       do i = 1, n_generators
          do j = -1, 1, 2
             next_x = modulo(queue(queue_head, 1) + j * generators(i, 1), modulus)
             next_y = modulo(queue(queue_head, 2) + j * generators(i, 2), modulus)
             if (seen(next_x, next_y)) cycle
             seen(next_x, next_y) = .true.
             queue_tail = queue_tail + 1
             queue(queue_tail, 1) = next_x
             queue(queue_tail, 2) = next_y
             candidate = 0._real32
             candidate(1) = real(next_x, real32) / real(modulus, real32)
             candidate(2) = real(next_y, real32) / real(modulus, real32)
             call append_unique_translation(candidates, n_candidates, candidate, &
                  frac_tol)
          end do
       end do
       queue_head = queue_head + 1
    end do

  end subroutine build_relative_candidates_from_bases
!###############################################################################


!###############################################################################
  subroutine integerise_translation_generators(&
       lower_basis, n_lower_basis, upper_basis, n_upper_basis, frac_tol, &
       modulus, generators, n_generators, success)
    !! Convert analytical fractional generators to an exact integer subgroup.
    implicit none

    real(real32), dimension(:,:), intent(in) :: lower_basis, upper_basis
    integer, intent(in) :: n_lower_basis, n_upper_basis
    real(real32), intent(in) :: frac_tol
    integer, intent(out) :: modulus, n_generators
    integer, allocatable, dimension(:,:), intent(out) :: generators
    logical, intent(out) :: success

    integer :: i, denominator


    modulus = 1
    success = .false.
    n_generators = 0
    allocate(generators(max(1, n_lower_basis + n_upper_basis), 2))
    generators = 0

    do i = 1, n_lower_basis
       call update_integer_modulus(lower_basis(i, :), frac_tol, modulus, success)
       if (.not.success) return
    end do
    do i = 1, n_upper_basis
       call update_integer_modulus(upper_basis(i, :), frac_tol, modulus, success)
       if (.not.success) return
    end do

    do i = 1, n_lower_basis
       call append_integer_generator(generators, n_generators, lower_basis(i, :), &
            modulus)
    end do
    do i = 1, n_upper_basis
       call append_integer_generator(generators, n_generators, upper_basis(i, :), &
            modulus)
    end do

    if (n_generators.eq.0) then
       success = .false.
       return
    end if

    denominator = modulus
    if (denominator.le.0) then
       success = .false.
       return
    end if
    success = .true.

  end subroutine integerise_translation_generators
!###############################################################################


!###############################################################################
  subroutine update_integer_modulus(translation, frac_tol, modulus, success)
    !! Update the common integer denominator for one analytical generator.
    implicit none

    real(real32), dimension(3), intent(in) :: translation
    real(real32), intent(in) :: frac_tol
    integer, intent(inout) :: modulus
    logical, intent(out) :: success

    integer :: axis, denominator
    real(real32) :: value


    success = .true.
    do axis = 1, 2
       value = abs(translation(axis))
       if (value.le.frac_tol) cycle
       denominator = get_frac_denom(value)
       if (denominator.le.0) then
          success = .false.
          return
       end if
       modulus = lcm(modulus, denominator)
    end do

  end subroutine update_integer_modulus
!###############################################################################


!###############################################################################
  subroutine append_integer_generator(generators, n_generators, translation, modulus)
    !! Add one integerised generator if it is non-zero and unique.
    implicit none

    integer, dimension(:,:), intent(inout) :: generators
    integer, intent(inout) :: n_generators
    real(real32), dimension(3), intent(in) :: translation
    integer, intent(in) :: modulus

    integer :: i, gx, gy


    gx = modulo(nint(translation(1) * real(modulus, real32)), modulus)
    gy = modulo(nint(translation(2) * real(modulus, real32)), modulus)
    if (gx.eq.0 .and. gy.eq.0) return

    do i = 1, n_generators
       if (generators(i, 1).eq.gx .and. generators(i, 2).eq.gy) return
    end do

    n_generators = n_generators + 1
    generators(n_generators, 1) = gx
    generators(n_generators, 2) = gy

  end subroutine append_integer_generator
!###############################################################################


!###############################################################################
  subroutine collect_pair_candidates(frac_atoms, spec_num, candidates, n_candidates, &
       frac_tol)
    !! Collect projected in-plane candidates from same-species atom pairs.
    implicit none

    real(real32), dimension(:,:,:), intent(in) :: frac_atoms
    integer, dimension(:), intent(in) :: spec_num
    real(real32), dimension(:,:), intent(inout) :: candidates
    integer, intent(inout) :: n_candidates
    real(real32), intent(in) :: frac_tol

    integer :: is, ia, ja
    real(real32), dimension(3) :: candidate


    do is = 1, size(spec_num)
       do ia = 1, spec_num(is)
          do ja = 1, spec_num(is)
             candidate(:) = frac_atoms(is, ia, :) - frac_atoms(is, ja, :)
             candidate(3) = 0._real32
             call append_unique_translation(candidates, n_candidates, candidate, &
                  frac_tol)
          end do
       end do
    end do

  end subroutine collect_pair_candidates
!###############################################################################


!###############################################################################
  subroutine collect_valid_translations(&
       frac_atoms, spec_num, lat, &
       candidates, n_candidates, &
       valid_vectors, n_valid, &
       frac_tol, cart_tol &
  )
    !! Validate candidate translations against the full structure.
    implicit none

    real(real32), dimension(:,:,:), intent(in) :: frac_atoms
    integer, dimension(:), intent(in) :: spec_num
    real(real32), dimension(3,3), intent(in) :: lat
    real(real32), dimension(:,:), intent(in) :: candidates
    integer, intent(in) :: n_candidates
    real(real32), dimension(:,:), intent(inout) :: valid_vectors
    integer, intent(inout) :: n_valid
    real(real32), intent(in) :: frac_tol, cart_tol

    integer :: i


    do i = 1, n_candidates
       if (.not. is_valid_translation(frac_atoms, spec_num, lat, candidates(i, :), &
            cart_tol)) cycle
       call append_unique_translation(valid_vectors, n_valid, candidates(i, :), &
            frac_tol)
    end do

  end subroutine collect_valid_translations
!###############################################################################


!###############################################################################
  subroutine build_translation_group(valid_vectors, n_valid, group, n_group, frac_tol)
    !! Build the full translation group generated by the valid slab translations.
    implicit none

    real(real32), dimension(:,:), intent(in) :: valid_vectors
    integer, intent(in) :: n_valid
    real(real32), dimension(:,:), intent(inout) :: group
    integer, intent(out) :: n_group
    real(real32), intent(in) :: frac_tol

    integer :: i, j
    real(real32), dimension(3) :: candidate


    n_group = 1
    group = 0._real32

    i = 1
    do while (i.le.n_group)
       do j = 1, n_valid
          candidate(:) = group(i, :) + valid_vectors(j, :)
          call canonicalise_translation(candidate, frac_tol, candidate)
          if (norm2(candidate(1:2)).le.frac_tol) candidate = 0._real32
          if (.not. vector_in_set(group, n_group, candidate, frac_tol)) then
             if (n_group.ge.size(group, dim=1)) &
                  error stop "translation group storage exhausted"
             n_group = n_group + 1
             group(n_group, :) = candidate(:)
          end if

          candidate(:) = group(i, :) - valid_vectors(j, :)
          call canonicalise_translation(candidate, frac_tol, candidate)
          if (norm2(candidate(1:2)).le.frac_tol) candidate = 0._real32
          if (.not. vector_in_set(group, n_group, candidate, frac_tol)) then
             if (n_group.ge.size(group, dim=1)) &
                  error stop "translation group storage exhausted"
             n_group = n_group + 1
             group(n_group, :) = candidate(:)
          end if
       end do
       i = i + 1
    end do

  end subroutine build_translation_group
!###############################################################################


!###############################################################################
  subroutine collect_relative_candidates(&
       lower_group, n_lower_group, upper_group, n_upper_group, &
       candidates, n_candidates, frac_tol &
  )
    !! Collect candidate relative shifts from upper and lower slab translations.
    implicit none

    real(real32), dimension(:,:), intent(in) :: lower_group, upper_group
    integer, intent(in) :: n_lower_group, n_upper_group
    real(real32), dimension(:,:), intent(inout) :: candidates
    integer, intent(out) :: n_candidates
    real(real32), intent(in) :: frac_tol

    integer :: i, j
    real(real32), dimension(3) :: candidate


    n_candidates = 0
    do i = 1, n_upper_group
       do j = 1, n_lower_group
          candidate(:) = upper_group(i, :) - lower_group(j, :)
          call append_unique_translation(candidates, n_candidates, candidate, &
               frac_tol)
       end do
    end do

  end subroutine collect_relative_candidates
!###############################################################################


!###############################################################################
  logical function is_valid_translation(frac_atoms, spec_num, lat, translation, &
       cart_tol)
    !! Check whether an in-plane translation maps each species onto itself.
    implicit none

    real(real32), dimension(:,:,:), intent(in) :: frac_atoms
    integer, dimension(:), intent(in) :: spec_num
    real(real32), dimension(3,3), intent(in) :: lat
    real(real32), dimension(3), intent(in) :: translation
    real(real32), intent(in) :: cart_tol

    integer :: is, ia, ja, best_idx, n_atom
    logical, allocatable, dimension(:) :: matched
    real(real32) :: best_dist, dist
    real(real32), dimension(3) :: shifted, diff, cart, shift


    shift(:) = translation(:)
    shift(3) = 0._real32
    is_valid_translation = .true.

    do is = 1, size(spec_num)
       n_atom = spec_num(is)
       if (n_atom.eq.0) cycle
       allocate(matched(n_atom))
       matched = .false.

       do ia = 1, n_atom
          shifted(:) = frac_atoms(is, ia, :) + shift(:)
          call wrap_zero_one(shifted)
          best_idx = 0
          best_dist = huge(0._real32)

          do ja = 1, n_atom
             if (matched(ja)) cycle
             diff(:) = frac_atoms(is, ja, :) - shifted(:)
             call wrap_half(diff)
             cart(:) = matmul(diff, lat)
             dist = norm2(cart)
             if (dist.lt.best_dist) then
                best_dist = dist
                best_idx = ja
             end if
          end do

          if (best_idx.eq.0 .or. best_dist.gt.cart_tol) then
             is_valid_translation = .false.
             deallocate(matched)
             return
          end if
          matched(best_idx) = .true.
       end do

       deallocate(matched)
    end do

  end function is_valid_translation
!###############################################################################


!###############################################################################
  subroutine append_unique_translation(storage, nstored, translation, frac_tol)
    !! Canonicalise and deduplicate a translation candidate.
    implicit none

    real(real32), dimension(:,:), intent(inout) :: storage
    integer, intent(inout) :: nstored
    real(real32), dimension(3), intent(in) :: translation
    real(real32), intent(in) :: frac_tol

    integer :: i
    real(real32), dimension(3) :: candidate


    call canonicalise_translation(translation, frac_tol, candidate)
    if (norm2(candidate(1:2)).le.frac_tol) return

    do i = 1, nstored
       if (all(abs(storage(i, 1:2) - candidate(1:2)).le.frac_tol)) return
    end do

    if (nstored.ge.size(storage, dim=1)) &
         error stop "interface translation storage exhausted"

    nstored = nstored + 1
    storage(nstored, :) = candidate(:)

  end subroutine append_unique_translation
!###############################################################################


!###############################################################################
  subroutine canonicalise_translation(translation, frac_tol, candidate)
    !! Reduce an in-plane translation to a deterministic representative.
    implicit none

    real(real32), dimension(3), intent(in) :: translation
    real(real32), intent(in) :: frac_tol
    real(real32), dimension(3), intent(out) :: candidate


    candidate(:) = 0._real32
    candidate(1:2) = translation(1:2)
    call wrap_half(candidate)
    candidate(3) = 0._real32

    if (candidate(1).lt.-frac_tol .or. &
         (abs(candidate(1)).le.frac_tol .and. candidate(2).lt.-frac_tol)) then
       candidate(:) = -candidate(:)
       candidate(3) = 0._real32
    end if

    if (abs(candidate(1)).le.frac_tol) candidate(1) = 0._real32
    if (abs(candidate(2)).le.frac_tol) candidate(2) = 0._real32

  end subroutine canonicalise_translation
!###############################################################################


!###############################################################################
  subroutine expand_valid_translations(frac_atoms, spec_num, lat, valid_vectors, &
       n_valid, frac_tol, cart_tol)
    !! Add equivalent short translations from sums and differences of valid vectors.
    implicit none

    real(real32), dimension(:,:,:), intent(in) :: frac_atoms
    integer, dimension(:), intent(in) :: spec_num
    real(real32), dimension(3,3), intent(in) :: lat
    real(real32), dimension(:,:), intent(inout) :: valid_vectors
    integer, intent(inout) :: n_valid
    real(real32), intent(in) :: frac_tol, cart_tol

    integer :: i, j, n_seed, n_expanded, expanded_size
    real(real32), allocatable, dimension(:,:) :: seeds, expanded
    real(real32), dimension(3) :: candidate


    n_seed = n_valid + 2
    expanded_size = max(1, 3 * n_seed * max(1, n_seed - 1) / 2)
    allocate(seeds(n_seed, 3))
    allocate(expanded(expanded_size, 3))
    seeds = 0._real32
    expanded = 0._real32
    n_expanded = 0

    if (n_valid.gt.0) seeds(1:n_valid, :) = valid_vectors(1:n_valid, :)
    seeds(n_valid + 1, :) = [ 1.0_real32, 0.0_real32, 0.0_real32 ]
    seeds(n_valid + 2, :) = [ 0.0_real32, 1.0_real32, 0.0_real32 ]

    do i = 1, n_seed - 1
       do j = i + 1, n_seed
          candidate(:) = seeds(i, :) + seeds(j, :)
          call canonicalise_translation(candidate, frac_tol, candidate)
          if (is_valid_translation(frac_atoms, spec_num, lat, candidate, &
               cart_tol)) then
             call append_unique_translation(expanded, n_expanded, candidate, frac_tol)
          end if

          candidate(:) = seeds(i, :) - seeds(j, :)
          call canonicalise_translation(candidate, frac_tol, candidate)
          if (is_valid_translation(frac_atoms, spec_num, lat, candidate, &
               cart_tol)) then
             call append_unique_translation(expanded, n_expanded, candidate, frac_tol)
          end if

          candidate(:) = seeds(j, :) - seeds(i, :)
          call canonicalise_translation(candidate, frac_tol, candidate)
          if (is_valid_translation(frac_atoms, spec_num, lat, candidate, &
               cart_tol)) then
             call append_unique_translation(expanded, n_expanded, candidate, frac_tol)
          end if
       end do
    end do

    do i = 1, n_expanded
       call append_unique_translation(valid_vectors, n_valid, expanded(i, :), &
            frac_tol)
    end do

    deallocate(seeds, expanded)

  end subroutine expand_valid_translations
!###############################################################################


!###############################################################################
  subroutine select_translation_pair_from_group(lat, valid_vectors, n_valid, &
       frac_tol, t1, t2)
    !! Select the primitive in-plane basis directly from a full relative group.
    implicit none

    real(real32), dimension(3,3), intent(in) :: lat
    real(real32), dimension(:,:), intent(in) :: valid_vectors
    integer, intent(in) :: n_valid
    real(real32), intent(in) :: frac_tol
    real(real32), dimension(3), intent(out) :: t1, t2

    integer :: i, j
    integer :: area_rank, best_area_rank
    integer :: support, best_support
    integer, dimension(2) :: key_left, key_right, best_key_left, best_key_right
    logical :: found
    real(real32) :: area, cell_area, area_resolution
    real(real32) :: len_left, len_right, sum_sq, best_sum_sq
    real(real32) :: frac_size, best_frac_size
    real(real32) :: orthogonality, best_orthogonality
    real(real32) :: big_len, best_big_len, best_area
    real(real32), dimension(3) :: left, right, cart_left, cart_right, tmp


    found = .false.
    cell_area = norm2(cross(lat(1, :), lat(2, :)))
    area_resolution = max(cell_area * frac_tol, 1.0e-8_real32)
    best_area_rank = huge(0)
    best_support = huge(0)
    best_frac_size = huge(0._real32)
    best_sum_sq = huge(0._real32)
    best_orthogonality = huge(0._real32)
    best_big_len = huge(0._real32)
    best_area = huge(0._real32)
    best_key_left = 0
    best_key_right = 0
    t1 = [ 1.0_real32, 0.0_real32, 0.0_real32 ]
    t2 = [ 0.0_real32, 1.0_real32, 0.0_real32 ]

    do i = 1, n_valid - 1
       do j = i + 1, n_valid
          left(:) = valid_vectors(i, :)
          right(:) = valid_vectors(j, :)
          if (translation_sort_less(lat, right, left)) then
             tmp(:) = left(:)
             left(:) = right(:)
             right(:) = tmp(:)
          end if

          cart_left(:) = matmul(left, lat)
          cart_right(:) = matmul(right, lat)
          len_left = norm2(cart_left)
          len_right = norm2(cart_right)
          area = norm2(cross(cart_left, cart_right))
          if (area.le.frac_tol) cycle

          area_rank = nint(area / area_resolution)
          support = pair_support(left, right, frac_tol)
          frac_size = pair_fractional_size(left, right)
          sum_sq = len_left**2 + len_right**2
          orthogonality = pair_cosine(cart_left, cart_right)
          big_len = len_right
          call translation_key(left, frac_tol, key_left)
          call translation_key(right, frac_tol, key_right)

          if (.not.found) then
             found = .true.
          elseif (.not. better_pair(&
               area_rank, support, frac_size, sum_sq, orthogonality, big_len, area, &
               key_left, key_right, &
               best_area_rank, best_support, best_frac_size, best_sum_sq, &
               best_orthogonality, best_big_len, best_area, best_key_left, &
               best_key_right &
          )) then
             cycle
          end if

          best_area_rank = area_rank
          best_support = support
          best_frac_size = frac_size
          best_sum_sq = sum_sq
          best_orthogonality = orthogonality
          best_big_len = big_len
          best_area = area
          best_key_left = key_left
          best_key_right = key_right
          t1(:) = left(:)
          t2(:) = right(:)
       end do
    end do

  end subroutine select_translation_pair_from_group
!###############################################################################


!###############################################################################
  subroutine select_translation_pair(lat, valid_vectors, n_valid, frac_tol, t1, t2)
    !! Select the primitive in-plane basis from the valid relative-shift group.
    implicit none

    real(real32), dimension(3,3), intent(in) :: lat
    real(real32), dimension(:,:), intent(in) :: valid_vectors
    integer, intent(in) :: n_valid
    real(real32), intent(in) :: frac_tol
    real(real32), dimension(3), intent(out) :: t1, t2

    integer :: i, j, n_span
    integer :: area_rank, best_area_rank
    integer, dimension(2) :: key_left, key_right, best_key_left, best_key_right
    logical :: found
    real(real32) :: area, cell_area, area_resolution
    integer :: support, best_support
    real(real32) :: len_left, len_right, sum_sq, best_sum_sq
    real(real32) :: frac_size, best_frac_size
    real(real32) :: orthogonality, best_orthogonality
    real(real32) :: big_len, best_big_len, best_area
    real(real32), allocatable, dimension(:,:) :: span_vectors
    real(real32), dimension(3) :: left, right, cart_left, cart_right


    if (n_valid.ge.2) then
       n_span = n_valid
    else
       n_span = 2
    end if
    allocate(span_vectors(n_span, 3))
    span_vectors = 0._real32
    if (n_valid.ge.2) then
       span_vectors(1:n_valid, :) = valid_vectors(1:n_valid, :)
    else
       span_vectors(1, :) = [ 1.0_real32, 0.0_real32, 0.0_real32 ]
       span_vectors(2, :) = [ 0.0_real32, 1.0_real32, 0.0_real32 ]
    end if

    cell_area = norm2(cross(lat(1, :), lat(2, :)))
    area_resolution = max(cell_area * frac_tol, 1.0e-8_real32)
    found = .false.
    best_area_rank = huge(0)
    best_support = huge(0)
    best_frac_size = huge(0._real32)
    best_sum_sq = huge(0._real32)
    best_orthogonality = huge(0._real32)
    best_big_len = huge(0._real32)
    best_area = huge(0._real32)
    best_key_left = 0
    best_key_right = 0
    t1 = [ 1.0_real32, 0.0_real32, 0.0_real32 ]
    t2 = [ 0.0_real32, 1.0_real32, 0.0_real32 ]

    do i = 1, n_span - 1
       do j = i + 1, n_span
          left(:) = span_vectors(i, :)
          right(:) = span_vectors(j, :)
          call reduce_translation_pair(lat, span_vectors, n_span, frac_tol, left, right)
          cart_left(:) = matmul(left, lat)
          cart_right(:) = matmul(right, lat)
          len_left = norm2(cart_left)
          len_right = norm2(cart_right)

          area = norm2(cross(cart_left, cart_right))
          if (area.le.frac_tol) cycle
          area_rank = nint(area / area_resolution)
          support = pair_support(left, right, frac_tol)
          frac_size = pair_fractional_size(left, right)
          sum_sq = len_left**2 + len_right**2
          orthogonality = pair_cosine(cart_left, cart_right)
          big_len = len_right
          call translation_key(left, frac_tol, key_left)
          call translation_key(right, frac_tol, key_right)

          if (.not.found) then
             found = .true.
          elseif (.not. better_pair(&
               area_rank, support, frac_size, sum_sq, orthogonality, big_len, area, &
               key_left, key_right, &
               best_area_rank, best_support, best_frac_size, best_sum_sq, &
               best_orthogonality, best_big_len, best_area, best_key_left, &
               best_key_right &
          )) then
             cycle
          end if

          best_area_rank = area_rank
          best_support = support
          best_frac_size = frac_size
          best_sum_sq = sum_sq
          best_orthogonality = orthogonality
          best_big_len = big_len
          best_area = area
          best_key_left = key_left
          best_key_right = key_right
          t1(:) = left(:)
          t2(:) = right(:)
       end do
    end do

    deallocate(span_vectors)

  end subroutine select_translation_pair
!###############################################################################


!###############################################################################
  subroutine reduce_translation_pair(lat, valid_vectors, n_valid, frac_tol, left, &
       right)
    !! Reduce a primitive pair through equivalent vectors that stay in the valid group.
    implicit none

    real(real32), dimension(3,3), intent(in) :: lat
    real(real32), dimension(:,:), intent(in) :: valid_vectors
    integer, intent(in) :: n_valid
    real(real32), intent(in) :: frac_tol
    real(real32), dimension(3), intent(inout) :: left, right

    integer :: iter, best_support, cand_support
    real(real32) :: best_frac_size, cand_frac_size
    real(real32) :: best_sum_sq, cand_sum_sq
    real(real32) :: best_orthogonality, cand_orthogonality
    real(real32) :: best_big_len, cand_big_len
    real(real32), dimension(3) :: best_left, best_right, cand_left, &
         cand_right, tmp
    real(real32), dimension(3) :: cart_left, cart_right, cand_cart_left, &
         cand_cart_right
    integer, dimension(2) :: best_key_left, best_key_right, cand_key_left, &
         cand_key_right
    logical :: improved


    do iter = 1, 12
       call canonicalise_translation(left, frac_tol, left)
       call canonicalise_translation(right, frac_tol, right)

       if (translation_sort_less(lat, right, left)) then
          tmp(:) = left(:)
          left(:) = right(:)
          right(:) = tmp(:)
       end if

       cart_left(:) = matmul(left, lat)
       cart_right(:) = matmul(right, lat)
       best_left(:) = left(:)
       best_right(:) = right(:)
       best_support = pair_support(left, right, frac_tol)
       best_frac_size = pair_fractional_size(left, right)
       best_sum_sq = norm2(cart_left)**2 + norm2(cart_right)**2
       best_orthogonality = pair_cosine(cart_left, cart_right)
       best_big_len = norm2(cart_right)
       call translation_key(left, frac_tol, best_key_left)
       call translation_key(right, frac_tol, best_key_right)
       improved = .false.

       call try_reduced_pair(lat, valid_vectors, n_valid, frac_tol, &
            left, right - left, &
            best_left, best_right, best_support, best_frac_size, best_sum_sq, &
            best_orthogonality, best_big_len, best_key_left, best_key_right, improved)
       call try_reduced_pair(lat, valid_vectors, n_valid, frac_tol, &
            left, right + left, &
            best_left, best_right, best_support, best_frac_size, best_sum_sq, &
            best_orthogonality, best_big_len, best_key_left, best_key_right, improved)
       call try_reduced_pair(lat, valid_vectors, n_valid, frac_tol, &
            left - right, right, &
            best_left, best_right, best_support, best_frac_size, best_sum_sq, &
            best_orthogonality, best_big_len, best_key_left, best_key_right, improved)
       call try_reduced_pair(lat, valid_vectors, n_valid, frac_tol, &
            left + right, right, &
            best_left, best_right, best_support, best_frac_size, best_sum_sq, &
            best_orthogonality, best_big_len, best_key_left, best_key_right, improved)

       if (.not.improved) exit
       left(:) = best_left(:)
       right(:) = best_right(:)
    end do

    call canonicalise_translation(left, frac_tol, left)
    call canonicalise_translation(right, frac_tol, right)
    if (translation_sort_less(lat, right, left)) then
       tmp(:) = left(:)
       left(:) = right(:)
       right(:) = tmp(:)
    end if

  end subroutine reduce_translation_pair
!###############################################################################


!###############################################################################
  subroutine try_reduced_pair(lat, valid_vectors, n_valid, frac_tol, raw_left, &
       raw_right, &
       best_left, best_right, best_support, best_frac_size, best_sum_sq, &
       best_orthogonality, best_big_len, best_key_left, best_key_right, improved)
    !! Try one equivalent pair transformation that remains inside the valid group.
    implicit none

    real(real32), dimension(3,3), intent(in) :: lat
    real(real32), dimension(:,:), intent(in) :: valid_vectors
    integer, intent(in) :: n_valid
    real(real32), intent(in) :: frac_tol
    real(real32), dimension(3), intent(in) :: raw_left, raw_right
    real(real32), dimension(3), intent(inout) :: best_left, best_right
    integer, intent(inout) :: best_support
    real(real32), intent(inout) :: best_frac_size, best_sum_sq, best_orthogonality, &
         best_big_len
    integer, dimension(2), intent(inout) :: best_key_left, best_key_right
    logical, intent(inout) :: improved

    integer :: cand_support
    integer, dimension(2) :: cand_key_left, cand_key_right
    real(real32) :: cand_frac_size, cand_sum_sq, cand_orthogonality, cand_big_len
    real(real32), dimension(3) :: cand_left, cand_right, tmp, cand_cart_left, &
         cand_cart_right


    cand_left(:) = raw_left(:)
    cand_right(:) = raw_right(:)
    call canonicalise_translation(cand_left, frac_tol, cand_left)
    call canonicalise_translation(cand_right, frac_tol, cand_right)
    if (norm2(cand_left(1:2)).le.frac_tol) return
    if (norm2(cand_right(1:2)).le.frac_tol) return
    if (.not.vector_in_set(valid_vectors, n_valid, cand_left, frac_tol)) return
    if (.not.vector_in_set(valid_vectors, n_valid, cand_right, frac_tol)) return

    if (translation_sort_less(lat, cand_right, cand_left)) then
       tmp(:) = cand_left(:)
       cand_left(:) = cand_right(:)
       cand_right(:) = tmp(:)
    end if

    cand_cart_left(:) = matmul(cand_left, lat)
    cand_cart_right(:) = matmul(cand_right, lat)
    if (norm2(cross(cand_cart_left, cand_cart_right)).le.frac_tol) return

    cand_support = pair_support(cand_left, cand_right, frac_tol)
    cand_frac_size = pair_fractional_size(cand_left, cand_right)
    cand_sum_sq = norm2(cand_cart_left)**2 + norm2(cand_cart_right)**2
    cand_orthogonality = pair_cosine(cand_cart_left, cand_cart_right)
    cand_big_len = norm2(cand_cart_right)
    call translation_key(cand_left, frac_tol, cand_key_left)
    call translation_key(cand_right, frac_tol, cand_key_right)

    if (.not. better_reduced_pair(&
         cand_support, cand_frac_size, cand_sum_sq, cand_orthogonality, cand_big_len, &
         cand_key_left, cand_key_right, &
         best_support, best_frac_size, best_sum_sq, best_orthogonality, best_big_len, &
         best_key_left, best_key_right)) return

    best_left(:) = cand_left(:)
    best_right(:) = cand_right(:)
    best_support = cand_support
    best_frac_size = cand_frac_size
    best_sum_sq = cand_sum_sq
    best_orthogonality = cand_orthogonality
    best_big_len = cand_big_len
    best_key_left = cand_key_left
    best_key_right = cand_key_right
    improved = .true.

  end subroutine try_reduced_pair
!###############################################################################


!###############################################################################
  logical function better_pair(&
       area_rank, support, frac_size, sum_sq, orthogonality, big_len, area, &
       key_left, key_right, &
       best_area_rank, best_support, best_frac_size, best_sum_sq, best_orthogonality, &
       best_big_len, best_area, best_key_left, best_key_right &
  )
    !! Lexicographic comparison of primitive-basis candidates.
    implicit none

    integer, intent(in) :: area_rank, support, best_area_rank, best_support
    integer, dimension(2), intent(in) :: key_left, key_right, best_key_left, &
         best_key_right
    real(real32), intent(in) :: frac_size, sum_sq, orthogonality, big_len, area
    real(real32), intent(in) :: best_frac_size, best_sum_sq, best_orthogonality, &
         best_big_len, best_area
    real(real32), parameter :: eps = 1.0e-6_real32


    better_pair = .false.
    if (area_rank.lt.best_area_rank) then
       better_pair = .true.
       return
    elseif (area_rank.gt.best_area_rank) then
       return
    end if

    if (support.lt.best_support) then
       better_pair = .true.
       return
    elseif (support.gt.best_support) then
       return
    end if

    if (frac_size.lt.best_frac_size - eps) then
       better_pair = .true.
       return
    elseif (frac_size.gt.best_frac_size + eps) then
       return
    end if

    if (sum_sq.lt.best_sum_sq - eps) then
       better_pair = .true.
       return
    elseif (sum_sq.gt.best_sum_sq + eps) then
       return
    end if

    if (orthogonality.lt.best_orthogonality - eps) then
       better_pair = .true.
       return
    elseif (orthogonality.gt.best_orthogonality + eps) then
       return
    end if

    if (big_len.lt.best_big_len - eps) then
       better_pair = .true.
       return
    elseif (big_len.gt.best_big_len + eps) then
       return
    end if

    if (area.lt.best_area - eps) then
       better_pair = .true.
       return
    elseif (area.gt.best_area + eps) then
       return
    end if

    if (translation_key_less(key_left, best_key_left)) then
       better_pair = .true.
       return
    elseif (translation_key_less(best_key_left, key_left)) then
       return
    end if

    if (translation_key_less(key_right, best_key_right)) better_pair = .true.

  end function better_pair
!###############################################################################


!###############################################################################
  logical function better_reduced_pair(&
       support, frac_size, sum_sq, orthogonality, big_len, key_left, key_right, &
       best_support, best_frac_size, best_sum_sq, best_orthogonality, best_big_len, &
       best_key_left, best_key_right)
    !! Lexicographic comparison for equivalent pair reductions inside the valid group.
    implicit none

    integer, intent(in) :: support, best_support
    integer, dimension(2), intent(in) :: key_left, key_right, best_key_left, &
         best_key_right
    real(real32), intent(in) :: frac_size, sum_sq, orthogonality, big_len
    real(real32), intent(in) :: best_frac_size, best_sum_sq, best_orthogonality, &
         best_big_len
    real(real32), parameter :: eps = 1.0e-6_real32


    better_reduced_pair = .false.
    if (support.lt.best_support) then
       better_reduced_pair = .true.
       return
    elseif (support.gt.best_support) then
       return
    end if

    if (frac_size.lt.best_frac_size - eps) then
       better_reduced_pair = .true.
       return
    elseif (frac_size.gt.best_frac_size + eps) then
       return
    end if

    if (sum_sq.lt.best_sum_sq - eps) then
       better_reduced_pair = .true.
       return
    elseif (sum_sq.gt.best_sum_sq + eps) then
       return
    end if

    if (orthogonality.lt.best_orthogonality - eps) then
       better_reduced_pair = .true.
       return
    elseif (orthogonality.gt.best_orthogonality + eps) then
       return
    end if

    if (big_len.lt.best_big_len - eps) then
       better_reduced_pair = .true.
       return
    elseif (big_len.gt.best_big_len + eps) then
       return
    end if

    if (translation_key_less(key_left, best_key_left)) then
       better_reduced_pair = .true.
       return
    elseif (translation_key_less(best_key_left, key_left)) then
       return
    end if

    if (translation_key_less(key_right, best_key_right)) better_reduced_pair = .true.

  end function better_reduced_pair
!###############################################################################


!###############################################################################
  real(real32) function pair_cosine(lhs, rhs) result(value)
  !! Return the absolute cosine between two Cartesian vectors.
  implicit none

  real(real32), dimension(3), intent(in) :: lhs, rhs
  real(real32) :: lhs_norm, rhs_norm


  lhs_norm = norm2(lhs)
  rhs_norm = norm2(rhs)
  if (lhs_norm.le.1.0e-16_real32 .or. rhs_norm.le.1.0e-16_real32) then
     value = 1.0_real32
     return
  end if

  value = abs(dot_product(lhs, rhs)) / (lhs_norm * rhs_norm)

  end function pair_cosine
!###############################################################################


!###############################################################################
  integer function pair_support(lhs, rhs, frac_tol) result(value)
    !! Count the number of non-zero in-plane components across a pair.
    implicit none

    real(real32), dimension(3), intent(in) :: lhs, rhs
    real(real32), intent(in) :: frac_tol


    value = 0
    if (abs(lhs(1)).gt.frac_tol) value = value + 1
    if (abs(lhs(2)).gt.frac_tol) value = value + 1
    if (abs(rhs(1)).gt.frac_tol) value = value + 1
    if (abs(rhs(2)).gt.frac_tol) value = value + 1

  end function pair_support
!###############################################################################


!###############################################################################
  real(real32) function pair_fractional_size(lhs, rhs) result(value)
  !! Return the total absolute in-plane fractional size of a pair.
  implicit none

  real(real32), dimension(3), intent(in) :: lhs, rhs


  value = abs(lhs(1)) + abs(lhs(2)) + abs(rhs(1)) + abs(rhs(2))

  end function pair_fractional_size
!###############################################################################


!###############################################################################
  subroutine translation_key(translation, frac_tol, key)
    !! Convert a translation into a deterministic integer key.
    implicit none

    real(real32), dimension(3), intent(in) :: translation
    real(real32), intent(in) :: frac_tol
    integer, dimension(2), intent(out) :: key


    key(1) = nint(translation(1) / frac_tol)
    key(2) = nint(translation(2) / frac_tol)

  end subroutine translation_key
!###############################################################################


!###############################################################################
  logical function translation_key_less(lhs, rhs)
    !! Lexicographic ordering helper for translation keys.
    implicit none

    integer, dimension(2), intent(in) :: lhs, rhs


    translation_key_less = .false.
    if (lhs(1).lt.rhs(1)) then
       translation_key_less = .true.
    elseif (lhs(1).eq.rhs(1) .and. lhs(2).lt.rhs(2)) then
       translation_key_less = .true.
    end if

  end function translation_key_less
!###############################################################################


!###############################################################################
  logical function vector_in_set(storage, nstored, candidate, frac_tol)
    !! Return whether a translation already exists in a storage array.
    implicit none

    real(real32), dimension(:,:), intent(in) :: storage
    integer, intent(in) :: nstored
    real(real32), dimension(3), intent(in) :: candidate
    real(real32), intent(in) :: frac_tol

    integer :: i


    vector_in_set = .false.
    do i = 1, nstored
       if (all(abs(storage(i, :) - candidate(:)).le.frac_tol)) then
          vector_in_set = .true.
          return
       end if
    end do

  end function vector_in_set
!###############################################################################


!###############################################################################
  logical function translation_sort_less(lat, lhs, rhs)
    !! Compare translations using length first, then absolute fractional components.
    implicit none

    real(real32), dimension(3,3), intent(in) :: lat
    real(real32), dimension(3), intent(in) :: lhs, rhs
    real(real32), parameter :: eps = 1.0e-6_real32
    real(real32) :: lhs_len, rhs_len


    translation_sort_less = .false.
    lhs_len = norm2(matmul(lhs, lat))
    rhs_len = norm2(matmul(rhs, lat))

    if (lhs_len.lt.rhs_len - eps) then
       translation_sort_less = .true.
       return
    elseif (lhs_len.gt.rhs_len + eps) then
       return
    end if

    if (abs(lhs(1)).lt.abs(rhs(1)) - eps) then
       translation_sort_less = .true.
    elseif (abs(lhs(1)).gt.abs(rhs(1)) + eps) then
       return
    elseif (abs(lhs(2)).lt.abs(rhs(2)) - eps) then
       translation_sort_less = .true.
    elseif (abs(lhs(2)).gt.abs(rhs(2)) + eps) then
       return
    elseif (abs(lhs(3)).lt.abs(rhs(3)) - eps) then
       translation_sort_less = .true.
    end if

  end function translation_sort_less
!###############################################################################


!###############################################################################
  subroutine wrap_zero_one(vector)
    !! Wrap a fractional vector into the [0, 1) interval.
    implicit none

    real(real32), dimension(3), intent(inout) :: vector
    integer :: i


    do i = 1, 3
       vector(i) = vector(i) - floor(vector(i))
       if (vector(i).ge.1.0_real32) vector(i) = vector(i) - 1.0_real32
    end do

  end subroutine wrap_zero_one
!###############################################################################


!###############################################################################
  subroutine wrap_half(vector)
    !! Wrap a fractional vector into the [-1/2, 1/2) interval.
    implicit none

    real(real32), dimension(3), intent(inout) :: vector
    integer :: i


    do i = 1, 3
       vector(i) = vector(i) - floor(vector(i) + 0.5_real32)
    end do

  end subroutine wrap_half
!###############################################################################


!###############################################################################
  subroutine unpack_translation(work_translation, work_axes, translation)
    !! Map a working translation basis back onto the original axis ordering.
    implicit none

    real(real32), dimension(3), intent(in) :: work_translation
    integer, dimension(3), intent(in) :: work_axes
    real(real32), dimension(3), intent(out) :: translation


    translation = 0._real32
    translation(work_axes(1)) = work_translation(1)
    translation(work_axes(2)) = work_translation(2)
    translation(work_axes(3)) = 0._real32

  end subroutine unpack_translation
!###############################################################################

end module artemis__interface_translations
