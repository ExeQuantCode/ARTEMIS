module artemis__crystal_compare
  use coreutils__kind, only: real32
  use coreutils__linalg, only: inverse_3x3
  use atomstruc, only: basis_type

  implicit none
  private

  public :: crystals_equivalent

  real(real32), parameter :: sort_eps = 64.0_real32 * epsilon(1.0_real32)

  type :: sorted_site_set_type
     integer :: nsite = 0
     integer :: nspecies = 0
     integer, allocatable :: species(:)
     real(real32), allocatable :: coords(:, :)
     integer, allocatable :: species_value(:)
     integer, allocatable :: species_start(:)
     integer, allocatable :: species_count(:)
  end type sorted_site_set_type

contains

  logical function crystals_equivalent( &
       basis_a_in, basis_b_in, tol, exact, allow_translation, plane_normal &
  ) result(is_equivalent)
    class(basis_type), intent(in) :: basis_a_in, basis_b_in
    real(real32), intent(in) :: tol
    logical, intent(in), optional :: exact, allow_translation
    real(real32), intent(in), optional :: plane_normal(3)

    type(basis_type) :: basis_a, basis_b
    real(real32), allocatable :: coords_a(:, :), coords_b(:, :)
    integer, allocatable :: species_a(:), species_b(:), pivot_a(:), pivot_b(:)
    logical :: compare_exact, compare_translation, has_plane
    real(real32) :: plane_unit(3)

    is_equivalent = .false.
    if (tol < 0.0_real32) return

    compare_exact = .false.
    if (present(exact)) compare_exact = exact

    compare_translation = .false.
    if (present(allow_translation)) compare_translation = allow_translation

    has_plane = present(plane_normal)
    if (has_plane) then
       compare_translation = .true.
       if (.not. normalise_vector(plane_normal, plane_unit)) return
    else
       plane_unit = 0.0_real32
    end if

    if (basis_a_in%natom /= basis_b_in%natom) return
    if (basis_a_in%nspec /= basis_b_in%nspec) return
    if (any(basis_a_in%pbc .neqv. basis_b_in%pbc)) return

    call prepare_fractional_basis(basis_a, basis_a_in)
    call prepare_fractional_basis(basis_b, basis_b_in)
    call flatten_basis(basis_a, coords_a, species_a)
    call flatten_basis(basis_b, coords_b, species_b)
    call analyse_species(species_a, species_b, compare_translation, pivot_a, pivot_b)
    if (.not. allocated(pivot_a)) return
    if (.not. allocated(pivot_b)) return

    is_equivalent = compare_in_basis_a_frame(&
         basis_a%lat, &
         coords_a, &
         species_a, &
         basis_b%lat, &
         coords_b, &
         species_b, &
         tol, &
         compare_exact, &
         compare_translation, &
         has_plane, &
         plane_unit, &
         pivot_a, &
         pivot_b &
    )
  end function crystals_equivalent

  subroutine prepare_fractional_basis(output_basis, input_basis)
    type(basis_type), intent(out) :: output_basis
    class(basis_type), intent(in) :: input_basis

    call output_basis%copy(input_basis)
    if (output_basis%lcart) call output_basis%convert()
    call output_basis%normalise()
  end subroutine prepare_fractional_basis

  subroutine flatten_basis(basis, coords, species)
    class(basis_type), intent(in) :: basis
    real(real32), allocatable, intent(out) :: coords(:, :)
    integer, allocatable, intent(out) :: species(:)

    integer :: atom_index, species_index, local_index
    integer :: species_code

    allocate(coords(basis%natom, 3))
    allocate(species(basis%natom))
    atom_index = 0

    do species_index = 1, basis%nspec
       species_code = encode_species_name(basis%spec(species_index)%name)
       do local_index = 1, basis%spec(species_index)%num
          atom_index = atom_index + 1
          coords(atom_index, :) = &
               real(basis%spec(species_index)%atom(local_index, 1:3), real32)
          species(atom_index) = species_code
       end do
    end do
  end subroutine flatten_basis

  subroutine analyse_species(species_a, species_b, allow_translation, pivot_a, pivot_b)
    integer, intent(in) :: species_a(:), species_b(:)
    logical, intent(in) :: allow_translation
    integer, allocatable, intent(out) :: pivot_a(:), pivot_b(:)

    integer, allocatable :: unique_species(:), count_a(:), count_b(:)
    integer :: nsite, unique_count, i, j, slot, best_species, best_count
    logical :: found

    nsite = size(species_a)
    if (nsite /= size(species_b)) return
    if (nsite == 0) then
       allocate(pivot_a(0), pivot_b(0))
       return
    end if

    allocate(unique_species(nsite), count_a(nsite), count_b(nsite))
    unique_species = 0
    count_a = 0
    count_b = 0
    unique_count = 0

    do i = 1, nsite
       slot = 0
       do j = 1, unique_count
          if (unique_species(j) == species_a(i)) then
             slot = j
             exit
          end if
       end do
       if (slot == 0) then
          unique_count = unique_count + 1
          unique_species(unique_count) = species_a(i)
          slot = unique_count
       end if
       count_a(slot) = count_a(slot) + 1
    end do

    do i = 1, nsite
       found = .false.
       do j = 1, unique_count
          if (unique_species(j) == species_b(i)) then
             count_b(j) = count_b(j) + 1
             found = .true.
             exit
          end if
       end do
       if (.not. found) return
    end do

    do i = 1, unique_count
       if (count_a(i) /= count_b(i)) return
    end do

    best_species = unique_species(1)
    best_count = count_a(1)
    if (allow_translation) then
       do i = 2, unique_count
          if (count_a(i) < best_count) then
             best_species = unique_species(i)
             best_count = count_a(i)
          end if
       end do
    end if

    allocate(pivot_a(best_count), pivot_b(best_count))
    slot = 0
    do i = 1, nsite
       if (species_a(i) == best_species) then
          slot = slot + 1
          pivot_a(slot) = i
       end if
    end do

    slot = 0
    do i = 1, nsite
       if (species_b(i) == best_species) then
          slot = slot + 1
          pivot_b(slot) = i
       end if
    end do
  end subroutine analyse_species

  logical function compare_in_basis_a_frame( &
       lat_a, coords_a, species_a, lat_b, coords_b, species_b, &
       tol, exact, allow_translation, has_plane, plane_unit, pivot_a, pivot_b &
  ) result(is_equivalent)
    real(real32), intent(in) :: lat_a(3, 3), coords_a(:, :), lat_b(3, 3), coords_b(:, :)
    integer, intent(in) :: species_a(:), species_b(:), pivot_a(:), pivot_b(:)
    real(real32), intent(in) :: tol, plane_unit(3)
    logical, intent(in) :: exact, allow_translation, has_plane

    real(real32) :: inv_a(3, 3), map_real(3, 3), map_int_real(3, 3), recon(3, 3)
    real(real32) :: gram(3, 3), frac_bound(3), lattice_tol, roundoff_tol, scale, tol2
    real(real32), allocatable :: &
         coords_b_in_a(:, :), ref_coords(:, :), candidate_coords(:, :)
    integer, allocatable :: candidate_order(:), sort_work(:)
    logical, allocatable :: used(:)
    type(sorted_site_set_type) :: reference_sites
    integer :: map_int(3, 3), nsite, pivot_index, site_index
    real(real32) :: shift(3), shift_cart(3)

    is_equivalent = .false.
    nsite = size(species_a)
    if (nsite /= size(species_b)) return
    if (nsite /= size(coords_a, 1)) return
    if (nsite /= size(coords_b, 1)) return

    scale = max(1.0_real32, maxval(abs(lat_a)), maxval(abs(lat_b)))
    roundoff_tol = 256.0_real32 * epsilon(scale) * scale
    lattice_tol = roundoff_tol
    if (.not. exact) lattice_tol = max(tol, roundoff_tol)

    if (.not. lattice_is_valid(lat_a, lattice_tol)) return
    if (.not. lattice_is_valid(lat_b, lattice_tol)) return

    inv_a = inverse_3x3(lat_a)
    map_real = matmul(lat_b, inv_a)
    map_int = nint(map_real)
    map_int_real = real(map_int, real32)
    if (abs(det3_int(map_int)) /= 1) return

    recon = matmul(map_int_real, lat_a)
    if (any(abs(recon - lat_b) > lattice_tol)) return

    gram = matmul(lat_a, transpose(lat_a))
    call get_fractional_bounds(inv_a, lattice_tol, roundoff_tol, frac_bound)
    tol2 = max(tol, roundoff_tol)**2

    allocate(coords_b_in_a(nsite, 3))
    coords_b_in_a = matmul(coords_b, map_int_real)
    do site_index = 1, nsite
       call wrap_zero_one_vector(coords_b_in_a(site_index, :))
    end do

    allocate( &
         candidate_coords(nsite, 3), &
         candidate_order(nsite), &
         sort_work(nsite), &
         used(nsite) &
    )

    if (.not. allow_translation) then
       call build_sorted_site_set(coords_a, species_a, reference_sites)
       candidate_coords = coords_b_in_a
       is_equivalent = candidate_matches_reference(&
            reference_sites, &
            candidate_coords, &
            species_b, &
            gram, &
            frac_bound, &
            tol2, &
            candidate_order, &
            sort_work, &
            used &
       )
       return
    end if

    allocate(ref_coords(nsite, 3))
    call build_relative_coordinates(coords_a, pivot_a(1), ref_coords)
    call build_sorted_site_set(ref_coords, species_a, reference_sites)

    do pivot_index = 1, size(pivot_b)
       shift = coords_a(pivot_a(1), :) - coords_b_in_a(pivot_b(pivot_index), :)
       call wrap_half_vector(shift)
       if (has_plane) then
          shift_cart = matmul(shift, lat_a)
          if (abs(dot_product(shift_cart, plane_unit)) > lattice_tol) cycle
       end if

       call build_relative_coordinates( &
            coords_b_in_a, &
            pivot_b(pivot_index), &
            candidate_coords &
       )
       if (candidate_matches_reference(&
            reference_sites, &
            candidate_coords, &
            species_b, &
            gram, &
            frac_bound, &
            tol2, &
            candidate_order, &
            sort_work, &
            used &
       )) then
          is_equivalent = .true.
          return
       end if
    end do
  end function compare_in_basis_a_frame

  subroutine build_relative_coordinates(coords, pivot_index, relative_coords)
    real(real32), intent(in) :: coords(:, :)
    integer, intent(in) :: pivot_index
    real(real32), intent(out) :: relative_coords(:, :)

    integer :: site_index

    do site_index = 1, size(coords, 1)
       relative_coords(site_index, :) = coords(site_index, :) - coords(pivot_index, :)
       call wrap_half_vector(relative_coords(site_index, :))
       call wrap_zero_one_vector(relative_coords(site_index, :))
    end do
    relative_coords(pivot_index, :) = 0.0_real32
  end subroutine build_relative_coordinates

  subroutine build_sorted_site_set(coords, species, output_sites)
    real(real32), intent(in) :: coords(:, :)
    integer, intent(in) :: species(:)
    type(sorted_site_set_type), intent(out) :: output_sites

    integer, allocatable :: order(:), work(:)
    integer :: site_index, species_index, current_species

    output_sites%nsite = size(species)
    allocate(output_sites%coords(output_sites%nsite, 3))
    allocate(output_sites%species(output_sites%nsite))
    allocate(order(output_sites%nsite), work(max(1, output_sites%nsite)))

    do site_index = 1, output_sites%nsite
       order(site_index) = site_index
    end do
    call sort_site_order(order, species, coords, work)

    do site_index = 1, output_sites%nsite
       output_sites%species(site_index) = species(order(site_index))
       output_sites%coords(site_index, :) = coords(order(site_index), :)
    end do

    output_sites%nspecies = 0
    current_species = huge(0)
    do site_index = 1, output_sites%nsite
       if (output_sites%species(site_index) /= current_species) then
          output_sites%nspecies = output_sites%nspecies + 1
          current_species = output_sites%species(site_index)
       end if
    end do

    allocate(output_sites%species_value(output_sites%nspecies))
    allocate(output_sites%species_start(output_sites%nspecies))
    allocate(output_sites%species_count(output_sites%nspecies))

    species_index = 0
    do site_index = 1, output_sites%nsite
       if( &
            site_index == 1 .or. &
            output_sites%species(site_index) /= output_sites%species(site_index - 1) &
       ) then
          species_index = species_index + 1
          output_sites%species_value(species_index) = output_sites%species(site_index)
          output_sites%species_start(species_index) = site_index
          output_sites%species_count(species_index) = 1
       else
          output_sites%species_count(species_index) = &
               output_sites%species_count(species_index) + 1
       end if
    end do
  end subroutine build_sorted_site_set

  logical function candidate_matches_reference( &
       reference_sites, candidate_coords, candidate_species, gram, &
       frac_bound, tol2, candidate_order, sort_work, used &
  ) result(is_match)
    type(sorted_site_set_type), intent(in) :: reference_sites
    real(real32), intent(in) :: candidate_coords(:, :), gram(3, 3), frac_bound(3), &
         tol2
    integer, intent(in) :: candidate_species(:)
    integer, intent(inout) :: candidate_order(:), sort_work(:)
    logical, intent(inout) :: used(:)

    integer :: site_index, start_index, stop_index, best_index
    real(real32), allocatable :: sorted_candidate_coords(:, :)
    integer, allocatable :: sorted_candidate_species(:)

    is_match = .false.
    if (reference_sites%nsite /= size(candidate_species)) return

    allocate(sorted_candidate_coords(reference_sites%nsite, 3))
    allocate(sorted_candidate_species(reference_sites%nsite))
    used = .false.

    do site_index = 1, reference_sites%nsite
       candidate_order(site_index) = site_index
    end do
    call sort_site_order( &
         candidate_order, candidate_species, candidate_coords, sort_work &
    )

    do site_index = 1, reference_sites%nsite
       sorted_candidate_species(site_index) = &
            candidate_species(candidate_order(site_index))
       sorted_candidate_coords(site_index, :) = &
            candidate_coords(candidate_order(site_index), :)
    end do

    do site_index = 1, reference_sites%nsite
       call get_species_range( &
            reference_sites, sorted_candidate_species(site_index), &
            start_index, stop_index &
       )
       if (start_index == 0) return
       best_index = find_best_match(&
            reference_sites%coords, &
            sorted_candidate_coords(site_index, :), &
            used, &
            start_index, &
            stop_index, &
            frac_bound, &
            gram, &
            tol2 &
       )
       if (best_index == 0) return
       used(best_index) = .true.
    end do

    is_match = .true.
  end function candidate_matches_reference

  subroutine get_species_range(reference_sites, species_value, start_index, stop_index)
    type(sorted_site_set_type), intent(in) :: reference_sites
    integer, intent(in) :: species_value
    integer, intent(out) :: start_index, stop_index

    integer :: species_index

    start_index = 0
    stop_index = -1
    do species_index = 1, reference_sites%nspecies
       if (reference_sites%species_value(species_index) == species_value) then
          start_index = reference_sites%species_start(species_index)
          stop_index = start_index + reference_sites%species_count(species_index) - 1
          return
       end if
    end do
  end subroutine get_species_range

  integer function find_best_match( &
       reference_coords, candidate_coord, used, &
       start_index, stop_index, frac_bound, gram, tol2 &
  ) result(best_index)
    real(real32), intent(in) :: reference_coords(:, :), &
         candidate_coord(3), frac_bound(3), gram(3, 3), tol2
    logical, intent(in) :: used(:)
    integer, intent(in) :: start_index, stop_index

    integer :: site_index
    real(real32) :: diff(3), dist2, best_d2

    best_index = 0
    best_d2 = huge(1.0_real32)

    do site_index = start_index, stop_index
       if (used(site_index)) cycle
       diff = candidate_coord - reference_coords(site_index, :)
       call wrap_half_vector(diff)
       if (abs(diff(1)) > frac_bound(1) + sort_eps) cycle
       if (abs(diff(2)) > frac_bound(2) + sort_eps) cycle
       if (abs(diff(3)) > frac_bound(3) + sort_eps) cycle
       dist2 = metric_distance2(diff, gram)
       if (dist2 <= tol2 .and. dist2 + tiny(1.0_real32) < best_d2) then
          best_d2 = dist2
          best_index = site_index
       end if
    end do
  end function find_best_match

  subroutine sort_site_order(order, species, coords, work)
    integer, intent(inout) :: order(:)
    integer, intent(in) :: species(:)
    real(real32), intent(in) :: coords(:, :)
    integer, intent(inout) :: work(:)

    integer :: width, left, middle, right, left_cursor, right_cursor, target
    integer :: nsite

    nsite = size(order)
    if (nsite <= 1) return

    width = 1
    do while (width < nsite)
       left = 1
       do while (left <= nsite)
          middle = min(left + width - 1, nsite)
          right = min(left + 2 * width - 1, nsite)
          if (middle < right) then
             left_cursor = left
             right_cursor = middle + 1
             target = left
             do while (left_cursor <= middle .and. right_cursor <= right)
                if( site_less( &
                     order(left_cursor), order(right_cursor), species, coords &
                ) )then
                   work(target) = order(left_cursor)
                   left_cursor = left_cursor + 1
                else
                   work(target) = order(right_cursor)
                   right_cursor = right_cursor + 1
                end if
                target = target + 1
             end do
             do while (left_cursor <= middle)
                work(target) = order(left_cursor)
                left_cursor = left_cursor + 1
                target = target + 1
             end do
             do while (right_cursor <= right)
                work(target) = order(right_cursor)
                right_cursor = right_cursor + 1
                target = target + 1
             end do
             order(left:right) = work(left:right)
          end if
          left = left + 2 * width
       end do
       width = width * 2
    end do
  end subroutine sort_site_order

  logical function site_less(lhs, rhs, species, coords) result(is_less)
    integer, intent(in) :: lhs, rhs, species(:)
    real(real32), intent(in) :: coords(:, :)

    integer :: axis_index

    if (species(lhs) /= species(rhs)) then
       is_less = species(lhs) < species(rhs)
       return
    end if

    do axis_index = 1, 3
       if (coords(lhs, axis_index) < coords(rhs, axis_index) - sort_eps) then
          is_less = .true.
          return
       end if
       if (coords(lhs, axis_index) > coords(rhs, axis_index) + sort_eps) then
          is_less = .false.
          return
       end if
    end do

    is_less = lhs < rhs
  end function site_less

  subroutine get_fractional_bounds(inv_lat, cart_tol, roundoff_tol, frac_bound)
    real(real32), intent(in) :: inv_lat(3, 3), cart_tol, roundoff_tol
    real(real32), intent(out) :: frac_bound(3)

    integer :: axis_index

    do axis_index = 1, 3
       frac_bound(axis_index) = cart_tol * norm2(inv_lat(:, axis_index)) + roundoff_tol
    end do
  end subroutine get_fractional_bounds

  logical function lattice_is_valid(lat, lattice_tol) result(is_valid)
    real(real32), intent(in) :: lat(3, 3), lattice_tol
    real(real32) :: scale

    scale = max(1.0_real32, maxval(abs(lat)))
    is_valid = &
         abs(det3_real(lat)) > &
         max(lattice_tol, 256.0_real32 * epsilon(scale) * scale**3)
  end function lattice_is_valid

  logical function normalise_vector(vector_in, vector_out) result(is_valid)
    real(real32), intent(in) :: vector_in(3)
    real(real32), intent(out) :: vector_out(3)
    real(real32) :: length

    length = norm2(vector_in)
    if (length <= epsilon(length)) then
       vector_out = 0.0_real32
       is_valid = .false.
       return
    end if

    vector_out = vector_in / length
    is_valid = .true.
  end function normalise_vector

  integer function encode_species_name(name) result(code)
    character(len=*), intent(in) :: name
    integer :: index

    code = 0
    do index = 1, len(name)
       code = code + iachar(name(index:index)) * 256**(index - 1)
    end do
  end function encode_species_name

  real(real32) function det3_real(mat) result(det)
  real(real32), intent(in) :: mat(3, 3)

  det = mat(1, 1) * (mat(2, 2) * mat(3, 3) - mat(2, 3) * mat(3, 2)) - &
       mat(1, 2) * (mat(2, 1) * mat(3, 3) - mat(2, 3) * mat(3, 1)) + &
       mat(1, 3) * (mat(2, 1) * mat(3, 2) - mat(2, 2) * mat(3, 1))
  end function det3_real

  integer function det3_int(mat) result(det)
    integer, intent(in) :: mat(3, 3)

    det = mat(1, 1) * (mat(2, 2) * mat(3, 3) - mat(2, 3) * mat(3, 2)) - &
         mat(1, 2) * (mat(2, 1) * mat(3, 3) - mat(2, 3) * mat(3, 1)) + &
         mat(1, 3) * (mat(2, 1) * mat(3, 2) - mat(2, 2) * mat(3, 1))
  end function det3_int

  subroutine wrap_zero_one_vector(vector)
    real(real32), intent(inout) :: vector(3)
    integer :: axis_index

    do axis_index = 1, 3
       vector(axis_index) = vector(axis_index) - floor(vector(axis_index))
       if (vector(axis_index) < 0.0_real32) &
            vector(axis_index) = vector(axis_index) + 1.0_real32
       if (vector(axis_index) >= 1.0_real32 - sort_eps) &
            vector(axis_index) = 0.0_real32
    end do
  end subroutine wrap_zero_one_vector

  subroutine wrap_half_vector(vector)
    real(real32), intent(inout) :: vector(3)
    integer :: axis_index

    do axis_index = 1, 3
       vector(axis_index) = wrap_half_scalar(vector(axis_index))
    end do
  end subroutine wrap_half_vector

  real(real32) function wrap_half_scalar(value) result(wrapped)
  real(real32), intent(in) :: value

  wrapped = value - anint(value)
  if (wrapped <= -0.5_real32 + sort_eps) wrapped = wrapped + 1.0_real32
  if (wrapped > 0.5_real32 + sort_eps) wrapped = wrapped - 1.0_real32
  end function wrap_half_scalar

  real(real32) function metric_distance2(diff, gram) result(distance2)
  real(real32), intent(in) :: diff(3), gram(3, 3)

  distance2 = dot_product(matmul(diff, gram), diff)
  end function metric_distance2

end module artemis__crystal_compare
