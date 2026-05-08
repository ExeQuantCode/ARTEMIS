from __future__ import annotations

from itertools import combinations
from typing import Iterable

import numpy as np
import spglib
from ase import Atoms


DEFAULT_CART_TOL = 7.5e-2
DEFAULT_SPGLIB_SYMPREC = 1.0e-2


def as_ase_atoms(structure: Atoms | object) -> Atoms:
    """Return an ASE atoms object from either ASE or ARTEMIS basis input."""
    if isinstance(structure, Atoms):
        return structure.copy()

    if hasattr(structure, "toase"):
        atoms = structure.toase()
        if isinstance(atoms, Atoms):
            return atoms

    raise TypeError("structure must be an ASE Atoms object or an ARTEMIS basis")


def is_valid_interface_translation(
    structure: Atoms | object,
    translation: Iterable[float],
    *,
    axis: int = 2,
    cart_tol: float = DEFAULT_CART_TOL,
) -> bool:
    """Return whether an in-plane translation maps the structure onto itself."""
    atoms = as_ase_atoms(structure)
    shift = _normalise_input_translation(translation, axis)
    scaled = atoms.get_scaled_positions(wrap=True)
    cell = np.asarray(atoms.cell.array, dtype=float)
    numbers = np.asarray(atoms.numbers, dtype=int)

    for number in np.unique(numbers):
        positions = scaled[numbers == number]
        shifted = positions + shift
        shifted -= np.floor(shifted)
        if not _species_mapping_succeeds(positions, shifted, cell, cart_tol):
            return False

    return True


def get_interface_translations(
    structure: Atoms | object,
    *,
    axis: int | None = None,
    cart_tol: float = DEFAULT_CART_TOL,
    spglib_symprec: float = DEFAULT_SPGLIB_SYMPREC,
) -> tuple[np.ndarray, np.ndarray]:
    """Return the two smallest independent relative interface shifts.

    The returned vectors are expressed in fractional coordinates of the input
    cell, as combinations of the input in-plane lattice vectors. The lower
    slab is held fixed while the upper slab is shifted within the interface
    plane.
    """
    atoms = as_ase_atoms(structure)
    bounds, axis = _get_interface_definition(atoms, axis)
    _validate_interface_cell(atoms, axis)
    frac_tol = _fractional_tolerance(atoms.cell.array, axis, cart_tol)

    upper_atoms, lower_atoms = _split_interface_materials(atoms, bounds, axis)

    lower_valid = _collect_valid_translations(
        lower_atoms,
        axis=axis,
        cart_tol=cart_tol,
        spglib_symprec=spglib_symprec,
        frac_tol=frac_tol,
    )
    upper_valid = _collect_valid_translations(
        upper_atoms,
        axis=axis,
        cart_tol=cart_tol,
        spglib_symprec=spglib_symprec,
        frac_tol=frac_tol,
    )

    lower_group = _translation_group(
        lower_valid,
        cell=np.asarray(lower_atoms.cell.array, dtype=float),
        axis=axis,
        frac_tol=frac_tol,
    )
    upper_group = _translation_group(
        upper_valid,
        cell=np.asarray(upper_atoms.cell.array, dtype=float),
        axis=axis,
        frac_tol=frac_tol,
    )
    candidate_vectors = _generate_relative_candidate_vectors(
        lower_group,
        upper_group,
        axis=axis,
        frac_tol=frac_tol,
    )
    valid_vectors = candidate_vectors.copy()

    t1, t2 = _select_primitive_pair(
        atoms.cell.array,
        valid_vectors,
        axis=axis,
        frac_tol=frac_tol,
    )
    return t1, t2


def _validate_interface_cell(atoms: Atoms, axis: int) -> None:
    if axis not in (0, 1, 2):
        raise ValueError("axis must be 0, 1, or 2")

    if not atoms.pbc[axis]:
        raise ValueError("interface normal axis must be periodic")

    plane_axes = _plane_axes(axis)
    if not bool(atoms.pbc[plane_axes[0]]) or not bool(atoms.pbc[plane_axes[1]]):
        raise ValueError("in-plane lattice directions must be periodic")


def _get_interface_definition(atoms: Atoms, axis: int | None) -> tuple[np.ndarray, int]:
    from artemis import generator as artemis_generator_module

    generator = artemis_generator_module.artemis_generator()
    axis_hint = None if axis is None else axis + 1
    bounds, detected_axis = generator.get_interface_location(
        atoms,
        axis=axis_hint,
        return_fractional=True,
    )
    return np.asarray(bounds, dtype=float), int(detected_axis) - 1


def _split_interface_materials(
    atoms: Atoms,
    bounds: np.ndarray,
    axis: int,
) -> tuple[Atoms, Atoms]:
    scaled = atoms.get_scaled_positions(wrap=True)
    lower_bound = float(bounds[0])
    upper_bound = float(bounds[1])
    if upper_bound < lower_bound:
        upper_bound += 1.0

    region_positions = scaled[:, axis] - np.floor(scaled[:, axis] - lower_bound)
    upper_mask = region_positions < upper_bound
    upper_atoms = atoms[upper_mask]
    lower_atoms = atoms[~upper_mask]

    if len(upper_atoms) == 0 or len(lower_atoms) == 0:
        raise RuntimeError("unable to partition interface into upper and lower materials")

    return upper_atoms, lower_atoms


def _plane_axes(axis: int) -> tuple[int, int]:
    if axis == 0:
        return (1, 2)
    if axis == 1:
        return (0, 2)
    return (0, 1)


def _fractional_tolerance(cell: np.ndarray, axis: int, cart_tol: float) -> float:
    plane_axes = _plane_axes(axis)
    plane_lengths = [np.linalg.norm(cell[idx]) for idx in plane_axes]
    return max(cart_tol / max(min(plane_lengths), 1.0e-12), 1.0e-6)


def _normalise_input_translation(translation: Iterable[float], axis: int) -> np.ndarray:
    shift = np.asarray(list(translation), dtype=float)
    if shift.shape == (2,):
        output = np.zeros(3, dtype=float)
        plane_axes = _plane_axes(axis)
        output[plane_axes[0]] = shift[0]
        output[plane_axes[1]] = shift[1]
        return output

    if shift.shape != (3,):
        raise ValueError("translation must have either two or three components")

    output = shift.copy()
    output[axis] = 0.0
    return output


def _wrap_to_half(values: np.ndarray) -> np.ndarray:
    return values - np.floor(values + 0.5)


def _canonicalise_translation(translation: np.ndarray, axis: int, frac_tol: float) -> np.ndarray:
    plane_axes = _plane_axes(axis)
    plane_idx = list(plane_axes)
    output = np.zeros(3, dtype=float)
    output[plane_idx] = _wrap_to_half(np.asarray(translation, dtype=float)[plane_idx])

    first = output[plane_axes[0]]
    second = output[plane_axes[1]]
    if first < -frac_tol or (abs(first) <= frac_tol and second < -frac_tol):
        output *= -1.0

    if abs(output[plane_axes[0]]) <= frac_tol:
        output[plane_axes[0]] = 0.0
    if abs(output[plane_axes[1]]) <= frac_tol:
        output[plane_axes[1]] = 0.0
    return output


def _vector_key(vector: np.ndarray, axis: int, frac_tol: float) -> tuple[int, int]:
    plane_axes = _plane_axes(axis)
    scaled = np.rint(vector[list(plane_axes)] / frac_tol).astype(int)
    return int(scaled[0]), int(scaled[1])


def _append_unique_vector(
    storage: list[np.ndarray],
    translation: np.ndarray,
    *,
    axis: int,
    frac_tol: float,
) -> None:
    candidate = _canonicalise_translation(translation, axis, frac_tol)
    plane_axes = _plane_axes(axis)
    if np.linalg.norm(candidate[list(plane_axes)]) <= frac_tol:
        return

    key = _vector_key(candidate, axis, frac_tol)
    for existing in storage:
        if key == _vector_key(existing, axis, frac_tol) and np.allclose(existing, candidate, atol=frac_tol):
            return

    storage.append(candidate)


def _collect_valid_translations(
    atoms: Atoms,
    *,
    axis: int,
    cart_tol: float,
    spglib_symprec: float,
    frac_tol: float,
) -> list[np.ndarray]:
    candidate_vectors = _generate_candidate_vectors(
        atoms,
        axis=axis,
        cart_tol=cart_tol,
        spglib_symprec=spglib_symprec,
        frac_tol=frac_tol,
    )

    valid_vectors: list[np.ndarray] = []
    for candidate in candidate_vectors:
        if is_valid_interface_translation(atoms, candidate, axis=axis, cart_tol=cart_tol):
            _append_unique_vector(valid_vectors, candidate, axis=axis, frac_tol=frac_tol)

    _expand_valid_vectors(
        atoms,
        valid_vectors,
        axis=axis,
        cart_tol=cart_tol,
        frac_tol=frac_tol,
    )
    return valid_vectors


def _translation_group(
    valid_vectors: list[np.ndarray],
    *,
    cell: np.ndarray,
    axis: int,
    frac_tol: float,
) -> list[np.ndarray]:
    generators = [
        _canonicalise_translation(candidate, axis, frac_tol)
        for candidate in valid_vectors
        if np.linalg.norm(_canonicalise_translation(candidate, axis, frac_tol)[list(_plane_axes(axis))]) > frac_tol
    ]
    if len(generators) >= 2:
        generators = list(_select_primitive_pair(cell, generators, axis=axis, frac_tol=frac_tol))

    group = [np.zeros(3, dtype=float)]
    frontier = [np.zeros(3, dtype=float)]
    max_group_size = max(64, 16 * max(1, len(generators)))

    while frontier:
        current = frontier.pop()
        for generator in generators:
            for candidate in (current + generator, current - generator):
                wrapped = _canonicalise_translation(candidate, axis, frac_tol)
                if np.linalg.norm(wrapped[list(_plane_axes(axis))]) <= frac_tol:
                    wrapped = np.zeros(3, dtype=float)
                if _contains_vector(group, wrapped, axis=axis, frac_tol=frac_tol):
                    continue
                group.append(wrapped)
                frontier.append(wrapped)
                if len(group) > max_group_size:
                    raise RuntimeError("translation group expansion did not converge")

    return group


def _generate_relative_candidate_vectors(
    lower_group: list[np.ndarray],
    upper_group: list[np.ndarray],
    *,
    axis: int,
    frac_tol: float,
) -> list[np.ndarray]:
    candidates: list[np.ndarray] = []

    for upper_vector in upper_group:
        for lower_vector in lower_group:
            _append_unique_vector(
                candidates,
                upper_vector - lower_vector,
                axis=axis,
                frac_tol=frac_tol,
            )

    return candidates


def _contains_vector(
    storage: list[np.ndarray],
    candidate: np.ndarray,
    *,
    axis: int,
    frac_tol: float,
) -> bool:
    key = _vector_key(candidate, axis, frac_tol)
    for existing in storage:
        if key == _vector_key(existing, axis, frac_tol) and np.allclose(existing, candidate, atol=frac_tol):
            return True
    return False


def _is_valid_relative_shift(
    lower_atoms: Atoms,
    upper_atoms: Atoms,
    translation: np.ndarray,
    *,
    axis: int,
    cart_tol: float,
    frac_tol: float,
    lower_group: list[np.ndarray],
) -> bool:
    candidate = _canonicalise_translation(translation, axis, frac_tol)
    plane_axes = _plane_axes(axis)
    if np.linalg.norm(candidate[list(plane_axes)]) <= frac_tol:
        return False

    for global_shift in lower_group:
        if is_valid_interface_translation(
            upper_atoms,
            candidate + global_shift,
            axis=axis,
            cart_tol=cart_tol,
        ):
            return True

    return False


def _generate_candidate_vectors(
    atoms: Atoms,
    *,
    axis: int,
    cart_tol: float,
    spglib_symprec: float,
    frac_tol: float,
) -> list[np.ndarray]:
    scaled = atoms.get_scaled_positions(wrap=True)
    numbers = np.asarray(atoms.numbers, dtype=int)
    candidates: list[np.ndarray] = []

    _collect_spglib_candidates(
        atoms,
        candidates,
        axis=axis,
        frac_tol=frac_tol,
        spglib_symprec=spglib_symprec,
    )

    plane_axes = _plane_axes(axis)
    for number in np.unique(numbers):
        positions = scaled[numbers == number]
        diffs = positions[:, None, :] - positions[None, :, :]
        for candidate in diffs.reshape(-1, 3):
            projected = np.zeros(3, dtype=float)
            projected[list(plane_axes)] = candidate[list(plane_axes)]
            _append_unique_vector(candidates, projected, axis=axis, frac_tol=frac_tol)

    basis_a, basis_b = _basis_vectors(axis)
    _append_unique_vector(candidates, basis_a * 0.5, axis=axis, frac_tol=frac_tol)
    _append_unique_vector(candidates, basis_b * 0.5, axis=axis, frac_tol=frac_tol)
    _append_unique_vector(candidates, basis_a + basis_b, axis=axis, frac_tol=frac_tol)
    _append_unique_vector(candidates, basis_a - basis_b, axis=axis, frac_tol=frac_tol)

    return candidates


def _collect_spglib_candidates(
    atoms: Atoms,
    storage: list[np.ndarray],
    *,
    axis: int,
    frac_tol: float,
    spglib_symprec: float,
) -> None:
    symmetry = spglib.get_symmetry(
        (
            np.asarray(atoms.cell.array, dtype=float),
            np.asarray(atoms.get_scaled_positions(wrap=True), dtype=float),
            np.asarray(atoms.numbers, dtype=int),
        ),
        symprec=spglib_symprec,
    )
    if symmetry is None:
        return

    ident = np.eye(3, dtype=int)
    for rotation, translation in zip(symmetry["rotations"], symmetry["translations"]):
        if not np.array_equal(rotation, ident):
            continue
        _append_unique_vector(storage, np.asarray(translation, dtype=float), axis=axis, frac_tol=frac_tol)


def _species_mapping_succeeds(
    reference_positions: np.ndarray,
    shifted_positions: np.ndarray,
    cell: np.ndarray,
    cart_tol: float,
) -> bool:
    if len(reference_positions) != len(shifted_positions):
        return False

    ref_order = np.lexsort(reference_positions.T[::-1])
    shifted_order = np.lexsort(shifted_positions.T[::-1])
    ref_sorted = reference_positions[ref_order]
    shifted_sorted = shifted_positions[shifted_order]
    unmatched = np.ones(len(ref_sorted), dtype=bool)

    for point in shifted_sorted:
        available = np.flatnonzero(unmatched)
        diffs = ref_sorted[available] - point
        diffs = _wrap_to_half(diffs)
        distances = np.linalg.norm(diffs @ cell, axis=1)
        best = int(np.argmin(distances))
        if distances[best] > cart_tol:
            return False
        unmatched[available[best]] = False

    return True


def _basis_vectors(axis: int) -> tuple[np.ndarray, np.ndarray]:
    plane_axes = _plane_axes(axis)
    basis_a = np.zeros(3, dtype=float)
    basis_b = np.zeros(3, dtype=float)
    basis_a[plane_axes[0]] = 1.0
    basis_b[plane_axes[1]] = 1.0
    return basis_a, basis_b


def _expand_valid_vectors(
    atoms: Atoms,
    valid_vectors: list[np.ndarray],
    *,
    axis: int,
    cart_tol: float,
    frac_tol: float,
) -> None:
    basis_a, basis_b = _basis_vectors(axis)
    seeds = [*valid_vectors, basis_a, basis_b]
    expansions: list[np.ndarray] = []

    for left, right in combinations(seeds, 2):
        for candidate in (left + right, left - right, right - left):
            wrapped = _canonicalise_translation(candidate, axis, frac_tol)
            if not is_valid_interface_translation(atoms, wrapped, axis=axis, cart_tol=cart_tol):
                continue
            _append_unique_vector(expansions, wrapped, axis=axis, frac_tol=frac_tol)

    for candidate in expansions:
        _append_unique_vector(valid_vectors, candidate, axis=axis, frac_tol=frac_tol)


def _expand_relative_vectors(
    lower_atoms: Atoms,
    upper_atoms: Atoms,
    valid_vectors: list[np.ndarray],
    *,
    axis: int,
    cart_tol: float,
    frac_tol: float,
    lower_group: list[np.ndarray],
) -> None:
    basis_a, basis_b = _basis_vectors(axis)
    seeds = [*valid_vectors, basis_a, basis_b]
    expansions: list[np.ndarray] = []

    for left, right in combinations(seeds, 2):
        for candidate in (left + right, left - right, right - left):
            if not _is_valid_relative_shift(
                lower_atoms,
                upper_atoms,
                candidate,
                axis=axis,
                cart_tol=cart_tol,
                frac_tol=frac_tol,
                lower_group=lower_group,
            ):
                continue
            _append_unique_vector(expansions, candidate, axis=axis, frac_tol=frac_tol)

    for candidate in expansions:
        _append_unique_vector(valid_vectors, candidate, axis=axis, frac_tol=frac_tol)


def _select_primitive_pair(
    cell: np.ndarray,
    valid_vectors: list[np.ndarray],
    *,
    axis: int,
    frac_tol: float,
) -> tuple[np.ndarray, np.ndarray]:
    basis_a, basis_b = _basis_vectors(axis)
    span_vectors = valid_vectors if len(valid_vectors) >= 2 else [basis_a, basis_b]
    plane_axes = _plane_axes(axis)
    cell_area = np.linalg.norm(np.cross(cell[plane_axes[0]], cell[plane_axes[1]]))
    area_resolution = max(cell_area * frac_tol, 1.0e-8)
    best_score = None
    best_pair = None

    for left, right in combinations(span_vectors, 2):
        pair = _reduce_translation_pair_valid(
            cell,
            left,
            right,
            valid_vectors=span_vectors,
            axis=axis,
            frac_tol=frac_tol,
        )
        cart_left = pair[0] @ cell
        cart_right = pair[1] @ cell
        area = np.linalg.norm(np.cross(cart_left, cart_right))
        if area <= frac_tol:
            continue

        score = (
            int(np.rint(area / area_resolution)),
            _pair_support(pair, axis=axis, frac_tol=frac_tol),
            round(_pair_fractional_size(pair, axis=axis), 12),
            round(np.linalg.norm(cart_left) ** 2 + np.linalg.norm(cart_right) ** 2, 12),
            round(_pair_cosine(cart_left, cart_right), 12),
            round(np.linalg.norm(cart_right), 12),
            round(area, 12),
            _vector_key(pair[0], axis, frac_tol),
            _vector_key(pair[1], axis, frac_tol),
        )
        if best_score is None or score < best_score:
            best_score = score
            best_pair = (pair[0], pair[1])

    if best_pair is None:
        raise RuntimeError("unable to determine a valid in-plane translation basis")

    return best_pair


def _reduce_translation_pair_valid(
    cell: np.ndarray,
    left: np.ndarray,
    right: np.ndarray,
    *,
    valid_vectors: list[np.ndarray],
    axis: int,
    frac_tol: float,
) -> tuple[np.ndarray, np.ndarray]:
    sort_key = _translation_sort_key(cell)
    current_left, current_right = sorted(
        (
            _canonicalise_translation(left, axis, frac_tol),
            _canonicalise_translation(right, axis, frac_tol),
        ),
        key=sort_key,
    )

    for _ in range(12):
        current_cart_left = current_left @ cell
        current_cart_right = current_right @ cell
        best_pair = (current_left, current_right)
        best_score = (
            _pair_support(best_pair, axis=axis, frac_tol=frac_tol),
            round(_pair_fractional_size(best_pair, axis=axis), 12),
            round(np.linalg.norm(current_cart_left) ** 2 + np.linalg.norm(current_cart_right) ** 2, 12),
            round(_pair_cosine(current_cart_left, current_cart_right), 12),
            round(np.linalg.norm(current_cart_right), 12),
            _vector_key(current_left, axis, frac_tol),
            _vector_key(current_right, axis, frac_tol),
        )

        for next_left, next_right in (
            (current_left, current_right - current_left),
            (current_left, current_right + current_left),
            (current_left - current_right, current_right),
            (current_left + current_right, current_right),
        ):
            next_left = _canonicalise_translation(next_left, axis, frac_tol)
            next_right = _canonicalise_translation(next_right, axis, frac_tol)
            if np.linalg.norm(next_left[list(_plane_axes(axis))]) <= frac_tol:
                continue
            if np.linalg.norm(next_right[list(_plane_axes(axis))]) <= frac_tol:
                continue
            if not _contains_vector(valid_vectors, next_left, axis=axis, frac_tol=frac_tol):
                continue
            if not _contains_vector(valid_vectors, next_right, axis=axis, frac_tol=frac_tol):
                continue
            ordered = tuple(sorted((next_left, next_right), key=sort_key))
            ordered_cart_left = ordered[0] @ cell
            ordered_cart_right = ordered[1] @ cell
            if np.linalg.norm(np.cross(ordered_cart_left, ordered_cart_right)) <= frac_tol:
                continue
            score = (
                _pair_support(ordered, axis=axis, frac_tol=frac_tol),
                round(_pair_fractional_size(ordered, axis=axis), 12),
                round(np.linalg.norm(ordered_cart_left) ** 2 + np.linalg.norm(ordered_cart_right) ** 2, 12),
                round(_pair_cosine(ordered_cart_left, ordered_cart_right), 12),
                round(np.linalg.norm(ordered_cart_right), 12),
                _vector_key(ordered[0], axis, frac_tol),
                _vector_key(ordered[1], axis, frac_tol),
            )
            if score < best_score:
                best_pair = ordered
                best_score = score

        if np.allclose(best_pair[0], current_left, atol=frac_tol) and np.allclose(best_pair[1], current_right, atol=frac_tol):
            break
        current_left, current_right = best_pair

    return current_left, current_right


def _pair_cosine(left: np.ndarray, right: np.ndarray) -> float:
    left_norm = np.linalg.norm(left)
    right_norm = np.linalg.norm(right)
    if left_norm <= 1.0e-16 or right_norm <= 1.0e-16:
        return 1.0
    return float(abs(np.dot(left, right)) / (left_norm * right_norm))


def _pair_support(
    pair: tuple[np.ndarray, np.ndarray],
    *,
    axis: int,
    frac_tol: float,
) -> int:
    plane_axes = _plane_axes(axis)
    return sum(
        int(abs(vector[plane_axis]) > frac_tol)
        for vector in pair
        for plane_axis in plane_axes
    )


def _pair_fractional_size(
    pair: tuple[np.ndarray, np.ndarray],
    *,
    axis: int,
) -> float:
    plane_axes = _plane_axes(axis)
    return float(
        sum(abs(vector[plane_axis]) for vector in pair for plane_axis in plane_axes)
    )


def _translation_sort_key(cell: np.ndarray):
    def sort_key(vector: np.ndarray) -> tuple[float, float, float, float]:
        cart = vector @ cell
        return (
            round(np.linalg.norm(cart), 12),
            round(abs(vector[0]), 12),
            round(abs(vector[1]), 12),
            round(abs(vector[2]), 12),
        )

    return sort_key
