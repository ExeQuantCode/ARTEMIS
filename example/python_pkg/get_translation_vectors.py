from __future__ import annotations

import argparse
import importlib.util
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.io import read
from ase.visualize import view

from artemis import generator as artemis_generator_module


REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_FILE = Path(__file__).with_name("MoS2-Ag_0.xyz")
PYTHON_REFERENCE = REPO_ROOT / "src" / "artemis" / "interface_translations.py"


def load_python_reference_module():
	spec = importlib.util.spec_from_file_location("artemis_interface_translations", PYTHON_REFERENCE)
	if spec is None or spec.loader is None:
		raise RuntimeError(f"Unable to load Python reference module from {PYTHON_REFERENCE}")
	module = importlib.util.module_from_spec(spec)
	spec.loader.exec_module(module)
	return module


def array_to_string(values: np.ndarray) -> str:
	return np.array2string(np.asarray(values, dtype=float), precision=6, suppress_small=False)


def fractional_to_cartesian(atoms: Atoms, vector: np.ndarray) -> np.ndarray:
	return np.asarray(vector, dtype=float) @ np.asarray(atoms.cell.array, dtype=float)


def canonicalise_display_vector(vector: np.ndarray, axis: int) -> np.ndarray:
	output = np.asarray(vector, dtype=float).copy()
	output[axis] = 0.0
	plane_axes = [idx for idx in range(3) if idx != axis]
	first = output[plane_axes[0]]
	second = output[plane_axes[1]]
	if first < -1.0e-8 or (abs(first) <= 1.0e-8 and second < -1.0e-8):
		output *= -1.0
	return output


def get_interface_plane(atoms: Atoms) -> tuple[int, np.ndarray, float]:
	generator = artemis_generator_module.artemis_generator()
	bounds, axis_1based = generator.get_interface_location(atoms, return_fractional=True)
	axis = axis_1based - 1
	bounds = np.asarray(bounds, dtype=float)
	upper_bound = bounds[1] if bounds[1] >= bounds[0] else bounds[1] + 1.0
	midpoint = float(0.5 * (bounds[0] + upper_bound)) % 1.0
	return axis, bounds, midpoint


def build_vector_overlay(
	atoms: Atoms,
	t1: np.ndarray,
	t2: np.ndarray,
	axis: int,
	interface_bounds: np.ndarray,
) -> Atoms:
	marker_symbols: list[str] = []
	marker_positions: list[np.ndarray] = []

	plane_axes = [idx for idx in range(3) if idx != axis]
	labelled_vectors = (("Ar", "He", canonicalise_display_vector(t1, axis)), ("Kr", "Ne", canonicalise_display_vector(t2, axis)))
	interface_heights = np.asarray(interface_bounds, dtype=float) % 1.0

	for interface_height in interface_heights:
		for start_symbol, line_symbol, vector in labelled_vectors:
			start_frac = np.zeros(3, dtype=float)
			start_frac[axis] = interface_height
			for plane_axis in plane_axes:
				if vector[plane_axis] < 0.0:
					start_frac[plane_axis] = 1.0

			end_frac = start_frac + vector
			start_cart = fractional_to_cartesian(atoms, start_frac)
			end_cart = fractional_to_cartesian(atoms, end_frac)
			cart_vector = end_cart - start_cart
			vector_length = np.linalg.norm(cart_vector)
			n_markers = max(8, int(np.ceil(vector_length / 0.75)))

			marker_symbols.append(start_symbol)
			marker_positions.append(start_cart)

			for weight in np.linspace(0.0, 1.0, n_markers, endpoint=True)[1:]:
				marker_symbols.append(line_symbol)
				marker_positions.append(start_cart + weight * cart_vector)

	markers = Atoms(symbols=marker_symbols, positions=marker_positions, cell=atoms.cell, pbc=atoms.pbc)
	annotated = atoms.copy()
	annotated += markers
	return annotated


def main() -> None:
	parser = argparse.ArgumentParser(
		description="Print ARTEMIS in-plane translation vectors and view them in ASE."
	)
	parser.add_argument(
		"--path",
		type=Path,
		default=DATA_FILE,
		help="Path to the interface structure file.",
	)
	parser.add_argument(
		"--no-view",
		action="store_true",
		help="Print the vectors without opening the ASE viewer.",
	)
	args = parser.parse_args()

	structure_path = args.path.expanduser().resolve()
	atoms = read(structure_path)
	python_reference = load_python_reference_module()
	t1, t2 = python_reference.get_interface_translations(atoms)
	axis, bounds_frac, midpoint_frac = get_interface_plane(atoms)
	bounds_cart = bounds_frac * np.linalg.norm(np.asarray(atoms.cell.array, dtype=float)[axis])

	print(f"Loaded structure: {structure_path}")
	print(f"Number of atoms: {len(atoms)}")
	print(f"Interface normal axis: {axis + 1}")
	print(f"Interface bounds (fractional): {array_to_string(bounds_frac)}")
	print(f"Interface bounds (angstrom):   {array_to_string(bounds_cart)}")
	print(f"Interface midpoint (fractional): {midpoint_frac:.6f}")
	print(f"t1 fractional: {array_to_string(t1)}")
	print(f"t1 cartesian:  {array_to_string(fractional_to_cartesian(atoms, t1))} A")
	print(f"t2 fractional: {array_to_string(t2)}")
	print(f"t2 cartesian:  {array_to_string(fractional_to_cartesian(atoms, t2))} A")

	if args.no_view:
		return

	annotated = build_vector_overlay(atoms, t1, t2, axis, bounds_frac)
	print("Opening ASE viewer with vector markers on the interface plane.")
	print("Marker legend: Ar/He = t1 start/line, Kr/Ne = t2 start/line.")
	view(annotated)


if __name__ == "__main__":
	main()
