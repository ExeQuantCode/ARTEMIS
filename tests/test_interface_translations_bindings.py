from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from ase.io import read


ROOT = Path(__file__).resolve().parents[1]
artemis = pytest.importorskip("artemis")


REAL_INTERFACE = ROOT / "example" / "python_pkg" / "MoS2-Ag_0.xyz"


def load_real_interface():
    return read(REAL_INTERFACE)


def make_generator():
    return artemis.generator.artemis_generator()


def canonicalise_pair(v1: np.ndarray, v2: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    vectors = [np.array(v1, dtype=float), np.array(v2, dtype=float)]
    output = []
    for vector in vectors:
        if vector[0] < -1.0e-8 or (abs(vector[0]) <= 1.0e-8 and vector[1] < -1.0e-8):
            vector *= -1.0
        output.append(vector)
    output.sort(key=lambda vec: (round(np.linalg.norm(vec), 8), round(vec[0], 8), round(vec[1], 8), round(vec[2], 8)))
    return output[0], output[1]


def assert_pair_close(lhs: tuple[np.ndarray, np.ndarray], rhs: tuple[np.ndarray, np.ndarray], atol: float = 1.0e-2) -> None:
    left = canonicalise_pair(*lhs)
    right = canonicalise_pair(*rhs)
    assert np.allclose(left[0], right[0], atol=atol)
    assert np.allclose(left[1], right[1], atol=atol)


def test_generator_methods_match_basis_input_on_real_interface() -> None:
    atoms = load_real_interface()
    basis = artemis.geom.basis(atoms=atoms)
    generator = make_generator()
    atom_pair = generator.get_interface_translations(atoms)
    basis_pair = generator.get_interface_translations(basis)

    assert_pair_close(atom_pair, basis_pair)


def test_generator_methods_match_basis_input_on_single_axis_supercell() -> None:
    atoms = load_real_interface().repeat((2, 1, 1))
    basis = artemis.geom.basis(atoms=atoms)
    generator = make_generator()
    atom_pair = generator.get_interface_translations(atoms)
    basis_pair = generator.get_interface_translations(basis)

    assert_pair_close(atom_pair, basis_pair)


def test_generator_method_returns_expected_relative_shifts() -> None:
    atoms = load_real_interface()
    generator = make_generator()
    binding_pair = generator.get_interface_translations(atoms)

    assert_pair_close(
        binding_pair,
        (np.array([1.0 / 12.0, 0.0, 0.0]), np.array([0.0, 1.0 / 12.0, 0.0])),
    )
