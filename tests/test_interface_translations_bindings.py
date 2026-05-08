from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pytest
from ase.io import read


ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "src" / "artemis" / "interface_translations.py"
SPEC = importlib.util.spec_from_file_location("artemis_interface_translations", MODULE_PATH)
REFERENCE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(REFERENCE)

artemis = pytest.importorskip("artemis")


REAL_INTERFACE = ROOT / "example" / "python_pkg" / "MoS2-Ag_0.xyz"


def load_real_interface():
    return read(REAL_INTERFACE)


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


def test_python_binding_matches_reference_on_real_interface() -> None:
    atoms = load_real_interface()
    ref_pair = REFERENCE.get_interface_translations(atoms)
    binding_pair = artemis.get_interface_translations(atoms)

    assert_pair_close(binding_pair, ref_pair)


def test_python_binding_matches_reference_on_single_axis_supercell() -> None:
    atoms = load_real_interface().repeat((2, 1, 1))
    ref_pair = REFERENCE.get_interface_translations(atoms)
    binding_pair = artemis.get_interface_translations(atoms)

    assert_pair_close(binding_pair, ref_pair)
