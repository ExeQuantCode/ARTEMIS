from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
from ase.io import read


ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "src" / "artemis" / "interface_translations.py"
SPEC = importlib.util.spec_from_file_location("artemis_interface_translations", MODULE_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


REAL_INTERFACE = ROOT / "example" / "python_pkg" / "MoS2-Ag_0.xyz"


def load_real_interface():
    return read(REAL_INTERFACE)


def canonical_pair(v1: np.ndarray, v2: np.ndarray) -> tuple[tuple[float, float, float], tuple[float, float, float]]:
    vectors = [np.array(v1, dtype=float), np.array(v2, dtype=float)]
    output = []
    for vector in vectors:
        if vector[0] < -1.0e-8 or (abs(vector[0]) <= 1.0e-8 and vector[1] < -1.0e-8):
            vector *= -1.0
        output.append(tuple(np.round(vector, 8)))
    output.sort()
    return output[0], output[1]


def assert_pair_close(
    actual: tuple[np.ndarray, np.ndarray],
    expected: tuple[np.ndarray, np.ndarray],
    *,
    atol: float = 1.0e-8,
) -> None:
    actual_vectors = sorted(
        (np.array(actual[0], dtype=float), np.array(actual[1], dtype=float)),
        key=lambda vec: (round(np.linalg.norm(vec), 8), round(vec[0], 8), round(vec[1], 8), round(vec[2], 8)),
    )
    expected_vectors = sorted(
        (np.array(expected[0], dtype=float), np.array(expected[1], dtype=float)),
        key=lambda vec: (round(np.linalg.norm(vec), 8), round(vec[0], 8), round(vec[1], 8), round(vec[2], 8)),
    )

    for actual_vec, expected_vec in zip(actual_vectors, expected_vectors, strict=True):
        if actual_vec[0] < -atol or (abs(actual_vec[0]) <= atol and actual_vec[1] < -atol):
            actual_vec *= -1.0
        if expected_vec[0] < -atol or (abs(expected_vec[0]) <= atol and expected_vec[1] < -atol):
            expected_vec *= -1.0
        assert np.allclose(actual_vec, expected_vec, atol=atol)


def test_real_interface_returns_expected_relative_shifts() -> None:
    atoms = load_real_interface()
    t1, t2 = MODULE.get_interface_translations(atoms)

    assert_pair_close(
        (t1, t2),
        (np.array([1.0 / 12.0, 0.0, 0.0]), np.array([0.0, 1.0 / 12.0, 0.0])),
    )
    assert not MODULE.is_valid_interface_translation(atoms, t1)
    assert not MODULE.is_valid_interface_translation(atoms, t2)


def test_rotated_real_interface_keeps_fractional_shifts() -> None:
    atoms = load_real_interface()
    atoms.rotate(27.0, "z", rotate_cell=True)
    t1, t2 = MODULE.get_interface_translations(atoms)

    assert_pair_close(
        (t1, t2),
        (np.array([1.0 / 12.0, 0.0, 0.0]), np.array([0.0, 1.0 / 12.0, 0.0])),
    )


def test_single_axis_supercell_reduces_relative_shifts() -> None:
    atoms = load_real_interface().repeat((2, 1, 1))
    t1, t2 = MODULE.get_interface_translations(atoms)

    assert_pair_close(
        (t1, t2),
        (np.array([1.0 / 24.0, 0.0, 0.0]), np.array([0.0, 1.0 / 12.0, 0.0])),
    )


def test_real_interface_is_deterministic() -> None:
    atoms = load_real_interface()
    first = MODULE.get_interface_translations(atoms)
    second = MODULE.get_interface_translations(atoms)

    assert canonical_pair(*first) == canonical_pair(*second)
