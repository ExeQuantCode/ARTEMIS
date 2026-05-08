from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from ase.io import read


ROOT = Path(__file__).resolve().parents[1]
artemis = pytest.importorskip("artemis")


def load_real_interface():
    return read(ROOT / "example" / "python_pkg" / "MoS2-Ag_0.xyz")


def test_basis_method_accepts_atoms_and_basis_inputs() -> None:
    atoms = load_real_interface()
    basis = artemis.geom.basis(atoms=atoms)

    assert basis.crystals_equivalent(atoms, 1.0e-5)
    assert basis.crystals_equivalent(artemis.geom.basis(atoms=atoms), 1.0e-5)


def test_basis_method_can_ignore_global_translation() -> None:
    atoms = load_real_interface()
    shifted = atoms.copy()
    scaled = shifted.get_scaled_positions(wrap=True)
    scaled += np.array([0.137, -0.241, 0.083])
    scaled -= np.floor(scaled)
    shifted.set_scaled_positions(scaled)

    basis = artemis.geom.basis(atoms=atoms)

    assert not basis.crystals_equivalent(shifted, 1.0e-4, allow_translation=False)
    assert basis.crystals_equivalent(shifted, 1.0e-4, allow_translation=True)
