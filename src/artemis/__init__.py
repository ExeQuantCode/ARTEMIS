"""
artemis package

Python interface to the ARTEMIS Fortran library for generating
interface lattice matches between crystal structures.
"""

from importlib.metadata import PackageNotFoundError, version
try:
    __version__ = version("artemis-materials")
except PackageNotFoundError:
    __version__ = "unknown"
del PackageNotFoundError, version

import sys
import types

from .artemis import generator as _generator_class
from .artemis import geom_rw as _geom_rw_class
from . import artemis as _artemis_module

# Create simulated 'generator' and 'geom' submodules
generator = types.ModuleType("artemis.generator")
generator.__package__ = __name__
generator.artemis_generator = _generator_class.artemis_generator

geom = types.ModuleType("artemis.geom")
geom.__package__ = __name__
geom.basis_array = _geom_rw_class.basis_array
geom.basis = _geom_rw_class.basis

sys.modules["artemis.generator"] = generator
sys.modules["artemis.geom"] = geom

del _generator_class, _geom_rw_class, types


def get_suppress_warnings():
    """Get the current state of warning suppression."""
    return _artemis_module.artemis.get_suppress_warnings()


def set_suppress_warnings(value):
    """
    Set whether to suppress warnings.

    Parameters
    ----------
    value : bool
        If True, suppress warnings. If False, show warnings.
    """
    _artemis_module.artemis.set_suppress_warnings(value)


__all__ = [
    "__version__",
    "generator",
    "geom",
    "get_suppress_warnings",
    "set_suppress_warnings",
    "suppress_warnings",
]


def __getattr__(name):
    if name == "suppress_warnings":
        return _artemis_module.artemis.get_suppress_warnings()
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __setattr__(name, value):
    if name == "suppress_warnings":
        _artemis_module.artemis.set_suppress_warnings(value)
    else:
        object.__setattr__(sys.modules[__name__], name, value)