"""
artemis package

This package provides functionality to interface with a Fortran library,
including a Python wrapper around the Fortran code.
"""

from importlib.metadata import PackageNotFoundError, version
try:
    __version__ = version(__name__)
except PackageNotFoundError:
    __version__ = "unknown"

from .artemis import generator as _generator_class
from .artemis import geom_rw as _geom_rw_class
from . import artemis as _artemis_module
# from .artemis import generator


# Use the 'types' module to create simulated 'generator' and 'geom submodules
import types
generator = types.ModuleType('generator')
geom = types.ModuleType('geom')

# Assign the respective class to the simulated 'generator' and 'geom' modules
generator.artemis_generator = _generator_class.artemis_generator

# Assign the class to the simulated 'geom' module
geom.basis_array = _geom_rw_class.basis_array
geom.basis = _geom_rw_class.basis


# Add the simulated 'generator' and 'geom' module to the current package
import sys
sys.modules['artemis.generator'] = generator
sys.modules['artemis.geom'] = geom

# Expose suppress_warnings functions at package level
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

# Clean up internal imports (remove access to the direct classes)
del _generator_class
del _geom_rw_class
del PackageNotFoundError
del version
del sys
del types

__all__ = ['__version__', 'generator', 'geom', 'get_suppress_warnings', 'set_suppress_warnings', 'suppress_warnings']

def __getattr__(name):
    if name == "generator":
        return generator
    elif name == "geom":
        return geom
    elif name == "suppress_warnings":
        return _artemis_module.artemis.get_suppress_warnings()
    raise AttributeError(f"module {__name__} has no attribute {name}")

def __setattr__(name, value):
    if name == "suppress_warnings":
        _artemis_module.artemis.set_suppress_warnings(value)
    else:
        object.__setattr__(sys.modules[__name__], name, value)