from __future__ import print_function, absolute_import, division
import artemis._artemis as _artemis
import f90wrap.runtime
import logging
import numpy
from ase import Atoms

class Geom_Rw(f90wrap.runtime.FortranModule):
    """
    Code for handling geometry read/write operations.

    This module provides the necessary functionality to read, write, and
    store atomic geometries.
    In this module, and all of the codebase, element and species are used
    interchangeably.

    Defined in ../src/lib/mod_geom_rw.f90

    .. note::
        It is recommended not to use this module directly, but to handle
        atom objects through the ASE interface.
        This is provided mostly for compatibility with the existing codebase
        and Fortran code.
    """
    @f90wrap.runtime.register_class("artemis.species_type")
    class species_type(f90wrap.runtime.FortranDerivedType):
        def __init__(self, handle=None):
            """
            Create a ``species_type`` object.

            Returns:
                species (species_type):
                    Object to be constructed
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _artemis.f90wrap_geom_rw__species_type_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result

        def __del__(self):
            """
            Destructor for class species_type


            Defined at ../src/lib/mod_geom_rw.f90 lines \
                26-32

            Parameters
            ----------
            this : species_type
            	Object to be destructed


            Automatically generated destructor for species_type
            """
            if self._alloc:
                _artemis.f90wrap_geom_rw__species_type_finalise(this=self._handle)

        @property
        def atom(self):
            """
            Derived type containing the atomic information of a crystal.
            """
            array_ndim, array_type, array_shape, array_handle = \
                _artemis.f90wrap_species_type__array__atom(self._handle)
            if array_handle in self._arrays:
                atom = self._arrays[array_handle]
            else:
                atom = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _artemis.f90wrap_species_type__array__atom)
                self._arrays[array_handle] = atom
            return atom

        @atom.setter
        def atom(self, atom):
            self.atom[...] = atom

        @property
        def mass(self):
            """
            The mass of the element.
            """
            return _artemis.f90wrap_species_type__get__mass(self._handle)

        @mass.setter
        def mass(self, mass):
            _artemis.f90wrap_species_type__set__mass(self._handle, mass)

        @property
        def charge(self):
            """
            The charge of the element.
            """
            return _artemis.f90wrap_species_type__get__charge(self._handle)

        @property
        def radius(self):
            """
            The radius of the element.
            """
            return _artemis.f90wrap_species_type__get__radius(self._handle)

        @radius.setter
        def radius(self, radius):
            _artemis.f90wrap_species_type__set__radius(self._handle, radius)

        @charge.setter
        def charge(self, charge):
            _artemis.f90wrap_species_type__set__charge(self._handle, charge)

        @property
        def name(self):
            """
            The symbol of the element.
            """
            return _artemis.f90wrap_species_type__get__name(self._handle)

        @name.setter
        def name(self, name):
            _artemis.f90wrap_species_type__set__name(self._handle, name)

        @property
        def num(self):
            """
            The number of atoms of this species/element.
            """
            return _artemis.f90wrap_species_type__get__num(self._handle)

        @num.setter
        def num(self, num):
            _artemis.f90wrap_species_type__set__num(self._handle, num)

        def __str__(self):
            ret = ['<species_type>{\n']
            ret.append('    atom : ')
            ret.append(repr(self.atom))
            ret.append(',\n    mass : ')
            ret.append(repr(self.mass))
            ret.append(',\n    charge : ')
            ret.append(repr(self.charge))
            ret.append(',\n    name : ')
            ret.append(repr(self.name))
            ret.append(',\n    num : ')
            ret.append(repr(self.num))
            ret.append('}')
            return ''.join(ret)

        _dt_array_initialisers = []


    @f90wrap.runtime.register_class("artemis.basis")
    class basis(f90wrap.runtime.FortranDerivedType):
        def __init__(self, atoms=None, handle=None):
            """
            Create a ``basis`` object.

            This object is used to store the atomic information of a crystal,
            including lattice and basis information.
            This is confusingly named as a crystal = lattice + basis.

            Returns:
                basis (basis):
                    Object to be constructed
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _artemis.f90wrap_geom_rw__basis_type_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result

            if atoms is not None:
                self.fromase(atoms)

        def __del__(self):
            """
            Destructor for class basis


            Defined at ../src/lib/mod_geom_rw.f90 lines \
                34-42

            Parameters
            ----------
            this : basis
            	Object to be destructed


            Automatically generated destructor for basis
            """
            if self._alloc:
                _artemis.f90wrap_geom_rw__basis_type_finalise(this=self._handle)

        def allocate_species(self, num_species=None, species_symbols=None, species_count=None, \
            positions=None):
            """
            Allocate memory for the species list.

            Parameters:
                num_species (int):
                    Number of species
                species_symbols (list of str):
                    List of species symbols
                species_count (list of int):
                    List of species counts
                atoms (list of float):
                    List of atomic positions
            """
            _artemis.f90wrap_geom_rw__allocate_species__binding__basis_type(this=self._handle, \
                num_species=num_species, species_symbols=species_symbols, species_count=species_count, \
                atoms=positions)

        def _init_array_spec(self):
            """
            Initialise the species array.
            """
            self.spec = f90wrap.runtime.FortranDerivedTypeArray(self,
                                            _artemis.f90wrap_basis_type__array_getitem__spec,
                                            _artemis.f90wrap_basis_type__array_setitem__spec,
                                            _artemis.f90wrap_basis_type__array_len__spec,
                                            """
            Element spec ftype=type(species_type) pytype=species_type


            Defined at ../src/lib/mod_geom_rw.f90 line 35

            """, Geom_Rw.species_type)
            return self.spec

        def toase(self, calculator=None):
            """
            Convert the basis object to an ASE Atoms object.

            Parameters:
                calculator (ASE Calculator):
                    ASE calculator object to be assigned to the Atoms object.
            """
            from ase import Atoms

            # Set the species list
            positions = []
            species_string = ""
            for i in range(self.nspec):
                for j in range(self.spec[i].num):
                    species_string += str(self.spec[i].name.decode()).strip()
                    positions.append(self.spec[i].atom[j])

            # Set the atoms
            if(self.lcart):
                atoms = Atoms(species_string, positions=positions, cell=self.lat, pbc=self.pbc)
            else:
                atoms = Atoms(species_string, scaled_positions=positions, cell=self.lat, pbc=self.pbc)

            if calculator is not None:
                atoms.calc = calculator
            return atoms

        def fromase(self, atoms, verbose=False):
            """
            Convert the ASE Atoms object to a basis object.

            Parameters:
                atoms (ASE Atoms):
                    ASE Atoms object to be converted.
                verbose (bool):
                    Boolean whether to print warnings.
            """
            from ase.calculators.singlepoint import SinglePointCalculator

            # Get the species symbols
            species_symbols = atoms.get_chemical_symbols()
            species_symbols_unique = sorted(set(species_symbols))

            # Set the number of species
            self.nspec = len(species_symbols_unique)

            # Set the number of atoms
            self.natom = len(atoms)

            # check if calculator is present
            if atoms.calc is None:
                if verbose:
                    print("WARNING: No calculator present, setting energy to 0.0")
                atoms.calc = SinglePointCalculator(atoms, energy=0.0)
            self.energy = atoms.get_potential_energy()

            # # Set the lattice vectors
            self.lat = numpy.reshape(atoms.get_cell().flatten(), [3,3], order='A')
            self.pbc = atoms.pbc

            # Set the system name
            self.sysname = atoms.get_chemical_formula()

            # Set the species list
            species_count = []
            atom_positions = []
            positions = atoms.get_scaled_positions()
            for species in species_symbols_unique:
                species_count.append(sum([1 for symbol in species_symbols if symbol == species]))
                for j, symbol in enumerate(species_symbols):
                    if symbol == species:
                        atom_positions.append(positions[j])

            # Allocate memory for the atom list
            self.lcart = False
            self.allocate_species(species_symbols=species_symbols_unique, species_count=species_count, positions=atom_positions)

        @property
        def nspec(self):
            """
            The number of species in the basis.
            """
            return _artemis.f90wrap_basis_type__get__nspec(self._handle)

        @nspec.setter
        def nspec(self, nspec):
            _artemis.f90wrap_basis_type__set__nspec(self._handle, nspec)

        @property
        def natom(self):
            """
            The number of atoms in the basis.
            """
            return _artemis.f90wrap_basis_type__get__natom(self._handle)

        @natom.setter
        def natom(self, natom):
            _artemis.f90wrap_basis_type__set__natom(self._handle, natom)

        @property
        def energy(self):
            """
            The energy associated with the basis (or crystal).
            """
            return _artemis.f90wrap_basis_type__get__energy(self._handle)

        @energy.setter
        def energy(self, energy):
            _artemis.f90wrap_basis_type__set__energy(self._handle, energy)

        @property
        def lat(self):
            """
            The lattice vectors of the basis.
            """
            array_ndim, array_type, array_shape, array_handle = \
                _artemis.f90wrap_basis_type__array__lat(self._handle)
            if array_handle in self._arrays:
                lat = self._arrays[array_handle]
            else:
                lat = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _artemis.f90wrap_basis_type__array__lat)
                self._arrays[array_handle] = lat
            return lat

        @lat.setter
        def lat(self, lat):
            self.lat[...] = lat

        @property
        def lcart(self):
            """
            Boolean whether the atomic positions are in cartesian coordinates.
            """
            return _artemis.f90wrap_basis_type__get__lcart(self._handle)

        @lcart.setter
        def lcart(self, lcart):
            _artemis.f90wrap_basis_type__set__lcart(self._handle, lcart)

        @property
        def pbc(self):
            """
            Boolean array indicating the periodic boundary conditions.
            """
            array_ndim, array_type, array_shape, array_handle = \
                _artemis.f90wrap_basis_type__array__pbc(self._handle)
            if array_handle in self._arrays:
                pbc = self._arrays[array_handle]
            else:
                pbc = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _artemis.f90wrap_basis_type__array__pbc)
                self._arrays[array_handle] = pbc
            return pbc

        @pbc.setter
        def pbc(self, pbc):
            self.pbc[...] = pbc

        @property
        def sysname(self):
            """
            The name of the system.
            """
            return _artemis.f90wrap_basis_type__get__sysname(self._handle)

        @sysname.setter
        def sysname(self, sysname):
            _artemis.f90wrap_basis_type__set__sysname(self._handle, sysname)

        def __str__(self):
            ret = ['<basis>{\n']
            ret.append('    nspec : ')
            ret.append(repr(self.nspec))
            ret.append(',\n    natom : ')
            ret.append(repr(self.natom))
            ret.append(',\n    energy : ')
            ret.append(repr(self.energy))
            ret.append(',\n    lat : ')
            ret.append(repr(self.lat))
            ret.append(',\n    lcart : ')
            ret.append(repr(self.lcart))
            ret.append(',\n    pbc : ')
            ret.append(repr(self.pbc))
            ret.append(',\n    sysname : ')
            ret.append(repr(self.sysname))
            ret.append('}')
            return ''.join(ret)

        _dt_array_initialisers = [_init_array_spec]



    @f90wrap.runtime.register_class("artemis.basis_array")
    class basis_array(f90wrap.runtime.FortranDerivedType):
        def __init__(self, atoms=None, handle=None):
            """
            Create a ``basis_array`` object.


            Returns:
                basis_array (basis_array):
                    Object to be constructed
            """

            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _artemis.f90wrap_geom_rw__basis_type_xnum_array_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result


            # check if atoms is an ASE Atoms object or a list of ASE Atoms objects
            if atoms:
                from ase import Atoms
                if isinstance(atoms, Atoms):
                    self.allocate(1)
                    self.items[0].fromase(atoms)
                elif isinstance(atoms, list):
                    self.allocate(len(atoms))
                    for i, atom in enumerate(atoms):
                        self.items[i].fromase(atom)

        def __del__(self):
            """
            Destructor for class basis_array


            Defined at ../src/lib/mod_generator.f90 lines \
                19-21

            Parameters
            ----------
            this : basis_array
            	Object to be destructed


            Automatically generated destructor for basis_array
            """
            if self._alloc:
                _artemis.f90wrap_geom_rw__basis_type_xnum_array_finalise(this=self._handle)

        def _init_array_items(self):
            """
            Initialise the items array.
            """
            self.items = f90wrap.runtime.FortranDerivedTypeArray(self,
                                            _artemis.f90wrap_basis_type_xnum_array__array_getitem__items,
                                            _artemis.f90wrap_basis_type_xnum_array__array_setitem__items,
                                            _artemis.f90wrap_basis_type_xnum_array__array_len__items,
                                            """
            Element items ftype=type(basis_type) pytype=basis


            Defined at  line 0

            """, Geom_Rw.basis)
            return self.items

        def toase(self):
            """
            Convert the basis_array object to a list of ASE Atoms objects.
            """

            # Set the species list
            atoms = []
            for i in range(len(self.items)):
                atoms.append(self.items[i].toase())
            return atoms

        def allocate(self, size):
            """
            Allocate the items array with the given size.

            Parameters:
                size (int):
                    Size of the items array
            """
            _artemis.f90wrap_basis_type_xnum_array__array_alloc__items(self._handle, num=size)

        def deallocate(self):
            """
            Deallocate the items array
            """
            _artemis.f90wrap_basis_type_xnum_array__array_dealloc__items(self._handle)

        _dt_array_initialisers = [_init_array_items]

    _dt_array_initialisers = []


geom_rw = Geom_Rw()

# class Geom_Rw(f90wrap.runtime.FortranModule):
#     """
#     Module geom_rw
    
    
#     Defined at ../src/fortran/lib/mod_geom_rw.f90 \
#         lines 13-1907
    
#     """
#     @f90wrap.runtime.register_class("artemis.species_type")
#     class species_type(f90wrap.runtime.FortranDerivedType):
#         """
#         Type(name=species_type)
        
        
#         Defined at ../src/fortran/lib/mod_geom_rw.f90 \
#             lines 34-47
        
#         """
#         def __init__(self, handle=None):
#             """
#             self = Species_Type()
            
            
#             Defined at ../src/fortran/lib/mod_geom_rw.f90 \
#                 lines 34-47
            
            
#             Returns
#             -------
#             this : Species_Type
#             	Object to be constructed
            
            
#             Automatically generated constructor for species_type
#             """
#             f90wrap.runtime.FortranDerivedType.__init__(self)
#             result = _artemis.f90wrap_geom_rw__species_type_initialise()
#             self._handle = result[0] if isinstance(result, tuple) else result
        
#         def __del__(self):
#             """
#             Destructor for class Species_Type
            
            
#             Defined at ../src/fortran/lib/mod_geom_rw.f90 \
#                 lines 34-47
            
#             Parameters
#             ----------
#             this : Species_Type
#             	Object to be destructed
            
            
#             Automatically generated destructor for species_type
#             """
#             if self._alloc:
#                 _artemis.f90wrap_geom_rw__species_type_finalise(this=self._handle)
        
#         @property
#         def atom(self):
#             """
#             Element atom ftype=real(real32) pytype=float
            
            
#             Defined at ../src/fortran/lib/mod_geom_rw.f90 \
#                 line 36
            
#             """
#             array_ndim, array_type, array_shape, array_handle = \
#                 _artemis.f90wrap_species_type__array__atom(self._handle)
#             if array_handle in self._arrays:
#                 atom = self._arrays[array_handle]
#             else:
#                 atom = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
#                                         self._handle,
#                                         _artemis.f90wrap_species_type__array__atom)
#                 self._arrays[array_handle] = atom
#             return atom
        
#         @atom.setter
#         def atom(self, atom):
#             self.atom[...] = atom
        
#         @property
#         def mass(self):
#             """
#             Element mass ftype=real(real32) pytype=float
            
            
#             Defined at ../src/fortran/lib/mod_geom_rw.f90 \
#                 line 38
            
#             """
#             return _artemis.f90wrap_species_type__get__mass(self._handle)
        
#         @mass.setter
#         def mass(self, mass):
#             _artemis.f90wrap_species_type__set__mass(self._handle, mass)
        
#         @property
#         def charge(self):
#             """
#             Element charge ftype=real(real32) pytype=float
            
            
#             Defined at ../src/fortran/lib/mod_geom_rw.f90 \
#                 line 40
            
#             """
#             return _artemis.f90wrap_species_type__get__charge(self._handle)
        
#         @charge.setter
#         def charge(self, charge):
#             _artemis.f90wrap_species_type__set__charge(self._handle, charge)
        
#         @property
#         def radius(self):
#             """
#             Element radius ftype=real(real32) pytype=float
            
            
#             Defined at ../src/fortran/lib/mod_geom_rw.f90 \
#                 line 42
            
#             """
#             return _artemis.f90wrap_species_type__get__radius(self._handle)
        
#         @radius.setter
#         def radius(self, radius):
#             _artemis.f90wrap_species_type__set__radius(self._handle, radius)
        
#         @property
#         def name(self):
#             """
#             Element name ftype=character(len=3) pytype=str
            
            
#             Defined at ../src/fortran/lib/mod_geom_rw.f90 \
#                 line 44
            
#             """
#             return _artemis.f90wrap_species_type__get__name(self._handle)
        
#         @name.setter
#         def name(self, name):
#             _artemis.f90wrap_species_type__set__name(self._handle, name)
        
#         @property
#         def num(self):
#             """
#             Element num ftype=integer  pytype=int
            
            
#             Defined at ../src/fortran/lib/mod_geom_rw.f90 \
#                 line 46
            
#             """
#             return _artemis.f90wrap_species_type__get__num(self._handle)
        
#         @num.setter
#         def num(self, num):
#             _artemis.f90wrap_species_type__set__num(self._handle, num)
        
#         def __str__(self):
#             ret = ['<species_type>{\n']
#             ret.append('    atom : ')
#             ret.append(repr(self.atom))
#             ret.append(',\n    mass : ')
#             ret.append(repr(self.mass))
#             ret.append(',\n    charge : ')
#             ret.append(repr(self.charge))
#             ret.append(',\n    radius : ')
#             ret.append(repr(self.radius))
#             ret.append(',\n    name : ')
#             ret.append(repr(self.name))
#             ret.append(',\n    num : ')
#             ret.append(repr(self.num))
#             ret.append('}')
#             return ''.join(ret)
        
#         _dt_array_initialisers = []
        
    
#     @f90wrap.runtime.register_class("artemis.basis_type")
#     class basis_type(f90wrap.runtime.FortranDerivedType):
#         """
#         Type(name=basis_type)
        
        
#         Defined at ../src/fortran/lib/mod_geom_rw.f90 \
#             lines 49-83
        
#         """
#         def __init__(self, handle=None):
#             """
#             self = Basis_Type()
            
            
#             Defined at ../src/fortran/lib/mod_geom_rw.f90 \
#                 lines 49-83
            
            
#             Returns
#             -------
#             this : Basis_Type
#             	Object to be constructed
            
            
#             Automatically generated constructor for basis_type
#             """
#             f90wrap.runtime.FortranDerivedType.__init__(self)
#             result = _artemis.f90wrap_geom_rw__basis_type_initialise()
#             self._handle = result[0] if isinstance(result, tuple) else result
        
#         def __del__(self):
#             """
#             Destructor for class Basis_Type
            
            
#             Defined at ../src/fortran/lib/mod_geom_rw.f90 \
#                 lines 49-83
            
#             Parameters
#             ----------
#             this : Basis_Type
#             	Object to be destructed
            
            
#             Automatically generated destructor for basis_type
#             """
#             if self._alloc:
#                 _artemis.f90wrap_geom_rw__basis_type_finalise(this=self._handle)
        
#         def allocate_species(self, num_species=None, species_symbols=None, \
#             species_count=None, atoms=None):
#             """
#             allocate_species__binding__basis_type(self[, num_species, species_symbols, \
#                 species_count, atoms])
            
            
#             Defined at ../src/fortran/lib/mod_geom_rw.f90 \
#                 lines 110-153
            
#             Parameters
#             ----------
#             this : Basis_Type
#             num_species : int
#             species_symbols : str array
#             species_count : int array
#             atoms : float array
            
#             """
#             _artemis.f90wrap_geom_rw__allocate_species__binding__basis_type(this=self._handle, \
#                 num_species=num_species, species_symbols=species_symbols, \
#                 species_count=species_count, atoms=atoms)
        
#         def convert(self):
#             """
#             convert__binding__basis_type(self)
            
            
#             Defined at ../src/fortran/lib/mod_geom_rw.f90 \
#                 lines 1035-1057
            
#             Parameters
#             ----------
#             this : Basis_Type
            
#             """
#             _artemis.f90wrap_geom_rw__convert__binding__basis_type(this=self._handle)
        
#         def change_lattice(self, lattice):
#             """
#             change_lattice__binding__basis_type(self, lattice)
            
            
#             Defined at ../src/fortran/lib/mod_geom_rw.f90 \
#                 lines 1061-1088
            
#             Parameters
#             ----------
#             this : Basis_Type
#             lattice : float array
            
#             """
#             _artemis.f90wrap_geom_rw__change_lattice__binding__basis_type(this=self._handle, \
#                 lattice=lattice)
        
#         def normalise(self, ceil_val=None, floor_coords=None, round_coords=None, \
#             zero_round=None):
#             """
#             normalise__binding__basis_type(self[, ceil_val, floor_coords, round_coords, \
#                 zero_round])
            
            
#             Defined at ../src/fortran/lib/mod_geom_rw.f90 \
#                 lines 1097-1147
            
#             Parameters
#             ----------
#             this : Basis_Type
#             ceil_val : float
#             floor_coords : bool
#             round_coords : bool
#             zero_round : float
            
#             """
#             _artemis.f90wrap_geom_rw__normalise__binding__basis_type(this=self._handle, \
#                 ceil_val=ceil_val, floor_coords=floor_coords, round_coords=round_coords, \
#                 zero_round=zero_round)
        
#         def copy(self, basis, length=None):
#             """
#             copy__binding__basis_type(self, basis[, length])
            
            
#             Defined at ../src/fortran/lib/mod_geom_rw.f90 \
#                 lines 1229-1290
            
#             Parameters
#             ----------
#             this : Basis_Type
#             basis : Basis_Type
#             length : int
            
#             ---------------------------------------------------------------------------
#              determines whether user wants output basis extra translational dimension
#             ---------------------------------------------------------------------------
#             """
#             _artemis.f90wrap_geom_rw__copy__binding__basis_type(this=self._handle, \
#                 basis=basis._handle, length=length)
        
#         def get_lattice_constants(self, radians=None):
#             """
#             output = get_lattice_constants__binding__basis_type(self[, radians])
            
            
#             Defined at ../src/fortran/lib/mod_geom_rw.f90 \
#                 lines 1210-1225
            
#             Parameters
#             ----------
#             this : Basis_Type
#             radians : bool
            
#             Returns
#             -------
#             output : float array
            
#             """
#             output = \
#                 _artemis.f90wrap_geom_rw__get_lattice_constants__binding__bc9a1(this=self._handle, \
#                 radians=radians)
#             return output
        
#         def remove_atom(self, ispec, iatom):
#             """
#             remove_atom__binding__basis_type(self, ispec, iatom)
            
            
#             Defined at ../src/fortran/lib/mod_geom_rw.f90 \
#                 lines 1294-1336
            
#             Parameters
#             ----------
#             this : Basis_Type
#             ispec : int
#             iatom : int
            
#             ---------------------------------------------------------------------------
#              remove atom from basis
#             ---------------------------------------------------------------------------
#             """
#             _artemis.f90wrap_geom_rw__remove_atom__binding__basis_type(this=self._handle, \
#                 ispec=ispec, iatom=iatom)
        
#         def remove_atoms(self, atoms):
#             """
#             remove_atoms__binding__basis_type(self, atoms)
            
            
#             Defined at ../src/fortran/lib/mod_geom_rw.f90 \
#                 lines 1340-1403
            
#             Parameters
#             ----------
#             this : Basis_Type
#             atoms : int array
            
#             ---------------------------------------------------------------------------
#              reorder atoms to remove
#             ---------------------------------------------------------------------------
#             """
#             _artemis.f90wrap_geom_rw__remove_atoms__binding__basis_type(this=self._handle, \
#                 atoms=atoms)
        
#         def init_array_spec(self):
#             self.spec = f90wrap.runtime.FortranDerivedTypeArray(self,
#                                             _artemis.f90wrap_basis_type__array_getitem__spec,
#                                             _artemis.f90wrap_basis_type__array_setitem__spec,
#                                             _artemis.f90wrap_basis_type__array_len__spec,
#                                             """
#             Element spec ftype=type(species_type) pytype=Species_Type
            
            
#             Defined at ../src/fortran/lib/mod_geom_rw.f90 \
#                 line 51
            
#             """, Geom_Rw.species_type)
#             return self.spec
        
#         @property
#         def nspec(self):
#             """
#             Element nspec ftype=integer  pytype=int
            
            
#             Defined at ../src/fortran/lib/mod_geom_rw.f90 \
#                 line 53
            
#             """
#             return _artemis.f90wrap_basis_type__get__nspec(self._handle)
        
#         @nspec.setter
#         def nspec(self, nspec):
#             _artemis.f90wrap_basis_type__set__nspec(self._handle, nspec)
        
#         @property
#         def natom(self):
#             """
#             Element natom ftype=integer  pytype=int
            
            
#             Defined at ../src/fortran/lib/mod_geom_rw.f90 \
#                 line 55
            
#             """
#             return _artemis.f90wrap_basis_type__get__natom(self._handle)
        
#         @natom.setter
#         def natom(self, natom):
#             _artemis.f90wrap_basis_type__set__natom(self._handle, natom)
        
#         @property
#         def energy(self):
#             """
#             Element energy ftype=real(real32) pytype=float
            
            
#             Defined at ../src/fortran/lib/mod_geom_rw.f90 \
#                 line 57
            
#             """
#             return _artemis.f90wrap_basis_type__get__energy(self._handle)
        
#         @energy.setter
#         def energy(self, energy):
#             _artemis.f90wrap_basis_type__set__energy(self._handle, energy)
        
#         @property
#         def lat(self):
#             """
#             Element lat ftype=real(real32) pytype=float
            
            
#             Defined at ../src/fortran/lib/mod_geom_rw.f90 \
#                 line 59
            
#             """
#             array_ndim, array_type, array_shape, array_handle = \
#                 _artemis.f90wrap_basis_type__array__lat(self._handle)
#             if array_handle in self._arrays:
#                 lat = self._arrays[array_handle]
#             else:
#                 lat = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
#                                         self._handle,
#                                         _artemis.f90wrap_basis_type__array__lat)
#                 self._arrays[array_handle] = lat
#             return lat
        
#         @lat.setter
#         def lat(self, lat):
#             self.lat[...] = lat
        
#         @property
#         def lcart(self):
#             """
#             Element lcart ftype=logical pytype=bool
            
            
#             Defined at ../src/fortran/lib/mod_geom_rw.f90 \
#                 line 61
            
#             """
#             return _artemis.f90wrap_basis_type__get__lcart(self._handle)
        
#         @lcart.setter
#         def lcart(self, lcart):
#             _artemis.f90wrap_basis_type__set__lcart(self._handle, lcart)
        
#         @property
#         def pbc(self):
#             """
#             Element pbc ftype=logical pytype=bool
            
            
#             Defined at ../src/fortran/lib/mod_geom_rw.f90 \
#                 line 63
            
#             """
#             array_ndim, array_type, array_shape, array_handle = \
#                 _artemis.f90wrap_basis_type__array__pbc(self._handle)
#             if array_handle in self._arrays:
#                 pbc = self._arrays[array_handle]
#             else:
#                 pbc = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
#                                         self._handle,
#                                         _artemis.f90wrap_basis_type__array__pbc)
#                 self._arrays[array_handle] = pbc
#             return pbc
        
#         @pbc.setter
#         def pbc(self, pbc):
#             self.pbc[...] = pbc
        
#         @property
#         def sysname(self):
#             """
#             Element sysname ftype=character(len=128) pytype=str
            
            
#             Defined at ../src/fortran/lib/mod_geom_rw.f90 \
#                 line 65
            
#             """
#             return _artemis.f90wrap_basis_type__get__sysname(self._handle)
        
#         @sysname.setter
#         def sysname(self, sysname):
#             _artemis.f90wrap_basis_type__set__sysname(self._handle, sysname)
        
#         def __str__(self):
#             ret = ['<basis_type>{\n']
#             ret.append('    nspec : ')
#             ret.append(repr(self.nspec))
#             ret.append(',\n    natom : ')
#             ret.append(repr(self.natom))
#             ret.append(',\n    energy : ')
#             ret.append(repr(self.energy))
#             ret.append(',\n    lat : ')
#             ret.append(repr(self.lat))
#             ret.append(',\n    lcart : ')
#             ret.append(repr(self.lcart))
#             ret.append(',\n    pbc : ')
#             ret.append(repr(self.pbc))
#             ret.append(',\n    sysname : ')
#             ret.append(repr(self.sysname))
#             ret.append('}')
#             return ''.join(ret)
        
#         _dt_array_initialisers = [init_array_spec]
        
    
#     @staticmethod
#     def geom_read(unit, length=None, iostat=None):
#         """
#         basis = geom_read(unit[, length, iostat])
        
        
#         Defined at ../src/fortran/lib/mod_geom_rw.f90 \
#             lines 157-212
        
#         Parameters
#         ----------
#         unit : int
#         length : int
#         iostat : int
        
#         Returns
#         -------
#         basis : Basis_Type
        
#         """
#         basis = _artemis.f90wrap_geom_rw__geom_read(unit=unit, length=length, \
#             iostat=iostat)
#         basis = f90wrap.runtime.lookup_class("artemis.basis_type").from_handle(basis, \
#             alloc=True)
#         return basis
    
#     @staticmethod
#     def geom_write(unit, basis):
#         """
#         geom_write(unit, basis)
        
        
#         Defined at ../src/fortran/lib/mod_geom_rw.f90 \
#             lines 216-240
        
#         Parameters
#         ----------
#         unit : int
#         basis : Basis_Type
        
#         """
#         _artemis.f90wrap_geom_rw__geom_write(unit=unit, basis=basis._handle)
    
#     @staticmethod
#     def get_element_properties(element, charge=None, mass=None, radius=None):
#         """
#         get_element_properties(element[, charge, mass, radius])
        
        
#         Defined at ../src/fortran/lib/mod_geom_rw.f90 \
#             lines 1407-1906
        
#         Parameters
#         ----------
#         element : str
#         charge : float
#         mass : float
#         radius : float
        
#         ---------------------------------------------------------------------------
#          Return the values
#         ---------------------------------------------------------------------------
#         """
#         _artemis.f90wrap_geom_rw__get_element_properties(element=element, \
#             charge=charge, mass=mass, radius=radius)
    
#     @property
#     def igeom_input(self):
#         """
#         Element igeom_input ftype=integer  pytype=int
        
        
#         Defined at ../src/fortran/lib/mod_geom_rw.f90 \
#             line 24
        
#         """
#         return _artemis.f90wrap_geom_rw__get__igeom_input()
    
#     @igeom_input.setter
#     def igeom_input(self, igeom_input):
#         _artemis.f90wrap_geom_rw__set__igeom_input(igeom_input)
    
#     @property
#     def igeom_output(self):
#         """
#         Element igeom_output ftype=integer  pytype=int
        
        
#         Defined at ../src/fortran/lib/mod_geom_rw.f90 \
#             line 32
        
#         """
#         return _artemis.f90wrap_geom_rw__get__igeom_output()
    
#     @igeom_output.setter
#     def igeom_output(self, igeom_output):
#         _artemis.f90wrap_geom_rw__set__igeom_output(igeom_output)
    
#     def __str__(self):
#         ret = ['<geom_rw>{\n']
#         ret.append('    igeom_input : ')
#         ret.append(repr(self.igeom_input))
#         ret.append(',\n    igeom_output : ')
#         ret.append(repr(self.igeom_output))
#         ret.append('}')
#         return ''.join(ret)
    
#     _dt_array_initialisers = []
    

# geom_rw = Geom_Rw()

class Termination_Generator(f90wrap.runtime.FortranModule):
    """
    Module artemis__termination_generator
    
    
    Defined at \
        ../src/fortran/lib/mod_term_generator.f90 \
        lines 7-202
    
    """
    @f90wrap.runtime.register_class("artemis.artemis_termination_generator")
    class artemis_termination_generator(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=artemis_termination_generator_type)
        
        
        Defined at \
            ../src/fortran/lib/mod_term_generator.f90 \
            lines 21-24
        
        """
        def __init__(self, handle=None):
            """
            self = Artemis_Termination_Generator_Type()
            
            
            Defined at \
                ../src/fortran/lib/mod_term_generator.f90 \
                lines 21-24
            
            
            Returns
            -------
            this : Artemis_Termination_Generator_Type
            	Object to be constructed
            
            
            Automatically generated constructor for artemis_termination_generator_type
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = \
                _artemis.f90wrap_term_gen__artemis_termination293d()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Artemis_Termination_Generator_Type
            
            
            Defined at \
                ../src/fortran/lib/mod_term_generator.f90 \
                lines 21-24
            
            Parameters
            ----------
            this : Artemis_Termination_Generator_Type
            	Object to be destructed
            
            
            Automatically generated destructor for artemis_termination_generator_type
            """
            if self._alloc:
                _artemis.f90wrap_term_gen__artemis_terminationdf16(this=self._handle)
        
        def generate(self, basis, miller_plane, axis, surface=None, num_layers=None, \
            thickness=None, orthogonalise=None, normalise=None, break_on_fail=None):
            """
            generate__binding__artemis_termination_generator_type(self, basis, miller_plane, \
                axis[, surface, num_layers, thickness, orthogonalise, normalise, \
                break_on_fail])
            
            
            Defined at \
                ../src/fortran/lib/mod_term_generator.f90 \
                lines 31-201
            
            Parameters
            ----------
            this : Artemis_Termination_Generator_Type
            basis : Basis_Type
            miller_plane : int array
            axis : int
            surface : int array
            num_layers : int
            thickness : float
            orthogonalise : bool
            normalise : bool
            break_on_fail : bool
            
            ---------------------------------------------------------------------------
             Finds smallest thickness of the slab and increases to ...
             ... user-defined thickness
            ---------------------------------------------------------------------------
            """

            # check if host is ase.Atoms object or a Fortran derived type basis_type
            if isinstance(basis, Atoms):
                basis = geom_rw.basis(atoms=basis)

            _artemis.f90wrap_term_gen__generate__binding__2af7(this=self._handle, \
                basis=basis._handle, miller_plane=miller_plane, axis=axis, surface=surface, \
                num_layers=num_layers, thickness=thickness, orthogonalise=orthogonalise, \
                normalise=normalise, break_on_fail=break_on_fail)
        
        @property
        def layer_separation_cutoff(self):
            """
            Element layer_separation_cutoff ftype=real(real32) pytype=float
            
            
            Defined at \
                ../src/fortran/lib/mod_term_generator.f90 \
                line 22
            
            """
            return \
                _artemis.f90wrap_artemis_termination_generator_type__get__layer_sepace78(self._handle)
        
        @layer_separation_cutoff.setter
        def layer_separation_cutoff(self, layer_separation_cutoff):
            _artemis.f90wrap_artemis_termination_generator_type__set__layer_sepae7ef(self._handle, \
                layer_separation_cutoff)
        
        def __str__(self):
            ret = ['<artemis_termination_generator_type>{\n']
            ret.append('    layer_separation_cutoff : ')
            ret.append(repr(self.layer_separation_cutoff))
            ret.append('}')
            return ''.join(ret)
        
        _dt_array_initialisers = []
        
    
    _dt_array_initialisers = []
    

termination_generator = Termination_Generator()

class Interface_Generator(f90wrap.runtime.FortranModule):
    """
    Module artemis__interface_generator
    
    
    Defined at \
        ../src/fortran/lib/mod_intf_generator.f90 \
        lines 7-1373
    
    """
    @f90wrap.runtime.register_class("artemis.artemis_interface_generator")
    class artemis_interface_generator(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=artemis_interface_generator_type)
        
        
        Defined at \
            ../src/fortran/lib/mod_intf_generator.f90 \
            lines 30-75
        
        """
        def __init__(self, handle=None):
            """
            self = Artemis_Interface_Generator_Type()
            
            
            Defined at \
                ../src/fortran/lib/mod_intf_generator.f90 \
                lines 30-75
            
            
            Returns
            -------
            this : Artemis_Interface_Generator_Type
            	Object to be constructed
            
            
            Automatically generated constructor for artemis_interface_generator_type
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = \
                _artemis.f90wrap_intf_gen__artemis_interface_gen0ea8()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Artemis_Interface_Generator_Type
            
            
            Defined at \
                ../src/fortran/lib/mod_intf_generator.f90 \
                lines 30-75
            
            Parameters
            ----------
            this : Artemis_Interface_Generator_Type
            	Object to be destructed
            
            
            Automatically generated destructor for artemis_interface_generator_type
            """
            if self._alloc:
                _artemis.f90wrap_intf_gen__artemis_interface_genbc51(this=self._handle)
        
        def set_tolerance(self, vector_mismatch=None, angle_mismatch=None, \
            area_mismatch=None, max_length=None, max_area=None, max_fit=None, \
            max_extension=None, angle_weight=None, area_weight=None):
            """
            set_tolerance__binding__artemis_interface_generator_type(self[, vector_mismatch, \
                angle_mismatch, area_mismatch, max_length, max_area, max_fit, max_extension, \
                angle_weight, area_weight])
            
            
            Defined at \
                ../src/fortran/lib/mod_intf_generator.f90 \
                lines 85-125
            
            Parameters
            ----------
            this : Artemis_Interface_Generator_Type
            vector_mismatch : float
            angle_mismatch : float
            area_mismatch : float
            max_length : float
            max_area : float
            max_fit : int
            max_extension : int
            angle_weight : float
            area_weight : float
            
            """
            _artemis.f90wrap_intf_gen__set_tolerance__bindinfd58(this=self._handle, \
                vector_mismatch=vector_mismatch, angle_mismatch=angle_mismatch, \
                area_mismatch=area_mismatch, max_length=max_length, max_area=max_area, \
                max_fit=max_fit, max_extension=max_extension, angle_weight=angle_weight, \
                area_weight=area_weight)
        
        def set_shift_method(self, method=None, num_shifts=None, shifts=None, \
            interface_depth=None, separation_scale=None, depth_method=None):
            """
            set_shift_method__binding__artemis_interface_generator_type(self[, method, \
                num_shifts, shifts, interface_depth, separation_scale, depth_method])
            
            
            Defined at \
                ../src/fortran/lib/mod_intf_generator.f90 \
                lines 133-196
            
            Parameters
            ----------
            this : Artemis_Interface_Generator_Type
            method : int
            num_shifts : int
            shifts : float array
            interface_depth : float
            separation_scale : float
            depth_method : int
            
            """
            _artemis.f90wrap_intf_gen__set_shift_method__bin4dc1(this=self._handle, \
                method=method, num_shifts=num_shifts, shifts=shifts, \
                interface_depth=interface_depth, separation_scale=separation_scale, \
                depth_method=depth_method)
        
        def generate(self, basis_lw, basis_up, miller_lw=None, miller_up=None, \
            surface_lw=None, surface_up=None, thickness_lw=None, thickness_up=None, \
            num_layers_lw=None, num_layers_up=None, use_pricel_lw=None, \
            use_pricel_up=None, is_layered_lw=None, is_layered_up=None, \
            elastic_constants_lw=None, elastic_constants_up=None, \
            print_lattice_match_info=None, print_termination_info=None, \
            print_shift_info=None, break_on_fail=None, icheck_match=None, \
            interface_idx=None, generate_structures=None, seed=None, 
            calc=None):
            """
            generate__binding__artemis_interface_generator_type(self, basis_lw, basis_up[, \
                miller_lw, miller_up, surface_lw, surface_up, thickness_lw, thickness_up, \
                num_layers_lw, num_layers_up, use_pricel_lw, use_pricel_up, is_layered_lw, \
                is_layered_up, elastic_constants_lw, elastic_constants_up, \
                print_lattice_match_info, print_termination_info, print_shift_info, \
                break_on_fail, icheck_match, interface_idx, generate_structures, seed])
            
            
            Defined at \
                ../src/fortran/lib/mod_intf_generator.f90 \
                lines 315-1111
            
            Parameters
            ----------
            this : Artemis_Interface_Generator_Type
            basis_lw : Basis_Type
            basis_up : Basis_Type
            miller_lw : int array
            miller_up : int array
            surface_lw : int array
            surface_up : int array
            thickness_lw : float
            thickness_up : float
            num_layers_lw : int
            num_layers_up : int
            use_pricel_lw : bool
            use_pricel_up : bool
            is_layered_lw : bool
            is_layered_up : bool
            elastic_constants_lw : float array
            elastic_constants_up : float array
            print_lattice_match_info : bool
            print_termination_info : bool
            print_shift_info : bool
            break_on_fail : bool
            icheck_match : int
            interface_idx : int
            generate_structures : bool
            seed : int
            
            """

            exit_code = 0
            structures = None

            # check if host is ase.Atoms object or a Fortran derived type basis_type
            if isinstance(basis_lw, Atoms):
                basis_lw = geom_rw.basis(atoms=basis_lw)

            if isinstance(basis_up, Atoms):
                basis_up = geom_rw.basis(atoms=basis_up)

            exit_code = _artemis.f90wrap_intf_gen__generate__binding__aigt(this=self._handle, \
                basis_lw=basis_lw._handle, basis_up=basis_up._handle, miller_lw=miller_lw, \
                miller_up=miller_up, surface_lw=surface_lw, surface_up=surface_up, \
                thickness_lw=thickness_lw, thickness_up=thickness_up, \
                num_layers_lw=num_layers_lw, num_layers_up=num_layers_up, \
                use_pricel_lw=use_pricel_lw, use_pricel_up=use_pricel_up, \
                is_layered_lw=is_layered_lw, is_layered_up=is_layered_up, \
                elastic_constants_lw=elastic_constants_lw, \
                elastic_constants_up=elastic_constants_up, \
                print_lattice_match_info=print_lattice_match_info, \
                print_termination_info=print_termination_info, \
                print_shift_info=print_shift_info, break_on_fail=break_on_fail, \
                icheck_match=icheck_match, interface_idx=interface_idx, \
                generate_structures=generate_structures, seed=seed \
            )
        
            structures = self.get_structures(calc)
            return structures, exit_code

        def restart(self, basis, interface_location=None, print_shift_info=None, \
            seed=None):
            """
            restart__binding__artemis_interface_generator_type(self, basis[, \
                interface_location, print_shift_info, seed])
            
            
            Defined at \
                ../src/fortran/lib/mod_intf_generator.f90 \
                lines 202-297
            
            Parameters
            ----------
            this : Artemis_Interface_Generator_Type
            basis : Basis_Type
            interface_location : float array
            print_shift_info : bool
            seed : int
            
            ---------------------------------------------------------------------------
             Set the random seed
            ---------------------------------------------------------------------------
            """
            _artemis.f90wrap_intf_gen__restart__binding__aigt(this=self._handle, \
                basis=basis._handle, interface_location=interface_location, \
                print_shift_info=print_shift_info, seed=seed)
        
        def get_structures(self, calculator=None):
            """
            Get the generated structures as a list of ASE Atoms objects.

            Parameters:
                calculator (ASE calculator):
                    The calculator to use for the generated structures.
            """
            atoms = []
            for structure in self.structures:
                atoms.append(structure.toase(calculator))
            return atoms

        @property
        def num_structures(self):
            """
            The number of generated structures currently stored in the generator.
            """
            return _artemis.f90wrap_artemis_intf_gen_type__get__num_structures(self._handle)

        @num_structures.setter
        def num_structures(self, num_structures):
            _raffle.f90wrap_artemis_intf_gen_type__set__num_structures(self._handle, \
                num_structures)

        @property
        def max_num_structures(self):
            """
            The maximum number of generated structures that can be stored in the generator.
            """
            return _artemis.f90wrap_artemis_intf_gen_type__get__num_structures(self._handle)

        @max_num_structures.setter
        def max_num_structures(self, max_num_structures):
            _raffle.f90wrap_artemis_intf_gen_type__set__max_num_structures(self._handle, \
                max_num_structures)

        @property
        def shift_method(self):
            """
            Element shift_method ftype=integer  pytype=int
            
            
            Defined at \
                ../src/fortran/lib/mod_intf_generator.f90 \
                line 31
            
            """
            return \
                _artemis.f90wrap_artemis_intf_gen_type__get__shift_method(self._handle)
        
        @shift_method.setter
        def shift_method(self, shift_method):
            _artemis.f90wrap_artemis_intf_gen_type__set__shift_method(self._handle, \
                shift_method)
        
        @property
        def num_shifts(self):
            """
            Element num_shifts ftype=integer  pytype=int
            
            
            Defined at \
                ../src/fortran/lib/mod_intf_generator.f90 \
                line 33
            
            """
            return \
                _artemis.f90wrap_artemis_intf_gen_type__get__num_shifts(self._handle)
        
        @num_shifts.setter
        def num_shifts(self, num_shifts):
            _artemis.f90wrap_artemis_intf_gen_type__set__num_shifts(self._handle, \
                num_shifts)
        
        @property
        def shifts(self):
            """
            Element shifts ftype=real(real32) pytype=float
            
            
            Defined at \
                ../src/fortran/lib/mod_intf_generator.f90 \
                line 35
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _artemis.f90wrap_artemis_intf_gen_type__array__shifts(self._handle)
            if array_handle in self._arrays:
                shifts = self._arrays[array_handle]
            else:
                shifts = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _artemis.f90wrap_artemis_intf_gen_type__array__shifts)
                self._arrays[array_handle] = shifts
            return shifts
        
        @shifts.setter
        def shifts(self, shifts):
            self.shifts[...] = shifts
        
        @property
        def interface_depth(self):
            """
            Element interface_depth ftype=real(real32) pytype=float
            
            
            Defined at \
                ../src/fortran/lib/mod_intf_generator.f90 \
                line 37
            
            """
            return \
                _artemis.f90wrap_artemis_intf_gen_type__get__interface_depth(self._handle)
        
        @interface_depth.setter
        def interface_depth(self, interface_depth):
            _artemis.f90wrap_artemis_intf_gen_type__set__interface_depth(self._handle, \
                interface_depth)
        
        @property
        def separation_scale(self):
            """
            Element separation_scale ftype=real(real32) pytype=float
            
            
            Defined at \
                ../src/fortran/lib/mod_intf_generator.f90 \
                line 39
            
            """
            return \
                _artemis.f90wrap_artemis_intf_gen_type__get__separation_scale(self._handle)
        
        @separation_scale.setter
        def separation_scale(self, separation_scale):
            _artemis.f90wrap_artemis_intf_gen_type__set__separation_scale(self._handle, \
                separation_scale)
        
        @property
        def depth_method(self):
            """
            Element depth_method ftype=integer  pytype=int
            
            
            Defined at \
                ../src/fortran/lib/mod_intf_generator.f90 \
                line 41
            
            """
            return \
                _artemis.f90wrap_artemis_intf_gen_type__get__depth_method(self._handle)
        
        @depth_method.setter
        def depth_method(self, depth_method):
            _artemis.f90wrap_artemis_intf_gen_type__set__depth_method(self._handle, \
                depth_method)
        
        @property
        def shift_data(self):
            """
            Element shift_data ftype=real(real32) pytype=float
            
            
            Defined at \
                ../src/fortran/lib/mod_intf_generator.f90 \
                line 43
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _artemis.f90wrap_artemis_intf_gen_type__array__shift_data(self._handle)
            if array_handle in self._arrays:
                shift_data = self._arrays[array_handle]
            else:
                shift_data = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _artemis.f90wrap_artemis_intf_gen_type__array__shift_data)
                self._arrays[array_handle] = shift_data
            return shift_data
        
        @shift_data.setter
        def shift_data(self, shift_data):
            self.shift_data[...] = shift_data
        
        @property
        def swap_method(self):
            """
            Element swap_method ftype=integer  pytype=int
            
            
            Defined at \
                ../src/fortran/lib/mod_intf_generator.f90 \
                line 45
            
            """
            return \
                _artemis.f90wrap_artemis_intf_gen_type__get__swap_method(self._handle)
        
        @swap_method.setter
        def swap_method(self, swap_method):
            _artemis.f90wrap_artemis_intf_gen_type__set__swap_method(self._handle, \
                swap_method)
        
        @property
        def num_swaps(self):
            """
            Element num_swaps ftype=integer  pytype=int
            
            
            Defined at \
                ../src/fortran/lib/mod_intf_generator.f90 \
                line 47
            
            """
            return \
                _artemis.f90wrap_artemis_intf_gen_type__get__num_swaps(self._handle)
        
        @num_swaps.setter
        def num_swaps(self, num_swaps):
            _artemis.f90wrap_artemis_intf_gen_type__set__num_swaps(self._handle, \
                num_swaps)
        
        @property
        def swap_density(self):
            """
            Element swap_density ftype=real(real32) pytype=float
            
            
            Defined at \
                ../src/fortran/lib/mod_intf_generator.f90 \
                line 49
            
            """
            return \
                _artemis.f90wrap_artemis_intf_gen_type__get__swap_density(self._handle)
        
        @swap_density.setter
        def swap_density(self, swap_density):
            _artemis.f90wrap_artemis_intf_gen_type__set__swap_density(self._handle, \
                swap_density)
        
        @property
        def swap_depth(self):
            """
            Element swap_depth ftype=real(real32) pytype=float
            
            
            Defined at \
                ../src/fortran/lib/mod_intf_generator.f90 \
                line 51
            
            """
            return \
                _artemis.f90wrap_artemis_intf_gen_type__get__swap_depth(self._handle)
        
        @swap_depth.setter
        def swap_depth(self, swap_depth):
            _artemis.f90wrap_artemis_intf_gen_type__set__swap_depth(self._handle, \
                swap_depth)
        
        @property
        def swap_sigma(self):
            """
            Element swap_sigma ftype=real(real32) pytype=float
            
            
            Defined at \
                ../src/fortran/lib/mod_intf_generator.f90 \
                line 53
            
            """
            return \
                _artemis.f90wrap_artemis_intf_gen_type__get__swap_sigma(self._handle)
        
        @swap_sigma.setter
        def swap_sigma(self, swap_sigma):
            _artemis.f90wrap_artemis_intf_gen_type__set__swap_sigma(self._handle, \
                swap_sigma)
        
        @property
        def require_mirror_swaps(self):
            """
            Element require_mirror_swaps ftype=logical pytype=bool
            
            
            Defined at \
                ../src/fortran/lib/mod_intf_generator.f90 \
                line 55
            
            """
            return \
                _artemis.f90wrap_artemis_intf_gen_type__get__require_mirr41cf(self._handle)
        
        @require_mirror_swaps.setter
        def require_mirror_swaps(self, require_mirror_swaps):
            _artemis.f90wrap_artemis_intf_gen_type__set__require_mirr3bfa(self._handle, \
                require_mirror_swaps)
        
        @property
        def match_method(self):
            """
            Element match_method ftype=integer  pytype=int
            
            
            Defined at \
                ../src/fortran/lib/mod_intf_generator.f90 \
                line 57
            
            """
            return \
                _artemis.f90wrap_artemis_intf_gen_type__get__match_method(self._handle)
        
        @match_method.setter
        def match_method(self, match_method):
            _artemis.f90wrap_artemis_intf_gen_type__set__match_method(self._handle, \
                match_method)
        
        @property
        def max_num_matches(self):
            """
            Element max_num_matches ftype=integer  pytype=int
            
            
            Defined at \
                ../src/fortran/lib/mod_intf_generator.f90 \
                line 58
            
            """
            return \
                _artemis.f90wrap_artemis_intf_gen_type__get__max_num_matches(self._handle)
        
        @max_num_matches.setter
        def max_num_matches(self, max_num_matches):
            _artemis.f90wrap_artemis_intf_gen_type__set__max_num_matches(self._handle, \
                max_num_matches)
        
        @property
        def max_num_terms(self):
            """
            Element max_num_terms ftype=integer  pytype=int
            
            
            Defined at \
                ../src/fortran/lib/mod_intf_generator.f90 \
                line 59
            
            """
            return \
                _artemis.f90wrap_artemis_intf_gen_type__get__max_num_terms(self._handle)
        
        @max_num_terms.setter
        def max_num_terms(self, max_num_terms):
            _artemis.f90wrap_artemis_intf_gen_type__set__max_num_terms(self._handle, \
                max_num_terms)
        
        @property
        def max_num_planes(self):
            """
            Element max_num_planes ftype=integer  pytype=int
            
            
            Defined at \
                ../src/fortran/lib/mod_intf_generator.f90 \
                line 60
            
            """
            return \
                _artemis.f90wrap_artemis_intf_gen_type__get__max_num_planes(self._handle)
        
        @max_num_planes.setter
        def max_num_planes(self, max_num_planes):
            _artemis.f90wrap_artemis_intf_gen_type__set__max_num_planes(self._handle, \
                max_num_planes)
        
        @property
        def fix_normal(self):
            """
            Element fix_normal ftype=logical pytype=bool
            
            
            Defined at \
                ../src/fortran/lib/mod_intf_generator.f90 \
                line 61
            
            """
            return \
                _artemis.f90wrap_artemis_intf_gen_type__get__fix_normal(self._handle)
        
        @fix_normal.setter
        def fix_normal(self, fix_normal):
            _artemis.f90wrap_artemis_intf_gen_type__set__fix_normal(self._handle, \
                fix_normal)
        
        @property
        def bondlength_cutoff(self):
            """
            Element bondlength_cutoff ftype=real(real32) pytype=float
            
            
            Defined at \
                ../src/fortran/lib/mod_intf_generator.f90 \
                line 65
            
            """
            return \
                _artemis.f90wrap_artemis_intf_gen_type__get__bondlength_c21a8(self._handle)
        
        @bondlength_cutoff.setter
        def bondlength_cutoff(self, bondlength_cutoff):
            _artemis.f90wrap_artemis_intf_gen_type__set__bondlength_cbd11(self._handle, \
                bondlength_cutoff)
        
        @property
        def layer_separation_cutoff(self):
            """
            Element layer_separation_cutoff ftype=real(real32) pytype=float
            
            
            Defined at \
                ../src/fortran/lib/mod_intf_generator.f90 \
                line 66
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _artemis.f90wrap_artemis_intf_gen_type__array__layer_sepa90a5(self._handle)
            if array_handle in self._arrays:
                layer_separation_cutoff = self._arrays[array_handle]
            else:
                layer_separation_cutoff = \
                    f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _artemis.f90wrap_artemis_intf_gen_type__array__layer_sepa90a5)
                self._arrays[array_handle] = layer_separation_cutoff
            return layer_separation_cutoff
        
        @layer_separation_cutoff.setter
        def layer_separation_cutoff(self, layer_separation_cutoff):
            self.layer_separation_cutoff[...] = layer_separation_cutoff
        
        def _init_array_structures(self):
            """
            Initialise the structures array.

            It is not recommended to use this function directly. Use the `structures` property instead.
            """
            self.structures = f90wrap.runtime.FortranDerivedTypeArray(self,
                                            _artemis.f90wrap_artemis_intf_gen_type__array_getitem__structures,
                                            _artemis.f90wrap_artemis_intf_gen_type__array_setitem__structures,
                                            _artemis.f90wrap_artemis_intf_gen_type__array_len__structures,
                                            """
            Element items ftype=type(basis_type) pytype=basis


            Defined at ../src/lib/mod_generator.f90 line \
                29

            """, Geom_Rw.basis)
            return self.structures

        def __str__(self):
            ret = ['<artemis_interface_generator_type>{\n']
            ret.append('    shift_method : ')
            ret.append(repr(self.shift_method))
            ret.append(',\n    num_shifts : ')
            ret.append(repr(self.num_shifts))
            ret.append(',\n    shifts : ')
            ret.append(repr(self.shifts))
            ret.append(',\n    interface_depth : ')
            ret.append(repr(self.interface_depth))
            ret.append(',\n    separation_scale : ')
            ret.append(repr(self.separation_scale))
            ret.append(',\n    depth_method : ')
            ret.append(repr(self.depth_method))
            ret.append(',\n    shift_data : ')
            ret.append(repr(self.shift_data))
            ret.append(',\n    swap_method : ')
            ret.append(repr(self.swap_method))
            ret.append(',\n    num_swaps : ')
            ret.append(repr(self.num_swaps))
            ret.append(',\n    swap_density : ')
            ret.append(repr(self.swap_density))
            ret.append(',\n    swap_depth : ')
            ret.append(repr(self.swap_depth))
            ret.append(',\n    swap_sigma : ')
            ret.append(repr(self.swap_sigma))
            ret.append(',\n    require_mirror_swaps : ')
            ret.append(repr(self.require_mirror_swaps))
            ret.append(',\n    match_method : ')
            ret.append(repr(self.match_method))
            ret.append(',\n    max_num_matches : ')
            ret.append(repr(self.max_num_matches))
            ret.append(',\n    max_num_terms : ')
            ret.append(repr(self.max_num_terms))
            ret.append(',\n    max_num_planes : ')
            ret.append(repr(self.max_num_planes))
            ret.append(',\n    fix_normal : ')
            ret.append(repr(self.fix_normal))
            ret.append(',\n    bondlength_cutoff : ')
            ret.append(repr(self.bondlength_cutoff))
            ret.append(',\n    layer_separation_cutoff : ')
            ret.append(repr(self.layer_separation_cutoff))
            ret.append(',\n    structures : ')
            ret.append(repr(self.structures))
            ret.append('}')
            return ''.join(ret)
        
        _dt_array_initialisers = [_init_array_structures]
        
    
    _dt_array_initialisers = []
    

interface_generator = Interface_Generator()

class Artemis(f90wrap.runtime.FortranModule):
    """
    Module artemis
    
    
    Defined at ../src/fortran/artemis.f90 lines \
        1-4
    
    """
    pass
    _dt_array_initialisers = []
    

artemis = Artemis()

