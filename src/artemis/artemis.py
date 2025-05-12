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
                    positions.append(self.spec[i].atom[j][:3])

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

        def toase(self, calculator=None):
            """
            Convert the basis_array object to a list of ASE Atoms objects.
            """

            # Set the species list
            atoms = []
            for i in range(len(self.items)):
                atoms.append(self.items[i].toase(calculator=calculator))
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




class Misc_Types(f90wrap.runtime.FortranModule):
    """
    Module artemis__misc_types
    
    
    Defined at \
        ../fortran/lib/mod_misc_types.f90 \
        lines 1-261
    
    """
    @f90wrap.runtime.register_class("artemis.struc_data_type")
    class struc_data_type(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=struc_data_type)
        
        
        Defined at \
            ../fortran/lib/mod_misc_types.f90 \
            lines 24-42
        
        """
        def __init__(self, handle=None):
            """
            self = Struc_Data_Type()
            
            
            Defined at \
                ../fortran/lib/mod_misc_types.f90 \
                lines 24-42
            
            
            Returns
            -------
            this : Struc_Data_Type
            	Object to be constructed
            
            
            Automatically generated constructor for struc_data_type
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _artemis.f90wrap_misc_types__struc_data_type_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Struc_Data_Type
            
            
            Defined at \
                ../fortran/lib/mod_misc_types.f90 \
                lines 24-42
            
            Parameters
            ----------
            this : Struc_Data_Type
            	Object to be destructed
            
            
            Automatically generated destructor for struc_data_type
            """
            if self._alloc:
                _artemis.f90wrap_misc_types__struc_data_type_finalise(this=self._handle)
        
        @property
        def match_idx(self):
            """
            Element match_idx ftype=integer  pytype=int
            
            
            Defined at \
                ../fortran/lib/mod_misc_types.f90 \
                line 25
            
            """
            return _artemis.f90wrap_struc_data_type__get__match_idx(self._handle)
        
        @match_idx.setter
        def match_idx(self, match_idx):
            _artemis.f90wrap_struc_data_type__set__match_idx(self._handle, match_idx)
        
        @property
        def shift_idx(self):
            """
            Element shift_idx ftype=integer  pytype=int
            
            
            Defined at \
                ../fortran/lib/mod_misc_types.f90 \
                line 26
            
            """
            return _artemis.f90wrap_struc_data_type__get__shift_idx(self._handle)
        
        @shift_idx.setter
        def shift_idx(self, shift_idx):
            _artemis.f90wrap_struc_data_type__set__shift_idx(self._handle, shift_idx)
        
        @property
        def swap_idx(self):
            """
            Element swap_idx ftype=integer  pytype=int
            
            
            Defined at \
                ../fortran/lib/mod_misc_types.f90 \
                line 27
            
            """
            return _artemis.f90wrap_struc_data_type__get__swap_idx(self._handle)
        
        @swap_idx.setter
        def swap_idx(self, swap_idx):
            _artemis.f90wrap_struc_data_type__set__swap_idx(self._handle, swap_idx)
        
        @property
        def from_pricel_lw(self):
            """
            Element from_pricel_lw ftype=logical pytype=bool
            
            
            Defined at \
                ../fortran/lib/mod_misc_types.f90 \
                line 28
            
            """
            return _artemis.f90wrap_struc_data_type__get__from_pricel_lw(self._handle)
        
        @from_pricel_lw.setter
        def from_pricel_lw(self, from_pricel_lw):
            _artemis.f90wrap_struc_data_type__set__from_pricel_lw(self._handle, \
                from_pricel_lw)
        
        @property
        def from_pricel_up(self):
            """
            Element from_pricel_up ftype=logical pytype=bool
            
            
            Defined at \
                ../fortran/lib/mod_misc_types.f90 \
                line 29
            
            """
            return _artemis.f90wrap_struc_data_type__get__from_pricel_up(self._handle)
        
        @from_pricel_up.setter
        def from_pricel_up(self, from_pricel_up):
            _artemis.f90wrap_struc_data_type__set__from_pricel_up(self._handle, \
                from_pricel_up)
        
        @property
        def term_lw_idx(self):
            """
            Element term_lw_idx ftype=integer pytype=int
            
            
            Defined at \
                ../fortran/lib/mod_misc_types.f90 \
                line 30
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _artemis.f90wrap_struc_data_type__array__term_lw_idx(self._handle)
            if array_handle in self._arrays:
                term_lw_idx = self._arrays[array_handle]
            else:
                term_lw_idx = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _artemis.f90wrap_struc_data_type__array__term_lw_idx)
                self._arrays[array_handle] = term_lw_idx
            return term_lw_idx
        
        @term_lw_idx.setter
        def term_lw_idx(self, term_lw_idx):
            self.term_lw_idx[...] = term_lw_idx
        
        @property
        def term_up_idx(self):
            """
            Element term_up_idx ftype=integer pytype=int
            
            
            Defined at \
                ../fortran/lib/mod_misc_types.f90 \
                line 31
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _artemis.f90wrap_struc_data_type__array__term_up_idx(self._handle)
            if array_handle in self._arrays:
                term_up_idx = self._arrays[array_handle]
            else:
                term_up_idx = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _artemis.f90wrap_struc_data_type__array__term_up_idx)
                self._arrays[array_handle] = term_up_idx
            return term_up_idx
        
        @term_up_idx.setter
        def term_up_idx(self, term_up_idx):
            self.term_up_idx[...] = term_up_idx
        
        @property
        def transform_lw(self):
            """
            Element transform_lw ftype=integer pytype=int
            
            
            Defined at \
                ../fortran/lib/mod_misc_types.f90 \
                line 32
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _artemis.f90wrap_struc_data_type__array__transform_lw(self._handle)
            if array_handle in self._arrays:
                transform_lw = self._arrays[array_handle]
            else:
                transform_lw = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _artemis.f90wrap_struc_data_type__array__transform_lw)
                self._arrays[array_handle] = transform_lw
            return transform_lw
        
        @transform_lw.setter
        def transform_lw(self, transform_lw):
            self.transform_lw[...] = transform_lw
        
        @property
        def transform_up(self):
            """
            Element transform_up ftype=integer pytype=int
            
            
            Defined at \
                ../fortran/lib/mod_misc_types.f90 \
                line 33
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _artemis.f90wrap_struc_data_type__array__transform_up(self._handle)
            if array_handle in self._arrays:
                transform_up = self._arrays[array_handle]
            else:
                transform_up = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _artemis.f90wrap_struc_data_type__array__transform_up)
                self._arrays[array_handle] = transform_up
            return transform_up
        
        @transform_up.setter
        def transform_up(self, transform_up):
            self.transform_up[...] = transform_up
        
        @property
        def approx_thickness_lw(self):
            """
            Element approx_thickness_lw ftype=real(real32) pytype=float
            
            
            Defined at \
                ../fortran/lib/mod_misc_types.f90 \
                line 34
            
            """
            return _artemis.f90wrap_struc_data_type__get__approx_thickness_lw(self._handle)
        
        @approx_thickness_lw.setter
        def approx_thickness_lw(self, approx_thickness_lw):
            _artemis.f90wrap_struc_data_type__set__approx_thickness_lw(self._handle, \
                approx_thickness_lw)
        
        @property
        def approx_thickness_up(self):
            """
            Element approx_thickness_up ftype=real(real32) pytype=float
            
            
            Defined at \
                ../fortran/lib/mod_misc_types.f90 \
                line 35
            
            """
            return _artemis.f90wrap_struc_data_type__get__approx_thickness_up(self._handle)
        
        @approx_thickness_up.setter
        def approx_thickness_up(self, approx_thickness_up):
            _artemis.f90wrap_struc_data_type__set__approx_thickness_up(self._handle, \
                approx_thickness_up)
        
        @property
        def mismatch(self):
            """
            Element mismatch ftype=real(real32) pytype=float
            
            
            Defined at \
                ../fortran/lib/mod_misc_types.f90 \
                line 36
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _artemis.f90wrap_struc_data_type__array__mismatch(self._handle)
            if array_handle in self._arrays:
                mismatch = self._arrays[array_handle]
            else:
                mismatch = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _artemis.f90wrap_struc_data_type__array__mismatch)
                self._arrays[array_handle] = mismatch
            return mismatch
        
        @mismatch.setter
        def mismatch(self, mismatch):
            self.mismatch[...] = mismatch
        
        @property
        def shift(self):
            """
            Element shift ftype=real(real32) pytype=float
            
            
            Defined at \
                ../fortran/lib/mod_misc_types.f90 \
                line 37
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _artemis.f90wrap_struc_data_type__array__shift(self._handle)
            if array_handle in self._arrays:
                shift = self._arrays[array_handle]
            else:
                shift = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _artemis.f90wrap_struc_data_type__array__shift)
                self._arrays[array_handle] = shift
            return shift
        
        @shift.setter
        def shift(self, shift):
            self.shift[...] = shift
        
        @property
        def swap_density(self):
            """
            Element swap_density ftype=real(real32) pytype=float
            
            
            Defined at \
                ../fortran/lib/mod_misc_types.f90 \
                line 39
            
            """
            return _artemis.f90wrap_struc_data_type__get__swap_density(self._handle)
        
        @swap_density.setter
        def swap_density(self, swap_density):
            _artemis.f90wrap_struc_data_type__set__swap_density(self._handle, swap_density)
        
        @property
        def approx_eff_swap_conc(self):
            """
            Element approx_eff_swap_conc ftype=real(real32) pytype=float
            
            
            Defined at \
                ../fortran/lib/mod_misc_types.f90 \
                line 40
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _artemis.f90wrap_struc_data_type__array__approx_eff_swap_conc(self._handle)
            if array_handle in self._arrays:
                approx_eff_swap_conc = self._arrays[array_handle]
            else:
                approx_eff_swap_conc = \
                    f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _artemis.f90wrap_struc_data_type__array__approx_eff_swap_conc)
                self._arrays[array_handle] = approx_eff_swap_conc
            return approx_eff_swap_conc
        
        @approx_eff_swap_conc.setter
        def approx_eff_swap_conc(self, approx_eff_swap_conc):
            self.approx_eff_swap_conc[...] = approx_eff_swap_conc
        
        def __str__(self):
            ret = ['<struc_data_type>{\n']
            ret.append('    match_idx : ')
            ret.append(repr(self.match_idx))
            ret.append(',\n    shift_idx : ')
            ret.append(repr(self.shift_idx))
            ret.append(',\n    swap_idx : ')
            ret.append(repr(self.swap_idx))
            ret.append(',\n    from_pricel_lw : ')
            ret.append(repr(self.from_pricel_lw))
            ret.append(',\n    from_pricel_up : ')
            ret.append(repr(self.from_pricel_up))
            ret.append(',\n    term_lw_idx : ')
            ret.append(repr(self.term_lw_idx))
            ret.append(',\n    term_up_idx : ')
            ret.append(repr(self.term_up_idx))
            ret.append(',\n    transform_lw : ')
            ret.append(repr(self.transform_lw))
            ret.append(',\n    transform_up : ')
            ret.append(repr(self.transform_up))
            ret.append(',\n    approx_thickness_lw : ')
            ret.append(repr(self.approx_thickness_lw))
            ret.append(',\n    approx_thickness_up : ')
            ret.append(repr(self.approx_thickness_up))
            ret.append(',\n    mismatch : ')
            ret.append(repr(self.mismatch))
            ret.append(',\n    shift : ')
            ret.append(repr(self.shift))
            ret.append(',\n    swap_density : ')
            ret.append(repr(self.swap_density))
            ret.append(',\n    approx_eff_swap_conc : ')
            ret.append(repr(self.approx_eff_swap_conc))
            ret.append('}')
            return ''.join(ret)
        
        _dt_array_initialisers = []


misc_types = Misc_Types()




class Generator(f90wrap.runtime.FortranModule):
    """
    Module artemis__generator
    
    
    Defined at \
        ../src/fortran/lib/mod_intf_generator.f90 \
        lines 7-1373
    
    """
    @f90wrap.runtime.register_class("artemis.artemis_generator")
    class artemis_generator(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=artemis_generator_type)
        
        
        Defined at \
            ../src/fortran/lib/mod_intf_generator.f90 \
            lines 30-75
        
        """
        def __init__(self, handle=None):
            """
            self = Artemis_generator_Type()
            
            
            Defined at \
                ../src/fortran/lib/mod_intf_generator.f90 \
                lines 30-75
            
            
            Returns
            -------
            this : Artemis_generator_Type
            	Object to be constructed
            
            
            Automatically generated constructor for artemis_generator_type
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = \
                _artemis.f90wrap_intf_gen__artemis_gen_type_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Artemis_generator_Type
            
            
            Defined at \
                ../src/fortran/lib/mod_intf_generator.f90 \
                lines 30-75
            
            Parameters
            ----------
            this : Artemis_generator_Type
            	Object to be destructed
            
            
            Automatically generated destructor for artemis_generator_type
            """
            if self._alloc:
                _artemis.f90wrap_intf_gen__artemis_gen_type_finalise(this=self._handle)
        
        def get_all_structures_data(self):
            """
            output = get_all_structure_data__binding__artemis_generator_type(self)
            
            
            Defined at \
                ../fortran/lib/mod_generator.f90 \
                lines 134-146
            
            Parameters
            ----------
            this : Artemis_Generator_Type
            
            Returns
            -------
            output : Struc_Data_Type array
            
            """
            output = []
            for i in range(self.num_structures):
                output.append(self.get_structure_data(i))

            # output = \
            #     _artemis.f90wrap_intf_gen__get_all_structures_data__binding_agt(this=self._handle)
            # output = \
            #     f90wrap.runtime.lookup_class("artemis.struc_data_type").from_handle(output, \
            #     alloc=True)
            return output
        
        def get_structure_data(self, idx):
            """
            output = get_structure_data__binding__artemis_generator_type(self, idx)
            
            
            Defined at \
                ../fortran/lib/mod_generator.f90 \
                lines 150-160
            
            Parameters
            ----------
            this : Artemis_Generator_Type
            idx : int
            
            Returns
            -------
            output : Struc_Data_Type
            
            """
            output = \
                _artemis.f90wrap_intf_gen__get_structure_data__binding_agt(this=self._handle, \
                idx=idx)
            output = \
                f90wrap.runtime.lookup_class("artemis.struc_data_type").from_handle(output, \
                alloc=True)
            
            output_dict = {
                'match_idx': output.match_idx,
                'shift_idx': output.shift_idx,
                'swap_idx': output.swap_idx,
                'from_pricel_lw': output.from_pricel_lw,
                'from_pricel_up': output.from_pricel_up,
                'term_lw_idx': output.term_lw_idx,
                'term_up_idx': output.term_up_idx,
                'transform_lw': output.transform_lw,
                'transform_up': output.transform_up,
                'approx_thickness_lw': output.approx_thickness_lw,
                'approx_thickness_up': output.approx_thickness_up,
                'mismatch': output.mismatch,
                'shift': output.shift,
                'swap_density': output.swap_density,
                'approx_eff_swap_conc': output.approx_eff_swap_conc
            }

            return output_dict
        
        def get_all_structures_mismatch(self):
            """
            output = get_all_structures_mismatch__binding__artemis_generator_type(self)
            
            
            Defined at \
                ../fortran/lib/mod_generator.f90 \
                lines 130-142
            
            Parameters
            ----------
            this : Artemis_Generator_Type
            
            Returns
            -------
            output : float array
            
            """
            output = \
                _artemis.f90wrap_intf_gen__get_all_structures_mismatch__binding_agt(this=self._handle)
            return output
        
        def get_structure_mismatch(self, idx):
            """
            output = get_structure_mismatch__binding__artemis_generator_type(self, idx)
            
            
            Defined at \
                ../fortran/lib/mod_generator.f90 \
                lines 146-156
            
            Parameters
            ----------
            this : Artemis_Generator_Type
            idx : int
            
            Returns
            -------
            output : float array
            
            """
            output = \
                _artemis.f90wrap_intf_gen__get_structure_mismatch__binding_agt(this=self._handle, \
                idx=idx)
            return output
        
        def get_all_structures_transform(self):
            """
            output = get_all_structures_transform__binding__artemis_generator_type(self)
            
            
            Defined at \
                ../fortran/lib/mod_generator.f90 \
                lines 160-174
            
            Parameters
            ----------
            this : Artemis_Generator_Type
            
            Returns
            -------
            output : int array
            
            """
            output = \
                _artemis.f90wrap_intf_gen__get_all_structures_transform__binding_agt(this=self._handle)
            return output
        
        def get_structure_transform(self, idx):
            """
            output = get_structure_transform__binding__artemis_generator_type(self, idx)
            
            
            Defined at \
                ../fortran/lib/mod_generator.f90 \
                lines 178-189
            
            Parameters
            ----------
            this : Artemis_Generator_Type
            idx : int
            
            Returns
            -------
            output : int array
            
            """
            output = \
                _artemis.f90wrap_intf_gen__get_structure_transform__binding_agt(this=self._handle, \
                idx=idx)
            return output
        
        def get_all_structures_shift(self):
            """
            output = get_all_structures_shifts__binding__artemis_generator_type(self)
            
            
            Defined at \
                ../fortran/lib/mod_generator.f90 \
                lines 193-205
            
            Parameters
            ----------
            this : Artemis_Generator_Type
            
            Returns
            -------
            output : float array
            
            """
            output = \
                _artemis.f90wrap_intf_gen__get_all_structures_shift__binding_agt(this=self._handle)
            return output
        
        def get_structure_shift(self, idx):
            """
            output = get_structure_shifts__binding__artemis_generator_type(self, idx)
            
            
            Defined at \
                ../fortran/lib/mod_generator.f90 \
                lines 209-219
            
            Parameters
            ----------
            this : Artemis_Generator_Type
            idx : int
            
            Returns
            -------
            output : float array
            
            """
            output = \
                _artemis.f90wrap_intf_gen__get_structure_shift__binding_agt(this=self._handle, \
                idx=idx)
            return output
        
        def set_tolerance(self, vector_mismatch=None, angle_mismatch=None, \
            area_mismatch=None, max_length=None, max_area=None, max_fit=None, \
            max_extension=None, angle_weight=None, area_weight=None):
            """
            set_tolerance__binding__artemis_gen_type(self[, vector_mismatch, \
                angle_mismatch, area_mismatch, max_length, max_area, max_fit, max_extension, \
                angle_weight, area_weight])
            
            
            Defined at \
                ../src/fortran/lib/mod_intf_generator.f90 \
                lines 85-125
            
            Parameters
            ----------
            this : Artemis_generator_Type
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
            _artemis.f90wrap_intf_gen__set_tolerance__bindind_agt(this=self._handle, \
                vector_mismatch=vector_mismatch, angle_mismatch=angle_mismatch, \
                area_mismatch=area_mismatch, max_length=max_length, max_area=max_area, \
                max_fit=max_fit, max_extension=max_extension, angle_weight=angle_weight, \
                area_weight=area_weight)
        
        def set_shift_method(self, method: int =None, num_shifts: int =None, shifts=None, \
            interface_depth=None, separation_scale=None, depth_method=None, \
            bondlength_cutoff=None):
            """
            set_shift_method__binding__artemis_generator_type(self[, method, num_shifts, \
                shifts, interface_depth, separation_scale, depth_method, bondlength_cutoff])
            
            
            Defined at \
                ../fortran/lib/mod_generator.f90 \
                lines 180-252
            
            Parameters
            ----------
            this : Artemis_Generator_Type
            method : int
            num_shifts : int
            shifts : float array
            interface_depth : float
            separation_scale : float
            depth_method : int
            bondlength_cutoff : float
            
            """
            
            if shifts is not None:
                # if shifts is a scalar, convert it to a 2D array, fortran order
                if isinstance(shifts, float) or isinstance(shifts, int):
                    shifts = numpy.array([[shifts]], order='F')
                # if shifts is a 1D array, convert it to a 2D array, fortran order
                elif isinstance(shifts, list):
                    shifts = numpy.array([shifts], order='F')
                elif len(shifts.shape) == 1:
                    shifts = numpy.array([shifts], order='F')
                # if shifts is a 2D array, convert it to a 2D array, fortran order
                elif len(shifts.shape) == 2:
                    shifts = numpy.array(shifts, order='F')

            _artemis.f90wrap_intf_gen__set_shift_method__binding__agt(this=self._handle, \
                method=method, num_shifts=num_shifts, shifts=shifts, \
                interface_depth=interface_depth, separation_scale=separation_scale, \
                depth_method=depth_method, bondlength_cutoff=bondlength_cutoff)
        
        def set_swap_method(self, method=None, num_swaps=None, swap_density=None, \
            swap_depth=None, swap_sigma=None, require_mirror_swaps=None):
            """
            set_swap_method__binding__artemis_generator_type(self[, method, num_swaps, \
                swap_density, swap_depth, swap_sigma, require_mirror_swaps])
            
            
            Defined at \
                ../fortran/lib/mod_generator.f90 \
                lines 259-283
            
            Parameters
            ----------
            this : Artemis_Generator_Type
            method : int
            num_swaps : int
            swap_density : float
            swap_depth : float
            swap_sigma : float
            require_mirror_swaps : bool
            
            """
            _artemis.f90wrap_intf_gen__set_swap_method__binding__agt(this=self._handle, \
                method=method, num_swaps=num_swaps, swap_density=swap_density, \
                swap_depth=swap_depth, swap_sigma=swap_sigma, \
                require_mirror_swaps=require_mirror_swaps)
        
        def set_match_method(self, method=None, max_num_matches=None, \
            max_num_terms=None, max_num_planes=None, compensate_normal=None):
            """
            set_match_method__binding__artemis_generator_type(self[, method, \
                max_num_matches, max_num_terms, max_num_planes, compensate_normal])
            
            
            Defined at \
                ../fortran/lib/mod_generator.f90 \
                lines 290-310
            
            Parameters
            ----------
            this : Artemis_Generator_Type
            method : int
            max_num_matches : int
            max_num_terms : int
            max_num_planes : int
            compensate_normal : bool
            
            """
            _artemis.f90wrap_intf_gen__set_match_method__binding__agt(this=self._handle, \
                method=method, max_num_matches=max_num_matches, max_num_terms=max_num_terms, \
                max_num_planes=max_num_planes, compensate_normal=compensate_normal)
        
        def set_materials(self,
                          structure_lw: Atoms | Geom_Rw.basis = None,
                          structure_up: Atoms | Geom_Rw.basis = None,
                          elastic_constants_lw=None,
                          elastic_constants_up=None,
                          use_pricel_lw=None,
                          use_pricel_up=None
        ):
            """
            set_materials__binding__artemis_gen_type(self, structure_lw, \
                structure_up[, elastic_constants_lw, elastic_constants_up, use_pricel_lw, \
                use_pricel_up])
            
            
            Defined at \
                ../fortran/lib/mod_intf_generator.f90 \
                lines 252-287
            
            Parameters
            ----------
            this : Artemis_generator_Type
            structure_lw : Basis_Type
            structure_up : Basis_Type
            elastic_constants_lw : float array
            elastic_constants_up : float array
            use_pricel_lw : bool
            use_pricel_up : bool
            
            ---------------------------------------------------------------------------
             Handle the elastic constants
            ---------------------------------------------------------------------------
            """

            # check if host is ase.Atoms object or a Fortran derived type basis_type
            if structure_lw is None:
                structure_lw_handle = None
            else:
                if isinstance(structure_lw, Atoms):
                    structure_lw = geom_rw.basis(atoms=structure_lw)
                structure_lw_handle = structure_lw._handle

            if structure_up is None:
                structure_up_handle = None
            else:
                if isinstance(structure_up, Atoms):
                    structure_up = geom_rw.basis(atoms=structure_up)
                structure_up_handle = structure_up._handle

            _artemis.f90wrap_intf_gen__set_materials__binding__agt(this=self._handle, \
                structure_lw=structure_lw_handle, structure_up=structure_up_handle, \
                elastic_constants_lw=elastic_constants_lw, \
                elastic_constants_up=elastic_constants_up, use_pricel_lw=use_pricel_lw, \
                use_pricel_up=use_pricel_up)
        
        def set_surface_properties(self, miller_lw=None, miller_up=None, \
            is_layered_lw=None, is_layered_up=None, \
            require_stoichiometry_lw=None, require_stoichiometry_up=None, \
            layer_separation_cutoff_lw=None, \
            layer_separation_cutoff_up=None, layer_separation_cutoff=None, \
            vacuum_gap=None):
            """
            set_surface_properties__binding__artemis_generator_type(self[, miller_lw, \
                miller_up, is_layered_lw, is_layered_up, layer_separation_cutoff_lw, \
                layer_separation_cutoff_up, layer_separation_cutoff, vacuum_gap])
            
            
            Defined at \
                ../fortran/lib/mod_generator.f90 \
                lines 364-435
            
            Parameters
            ----------
            this : Artemis_Generator_Type
            miller_lw : int array
            miller_up : int array
            is_layered_lw : bool
            is_layered_up : bool
            require_stoichiometry_lw : bool
            require_stoichiometry_up : bool
            layer_separation_cutoff_lw : float
            layer_separation_cutoff_up : float
            layer_separation_cutoff : float array
            vacuum_gap : float
            
            """
            _artemis.f90wrap_intf_gen__set_surface_properties__binding__agt(this=self._handle, \
                miller_lw=miller_lw, miller_up=miller_up, is_layered_lw=is_layered_lw, \
                is_layered_up=is_layered_up, \
                require_stoichiometry_lw=require_stoichiometry_lw, \
                require_stoichiometry_up=require_stoichiometry_up, \
                layer_separation_cutoff_lw=layer_separation_cutoff_lw, \
                layer_separation_cutoff_up=layer_separation_cutoff_up, \
                layer_separation_cutoff=layer_separation_cutoff, vacuum_gap=vacuum_gap)
        
        def reset_is_layered_lw(self):
            """
            reset_is_layered_lw__binding__artemis_gen_type(self)
            
            
            Defined at \
                ../fortran/lib/mod_intf_generator.f90 \
                lines 322-329
            
            Parameters
            ----------
            this : Artemis_generator_Type
            
            """
            _artemis.f90wrap_intf_gen__reset_is_layered_lw__binding__agt(this=self._handle)
        
        def reset_is_layered_up(self):
            """
            reset_is_layered_up__binding__artemis_gen_type(self)
            
            
            Defined at \
                ../fortran/lib/mod_intf_generator.f90 \
                lines 333-340
            
            Parameters
            ----------
            this : Artemis_generator_Type
            
            """
            _artemis.f90wrap_intf_gen__reset_is_layered_up__binding__agt(this=self._handle)
        
        def get_terminations_lw(self, miller=None, surface=None, num_layers=None, \
            thickness=None, orthogonalise=None, normalise=None, break_on_fail=None, 
            verbose=None, return_exit_code=False, calc=None):
            """
            
            Defined at \
                ../src/fortran/lib/mod_intf_generator.f90 \
            
            Parameters
            ----------
            """
            exit_code = 0
            structures = None

            exit_code, n_structs = _artemis.f90wrap_intf_gen__get_terminations__binding__agt(this=self._handle,
                identifier=1,
                miller=miller, surface=surface,
                num_layers=num_layers, thickness=thickness,
                orthogonalise=orthogonalise, normalise=normalise,
                break_on_fail=break_on_fail,
                verbose=verbose)
            if ( exit_code != 0 and exit_code != None ) and not return_exit_code:
                raise RuntimeError(f"Termination generation failed (exit code {exit_code})")

            # allocate the structures
            structures = geom_rw.basis_array() #.allocate(n_structs)
            structures.allocate(n_structs)
            _artemis.f90wrap_retrieve_last_generated_structures(structures._handle)
            structures = structures.toase()

            if return_exit_code:
                return structures, exit_code
            return structures
        
        def get_terminations_up(self, miller=None, surface=None, num_layers=None, \
            thickness=None, orthogonalise=None, normalise=None, break_on_fail=None, 
            verbose=None, return_exit_code=False, calc=None):
            """
            
            Defined at \
                ../src/fortran/lib/mod_intf_generator.f90 \
            
            Parameters
            ----------
            """
            exit_code = 0
            structures = None

            exit_code, n_structs = _artemis.f90wrap_intf_gen__get_terminations__binding__agt(this=self._handle,
                identifier=2,
                miller=miller, surface=surface,
                num_layers=num_layers, thickness=thickness,
                orthogonalise=orthogonalise, normalise=normalise,
                break_on_fail=break_on_fail,
                verbose=verbose)
            if ( exit_code != 0 and exit_code != None ) and not return_exit_code:
                raise RuntimeError(f"Termination generation failed (exit code {exit_code})")

            # allocate the structures
            structures = geom_rw.basis_array() #.allocate(n_structs)
            structures.allocate(n_structs)
            _artemis.f90wrap_retrieve_last_generated_structures(structures._handle)
            structures = structures.toase()

            if return_exit_code:
                return structures, exit_code
            return structures

        def get_interface_location(self, structure=None, axis=None, return_fractional=False):

            """
            get_interface_location__binding__artemis_gen_type(self, structure, axis)
            
            
            Defined at \
                ../src/fortran/lib/mod_intf_generator.f90 \
                lines 1112-1124
            
            Parameters:
                this : Artemis_generator_Type
                structure : Basis_Type
                axis : int
                return_fractional : bool
                    If True, return the location in fractional coordinates.
                    If False, return the location in angstroms.

            Returns:
                location : list of floats
                    The location of the interface in the structure (in Å).
                axis : int
                    The axis of the interface.
            """
            if isinstance(structure, Atoms):
                structure = geom_rw.basis(atoms=structure)

            ret_location, ret_axis = _artemis.f90wrap_intf_gen__get_interface_location__binding__agt(this=self._handle, \
                structure=structure._handle, axis=axis, return_fractional=return_fractional)
            
            if ret_axis != axis and axis is not None:
                raise RuntimeError(f"Interface location generation failed (axis {ret_axis} != {axis})")
            
            # convert the location from numpy array to list
            if isinstance(ret_location, numpy.ndarray):
                ret_location = ret_location.tolist()

            return ret_location, ret_axis
                

        def generate(self, surface_lw=None, surface_up=None, thickness_lw=None, \
            thickness_up=None, num_layers_lw=None, num_layers_up=None, \
            reduce_matches=None, \
            print_lattice_match_info=None, print_termination_info=None, \
            print_shift_info=None, break_on_fail=None, icheck_term_pair=None, \
            interface_idx=None, generate_structures=None, seed=None, verbose=None, \
            return_exit_code=False, calc=None):
            """
            generate__binding__artemis_gen_type(self[, surface_lw, \
                surface_up, thickness_lw, thickness_up, num_layers_lw, num_layers_up, \
                print_lattice_match_info, print_termination_info, print_shift_info, \
                break_on_fail, icheck_term_pair, interface_idx, generate_structures, seed, \
                verbose, exit_code])
            
            
            Defined at \
                ../src/fortran/lib/mod_intf_generator.f90 \
                lines 315-1111
            
            Parameters
            ----------
            this : Artemis_generator_Type
            surface_lw : int array
            surface_up : int array
            thickness_lw : float
            thickness_up : float
            num_layers_lw : int
            num_layers_up : int
            print_lattice_match_info : bool
            print_termination_info : bool
            print_shift_info : bool
            break_on_fail : bool
            icheck_term_pair : int
            interface_idx : int
            generate_structures : bool
            seed : int
            verbose : int
            exit_code : int
            
            """

            exit_code = 0
            structures = None

            exit_code = _artemis.f90wrap_intf_gen__generate__binding__agt(this=self._handle, \
                surface_lw=surface_lw, surface_up=surface_up,
                thickness_lw=thickness_lw, thickness_up=thickness_up,
                num_layers_lw=num_layers_lw, num_layers_up=num_layers_up, \
                reduce_matches=reduce_matches, \
                print_lattice_match_info=print_lattice_match_info, \
                print_termination_info=print_termination_info, \
                print_shift_info=print_shift_info, break_on_fail=break_on_fail, \
                icheck_term_pair=icheck_term_pair, interface_idx=interface_idx, \
                generate_structures=generate_structures, seed=seed, verbose=verbose
            )
            if ( exit_code != 0 and exit_code != None )  and not return_exit_code:
                raise RuntimeError(f"Interface generation failed (exit code {exit_code})")
        
            structures = self.get_structures(calc)
            if return_exit_code:
                return structures, exit_code
            return structures

        def restart(self, structure, interface_location=None, print_shift_info=None, \
            seed=None, verbose=None, return_exit_code=False, calc=None):
            """
            restart__binding__artemis_gen_type(self, basis[, \
                interface_location, print_shift_info, seed])
            
            
            Defined at \
                ../src/fortran/lib/mod_intf_generator.f90 \
                lines 202-297
            
            Parameters
            ----------
            this : Artemis_generator_Type
            basis : Basis_Type
            interface_location : float array
            print_shift_info : bool
            seed : int
            
            ---------------------------------------------------------------------------
             Set the random seed
            ---------------------------------------------------------------------------
            """
            exit_code = 0
            structures = None

            # check if host is ase.Atoms object or a Fortran derived type basis_type
            if isinstance(structure, Atoms):
                structure = geom_rw.basis(atoms=structure)

            exit_code = _artemis.f90wrap_intf_gen__restart__binding__agt(this=self._handle, \
                structure=structure._handle, interface_location=interface_location, \
                print_shift_info=print_shift_info, seed=seed, verbose=verbose)
            
            if ( exit_code != 0 and exit_code != None ) and not return_exit_code:
                raise RuntimeError(f"Interface generation failed (exit code {exit_code})")

            structures = self.get_structures(calc)
            if return_exit_code:
                return structures, exit_code
            return structures
        
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
            return _artemis.f90wrap_artemis_gen_type__get__num_structures(self._handle)

        @num_structures.setter
        def num_structures(self, num_structures):
            _artemis.f90wrap_artemis_gen_type__set__num_structures(self._handle, \
                num_structures)

        @property
        def max_num_structures(self):
            """
            The maximum number of generated structures that can be stored in the generator.
            """
            return _artemis.f90wrap_artemis_gen_type__get__num_structures(self._handle)

        @max_num_structures.setter
        def max_num_structures(self, max_num_structures):
            _artemis.f90wrap_artemis_gen_type__set__max_num_structures(self._handle, \
                max_num_structures)

        @property
        def structure_lw(self):
            """
            Element structure_lw ftype=type(basis_type) pytype=Basis_Type
            
            
            Defined at \
                ../fortran/lib/mod_intf_generator.f90 \
                line 32
            
            """
            structure_lw_handle = \
                _artemis.f90wrap_artemis_gen_type__get__structure_lw(self._handle)
            if tuple(structure_lw_handle) in self._objs:
                structure_lw = self._objs[tuple(structure_lw_handle)]
            else:
                structure_lw = geom_rw.basis.from_handle(structure_lw_handle)
                self._objs[tuple(structure_lw_handle)] = structure_lw
            return structure_lw
        
        @structure_lw.setter
        def structure_lw(self, structure_lw):
            structure_lw = structure_lw._handle
            _artemis.f90wrap_artemis_gen_type__set__structure_lw(self._handle, \
                structure_lw)
        
        @property
        def structure_up(self):
            """
            Element structure_up ftype=type(basis_type) pytype=Basis_Type
            
            
            Defined at \
                ../fortran/lib/mod_intf_generator.f90 \
                line 32
            
            """
            structure_up_handle = \
                _artemis.f90wrap_artemis_gen_type__get__structure_up(self._handle)
            if tuple(structure_up_handle) in self._objs:
                structure_up = self._objs[tuple(structure_up_handle)]
            else:
                structure_up = geom_rw.basis.from_handle(structure_up_handle)
                self._objs[tuple(structure_up_handle)] = structure_up
            return structure_up
        
        @structure_up.setter
        def structure_up(self, structure_up):
            structure_up = structure_up._handle
            _artemis.f90wrap_artemis_gen_type__set__structure_up(self._handle, \
                structure_up)
        
        @property
        def elastic_constants_lw(self):
            """
            Element elastic_constants_lw ftype=real(real32) pytype=float
            
            
            Defined at \
                ../fortran/lib/mod_intf_generator.f90 \
                line 34
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _artemis.f90wrap_artemis_gen_type__array__elastic_co4c3f(self._handle)

            if array_handle == 0:
                return None

            if array_handle in self._arrays:
                elastic_constants_lw = self._arrays[array_handle]
            else:
                elastic_constants_lw = \
                    f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _artemis.f90wrap_artemis_gen_type__array__elastic_co4c3f)
                self._arrays[array_handle] = elastic_constants_lw
            return elastic_constants_lw
        
        @elastic_constants_lw.setter
        def elastic_constants_lw(self, elastic_constants_lw):
            self.elastic_constants_lw[...] = elastic_constants_lw
        
        @property
        def elastic_constants_up(self):
            """
            Element elastic_constants_up ftype=real(real32) pytype=float
            
            
            Defined at \
                ../fortran/lib/mod_intf_generator.f90 \
                line 34
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _artemis.f90wrap_artemis_gen_type__array__elastic_coedb6(self._handle)

            if array_handle == 0:
                return None

            if array_handle in self._arrays:
                elastic_constants_up = self._arrays[array_handle]
            else:
                elastic_constants_up = \
                    f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _artemis.f90wrap_artemis_gen_type__array__elastic_coedb6)
                self._arrays[array_handle] = elastic_constants_up
            return elastic_constants_up
        
        @elastic_constants_up.setter
        def elastic_constants_up(self, elastic_constants_up):
            self.elastic_constants_up[...] = elastic_constants_up
        
        @property
        def use_pricel_lw(self):
            """
            Element use_pricel_lw ftype=logical pytype=bool
            
            
            Defined at \
                ../fortran/lib/mod_intf_generator.f90 \
                line 36
            
            """
            return \
                _artemis.f90wrap_artemis_gen_type__get__use_pricel_lw(self._handle)
        
        @use_pricel_lw.setter
        def use_pricel_lw(self, use_pricel_lw):
            _artemis.f90wrap_artemis_gen_type__set__use_pricel_lw(self._handle, \
                use_pricel_lw)
        
        @property
        def use_pricel_up(self):
            """
            Element use_pricel_up ftype=logical pytype=bool
            
            
            Defined at \
                ../fortran/lib/mod_intf_generator.f90 \
                line 36
            
            """
            return \
                _artemis.f90wrap_artemis_gen_type__get__use_pricel_up(self._handle)
        
        @use_pricel_up.setter
        def use_pricel_up(self, use_pricel_up):
            _artemis.f90wrap_artemis_gen_type__set__use_pricel_up(self._handle, \
                use_pricel_up)
        
        @property
        def miller_lw(self):
            """
            Element miller_lw ftype=integer pytype=int
            
            
            Defined at \
                ../fortran/lib/mod_intf_generator.f90 \
                line 38
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _artemis.f90wrap_artemis_gen_type__array__miller_lw(self._handle)
            if array_handle in self._arrays:
                miller_lw = self._arrays[array_handle]
            else:
                miller_lw = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _artemis.f90wrap_artemis_gen_type__array__miller_lw)
                self._arrays[array_handle] = miller_lw
            return miller_lw
        
        @miller_lw.setter
        def miller_lw(self, miller_lw):
            self.miller_lw[...] = miller_lw
        
        @property
        def miller_up(self):
            """
            Element miller_up ftype=integer pytype=int
            
            
            Defined at \
                ../fortran/lib/mod_intf_generator.f90 \
                line 38
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _artemis.f90wrap_artemis_gen_type__array__miller_up(self._handle)
            if array_handle in self._arrays:
                miller_up = self._arrays[array_handle]
            else:
                miller_up = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _artemis.f90wrap_artemis_gen_type__array__miller_up)
                self._arrays[array_handle] = miller_up
            return miller_up

        @miller_up.setter
        def miller_up(self, miller_up):
            self.miller_up[...] = miller_up
        
        @property
        def is_layered_lw(self):
            """
            Element is_layered_lw ftype=logical pytype=bool
            
            
            Defined at \
                ../fortran/lib/mod_intf_generator.f90 \
                line 40
            
            """
            return \
                _artemis.f90wrap_artemis_gen_type__get__is_layered_lw(self._handle)
        
        @is_layered_lw.setter
        def is_layered_lw(self, is_layered_lw):
            _artemis.f90wrap_artemis_gen_type__set__is_layered_lw(self._handle, \
                is_layered_lw)
        
        @property
        def is_layered_up(self):
            """
            Element is_layered_up ftype=logical pytype=bool
            
            
            Defined at \
                ../fortran/lib/mod_intf_generator.f90 \
                line 40
            
            """
            return \
                _artemis.f90wrap_artemis_gen_type__get__is_layered_up(self._handle)
        
        @is_layered_up.setter
        def is_layered_up(self, is_layered_up):
            _artemis.f90wrap_artemis_gen_type__set__is_layered_up(self._handle, \
                is_layered_up)
        
        @property
        def ludef_is_layered_lw(self):
            """
            Element ludef_is_layered_lw ftype=logical pytype=bool
            
            
            Defined at \
                ../fortran/lib/mod_intf_generator.f90 \
                line 42
            
            """
            return \
                _artemis.f90wrap_artemis_gen_type__get__ludef_is_lay4aa6(self._handle)
        
        @ludef_is_layered_lw.setter
        def ludef_is_layered_lw(self, ludef_is_layered_lw):
            _artemis.f90wrap_artemis_gen_type__set__ludef_is_lay87a5(self._handle, \
                ludef_is_layered_lw)
        
        @property
        def ludef_is_layered_up(self):
            """
            Element ludef_is_layered_up ftype=logical pytype=bool
            
            
            Defined at \
                ../fortran/lib/mod_intf_generator.f90 \
                line 42
            
            """
            return \
                _artemis.f90wrap_artemis_gen_type__get__ludef_is_lay60fd(self._handle)
        
        @ludef_is_layered_up.setter
        def ludef_is_layered_up(self, ludef_is_layered_up):
            _artemis.f90wrap_artemis_gen_type__set__ludef_is_laye6e4(self._handle, \
                ludef_is_layered_up)
        
        @property
        def shift_method(self):
            """
            Element shift_method ftype=integer  pytype=int
            
            
            Defined at \
                ../src/fortran/lib/mod_intf_generator.f90 \
                line 31
            
            """
            return \
                _artemis.f90wrap_artemis_gen_type__get__shift_method(self._handle)
        
        @shift_method.setter
        def shift_method(self, shift_method):
            _artemis.f90wrap_artemis_gen_type__set__shift_method(self._handle, \
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
                _artemis.f90wrap_artemis_gen_type__get__num_shifts(self._handle)
        
        @num_shifts.setter
        def num_shifts(self, num_shifts):
            _artemis.f90wrap_artemis_gen_type__set__num_shifts(self._handle, \
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
                _artemis.f90wrap_artemis_gen_type__array__shifts(self._handle)

            if array_handle == 0:
                return None

            if array_handle in self._arrays:
                shifts = self._arrays[array_handle]
            else:
                shifts = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _artemis.f90wrap_artemis_gen_type__array__shifts)
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
                _artemis.f90wrap_artemis_gen_type__get__interface_depth(self._handle)
        
        @interface_depth.setter
        def interface_depth(self, interface_depth):
            _artemis.f90wrap_artemis_gen_type__set__interface_depth(self._handle, \
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
                _artemis.f90wrap_artemis_gen_type__get__separation_scale(self._handle)
        
        @separation_scale.setter
        def separation_scale(self, separation_scale):
            _artemis.f90wrap_artemis_gen_type__set__separation_scale(self._handle, \
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
                _artemis.f90wrap_artemis_gen_type__get__depth_method(self._handle)
        
        @depth_method.setter
        def depth_method(self, depth_method):
            _artemis.f90wrap_artemis_gen_type__set__depth_method(self._handle, \
                depth_method)
        
        def init_array_structure_data(self):
            self.structure_data = f90wrap.runtime.FortranDerivedTypeArray(self,
                                            _artemis.f90wrap_artemis_gen_type__array_getitem__structure_data,
                                            _artemis.f90wrap_artemis_gen_type__array_setitem__structure_data,
                                            _artemis.f90wrap_artemis_gen_type__array_len__structure_data,
                                            """
            Element structure_data ftype=type(struc_data_type) pytype=Struc_Data_Type
            
            
            Defined at \
                ../fortran/lib/mod_generator.f90 line \
                56
            
            """, Misc_Types.struc_data_type)
            return self.structure_data
        
        @property
        def swap_method(self):
            """
            Element swap_method ftype=integer  pytype=int
            
            
            Defined at \
                ../src/fortran/lib/mod_intf_generator.f90 \
                line 45
            
            """
            return \
                _artemis.f90wrap_artemis_gen_type__get__swap_method(self._handle)
        
        @swap_method.setter
        def swap_method(self, swap_method):
            _artemis.f90wrap_artemis_gen_type__set__swap_method(self._handle, \
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
                _artemis.f90wrap_artemis_gen_type__get__num_swaps(self._handle)
        
        @num_swaps.setter
        def num_swaps(self, num_swaps):
            _artemis.f90wrap_artemis_gen_type__set__num_swaps(self._handle, \
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
                _artemis.f90wrap_artemis_gen_type__get__swap_density(self._handle)
        
        @swap_density.setter
        def swap_density(self, swap_density):
            _artemis.f90wrap_artemis_gen_type__set__swap_density(self._handle, \
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
                _artemis.f90wrap_artemis_gen_type__get__swap_depth(self._handle)
        
        @swap_depth.setter
        def swap_depth(self, swap_depth):
            _artemis.f90wrap_artemis_gen_type__set__swap_depth(self._handle, \
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
                _artemis.f90wrap_artemis_gen_type__get__swap_sigma(self._handle)
        
        @swap_sigma.setter
        def swap_sigma(self, swap_sigma):
            _artemis.f90wrap_artemis_gen_type__set__swap_sigma(self._handle, \
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
                _artemis.f90wrap_artemis_gen_type__get__require_mirror_swaps(self._handle)
        
        @require_mirror_swaps.setter
        def require_mirror_swaps(self, require_mirror_swaps):
            _artemis.f90wrap_artemis_gen_type__set__require_mirror_swaps(self._handle, \
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
                _artemis.f90wrap_artemis_gen_type__get__match_method(self._handle)
        
        @match_method.setter
        def match_method(self, match_method):
            _artemis.f90wrap_artemis_gen_type__set__match_method(self._handle, \
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
                _artemis.f90wrap_artemis_gen_type__get__max_num_matches(self._handle)
        
        @max_num_matches.setter
        def max_num_matches(self, max_num_matches):
            _artemis.f90wrap_artemis_gen_type__set__max_num_matches(self._handle, \
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
                _artemis.f90wrap_artemis_gen_type__get__max_num_terms(self._handle)
        
        @max_num_terms.setter
        def max_num_terms(self, max_num_terms):
            _artemis.f90wrap_artemis_gen_type__set__max_num_terms(self._handle, \
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
                _artemis.f90wrap_artemis_gen_type__get__max_num_planes(self._handle)
        
        @max_num_planes.setter
        def max_num_planes(self, max_num_planes):
            _artemis.f90wrap_artemis_gen_type__set__max_num_planes(self._handle, \
                max_num_planes)
        
        @property
        def compensate_normal(self):
            """
            Element compensate_normal ftype=logical pytype=bool
            
            
            Defined at \
                ../src/fortran/lib/mod_intf_generator.f90 \
                line 61
            
            """
            return \
                _artemis.f90wrap_artemis_gen_type__get__compensate_normal(self._handle)
        
        @compensate_normal.setter
        def compensate_normal(self, compensate_normal):
            _artemis.f90wrap_artemis_gen_type__set__compensate_normal(self._handle, \
                compensate_normal)
        
        @property
        def bondlength_cutoff(self):
            """
            Element bondlength_cutoff ftype=real(real32) pytype=float
            
            
            Defined at \
                ../src/fortran/lib/mod_intf_generator.f90 \
                line 65
            
            """
            return \
                _artemis.f90wrap_artemis_gen_type__get__bondlength_cutoff(self._handle)
        
        @bondlength_cutoff.setter
        def bondlength_cutoff(self, bondlength_cutoff):
            _artemis.f90wrap_artemis_gen_type__set__bondlength_cutoff(self._handle, \
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
                _artemis.f90wrap_artemis_gen_type__array__layer_separation_cutoff(self._handle)
            if array_handle in self._arrays:
                layer_separation_cutoff = self._arrays[array_handle]
            else:
                layer_separation_cutoff = \
                    f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _artemis.f90wrap_artemis_gen_type__array__layer_separation_cutoff)
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
                                            _artemis.f90wrap_artemis_gen_type__array_getitem__structures,
                                            _artemis.f90wrap_artemis_gen_type__array_setitem__structures,
                                            _artemis.f90wrap_artemis_gen_type__array_len__structures,
                                            """
            Element items ftype=type(basis_type) pytype=basis


            Defined at ../src/lib/mod_generator.f90 line \
                29

            """, Geom_Rw.basis)
            return self.structures

        def __str__(self):
            ret = ['<artemis_generator_type>{\n']
            ret.append('       num_structures : ')
            ret.append(repr(self.num_structures))
            ret.append(',\n    max_num_structures : ')
            ret.append(repr(self.max_num_structures))
            ret.append('\n     structure_lw : ')
            ret.append(repr(self.structure_lw))
            ret.append(',\n    structure_up : ')
            ret.append(repr(self.structure_up))
            ret.append(',\n    elastic_constants_lw : ')
            ret.append(repr(self.elastic_constants_lw))
            ret.append(',\n    elastic_constants_up : ')
            ret.append(repr(self.elastic_constants_up))
            ret.append(',\n    use_pricel_lw : ')
            ret.append(repr(self.use_pricel_lw))
            ret.append(',\n    use_pricel_up : ')
            ret.append(repr(self.use_pricel_up))
            ret.append(',\n    miller_lw : ')
            ret.append(repr(self.miller_lw))
            ret.append(',\n    miller_up : ')
            ret.append(repr(self.miller_up))
            ret.append(',\n    is_layered_lw : ')
            ret.append(repr(self.is_layered_lw))
            ret.append(',\n    is_layered_up : ')
            ret.append(repr(self.is_layered_up))
            ret.append(',\n    ludef_is_layered_lw : ')
            ret.append(repr(self.ludef_is_layered_lw))
            ret.append(',\n    ludef_is_layered_up : ')
            ret.append(repr(self.ludef_is_layered_up))
            ret.append('\n     shift_method : ')
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
            ret.append(',\n    compensate_normal : ')
            ret.append(repr(self.compensate_normal))
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
    

generator = Generator()

class Artemis(f90wrap.runtime.FortranModule):
    """
    Module artemis
    
    
    Defined at ../src/fortran/artemis.f90 lines \
        1-4
    
    """
    pass
    _dt_array_initialisers = []
    

artemis = Artemis()

