=======
ARTEMIS
=======

ARTEMIS (Ab Initio Restructuring Tool Enabling Modelling of Interface Structures) is a Python and Fortran package for generating lattice matched structured between materials.
ARTEMIS can be utilised as a Python package, a Fortran library, or a standalone Fortran executable.
The Python package provides a high-level interface to the Fortran library, which contains the core functionality.

The Python package interfaces seemlessly with `ASE (Atomic Simulation Environment) <https://wiki.fysik.dtu.dk/ase/>`_, allowing for easy reading, writing, and manipulation of atomic structures.
Although the package comes with a built-in atomic structure reader and writer, it is recommended to use ASE due to its greater functionality and wide-reaching support.

The code is provided freely available under the `GNU General Public License v3.0 <https://www.gnu.org/licenses/gpl-3.0.en.html>`_.

An example 

.. code-block:: python

    # A simple example of how to use ARTEMIS to generate lattice matches structures between silicon and germanium and write them to a single file.
    from ase import Atoms
    from ase.build import bulk
    from ase.io import write
    from artemis.generator import artemis_generator
    from mace.calculators import mace_mp
    from ase.calculators.singlepoint import SinglePointCalculator

    generator = artemis_generator()

    calc = mace_mp(model="medium", dispersion=False, default_dtype="float32", device='cpu')

    Si = bulk('Si', 'diamond', a=5.43, cubic=True)
    Ge = bulk('Ge', 'diamond', a=5.66, cubic=True)

    generator.set_materials(Si, Ge)

    generator.set_surface_properties(
        miller_lw = [ 1, 1, 0 ],
        miller_up = [ 1, 1, 0 ],
    )

    structures = generator.generate(calc=calc)

    for structure in structures:
        structure.calc = SinglePointCalculator(
            structure,
            energy=structure.get_potential_energy(),
            forces=structure.get_forces()
        )
 
    write('structures.traj', structures)

.. toctree::
   :maxdepth: 3
   :caption: Contents:

   about
   install
   tutorials/index
   faq
..    tutorials/index
..    Python API <modules>

.. Indices and tables
.. ==================

.. * :ref:`genindex`
.. * :ref:`modindex`
.. * :ref:`search`