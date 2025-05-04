.. parameters:

==================
Setting parameters
==================

This tutorial will detail how to initialise an ARTEMIS generator.
It will also explore the parameters associated with the lattice matching and interface alignment methods used by ARTEMIS, in addition to its surface termination identification parameters.

Initialisation
--------------
ARTEMIS is initialised by importing the generator object.
The object is the main interface for the user to interact with ARTEMIS.

.. code-block:: python

    # Initialise ARTEMIS generator
    from artemis.generator import artemis_generator

    generator = artemis_generator()
    
It is recommended to use the Atomic Simulation Environment (ASE)~\cite{ase-paper} for handling structure data.
Whilst ARTEMIS can handle its own atomic structure object, ASE is more widely used and has a more extensive feature set.


Constituent structures
----------------------

The first step in using ARTEMIS is to define the constituent structures.
The generator object has a method called ``set_materials`` which takes a list of ASE atoms objects.

.. code-block:: python

    from ase.build import bulk

    # Define the constituent structures
    Si = bulk('Si', 'diamond', a=5.43, cubic=True)
    Ge = bulk('Ge', 'diamond', a=5.66, cubic=True)

    generator.set_materials(Si, Ge)

The above code defines two bulk structures, silicon and germanium, and sets them as the constituent structures for the generator object.
This method can also be used to define the elastic constants of the constituent structures and define whether to identify and use the primitive cell for each structure.
These can be accessed by the following parameters:

.. code-block:: python

    # Set the elastic constants and primitive cell usage
    generator.set_materials(
        structure_lw=Si,
        structure_up=Ge,
        elastic_constants_lw=6,
        elastic_constants_up=12,
        use_pricel_lw=True,
        use_pricel_up=True
    )

The elastic constants are currently isotropic bulks moduli.
The elastic constants can be calculated using ASE or obtained from the literature, such as the Materials Project~\cite{materials-project}.
The primitive cell usage is a boolean value that indicates whether to use the primitive cell of the structure or not.


Surface properties
------------------

The next step is to define the surface properties of the interface.
The generator object has a method called ``set_surface_properties`` which takes the Miller indices of the surface planes to be used.
If no Miller indices are provided, the generator will search over the 10 lowest symmetry planes.

.. code-block:: python

    # Define the surface properties
    generator.set_surface_properties(
        miller_lw=[1, 1, 0],
        miller_up=[1, 1, 0]
    )

The above code sets the Miller indices of the surface planes to be used for the lower and upper structures.
The Miller indices are a set of three integers that describe the orientation of the surface planes in the crystal lattice.
Additional parameters can be set to define the surface properties, such as:

.. code-block:: python

    # Set additional surface properties
    generator.set_surface_properties(
        miller_lw=[1, 1, 0],
        miller_up=[1, 1, 0],
        is_layered_lw=True,
        is_layered_up=True,
        require_stoichiometry_lw=True,
        require_stoichiometry_up=True,
        layer_separation_cutoff_lw=0.5,
        layer_separation_cutoff_up=0.5,
    )

The above code sets the following additional parameters:
- ``is_layered_lw`` and ``is_layered_up``: boolean values that indicate whether the lower and upper structures are to be treated as layered or not.
- ``require_stoichiometry_lw`` and ``require_stoichiometry_up``: boolean values that indicate whether the generated lower and upper slabs should be stoichiometrically equivalent to their respective provided structures.
- ``layer_separation_cutoff_lw`` and ``layer_separation_cutoff_up``: float values that define the cutoff distance for the minimally accepted layer separation (in Angstroms) with which to define distinct planes of atoms.

