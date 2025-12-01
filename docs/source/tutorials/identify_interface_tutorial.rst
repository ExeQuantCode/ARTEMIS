.. identify_interface:

==================
Identify interface
==================

This tutorial demonstrates how to use the ARTEMIS library to return the interface location in an interface structure.


The following code snippet shows how to use ARTEMIS to identify the interface location in a structure.


.. code-block:: python

    # Import the necessary libraries
    from ase.io import read
    from artemis.generator import artemis_generator

    # Read the interface structure from a file
    atoms = read("interface.xyz")

    # Initialise the ARTEMIS generator
    generator = artemis_generator()

    # Get the interface location and axis using ARTEMIS
    location, axis = generator.get_interface_location(atoms, return_fractional=True)
    print("location", location)
    print("axis", axis)

The interface location is returned as a single value, which is the distance from the origin of the structure to the interface in the direction of the returned axis.
The axis is an integer specifying the direction of the interface in the structure (i.e. 0, 1, or 2 for a, b, or c respectively).

The `return_fractional` argument specifies whether to return the interface location in fractional coordinates (True) or in Cartesian coordinates (False).
The default value is False.

This can then be used in conjunction with `RAFFLE <https://raffle-fortran.readthedocs.io/>`_ to reconfigure atoms near to the interface to search for more stable configurations.
