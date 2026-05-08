.. identify_translations:

=======================================
Interface translation vectors
=======================================

Interface translation utilities are available through the Python generator API.
`get_interface_definition(...)` returns the detected interface bounds and axis,
`get_interface_translations(...)` returns the primitive in-plane translation
vectors for the unique interface-shift cell, and
`is_valid_interface_translation(...)` checks whether a trial in-plane shift is a
true slab self-translation rather than a distinct interface offset.

.. code-block:: python

    from ase.io import read
    from artemis.generator import artemis_generator

    generator = artemis_generator()
    atoms = read("example/python_pkg/MoS2-Ag_0.xyz")

    bounds, axis = generator.get_interface_definition(atoms)
    t1, t2 = generator.get_interface_translations(atoms, axis=axis, bounds=bounds)

    print(bounds, axis)
    print(t1, t2)
    print(generator.is_valid_interface_translation(atoms, t1))


The abrupt `method=4` shift search uses these primitive vectors directly, so the
in-plane search is confined to the unique primitive translation region rather
than a larger translationally redundant supercell.
