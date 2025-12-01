.. _about:

=====
About
=====


ARTEMIS (Ab Initio Restructuring Tool Enabling Modelling of Interface Structures) is a package for generating lattice matched interfaces between material.
ARTEMIS interfaces with the `Atomic Simulation Environment (ASE) <https://gitlab.com/ase/ase>`_.

ARTEMIS is both a Fortran and a Python library, with the option of a Fortran executable.
The code heavily relies on features of recent Fortran releases, so there is no backwards compatibility with Fortran95.

The library enables users to provide two crystal structures, from which it generates a set of lattice matched interfaces within user-defined tolerances.