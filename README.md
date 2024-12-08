<p align="center">
<img src="docs/artemis_logo_no_background.png" width="250"/>
</p>

[![License: CC BY-NC 3.0](https://img.shields.io/badge/License-CC_BY--NC_3.0-lightgrey.svg)](https://creativecommons.org/licenses/by-nc/3.0/ "View CC BY-NC 3.0 license")
[![Latest Release](https://img.shields.io/github/v/release/ExeQuantCode/ARTEMIS?sort=semver)](https://github.com/ExeQuantCode/ARTEMIS/releases "View on GitHub")
[![Paper](https://img.shields.io/badge/Paper-Comp_Phys_Comms-orange.svg)](https://doi.org/10.1016/j.cpc.2020.107515)
[![GCC compatibility](https://img.shields.io/badge/gcc-14.1.0-green)](https://gcc.gnu.org/gcc-14/ "View GCC")


Ab Initio Restructuring Tool Enabling Modelling of Interface Structures
=========================================================================
by Ned Thaddeus Taylor and Steven Paul Hepplestone, The ARTEMIS Research Group (Hepplestone Research Group)


## New Repository Location

This repository has been migrated from the University of Exeter GitLab to GitHub to facilitate community interaction and support. The latest version, updates, and collaboration now take place on this GitHub repository.

**GitLab Repository (Archived):** https://git.exeter.ac.uk/hepplestone/artemis

## Why the Migration?

It was decided that this project should be migrated to allow for better community support (i.e. allowing community users to raise issues).
All information has been ported over where possible.

---


ARTEMIS is a software package for the generation and modelling of interfaces between materials.

ARTEMIS is distributed with the following directories:

  docs/       Documentation  
  src/       Source code  
  tools/     Extra shell script tools  
  examples/  Example ARTEMIS files  

After ARTEMIS is compiled, the following directories may also exist:

  bin/       Contains binary executables  
  obj/       Contains object (built/indermetiate) files, which are compiled binary files that haven't been linked yet

For further information please see the User manual (docs/manual.pdf)



Setup
-----
Run the following command in the directory containing the Makefile:  
make

This should create a bin directory, in which the executable  
'artemis' can be found. This directory should be found in the
DARTEMIS directory.



How-to
------
ARTEMIS mainly works off of an input file, but can also perform some
actions via flags.

To get an example input file, run the following command:  
artemis -d example.in

This will generate the file 'example.in', with the structure of the
ARTEMIS input file.

To get descriptions of the tags within the input file, run either command:  
artemis --help [TAGNAME]  
artemis --search <STRING>


Further documentation on the workings of ARTEMIS can be found in the docs/
directory or on the wiki (linked below)



Websites
--------
Webpage: http:/www.artemis-materials.co.uk/  
Wiki:    http://www.artemis-materials.co.uk/HRG



Contact
-------
Please log issues, bug-reports, and feature requests on the issue tracker for this repository: https://github.com/ExeQuantCode/ARTEMIS/issues

For any serious or private concerns, please use the following email address:
support@artemis-materials.co.uk



Developers
------------
-Ned Thaddues Taylor  
-Francis Huw Davies  
-Isiah Edward Mikel Rudkin  
-Steven Paul Hepplestone  

Contributers
------------
-Conor Jason Price  
-Tsz Hin Chan  
-Joe Pitfield  
-Edward Allery Baker  
-Shane Graham Davies  

Advisors
------------
-Elizabeth L. Martin


License
------------
This work is licensed under a Creative Commons Attribution-NonCommercial 3.0 Unported (CC BY-NC 3.0) License.  
https://creativecommons.org/licenses/by-nc/3.0/


Source file descriptions
------------
src/main.f90           - main file that calls the functions and determines the task of the job  
src/inputs.f90         - handles input file and assigned default values to parameters  
src/interfaces.f90     - task 1 ARTEMIS job. Calls subroutines to generate interfaces  
src/aspect.f90         - task 0 ARTEMIS job. Calls subroutines to edit structure  
src/io.F90             - error handling file, help, search and startup printing  
src/mod_help.f90       - descriptions of all input tags  
src/default_infile.f90 - prints default/example input file of ARTEMIS  
src/mod_shifting.f90   - identifies and generates sets of interface shifts  
src/mod_swapping.f90   - generates sets of swaps (intermixing)  
src/mod_intf_identifier.f90 - identifies interface axis and location for pregen interface  
src/mod_lat_compare.f90     - performs lattice matching over a set of Miller planes  
src/mod_plane_matching.f90  - performs lattice matching over a single Miller plane  

src/lib/mod_constants.f90    - a set of global constants used in this code  
src/lib/mod_misc.f90         - miscellaneous functions and subroutines  
src/lib/mod_misc_maths.f90   - maths functions and subroutines  
src/lib/mod_misc_linalg.f90  - linear algebra functions and subroutines  
src/lib/mod_rw_geom.f90      - read and write structure (geometry) files  
src/lib/mod_edit_geom.f90    - tools to edit lattice and basis (geometry editing)  
src/lib/mod_sym.f90          - tools to apply and determine symmetries between bases  
src/lib/mod_tools_infile.f90 - tools to read input files  



Other files
------------
README.md         - a readme file with a brief description of the code and files  
Makefile          - the makefile used for compiling the code  
LICENSE           - license of ARTEMIS code  
CHANGE.LOG        - changelog for ARTEMIS  
artemis.ascii     - ARTEMIS logo in ascii form  
artemis_logo.pdf  - ARTEMIS logo  
docs/manual.pdf    - pdf of ARTEMIS manual/user guide  
docs/manual.tex    - tex file of ARTEMIS manual  
tools/compress.sh - script to compress ARTEMIS directory  
examples/generate_interface/param.in    - example input file  
examples/generate_interface/POSCAR_Si   - silicon 8 atom unit cell  
examples/generate_interface/POSCAR_Ge   - germanium 8 atom unit cell  
examples/generate_interface/DINTERFACES - directory containing example output interface structures  
examples/pregenerated_interface/param.in     - example input file  
examples/pregenerated_interface/POSCAR       - CaCu3Ti4O12|CuO interface structure  
examples/pregenerated_interface/DINTERFACES  - directory containing example output interface structures  
examples/identify_terminations/param.in      - example input file  
examples/identify_terminations/POSCAR        - silicon 2 atom primitive cell  
examples/identify_terminations/DTERMINATIONS - directory containing example output slab structures  
