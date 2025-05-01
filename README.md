<p align="center">
<img src="docs/artemis_logo_no_background.png" width="250"/>
</p>

[![GPLv3 workflow](https://img.shields.io/badge/License-GPLv3-yellow.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html "View GPLv3 license")
[![Latest Release](https://img.shields.io/github/v/release/ExeQuantCode/ARTEMIS?sort=semver)](https://github.com/ExeQuantCode/ARTEMIS/releases "View on GitHub")
[![Paper](https://img.shields.io/badge/Paper-Comp_Phys_Comms-orange.svg)](https://doi.org/10.1016/j.cpc.2020.107515)
[![Documentation Status](https://readthedocs.org/projects/artemis-materials/badge/?version=latest)](https://artemis-materials.readthedocs.io/en/latest/?badge=latest "ARTEMIS ReadTheDocs")
[![FPM](https://img.shields.io/badge/fpm-0.11.0-purple)](https://github.com/fortran-lang/fpm "View Fortran Package Manager")
[![CMAKE](https://img.shields.io/badge/cmake-3.27.7-red)](https://github.com/Kitware/CMake/releases/tag/v3.27.7 "View cmake")
[![GCC compatibility](https://img.shields.io/badge/gcc-14.2.0-green)](https://gcc.gnu.org/gcc-14/ "View GCC")


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

ARTEMIS is both a Fortran and a Python library, with the option of a Fortran executable.
The code relies on recent Fortran features, so has no backwards compatibility with Fortran95.


## Documentation


> **_NOTE_:**
> The Read*the*Docs is still under development.
> More guides will be added in the coming weeks and months.

Tutorials and documentation will be provided on the [docs](http://artemis-materials.readthedocs.io/) website.
The methodology is detailed in the [paper](https://doi.org/10.1016/j.cpc.2020.107515).

Refer to the [API Documentation section](#api-documentation) later in this document to see how to access the API-specific documentation.

The Fortran executable/app currently has the most extensive documentation.
This can be found in two forms:
1. [The PDF manual](docs/manual.pdf)
2. The executable help function (`--help` and `--search` flags)


## Requirements

- Fortran compiler supporting Fortran 2018 standard or later
- fpm or CMake (fpm works only for Fortran installation)

Python-specific installation:

- Python 3.11 or later (might work on earlier, have not tested)
- NumPy.f2py
- f90wrap
- cython
- scikit-build-core
- meson
- make or ninja
- CMake
- ASE

The library bas been developed and tested using the following Fortran compilers:
- gfortran -- gcc 11.4.0
- gfortran -- gcc 13.2.0
- gfortran -- gcc 14.1.0
- gfortran -- gcc 14.2.0

The library is known to not currently work with the intel Fortran compilers.

## Installation

For the Python library, the easiest method of installation is to install it directly from pip:

```
pip install artemis-materials
```

Once this is done, ARTEMIS is ready to be used.

Alternatively, to download development versions or, if, for some reason, the pip method does not work, then ARTEMIS can be installed from the source.
To do so, the source must be obtained from the git repository.
Use the following commands to get started:
```
 git clone https://github.com/ExeQuantCode/artemis.git
 cd artemis
```

Depending on what language will be used in, installation will vary from this point.

### Python

For Python, the easiest installation is through pip:
```
pip install .
```

Another option is installing it through cmake, which involves:
```
mkdir build
cd build
cmake ..
make install
```

Then, the path to the install directory (`${HOME}/.local/artemis`) needs to be added to the include path. NOTE: this method requires that the user manually installs the `ase`, `numpy` and `f90wrap` modules for Python.

### Fortran

For Fortran, either fpm or cmake are required.

#### fpm

fpm installation is as follows:

```
fpm build --profile release
```

This will install both the Fortran library and the Fortran application for ARTEMIS.
The library can then be called from other fpm-built Fortran programs through normal means (usually referencing the location of ARTEMIS in the program's own `fpm.toml` file).
The application can be run using
```
fpm run
```

#### cmake

cmake installation is as follows:
```
mkdir build
cd build
cmake [-DBUILD_PYTHON=Off] ..
make install
```
The optional filed (dentoted with `[...]`) can be used to turn off installation of the Python library.
This will build the library in the build/ directory.
All library files will then be found in:
```
${HOME}/.local/artemis
```
Inside this directory, the following files will be generated:
```
include/artemis.mod
lib/libartemis.a
```

How-to
------
Until recently, ARTEMIS has existed solely as a Fortran executable.
This version of the code is currently best documented, but this will change in the near future as the Python library is tested more.
To get an example input file, run the following command:  
```
artemis -d example.in
```

This will generate the file `example.in`, with the structure of the
ARTEMIS input file.

To get descriptions of the tags within the input file, run either command:

```
artemis --help [all|<TAGNAME>]  
artemis --search <STRING>
```



Websites
--------

Group webpage: http://www.artemis-materials.co.uk

Group wiki:    http://www.artemis-materials.co.uk/HRG

Guide and documentation: https://artemis-materials.readthedocs.io/en/latest/


API documentation
-----------------

> **_NOTE_:**
> API documentation is not yet set up.
> It is planned to be implemented in an upcoming release to work alongside the Read*the*Docs and Python library.

 <!-- API documentation can be generated using FORD (Fortran Documenter).
To do so, follow the installation guide on the [FORD website](https://forddocs.readthedocs.io/en/stable/) to ensure FORD is installed.
Once FORD is installed, run the following command in the root directory of the git repository:

```
  ford ford.md
``` -->

Contributing
------------

Please note that this project adheres to the [Contributing Guide](CONTRIBUTING.md).
If you want to contribute to this project, please first read through the guide.
If you have any questions, bug reports, or feature requests, please either discuss then in [issues](https://github.com/ExeQuantCode/artemis/issues).

For any serious or private concerns, please use the following email address:
support@artemis-materials.co.uk


License
------------
This work is licensed under a [GPL v3 license]([https://opensource.org/license/mit/](https://www.gnu.org/licenses/gpl-3.0.en.html)).

Developers
------------
- Ned Thaddues Taylor  
- Francis Huw Davies  
- Isiah Edward Mikel Rudkin  
- Steven Paul Hepplestone  

Contributers
------------
- Conor Jason Price  
- Tsz Hin Chan  
- Joe Pitfield  
- Edward Allery Baker  
- Shane Graham Davies  

Advisors
------------
- Elizabeth L. Martin



## References

If you use this code, please cite our paper:
```text
@article{Taylor2020ARTEMISAbInitio,
  title = {ARTEMIS: Ab initio restructuring tool enabling the modelling of interface structures},
  volume = {257},
  ISSN = {0010-4655},
  url = {http://dx.doi.org/10.1016/j.cpc.2020.107515},
  DOI = {10.1016/j.cpc.2020.107515},
  journal = {Computer Physics Communications},
  publisher = {Elsevier BV},
  author = {Taylor,  Ned Thaddeus and Davies,  Francis Huw and Rudkin,  Isiah Edward Mikel and Price,  Conor Jason and Chan,  Tsz Hin and Hepplestone,  Steven Paul},
  year = {2020},
  month = dec,
  pages = {107515}
}
```

This README has been copied from the [RAFFLE repository](https://github.com/ExeQuantCode/RAFFLE), with permission from the creator (Ned Taylor).