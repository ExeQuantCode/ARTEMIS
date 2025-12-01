.. _install:

============
Installation
============

For the Python library, the easiest method of installation is to install it directly from pip:

.. code-block:: bash

    pip install artemis-materials

Once this is done, ARTEMIS is ready to be used (both the Python library and the command line interface).

Alternatively, to install ARTEMIS from source, follow the instructions below.


ARTEMIS can be installed in one of three ways; as a Python package, as a Fortran library, or as a standalone Fortran executable.
All versions rely on the core Fortran code, with the Python package and standalone executable wrapping this code in a Python and Fortran interface, respectively.

The code is hosted on `GitHub <https://github.com/ExeQuantCode/artemis>`_.

This can be done by cloning the repository:

.. code-block:: bash

    git clone https://github.com/ExeQuantCode/artemis.git
    cd artemis

Depending on what language will be used in, installation will vary from this point.


Global requirements
===================

All installation methods require the following dependency:

- Fortran compiler (gfortran>=13.1, not compatible with intel compilers)

Python
======

Requirements
------------

- python (>=3.11)
- `pip <https://pip.pypa.io/en/stable/>`_
- `f90wrap <https://github.com/jameskermode/f90wrap>`_ (>=0.2.14)
- `numpy <https://numpy.org>`_ (>=1.26)
- `meson <https://mesonbuild.com>`_ (>=1.6)
- `cython <https://cython.org>`_ (>=3.0)
- `sckit-build-core <https://scikit-build-core.readthedocs.io/en/latest/>`_ (>=0.11)
- `cmake <https://cmake.org>`_ (>=3.17)
- `ninja <https://ninja-build.org>`_ (>=1.10) or `GNU Make <https://www.gnu.org/software/make/>`_
- `ASE <https://wiki.fysik.dtu.dk/ase/>`_ (>=3.23)


Installation using pip
-----------------------

The easiest way to install ARTEMIS is via pip.
The package is directly available via PyPI, so can be installed without downloading the repository. To do so, run:

.. code-block:: bash

    pip install artemis-materials

This will install the ARTEMIS package and all its dependencies in the default location.
This is the recommended method of installation, as it is the easiest and most straightforward way to get started with ARTEMIS.

Another option is to install ARTEMIS from the source code, which is recommended if you want to use the latest version of ARTEMIS or if you want to contribute to its development.
To do this, you will need to clone the repository from GitHub.

Once the library is cloned, navigate to the root directory of the repository and run:

.. code-block:: bash

    pip install .

Depending on your setup, this will install the Python package and all its dependencies in different places.
To find where this has been installed, you can run:

.. code-block:: bash

    pip show artemis-materials

This will show you the location of the installed package, in addition to other information about the package.

Installation using cmake
------------------------

Another option is installing it through cmake, which involves:
.. code-block:: bash

    mkdir build
    cd build
    cmake ..
    make install

Then, the path to the install directory (`${HOME}/.local/artemis`) needs to be added to the include path.
NOTE: this method requires that the user manually installs the `ase`, `numpy` and `f90wrap` modules for Python.

Fortran
=======

Requirements
------------

- `cmake <https://cmake.org>`_ (>=3.17) or `fpm <https://fpm.fortran-lang.org>`_ (>=0.9.0)
- `GNU Make <https://www.gnu.org/software/make/>`_ (if using cmake)


As mentioned, the Fortran library provides the same functionality as the Python package, but in Fortran instead.

To install the Fortran library or executable, the recommended method is to use the Fortran package manager (fpm).
Cmake is also supported.

Installation using fpm
----------------------

To install the Fortran library and the executable using fpm, navigate to the root directory of the repository and run:

.. code-block:: bash

    fpm build --profile release
    fpm install

This can also be set up as a dependency in your own fpm project by adding the following to your ``fpm.toml`` file:

.. code-block:: toml

    [dependencies]
    artemis = { git = "https://github.com/ExeQuantCode/ARTEMIS" }


Installation using cmake
------------------------

To install the Fortran library using cmake, navigate to the root directory of the repository and run:

.. code-block:: bash

    mkdir build
    cd build
    cmake -DBUILD_PYTHON=Off -DBUILD_EXECUTABLE=Off ..
    make
    make install

This will build the Fortran library and install it in the default location (``~/.local/artemis``).

To install the standalone executable, run:

.. code-block:: bash

    mkdir build
    cd build
    cmake -DBUILD_PYTHON=Off -DBUILD_EXECUTABLE=On ..
    make
    make install

This will build the Fortran library and install it in the default location (``~/.local/artemis``).


Installing on MacOS (Homebrew)
==============================

ARTEMIS is developed on Linux and MacOS, and should work on both.
However, there are likely some additional steps required to install ARTEMIS on MacOS.
This is because **it is not recommended to rely on the Mac system Python, or Fortran and C compilers**.

The recommended way to install Python, gfortran and gcc on MacOS is to use `Homebrew <https://brew.sh>`_.
First, install Homebrew by following the guide on their website.

Once Homebrew is installed, you can install the required dependencies by running:

.. code-block:: bash

    brew install python
    brew install gcc
    brew install gfortran
    export CC=$(brew --prefix gfortran)
    export FC=$(brew --prefix gcc)

Confirm a successful Python installation by running:

.. code-block:: bash

    python --version
    whereis python

This should show the correct Python version (3.11 or later) and path.

Next, if you are using ``pip``, then the following command is found to result in the least issues:

.. code-block:: bash

    python -m pip install --upgrade .

This ensures that the correct Python version is being called, and that the correct version of ``pip`` is being used.
