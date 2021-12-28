.. gyre-lc documentation master file, created by

===================================
Installation
===================================

This chapter discusses GYRE-lc installation in detail. If you just want to get up and running, have a look at the Quick Start chapter.

Prerequisites
-----------------------------------

A complete GYRE-lc workflow typically requires the use of additional software to produce the star and pulsation models that go into GYRE-lc as input for light curve synthesis. This includes:

- The `MESA Software Development Kit (SDK) <http://www.astro.wisc.edu/~townsend/static.php?ref=mesasdk>`_, which provides the compilers and supporting libraries needed to build GYRE-lc.
- `MESA <mesa.sourceforge.net>`_, which calculates the stellar models compatible with GYRE-lc.
- `GYRE <https://gyre.readthedocs.io/en/stable/>`_, which calculates the pulsation models compatible with GYRE-lc.
- `MSG <http://www.astro.wisc.edu/~townsend/resource/docs/msg/>`_, which rapidly interpolates stellar spectra from a multidimensional grid for GYRE-lc.

GYRE and MSG are currently officially compatible with Linux and MacOS platforms only- Windows at your own risk!

Most importantly, GYRE-lc requires Python 3.6+. 

To run GYRE-lc, you’ll need the following Python libraries installed:

- `<https://pypi.org/project/h5py/>`_, for HDF5 data management;
- `<https://pypi.org/project/f90nml/>`_, for namelist handling;
- `<https://pypi.org/project/scipy/>`_, for special math functions and operations;
- `<https://pypi.org/project/astropy/>`_, for MESA model handling; 

These components can be found via the PIP and Anaconda python package installers.


Setting up GYRE-lc
------------------------------------

Download
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Download the GYRE-lc source code, and unpack it from the command line using the tar utility:

``tar xf gyre-lc.tar.gz``

Set the GYRELC_DIR environment variable with the path to the newly created source directory; this can be achieved e.g. using the realpath command1:

``export GYRELC_DIR=$(realpath gyre-lc)``

You are ready to test and use GYRE-lc.

Test
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
To check that GYRE-lc functions as expected, you can run the calculation test suite via the command

``python $GYRELC_DIR/test.py``

The initial output from the tests should look something like this:

If things go awry, consult the troubleshooting chapter.

