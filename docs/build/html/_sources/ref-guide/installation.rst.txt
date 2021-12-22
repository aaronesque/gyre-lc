.. gyre-lc documentation master file, created by

===================================
Installation
===================================

This chapter discusses MSG installation in detail. If you just want to get up and running, have a look at the Quick Start chapter.

Prerequisites
-----------------------------------

To compile and run MSG, you’ll need the following software components:

- A modern (2003+) Fortran compiler
- The BLAS linear algebra library
- The LAPACK linear algebra library
- The LAPACK95 Fortran 95 interface to LAPACK
- The HDF5 data management library

On Linux and MacOS platforms, these components are bundled together in the MESA Software Development Kit (SDK), which can be downloaded from the MESA SDK homepage. Using this SDK is strongly recommended.

Building MSG
------------------------------------

Download
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Download the MSG source code, and unpack it from the command line using the tar utility:

``tar xf msg-main.tar.gz``

Set the MSG_DIR environment variable with the path to the newly created source directory; this can be achieved e.g. using the realpath command1:

``export MSG_DIR=$(realpath msg-main)``

Compile
Compile MSG using the make utility:

``make -j -C $MSG_DIR``
(the -j flags tells make to use multiple cores, speeding up the build).

Test
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
To check that MSG has compiled correctly and gives reasonable results, you can run the calculation test suite via the command

make -C $MSG_DIR test
The initial output from the tests should look something like this:

If things go awry, consult the troubleshooting chapter.

Custom Builds
Custom builds of MSG can be created by setting certain environment variables, and/or variables in the file $MSG_DIR/src/build/Makefile, to the value yes. The following variables are currently supported:

