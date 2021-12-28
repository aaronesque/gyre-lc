.. _python-walkthrough:

.. gyre-lc documentation master file, created by

==============================
Python Walkthrough
==============================

This chapter provides a walkthrough of using the GYRE-lc package to calculate a light curve for the eccentric ellipsoidal variable of :math:`{\iota}` Orionis.

The :math:`{\iota}` Ori Models
-------------------------------

The GitHub repository includes the model data necessary to create a light curve and test GYRE-lc's functionality.

Stellar models for Aa and Ab
Pulsation models for Aa and Ab


The GYRE-lc Module
-------------------------------

Install like this::

   # Import gyre-lc

   import sys
   import os
   sys.path.insert(0, os.path.join(os.environ['GYRELC_DIR'], 'lib'))
   import gyrelc

   # Import standard modules and configure them

   import numpy as np
   import matplotlib.pyplot as plt

   %matplotlib inline
   plt.rcParams.update({'font.size': 16})

.. make sure you include the build_spectrum script in the bundle

.. note:: This project is under active development.

