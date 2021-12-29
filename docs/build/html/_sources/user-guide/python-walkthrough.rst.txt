.. _python-walkthrough:

.. gyre-lc documentation master file, created by

==============================
Python Walkthrough
==============================

This chapter provides a walkthrough of using the GYRE-lc package to calculate a light curve for the eccentric ellipsoidal variable of :math:`{\iota}` Orionis.

The iota Orionis Models
-------------------------------

The GitHub repository includes the model data necessary to create a light curve and test GYRE-lc's functionality.

.. Stellar models for Aa and Ab
.. Pulsation models for Aa and Ab


The GYRE-lc Module
-------------------------------

Install like this::

    # Import standard modules

    import numpy as np
    import sys
    import os

    # Import pymsg

    sys.path.insert(0, os.path.join(os.environ['MSG_DIR'], 'lib'))
    import pymsg

    # Import gyrelc modules

    sys.path.insert(0, os.path.join(os.environ['GYRELC_DIR'], 'lib'))
    import atm_coeffs as ac
    import resp_coeffs as rc
    from star_class import Star
    from binary_class import Irradiation, Binary

.. make sure you include the build_spectrum script in the bundle

.. note:: This project is under active development.

