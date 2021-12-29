.. _python-walkthrough:

.. gyre-lc documentation master file, created by

==============================
Python Walkthrough
==============================

This chapter provides a walkthrough of using the GYRE-lc package to calculate a light curve for the eccentric ellipsoidal variable of :math:`{\iota}` Orionis.

The iota Orionis Models
-------------------------------

The GitHub repository includes the model data necessary to create a light curve and test GYRE-lc's functionality. You will be creating a light curve for the eccentric ellipsoidal variable of :math:`{\iota}` Orionis, which we'll refer to as simply :math:`{\iota}` Ori for brevity. 

The stellar models for each binary component, :math:`{\iota}` Ori Aa & Ab, were created with MESA using stellar parameters listed in :footcite:t:`2017:Pablo`. The inlist is included for reproducibility of results.

The pulsations models and their corresponding inlists are also included for each component. They are created with GYRE using the parameters listed in :footcite:t:`2017:Pablo`. These contain the amplitudes and frequencies for the first 100 normal modes of a star's tidally excited oscillations.

Lastly, model spectra (produced with `SYNSPEC <http://tlusty.oca.eu/Synspec49/synspec.html>`_) for each component are also included for testing purposes, but they are entirely optional. GYRE-lc works best with MSG-- for that, three spectral grids are included corresponding to the filters BRITE-R, BRITE-B, and Kepler.



The GYRE-lc Module
-------------------------------

To use GYRE-lc in Python, first make sure the ``GYRELC_DIR`` environment variable is set (see `Quick Start`). I use a Jupyter notebook for this walkthrough, but you may choose to write a Python script instead if you desire.

Next, copy and past the following imports::

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

.. footbibliography::
