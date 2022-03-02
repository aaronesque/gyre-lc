.. _python-walkthrough:

.. _iOri-Aa.mesa: https://github.com/aaronesque/gyre-lc/raw/master/model/iOri-Aa.mesa
.. _iOri-Ab.mesa: https://github.com/aaronesque/gyre-lc/raw/master/model/iOri-Ab.mesa
.. _iOri-Aa-response.h5: https://github.com/aaronesque/gyre-lc/raw/master/model/iOri-Aa-response.h5
.. _iOri-Ab-response.h5: https://github.com/aaronesque/gyre-lc/raw/master/model/iOri-Ab-response.h5

.. gyre-lc documentation master file, created by

#############################
Python Walkthrough
#############################

This chapter provides a walkthrough of using the GYRE-lc package to calculate a light curve for the eccentric ellipsoidal variable of :math:`{\iota}` Orionis.

*****************************
Seting up your inputs
*****************************

There are 3 sets of inputs to consider when producing a GYRE-lc light curve:

- 1-2 stellar models, depending on how many stars contribute to the overall light curve
- 1-2 tide models. one per stellar model
- the orbital parameters: :math:`a`, :math:`e`, and :math:`\Omega_{orb}`

The iota Orionis Models
=============================

The GitHub repository includes the model data necessary to create a light curve and test GYRE-lc's functionality. You will be creating a light curve for the eccentric ellipsoidal variable of :math:`{\iota}` Orionis, which we'll refer to as simply :math:`{\iota}` Ori for brevity. What follows is a list of required input files and descriptions thereof.  


`iOri-Aa.mesa`_ & `iOri-Ab.mesa`_
    The stellar models for each binary component, :math:`{\iota}` Ori Aa & Ab, were created with MESA using stellar parameters listed in :ads_citet:`Pablo:2017`. The MESA inlists are included for reproducibility of results.

`iOri-Aa-response.h5`_ & `iOri-Ab-response.h5`_
    The tide models and their corresponding GYRE inlists are also included for each component. They are created with GYRE using the parameters listed in :ads_citet:`Pablo:2017`. These contain the amplitudes and frequencies for the first 100 normal modes of a star's tidally excited oscillations.
   
:py:class:`pymsg.PhotGrid`
    Lastly, photometric data for each binary component are required. GYRE-lc works best with MSG, which rapidly interpolates desired spectra and photometry from a grid in :math:`log(g)-T_{eff}` space. For that, you will need to produce a :py:class:`pymsg.PhotGrid` object. Detailed instuctions can be found in the MSG documentation, however a brief walkthrough is included below.


******************************
The GYRE-lc Module
******************************

This walkthrough relies on `MSG <http://www.astro.wisc.edu/~townsend/resource/docs/msg/>`_ for rapid synthesis of photometric data. Download and install MSG, then set the :envvar:`MSG_DIR` environment variable as described in the `MSG Quick Start guide <http://www.astro.wisc.edu/~townsend/resource/docs/msg/user-guide/quick-start.html#quick-start>`_. 

To use GYRE-lc in Python, also make sure the :envvar:`GYRELC_DIR` environment variable is set (see `Quick Start`). I use a Jupyter notebook for this walkthrough, but you may later choose to write a Python script instead should it better suit your workflow.

First, create a new working directory and open a new Jupyter notebok there.

Copy and past the following imports::

    # Import standard modules

    import numpy as np
    import sys
    import os

    # Import pymsg

    sys.path.insert(0, os.path.join(os.environ['MSG_DIR'], 'lib'))
    import pymsg

    # Import gyrelc modules

    sys.path.insert(0, os.path.join(os.environ['GYRELC_DIR'], 'lib'))
    import gyrelc as lc

The :py:mod:`pymsg` and :py:mod:`gyrelc` modules both require :py:mod:`sys` and :py:mod:`os` to be imported, so we do that first. We also import the :py:mod:`numpy` module, which we use extensively.

Next, create a pair of :py:class:`gyrelc.Star` objects using the stellar and tide models provided::

    # Create Star objects
    Aa = lc.Star(mesa_model='iOri-Aa.mesa', gyre_model='iOri-Aa-response.h5')
    Ab = lc.Star(mesa_model='iOri-Ab.mesa', gyre_model='iOri-Ab-response.h5')

Use them, along with the corresponding orbital parameters, as inputs to create a :py:class:`gyrelc.Binary` object::

    # Create Binary object
    iOri = lc.Binary(Aa, Ab, a=132., e=0.764, omega_orb=0.03432)

Now create an ``Observer`` object::

    # Creat an Observer object
    obs = lc.Observer(iOri, 'BRITE-B')

The ``Binary`` object consists of two ``Star`` objects, an ``Irradiation`` object, as well as the various attributes and parameters required to provide the ``Observer`` object sufficient context to synthesize a light curve. The ``Observer`` object primarily contains functions for light curve synthesis and analysis thereof. The last parameter left to specify, the choice of passband, is left as an argument for the ``Observer`` class.

Finally, create a light curve::

    # Specify inclination and argument of periastron
    inc = 62.86
    omega = 122.2

    # Duration of 'observation' and number of points
    omega_orb = iOri.omega_orb
    t = np.linspace(0.5/omega_orb, 2.5/omega_orb, num=2000, endpoint=False)

    flux = obs.find_flux(inc, omega, t)

An important subtlety: the ``find_flux()`` function *requires* the observation time to be in units of the orbital period. Here, I'm simulating a BRITE-B passband observation of :math:`{\iota}` iOri that consists of 2000 data points over 2 orbital periods, begining at half a period past periastron. 

Using :py:mod:`matplotlib`, you may plot your results::

    # Plot

    fig, ax = plt.subplots(sharex=True, figsize=(8,4))

    legend_style = {'framealpha':1.0, 'handlelength':1.2, 'handletextpad':0.5, 'fontsize':'small'}

    ax.plot(t*omega_orb, flux, lw=1, label='BRITE-B')
    ax.legend(loc=1, **legend_style)

    ax.set_xlim(0.5,2.)

    ax.set_title(f'$\iota$Ori light curve, $\omega$={omega}')

    fig.text(0.01, 0.5, r'Mode Flux Perturbation', va='center', rotation='vertical')
    fig.text(0.5, 0.0, f'phase (P={1./omega_orb:4.4f} d)', ha='center')

The legend style and labels are entirely a matter of stylistic choice, but a plot with this *xlim* should look something like this:

.. image:: ./walkthrough-lightcurve.png

.. note:: This project is under active development.

.. rubric:: Footnotes:
