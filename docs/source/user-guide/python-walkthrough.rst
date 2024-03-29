.. _python-walkthrough:

.. _MSG: http://www.astro.wisc.edu/~townsend/resource/docs/msg/
.. _iOri-Aa.mesa: https://github.com/aaronesque/gyre-lc/raw/master/model/iOri-Aa.mesa
.. _iOri-Ab.mesa: https://github.com/aaronesque/gyre-lc/raw/master/model/iOri-Ab.mesa
.. _iOri-Aa-response.h5: https://github.com/aaronesque/gyre-lc/raw/master/model/iOri-Aa-response.h5
.. _iOri-Ab-response.h5: https://github.com/aaronesque/gyre-lc/raw/master/model/iOri-Ab-response.h5

.. gyre-lc documentation master file, created by

#############################
Python Walkthrough
#############################

This chapter provides a walkthrough of using the GYRE-lc package to calculate a light curve for the eccentric ellipsoidal variable of :math:`{\iota}` Orionis. Fig. 1 shows a visual summary of the user workflow for this, including the class objects you will pass into each other along with your inputs. 

.. figure:: ./walkthrough-flowchart.png

   Figure 1. A flowchart describing how a user can create a light curve using GYRE-lc. Notice how the :py:class:`gyrelc.Observer` object, which contains the :py:func:`find_flux()` method, takes :py:class:`gyrelc.Binary` as an input, which in turn takes 2 :py:class:`gyrelc.Star` as input. Although both :py:class:`gyrelc.Star` share the same :py:class:`msg.PhotGrid`, they each have their own corresponding tidal response and stellar models.

.. _python-walkthrough-inputs:

*****************************
Preparing Your Inputs
*****************************

There are several inputs to consider when producing a GYRE-lc light curve:

- 1-2 stellar models, depending on how many stars contribute to the overall light curve
- a tidal response model for each star
- a photometric model for each star
- the orbital parameters: :math:`a`, :math:`e`, and :math:`\Omega_{orb}` for the binary


The iota Orionis Models
=============================

The GitHub repository includes the model data necessary to create a light curve and test GYRE-lc's functionality. You will be creating a light curve for the eccentric ellipsoidal variable of :math:`{\iota}` Orionis, which we'll refer to as simply :math:`{\iota}` Ori for brevity. What follows is a list of required input files and descriptions thereof.  


`iOri-Aa.mesa`_ & `iOri-Ab.mesa`_
    The stellar models for each binary component, :math:`{\iota}` Ori Aa & Ab, were created with MESA using stellar parameters listed in :ads_citet:`Pablo:2017`. The MESA inlists are included for reproducibility of results.

`iOri-Aa-response.h5`_ & `iOri-Ab-response.h5`_
    The tide models and their corresponding GYRE inlists are also included for each component. They are created with GYRE using the parameters listed in :ads_citet:`Pablo:2017`. These contain the amplitudes and frequencies for the first 100 normal modes of a star's tidally excited oscillations.


The model filenames themselves don't have to be in any particular format, so feel free to name them whatever you find most convenient. Mind these files in `$GYRELC_DIR/models` and place them in your working directory.`.

The photometric grids
==============================

Lastly, photometric data for each binary component are required. GYRE-lc works best with MSG, which rapidly interpolates desired spectra and photometry from a grid in :math:`log(g)-T_{eff}` space. For that, you will need to produce a :py:class:`pymsg.PhotGrid` object using an MSG-produced spectral grid and a properly formatted passband file. Detailed instuctions can be found in the `MSG`_ documentation, however a brief walkthrough is included below.

Before starting Jupyter, download and place the following files in your working directory:

* `sg-demo.h5` from `MSG`_. This is a temperature-gravity grid of low-resolution intensity spectra (based on the solar-metallicity :ads_citet:`Castelli:2003` atmospheres). It's shipped with the MSG installation in `$MSG/data/grids/sg-demo.h5`.
* `pb-Generic-Johnson.B-Vega` from `MSG`_. This is a Johnson B passband file from the `table of tar archives <http://www.astro.wisc.edu/~townsend/resource/docs/msg/appendices/passband-files.html>`_ built using the `Spanish Virtual Observatory <https://svo.cab.inta-csic.es/main/index.php>`_ filter and calibration database.

******************************
Importing the GYRE-lc Module
******************************

This walkthrough relies on `MSG`_ for rapid synthesis of photometric data. Download and install MSG, then set the :envvar:`MSG_DIR` environment variable as described in the `MSG Quick Start guide <http://www.astro.wisc.edu/~townsend/resource/docs/msg/user-guide/quick-start.html#quick-start>`_. 

To use GYRE-lc in Python, also make sure the :envvar:`GYRELC_DIR` environment variable is set (see `Quick Start`). I use a Jupyter notebook for this walkthrough, but you may later choose to write a Python script instead should it better suit your workflow.

In a new Jupyter notebook, copy and past the following imports::

    # Import standard modules

    import numpy as np
    import sys
    import os

    # Import pymsg

    MSG_DIR = os.environ['MSG_DIR']
    sys.path.insert(0, os.path.join(MSG_DIR, 'python'))
    import pymsg

    # Import gyrelc modules

    sys.path.insert(0, os.path.join(os.environ['GYRELC_DIR'], 'src'))
    import gyrelc as lc

The :py:mod:`pymsg` and :py:mod:`gyrelc` modules both require :py:mod:`sys` and :py:mod:`os` to be imported, so we do that first. We also import the :py:mod:`numpy` module, which we use extensively.
At this point, you may also import and configure any plotting or visualization modules::

    # Import plotting module and configure
    
    import matplotlib.pyplot as plt
    %matplotlib inline
    plt.rcParams.update({'font.size': 16})

We must now create a photometric grid.

Creating a PhotGrid
=========================

With `sg-demo.h5` and `pb-Generic-Johnson.B-Vega.h5` in the current working directory::
    
    pg = pymsg.PhotGrid('sg-demo.h5', 'pb-Generic-Johnson.B-Vega.h5')

Modeling the "heartbeat"
=========================

Next, create a pair of :py:class:`gyrelc.Star` objects using the stellar and tide models provided::

    # Create Star objects
    Aa = lc.Star(mesa_model='iOri-Aa.mesa', gyre_model='iOri-Aa-response.h5', photgrid=pg)
    Ab = lc.Star(mesa_model='iOri-Ab.mesa', gyre_model='iOri-Ab-response.h5', photgrid=pg)

Use them, along with the corresponding orbital parameters, as inputs to create a :py:class:`gyrelc.Binary` object::

    # Create Binary object
    iOri = lc.Binary(Aa, Ab, a=132., e=0.764, omega_orb=0.03432)

Now create an ``Observer`` object::

    # Create an Observer object
    inc = 44.0
    omega = 112.5
    
    # Create an Observer object
    obs = lc.Observer(iOri, inc, omega)

The ``Binary`` object consists of two ``Star`` objects, an ``Irradiation`` object, as well as the various attributes and parameters required to provide the ``Observer`` object sufficient context to synthesize a light curve. The ``Observer`` object primarily contains functions for light curve synthesis and analysis thereof. The last parameter left to specify, the choice of passband, is left as an argument for the ``Observer`` class.

Finally, create a light curve::

    # Duration of 'observation' and number of points
    omega_orb = iOri.omega_orb
    t = np.linspace(0/omega_orb, 1/omega_orb, num=2000, endpoint=False)

    # Evaluate the flux
    flux = obs.find_flux(t, t_peri=0.25/omega_orb)

An important subtlety: the ``find_flux()`` function *requires* the observation time to be in units of the orbital period. Here, I'm simulating a BRITE-B passband observation of :math:`{\iota}` iOri that consists of 2000 data points over 2 orbital periods, begining at half a period past periastron. I chose the specific time of periastron for presentation purposes only. 

Using :py:mod:`matplotlib`, you may plot your results::

    # Plot

    fig, ax = plt.subplots(figsize=(8,4))

    legend_style = {'framealpha':1.0, 'handlelength':1.2, 'handletextpad':0.5, 'fontsize':'small'}

    ax.plot(t*omega_orb, flux, lw=1, label='B')
    ax.legend(loc=1, **legend_style)

    ax.set_xlim(0,1)

    ax.set_title(f'$\iota$Ori light curve, $\omega$={omega}')

    fig.text(0.01, 0.5, r'Mode Flux Perturbation', va='center', rotation='vertical')
    fig.text(0.5, 0.0, f'phase (P={1./omega_orb:4.4f} d)', ha='center')

The legend style and labels are entirely a matter of stylistic choice, but a plot with this *xlim* should look something like this:

.. image:: ./walkthrough-lightcurve.png

.. note:: This project is under active development.

.. rubric:: Footnotes:
