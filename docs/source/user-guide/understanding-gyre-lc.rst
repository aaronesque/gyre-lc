.. _understanding-gyre-lc:

.. gyre-lc documentation master file, created by

#############################
Understanding GYRE-lc
#############################

This chapter provides a deeper look into what the GYRE-lc package does and how it works. The primary function of GYRE-lc is the rapid forward modeling of light curves for tidal pulsators. It is the first light curve synthesizer:

- to be based on the semi-analytical formalism for modeling intensity variations detailed in :ads_citet:`Townsend:2003b`;
- to derive photospheric data from the spectral synthesis code MSG;
- to use the tidal response output from GYRE-tides :ads_citep:`Sun:2021`.

Companion irradiation is optionally modeled using first order approximations according to :ads_citet:`Burkart:2012`, but any flux variation due to eclipsing is ignored. The process for light curve modeling with GYRE-lc involves 3 major steps: 

- :ref:`step 1` for each component of a binary using GYRE-tides.
- :ref:`step 2` through GYRE-lc's deployment of the :ads_citet:`Townsend:2003b` formalism using MSG for photometric intensity moments. *Optional:* GYRE-lc deploys the :ads_citet:`Burkart:2012` formalism for irradiation contributions to the flux.
- :ref:`step 3` with GYRE-lc.

A short summary of each step follows.

.. _step 1:

******************************
Step 1: Find the Partial Tides
******************************

GYRE-tides models forced oscillations of a star in a binary due to its companion's gravitational field :ads_citep:`Sun:2021` as *partial tides*. As input for one such calculation, GYRE-tides takes a stellar model produced with `MESA <mesa.sourceforge.net>`_ and applies a forcing potential calculated via user-specified binary parameters (see :doc:`Python Walkthrough <python-walkthrough>`).

The forcing potential :math:`{\Phi_S}` can be written as an expansion of the gravitational potential at a point on the star's surface into spherical harmonics:

.. math:: 
   \Phi_S (\vec{r}; t) &= \frac{-q G M}{|\vec{r} - \vec{r}_S|} \\
   &= \sum^\infty_{\ell=0} \sum^l_{m=-\ell} \sum^\infty_{k=-\infty} \Phi_{r;\ell,m,k}(r) \; \Yml (\theta, \phi) \; e^{-i k \Omega_\textrm{orb} t}

Here, :math:`{\Phi_{r;\ell,m,k}}` is the radial component of the forcing potential amplitude, and :math:`{\Yml}` is the spherical harmonic of order :math:`m` and degree :math:`\ell`. The exponential term is the :math:`k`-th Fourier harmonic. 

Restricting ourselves to small amplitude tides allows us to write the response perturbation as a superposition of many different partial tides:

.. math::
   \boxed{ \xi_r(\vec{r}; t) = \sum_{\ell,m,k} \tilde{\xi}_{r; \ell,m,k}(r) \; \Yml (\theta, \phi) \; e^{-i k \Omega_\textrm{orb} t} }

It follows from :ads_citet:`Townsend:2003b` that we may also expand the radiative luminosity that way into partial surface luminosity variations:

.. math::
   \boxed{ \delta L(\vec{r};t)_\textrm{rad} = \widetilde{\delta L}_{\textrm{rad};\ell,m,k}(r) \; \Yml \; e^{-i k \Omega_\textrm{orb} t} }

It behooves us to consider the practical limitations of this approach. The net tidal force can be characterized by the tidal strength term

.. math::
   \epsilon_\mathrm{T} \equiv \left( \frac{R}{a} \right)^3 = \frac{R^3 \Omega_\textrm{orb}^2}{GM}\times \left( \frac{q}{1+q} \right).

For small amplitude tides, :math:`\epsilon_\mathrm{T} << 1`.

For wide binaries, this assertion easily holds as long as the primary's radius :math:`R` is much smaller than the semimajor axis :math:`a`. For some highly eccentric binaries on the other hand, such as eccentric ellipsoidal variables, a small mass ratio :math:`q=M_2/M` between the secondary and primary stars might be a good enough diagnostic.  We will probe the edge of where our 'weak tides' approach breaks down in a future work. For now, users should do their best to stay well-within the weak tide regime.

GYRE-tides calculates the tide model, i.e. the partial tide amplitudes :math:`\tilde{\xi}_{r;\ell,m,k}(R)` and surface luminosity variations :math:`\widetilde{\delta L}_{\textrm{rad};\ell,m,k}(R)`, and writes them to file. A corresponding tide model is then created for the companion's neighbor. Both tide models, along with their corresponding stellar models, are the 4 files required to build a single light curve using GYRE-lc.

.. _step 2:

*************************************
Step 2: Calculate Differential Fluxes
*************************************

In step 2, GYRE-lc represents flux variations due to tides as a linear combination of intensity moments. It does this according to the semi-analytical formalism for light variations described in :ads_citealt:`Townsend:2003b`, which applies to stellar perturbations that can be written as a superposition of partial perturbations–any well-converged GIRE-tides models.

In particular, it states that perturbations to the stellar flux :math:`\delta \FF_{x}` in some photometric passband :math:`x` can be expressed using the *differential flux functions* :math:`\{ \TT^m_{\lx}, \GG^m_{\lx}, \RR^m_{\lx} \}`, which depend on intensity moments :math:`\II_{\lx}`:

.. math::
   :nowrap:

   \begin{align}
   \II_{\lx} &= \int_0^1 \mu P_\ell(\mu)\II_x(\mu) d\mu \\
   \Aboxed{\RR^m_{\lx}(\theta_o,\phi_o) &\equiv \frac{(2+\ell)(1-\ell)}{\II_{0;x}} \II_{\lx} \Yml (\theta_o, \phi_o)} \\
   \Aboxed{\TT^m_{\lx}(\theta_o,\phi_o) &\equiv \frac{1}{\II_{0;x}} \frac{ \partial \II_{\lx}}{\partial \ln{ T_\eff}} \Yml (\theta_o, \phi_o)} \\
   \Aboxed{\GG^m_{\lx}(\theta_o,\phi_o) &\equiv\frac{1}{\II_{0;x}} \frac{ \partial \II_{\lx}}{\partial \ln{g}} \Yml (\theta_o, \phi_o)} \\
   \frac{\delta \FF_{\lx}}{\FF_{\lx}} (\theta_o, \phi_o; t) &= \mathrm{Re} \left[ \left\{ \Delta_R \RR^m_{\lx}(\theta_o, \phi_o) + \Delta_T \TT^m_{\lx}(\theta_o, \phi_o) + \Delta_g \GG^m_{\lx}(\theta_o, \phi_o) \right\} e^{-\ii \sigma t} \right]
   \end{align}

Here, :math:`\II_x(\mu)` is the specific intensity in passband :math:`x`, emergent from the stellar atmosphere at cosinus :math:`\mu` from the surface normal, and :math:`P_\ell(\mu)` is the Legendre polynomial of degree :math:`\ell`. :math:`\theta_o` and :math:`\phi_o` are calculated from the inclination to the observer and the argument of periastron, and :math:`\sigma=k\Omega_{orb}` is the  This math is handled by the :py:class:`Photosphere` class. On the other hand, the *perturbation coefficients* :math:`\Delta` are be retrieved from the GYRE-tides output through algebra in the :py:class:`Response` class:

.. math::
    \Delta_R &= \frac{\tilde{\xi}_r(R)}{R}\\
    \Delta_{T_\eff} &= \frac{1}{4} \left( \frac{\widetilde{\delta L}_\mathrm{rad}(R)}{L_\mathrm{rad}(R)} - 2 \frac{\tilde{\xi}_r(R)}{R} \right)\\
    \Delta_{g_\eff} &= (-\omega^2 - 2) \Delta_R

with :math:`\omega = -k\Omega_{orb} - m\Omega_{rot}` in the co-rotating frame.

Each :py:class:`Star` class object includes methods to evaluate the differential flux functions and perturbation coefficients on the fly by calling :py:class:`Response` and :py:class:`Photosphere` (see :ref:`Fig. 1 <class-diagram-star>`). The photospheric data that :py:class:`Photosphere` requires to return specific intensities is provided by MSG, which is why the :py:class:`pymsg.PhotGrid` must be passed to each :py:class:`Star` upon instantiation.

.. _class-diagram-star:

.. figure:: ./class-diagram-star.png
    :width: 50%
    :align: center

    Fig. 1. A :py:class:`Star` instance calls :py:class:`Response` and :py:class:`Photosphere` methods as needed.

.. .. math::
..    \frac{\delta R}{R} (\theta, \phi; t) &= \mathrm{Re} \left[ \Delta_R Y_l^m(\theta, \phi) e^{\ii \sigma t} \right] \\
..    \frac{\delta T_\eff }{T_\eff } (\theta, \phi; t) &= \mathrm{Re} \left[ \Delta_T Y_l^m(\theta, \phi) e^{\ii \sigma t} \right] \\
..    \frac{\delta g_\eff}{g_\eff} (\theta, \phi; t) &= \mathrm{Re} \left[ \Delta_g Y_l^m(\theta, \phi) e^{\ii \sigma t} \right] 

.. _optional:

***************
Irradiation
***************

Burkart's irradiation formalism describes the additional emergent flux from a stellar atmosphere that is due to radiative heating from an orbiting companion star. It applies to binaries within the current framework, and is therefore straightforward to implement (as in :ref:`Fig. 2<class-diagram-binary>`). However, Burkart's formalism assumes that all incident radiation is immediately reprocessed at the photosphere and emitted isotropically. This assumption only holds for stars of similar spectral types-- otherwise, we can expect a fraction of the incident radiation to be scattered. The decision to include flux contributions from irradiation (with or without a scatter coefficient) is left to the user.

Regardless, the :py:class:`Binary` class inherits methods from :py:class:`Irradiation` so the user may experiment with irradiation from either :py:class:`Star` freely.

.. _class-diagram-binary:

.. figure:: ./class-diagram-binary.png
    :width: 50%
    :align: center   

    Fig. 2. A :py:class:`Binary` instance contains 2 :py:class:`Star` instances and inherits methods from :py:class:`Irradiation`. 

.. _step 3:

*****************************
Step 3: Build the Light Curve
*****************************

Finally, you have everything you need to build the light curve. The :py:class:`Observer` class exists to facilitate user production of flux and other desirables (see :ref:`Fig. 3<class-diagram-observer>`), and must be passed a :py:class:`Binary`, the inclination, and the argument of periastron of the binary with respect to the observer. 

You may then pass a time or timeseries array to :py:func:`Observer.find_flux()`, which returns the sum of differential fluxes calculated from the intensity moments and perturbation coefficients provided by :py:class:`Star` and :py:class:`Irradiation` from within :py:class:`Binary`. 

.. _class-diagram-observer:

.. figure:: ./class-diagram-observer.png
    :width: 50%
    :align: center

    Figure 3. The :py:class:`Observer` class provides the user with methods for building the light curve, examining Fourier coefficients, and more.

Fig. 4 shows a class diagram representation of GYRE-lc's overall architecture, omitting some technical details like most private (mangled) methods and attributes. 

.. _class-diagram-architecture:

.. figure:: ./class-diagram-architecture.png
    :align: center

    Figure 4. GYRE-lc's architecture is roughly meant to imply a "zooming out" from the photosphere all the way out to the observer.

.. note:: This project is under active development.

.. rubric:: Footnote
