.. _how-gyre-lc-works:

.. gyre-lc documentation master file, created by

#############################
How GYRE-lc Works
#############################

Go :footcite:t:`2017:Pablo`. 

.. GYRE-tides models forced oscillations of a star in a binary due to its companion's gravitational field :footcite:t:`2021:Sun`. As input for one such calculation, GYRE-tides takes a stellar model produced with `MESA <mesa.sourceforge.net>`_ and applies a forcing potential calculated via user-specified binary parameters (see `inputs`).

.. Most hb stars still have amplitudes that are small. Even though they're dramatic, we can still say they're small.

The forcing potential :math:`{\Phi_S}` can be written as an expansion of the gravitational potential at a point on the star's surface into spherical harmonics:

.. math::    

    \Phi_S (\vec{r}; t) &= \frac{-q G M}{|\vec{r} - \vec{r}_S|} \\
    &= \sum^\infty_{l=0} \sum^l_{m=-l} \sum^\infty_{k=-\infty} \Phi_{r;l,m,k}(r) \; \Yml(\theta, \phi) \; e^{-i k \Omega_\textrm{orb} t}


.. Here, :math:`{\Phi_{r;l,m,k}}` is the radial component of the forcing potential amplitude, and :math:`{\Y^m_l}` is the spherical harmonic of order $m$ and degree $l$.  The exponential term is the $k$-th Fourier harmonic. Restricting ourselves to small amplitude tides allows us to write the response perturbation as a superposition of many different partial tides:

..    \xi_r(\vec{r}; t) = \sum_{l,m,k} \tilde{\xi}_{r; l,m,k}(r) \; Y^m_l (\theta, \phi) \; e^{-i k \Omega_\textrm{orb} t}

.. It follows from T03 (see \sref{formalism}) that we may also expand the radiative luminosity that way into surface luminosity variations:

    \delta L(\vec{r};t)_\textrm{rad} = \widetilde{\delta L}_{\textrm{rad};l,m,k}(r) \; Y^m_l \; e^{-i k \Omega_\textrm{orb} t }

.. It behooves us to probe the practical limitations of this approach. The net tidal force can be characterized by the tidal strength term
    \epsilon_\mathrm{T} \equiv \left( \frac{R}{a} \right)^3 = \frac{R^3 \Omega_\textrm{orb}^2}{GM}\times \left( \frac{q}{1+q} \right).
.. For small amplitude tides, $\epsilon_\mathrm{T} << 1$.


For wide binaries, this assertion easily holds as long as the primary's radius $R$ is much smaller than the semimajor axis $a$. For some highly eccentric binaries on the other hand, such as eccentric ellipsoidal variables, a small mass ratio $q=M_2/M$ between the secondary and primary stars might be a good enough diagnostic.  We will probe the edge of where our 'weak tides' approach breaks down in a future work.

.. GYRE-tides calculates the tide model, i.e. the partial tide amplitudes $\tilde{\xi}_{r;l,m,k}(R)$ and surface luminosity variations $\widetilde{\delta L}_{\textrm{rad};l,m,k}(R)$, and writes them to file. A corresponding tide model is then created for the companion's neighbor. Both tide models, along with their corresponding stellar models, are the 4 files required to build a single light curve using GYRE-LC.


.. \subsubsection{The semi-analytical formalism} \label{formalism}

.. The semi-analytical formalism for light variations due to tides extends earlier treatments of tides by \citet{Stamford_1981} and \citet{Watson_1988} to include the effects of the Coriolis force within the 'traditional approximation of rotation' (TAR; see, e.g. \citealt{Bildsten_1996}; \citealt{Lee_1997}; \citealt{Townsend_2003b}; and references therein). This is important because the Coriolis force can act as a waveguide confining oscillations to the equator. This phenomenon may significantly impact a star's observed variability, yet it has not been accounted for in previous studies of eccentric ellipsoidals.

.. Essentially, the semi-analytical formalism makes the statement that, for any stellar surface perturbation that can be written as a superposition of partial perturbations, we can write the resulting light variations in terms of intensity moments. 

.. Accordingly, we can express perturbations to stellar radius $R$, effective temperature $T_\mathrm{eff}$, and surface gravity $g_\eff$ in terms of spherical harmonics $Y_l^m(\theta, \phi)$ and perturbation coefficients:

    \Delta_R &= \frac{\tilde{\xi}_r(R)}{R}\\
    \Delta_{T_\eff} &= \frac{1}{4} \left( \frac{\widetilde{\delta L}_\mathrm{rad}(R)}{L_\mathrm{rad}(R)} - 2 \frac{\tilde{\xi}_r(R)}{R} \right)\\
    \omega &= -k\Omega_{orb} - m\Omega_{rot} \\
    \Delta_{g_\eff} &= (-\omega^2 - 2)\xi_{r_\mathrm{ref}}


    \frac{\delta R}{R} (\theta, \phi; t) &= \mathrm{Re} \left[ \Delta_R Y_l^m(\theta, \phi) e^{\ii \sigma t} \right] \\
    \frac{\delta T_\eff }{T_\eff } (\theta, \phi; t) &= \mathrm{Re} \left[ \Delta_T Y_l^m(\theta, \phi) e^{\ii \sigma t} \right] \\
    \frac{\delta g_\eff}{g_\eff} (\theta, \phi; t) &= \mathrm{Re} \left[ \Delta_g Y_l^m(\theta, \phi) e^{\ii \sigma t} \right] 

.. Therefore, perturbations $\delta \FF_{\lx}$ to the stellar flux $\FF_{\lx}$ in some photometric passband $x$ are modeled via the differential flux functions $\{ \TT^m_{\lx}, \GG^m_{\lx},
.. \RR^m_{\lx} \}$, which depend on intensity moments $\II_{\lx}$:

.. \frac{\delta \FF_{\lx}}{\FF_{\lx}} (\theta_o, \phi_o; t) &= \mathrm{Re} \left[ \left\{ \Delta_R \RR^m_{\lx}(\theta_o, \phi_o) + \Delta_T \TT^m_{\lx}(\theta_o, \phi_o) + \Delta_g \GG^m_{\lx}(\theta_o, \phi_o) \right\} e^{\ii \sigma t} \right] \\
.. \RR^m_{\lx}(\theta_o,\phi_o) &\equiv \frac{(2+\ell)(1-\ell)}{\II_{0;x}} \II_{\lx} Y^m_l (\theta_o, \phi_o) \\
.. \TT^m_{\lx}(\theta_o,\phi_o) &\equiv \frac{1}{\II_{0;x}} \frac{ \partial \II_{\lx}}{\partial \ln{ T_\eff}} Y^m_l (\theta_o, \phi_o) \\
.. \GG^m_{\lx}(\theta_o,\phi_o) &\equiv\frac{1}{\II_{0;x}} \frac{ \partial \II_{\lx}}{\partial \ln{g}} Y^m_l (\theta_o, \phi_o). \\
.. \II_{\lx} &= \int_0^1 \mu P_l(\mu)\II_x(\mu) d\mu

.. Here, $\II_x(\mu)$ is the specific intensity in passband $x$, emergent from the stellar atmosphere at cosinus $\mu$ from the surface normal, and $P_\ell(\mu)$ is the Legendre polynomial of degree $\ell$. The perturbation coefficients can be retrieved from the GYRE-tides output through algebra.

The photospheric data required to compute the specific intensities is provided by the spectral synthesis code for stars, MSG. A brief overview of its limitations and functionality follows.


.. note:: This project is under active development.

.. footbibliography::
