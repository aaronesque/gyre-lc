#!/usr/bin/env python
#
# Sketch for gyre-lc flow

import sys
import numpy as np

from star_func import *
import bin_func as bi
import obs_func as ob

if __name__ == '__main__':

    bin_iori = bi.binary('bin_params.in')
    obs = ob.observer(bin_iori)
    omega_orb = bin_iori.orbit['Omega_orb']

    # Evaluate fourier terms

    inc = 82.9
    omega = 122.2
    x = 'BRITE-B'

    t = np.linspace(-0.2/omega_orb, 2.2/omega_orb, 2001)

    flux_1 = obs.find_flux(inc, omega, 'BRITE-B', t)
    flux_2 = obs.find_flux(inc, omega, 'BRITE-R', t)
