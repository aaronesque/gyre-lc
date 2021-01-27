#!/usr/bin/env python
#
# wrapper for obs_func.py

import sys
import numpy as np
import f90nml

from obs_func import observer

####

def run_observer(bin_list):

    # Setup observer

    obs = observer(bin_list)

    inc = obs.nml['observer']['inc']
    omega = obs.nml['observer']['omega']
    x = obs.nml['observer']['bandpass']

    t_units = 1/obs.nml['orbit']['Omega_orb']

    t_start = obs.nml['observer']['t_start']*t_units
    t_end = obs.nml['observer']['t_end']*t_units
    n = obs.nml['observer']['n']

    t_ = np.linspace(t_start, t_end, n)

    # Evaluate fourier terms

    return obs.find_flux(inc, omega, x, t_), t_/t_units


####

if __name__ == '__main__':
    
    if len(sys.argv)==2:

        filename = sys.argv[1]
        flux, interval = run_observer(filename)

        print('flux [10^-3]:', flux*10e3)
        print('interval [Period]:', interval)
    
    else:
        raise Exception('Invalid number of inputs. Run as  "./gyre-lc.py [inlist]"')
