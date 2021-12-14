#!/usr/bin/env python

import grid as gr
import numpy as np
import sys
    

# Create function for testing

def func (Teff, logg):
    u = 0.4
    sigma = 5.670374E-5
    l = [ 0.5 - u/6, 1/3 - u/12, 1/8 + u/120, u/24, -1/48 + u/48 ]
    return ( (sigma*Teff**4)/(2*np.pi*(0.5 - u/6)) )*sum(l)
    #return Teff*logg**3 + logg**2

def bound_func (Teff, logg):
    return Teff**4/10**logg/1E15 < 1.

            
if __name__ == '__main__':

    # Set up axes

    Teff_axis = [2500., 5000., 7500., 10000., 15000., 20000., 25000., 30000., 40000., 50000.]
    logg_axis = [2.5, 3.0, 3.5, 4.0, 4.5]

    # Build the grid

    grid = gr.from_func(Teff_axis, logg_axis, func, bound_func=bound_func, debug=True)

    grid.show_topology()

