#!/usr/bin/env python

import grid as gr
import sys
            
if __name__ == '__main__':

    # Set up axes

    Teff_axis = [2500., 5000., 7500., 10000., 15000., 20000., 25000., 30000., 40000., 50000.]
    logg_axis = [2.5, 3.0, 3.5, 4.0, 4.5]

    # Create function for testing

    def func (Teff, logg):
        return Teff*logg**3 + logg**2

    def bound_func (Teff, logg):
        return Teff**4/10**logg/1E15 < 1.

    # Build the grid

    grid = gr.from_func(Teff_axis, logg_axis, func, bound_func=bound_func, debug=True)

    grid.show_topology()

