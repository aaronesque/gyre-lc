#!/usr/bin/env python

import grid as gr
import h5py
import numpy as np
import sys
            
def func (Teff, logg, l):
    u = 0.4
    sigma = 5.670374E-5
    P = [ 0.5 - u/6, 1/3 - u/12, 1/8 + u/120, u/24, -1/48 + u/48 ]
    return ( (sigma*Teff**4)/(2*np.pi*(0.5 - u/6)) )*P[l]
    #return Teff*logg**3 + logg**2

def df_dlnTeff (Teff, logg, l):
    u = 0.4
    sigma = 5.670374E-5
    P = [ 0.5 - u/6, 1/3 - u/12, 1/8 + u/120, u/24, -1/48 + u/48 ]
    return ( (4*sigma*Teff**3)/(2*np.pi*(0.5 - u/6)) )*P[l]/func(Teff, logg, l)

def df_dlng (Teff, logg, l):
    return 0


def bound_func (Teff, logg):
    return Teff**4/10**logg/1E15 < 1.


if __name__ == '__main__':

    # Set up axes

    #Teff_axis = [30890] 
    #logg_axis = [3.645] 
    Teff_axis = [29173] 
    logg_axis = [4.232] 
    #Teff_axis = [2500., 5000., 7500., 10000., 15000., 20000., 25000., 30000., 40000., 50000.]
    #logg_axis = [2.5, 3.0, 3.5, 4.0, 4.5]

    # Create intensity grid

    for Teff in Teff_axis:
        for logg in logg_axis:
            if bound_func(Teff, logg):
                
                fname = f't{int(Teff):05d}g{int(logg*100):02d}.h5'
                f = h5py.File(fname, "w")
                f.attrs['Teff'] = Teff
                f.attrs['logg'] = logg
                f.attrs['u'] = 0.4

                l_range = [0, 1, 2, 3, 4]
                
                I_tot = sum( [func(Teff, logg, l) for l in l_range] )
                dset = f.create_dataset("I_bol", data = [I_tot] )

                for l in l_range:
                    dset = f.create_dataset(f"I_bol_{l}", data = [func(Teff, logg, l)] )
                    dset = f.create_dataset(f"dlnTeff_bol_{l}", data = [df_dlnTeff(Teff, logg, l)] )
                    dset = f.create_dataset(f"dlng_bol_{l}", data = [df_dlng(Teff, logg, l)] )


                f.close()

    # Build the grid

