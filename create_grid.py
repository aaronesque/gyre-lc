#!/usr/bin/env python
#
#

import h5py
import random
import numpy as np

# Set up axes

Teffs = [2500., 5000., 7500., 10000., 15000., 20000., 25000., 30000., 40000., 50000.]
loggs = [2.5, 3.0, 3.5, 4.0, 4.5]

# Create function for testing

def data_function(Teff, logg):

    #return random.random()

    return (Teff**2)*(logg**2)

# Create nodes

for Teff in Teffs:
    for logg in loggs:

        # Decide whether to create a node

        Gamma = Teff**4/10**logg/1E15

        if Gamma < 1.:

            filename = 't{:05d}g{:03d}.h5'.format(int(Teff), int(logg*1e2))

            f = h5py.File(filename, 'w')

            f.attrs['Teff'] = Teff
            f.attrs['logg'] = logg

            data_list = [data_function(Teff,logg)]

            f.create_dataset('data', data=np.array(data_list))

            f.close()
            
        
        

        
