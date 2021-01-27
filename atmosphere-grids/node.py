# Node class

import h5py

### Class definition

class Node:

    def __init__ (self, Teff, logg, data):

        # Initialize values

        self.Teff = Teff
        self.logg = logg

        self.data = data

        self.ddata_dTeff = 0.
        self.ddata_dlogg = 0.
        self.ddata_cross = 0.

### Factory methods
        
def from_file (filename):

    # Read data from file

    f = h5py.File(filename, 'r')

    Teff = f.attrs['Teff']
    logg = f.attrs['logg']

    data = f['data'][...][0]

    f.close()

    # Return a new Node

    return Node(Teff, logg, data)

