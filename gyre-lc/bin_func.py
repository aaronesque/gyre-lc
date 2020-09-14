import numpy as np
import h5py
import f90nml as nml

### Class definition

class binary:
    
     def __init__ (self, bin_info):
        
        self.omega_orb = bin_info['omega_orb']
        self.a = bin_info['a']
        self.e = bin_info['e']
        
        self.L1 = bin_info['L1']
        self.R1 = bin_info['R1']
        self.M1 = bin_info['M1']
        
        self.L2 = bin_info['L2']
        self.R2 = bin_info['R2']
        self.M2 = bin_info['M2']
        
        #self.resp_data = resp_data

