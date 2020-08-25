import numpy as np
import h5py
from scipy.special import legendre, sph_harm
import f90nml as nml

### Class definition

class atm_coeffs:
    
    def __init__ (self, intensity_file, inlist_file):
        
        self.data = self.read_intensity(intensity_file)
        
        self.info = self.read_inlist(inlist_file)
        
        self.R_xl = {}
        self.T_xl = {}
        self.G_xl = {}
        
        for x in self.info:
            
            self.R_xl[x] = self.find_coeffs('R',x)
            self.T_xl[x] = self.find_coeffs('T',x)
            self.G_xl[x] = self.find_coeffs('G',x)
    
    
    def read_intensity(self, filename):
        
        return h5py.File(filename, 'r')
    
    
    def read_inlist(self, filename):
        
        inlist = nml.read(filename)
        
        info = {}
        
        if isinstance(inlist['color'], list):
            
            for color in inlist['color']:
                info[color['filter']] = {'l_min':color['l_min'], 'l_max':color['l_max'] }
            
        elif isinstance(inlist['color'], type(inlist)):
            
            color = inlist['color']
            info[color['filter']] = {'l_min':color['l_min'], 'l_max':color['l_max'] }
        
        else: raise Exception('ope')
        
        return info
    
            
    def find_coeffs(self, C, x, l=None):
        
        I = self.data
        
        def find_coeffs_C(I, C, x, l):
            
            if C=='R': 
                return (2 + l)*(1 - l)*I[f'I_{x}_{l}'][:] / I[f'I_{x}_0'][:]
        
            if C=='T':
                return I[f'dlnTeff_{x}_{l}'][:] / I[f'I_{x}_0'][:]
        
            if C=='G':
                return I[f'dlng_{x}_{l}'][:] / I[f'I_{x}_0'][:]
        
        if l==None:
        
            n_l = self.info[x]['l_max'] - self.info[x]['l_min']
            C_x = np.empty([n_l+1], dtype=complex)
        
            for l in np.arange(n_l+1,dtype=int):
                C_x[l] = find_coeffs_C(I,C,x,l)
            
            return C_x
        
        elif isinstance(l,int):
            return find_coeffs_C(I,C,x,l)
        
        else:
            raise Exception('TypeError: l must be int or int array')
            
