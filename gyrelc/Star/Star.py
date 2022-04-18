# Import standard modules
import numpy as np
import sys
import os

# Import special modules
import h5py
from scipy.special import lpmv, sph_harm
from scipy.integrate import quad as integrate
from scipy.optimize import fsolve
from astropy.io import ascii


### Class definitions

class resp_coeffs:
    
    def __init__ (self, response_file):
        
        self.data = self.read_response(response_file)
        
        # Sanity check
        
        if len(self.data) == 0:
            raise Exception('Empty response file')
        if not isinstance(self.data, dict):
            raise Exception('Response file type error')
            
        # Setup axes
        
        n_l = self.data['l_max']
        n_m = 2*n_l
        n_k = self.data['k_max']
        
        self.R_lmk = np.empty([n_l+1, n_m+1, n_k+1], dtype=complex)
        self.T_lmk = np.empty([n_l+1, n_m+1, n_k+1], dtype=complex)
        self.G_lmk = np.empty([n_l+1, n_m+1, n_k+1], dtype=complex)
        
        for l in range(2, n_l+1):
            for m in range(-l, l+1):
                for k in range(0, n_k+1):
                    
                    i_l = l
                    i_m = m + n_l
                    i_k = k

                    self.R_lmk[i_l,i_m,i_k] = self.find_Delta_R(i_l,i_m,i_k)
                    self.T_lmk[i_l,i_m,i_k] = self.find_Delta_T(i_l,i_m,i_k)
                    self.G_lmk[i_l,i_m,i_k] = self.find_Delta_G(i_l,i_m,i_k, m,k)
        
        
    def read_response(self, filename):
    
        # Read data from gyre_response
        
        if filename == '':
            
            # Fabricate dummy data
            
            xi_r_ref_re = np.zeros((0,0,0))
            xi_r_ref_im = np.zeros((0,0,0))
            
            lag_L_ref_re = np.zeros((0,0,0))
            lag_L_ref_im = np.zeros((0,0,0))
            
            l_max = 0
            k_max = 0
            
            Omega_rot = 0.
            Omega_orb = 0.
            
        else:
            
            f = h5py.File(filename, 'r')
            
            xi_r_ref_re = f['xi_r']['re'][...]
            xi_r_ref_im = f['xi_r']['im'][...]
            
            lag_L_ref_re = f['lag_L']['re'][...]
            lag_L_ref_im = f['lag_L']['im'][...]
            
            k_max = f.attrs['k_max']
            l_max = f.attrs['l_max']
            
            Omega_rot = f.attrs['Omega_rot']
            Omega_orb = f.attrs['Omega_orb']
            
            f.close()
        
        return {'xi_r_ref': xi_r_ref_re + 1j*xi_r_ref_im,
                'lag_L_ref': lag_L_ref_re + 1j*lag_L_ref_im,
                'k_max': k_max,
                'l_max': l_max,
                'Omega_rot': Omega_rot,
                'Omega_orb': Omega_orb}
    
    
    def find_Delta_R(self, i_l,i_m,i_k):
    
        return np.sqrt(4.*np.pi) * self.data['xi_r_ref'][i_k,i_m,i_l]
        
        
    def find_Delta_T(self, i_l,i_m,i_k):
            
        xi_r_ref = self.data['xi_r_ref'][i_k,i_m,i_l]
        
        lag_L_ref = self.data['lag_L_ref'][i_k,i_m,i_l]
        
        return np.sqrt(4.*np.pi)*(lag_L_ref - 2*xi_r_ref)/4.
      
    
    def find_Delta_G(self, i_l,i_m,i_k, m,k):
    
        xi_r_ref = self.data['xi_r_ref'][i_k,i_m,i_l]
        
        omega = -k*self.data['Omega_orb'] - m*self.data['Omega_rot']
        
        return np.sqrt(4*np.pi)*(-omega**2 - 2)*xi_r_ref

###

class atm_coeffs:
    
    def __init__ (self, intensity_file):
        
        # input to be replaced by Teff, logg
        
        self.data = self.read_intensity(intensity_file)
        
        self.info = self.read_info(intensity_file)
        
        self.R_xl = {}
        self.T_xl = {}
        self.G_xl = {}
        
        for x in self.info:
            
            self.R_xl[x] = self.find_coeffs('R',x)
            self.T_xl[x] = self.find_coeffs('T',x)
            self.G_xl[x] = self.find_coeffs('G',x)
    
    
    def read_intensity(self, filename):
        
        return h5py.File(filename, 'r')
    
    
    def read_info(self, filename):
        
        moments = self.read_intensity(filename)
        
        colors = []
        ells = []

        for moment in list(moments):
            words = moment.split('_')
            
            if len(words)==2:
                if words[-1] not in colors:
                    colors.append(words[-1])
                    
            if len(words)==3:
                if int(words[-1]) not in ells:
                    ells.append(int(words[-1]))
                    
        l_min, l_max = min(ells), max(ells)
        
        info = {}
        
        for color in colors:
            info[color] = {'l_min':l_min, 'l_max':l_max}
        
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
            
        C_x = {} 
        
        if l==None:

            n_l = self.info[x]['l_max'] # - self.info[x]['l_min'] 
            #this may break if l_min =/= 0
        
            for l in np.arange(n_l+1,dtype=int):
                C_x[l] = find_coeffs_C(I,C,x,l)[0]
            return C_x
        
        elif isinstance(l,int):
            C_x[l] = find_coeffs_C(I,C,x,l)[0]
            return C_x
        
        else:
            raise Exception('TypeError: l must be int or int array')

###
