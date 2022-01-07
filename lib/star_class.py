import numpy as np
import h5py
import f90nml as nml
from scipy.special import sph_harm
from astropy.io import ascii

import sys
import os
sys.path.insert(0, os.path.join(os.environ['MSG_DIR'], 'lib'))
import pymsg

import resp_coeffs as rc
import atm_coeffs as ac

### Class definitions

class Star:
    """Takes path to inlist specifying stellar parameters 
    necessary for plotting a light curve.
    
    Arguments
    ---------
    inlist_path : str
    comp_number : int
    
    Attributes
    ----------
    self.inlist_path : str
    self.params : dict
    self.resp_coeffs : dict
    
    Methods
    -------
    +read_mesa_params(comp_model_path) : dict
    +read_phot_coeffs_h5(filter_x, phot_file) : dict
    +read_phot_coeffs_msg(filter_x) : dict
    +read_phot_coeffs(filter_x, phot_file) : dict
    +D_moment(l, filter_x, deriv) : float
    +R_xl(filter_x, l) : float
    +T_xl(filter_x, l) : float
    +G_xl(filter_x, l) : float
    +disk_intg_factor(filter_x, l) : float
    +eval_fourier_moment(filter_x, theta, phi, l,m,k) : float
    +eval_fourier(filter_x, inc, omega, bol) : [array, array]
    """
    
    def __init__ (self, inlist_path, comp_number=1):
        
        # should check if 'inlist', this just checks 'str'
        
        if isinstance(inlist_path, str):
            self.inlist_path = inlist_path
            self.inlist = nml.read(inlist_path)[f'comp_{comp_number}']
        else: raise Exception('Inlist file error')
        
        if self.inlist['comp_model_type']=='MESA':
            
            self.params = self.read_mesa_params(self.inlist['comp_model_path'])
            self.resp_coeffs = rc.resp_coeffs(self.inlist['tide_model_path']) 
            self.phot_coeffs = {}
        
        elif self.inlist['comp_model_type']=='PT_MASS':
            
            self.params = self.read_pt_mass_params(self.inlist)
            self.resp_coeffs = rc.resp_coeffs('') #self.make_pt_mass_response()
            self.phot_coeffs = {}
            
        else: raise Exception("Invalid comp_model_type must be 'MESA' or 'PT_MASS'.")
        
    

    def read_pt_mass_params(self, comp_inlist):
        
        m = comp_inlist['mass']
        m_units = comp_inlist['mass_units']
         
        M_sol = 1.989e33
        
        if m_units=='CGS':
            m = m/M_sol
            
        return {'M': m, 
                'R': 0,
                'L': 0,
                'Teff': 0,
                'logg': 0}
           


    def make_pt_mass_response(self):

        return {'xi_r_ref': 0+0j,
                'lag_L_ref': 0+0j,
                'k_max': 0.,
                'l_max': 0.,
                'Omega_rot': 0.,
                'Omega_orb': 0.}



    def read_mesa_params(self, comp_model_path, units='SOLAR'):
        
        data = ascii.read(comp_model_path, data_start=0, data_end=1)
        
        # cgs constants
        
        R_sol = 6.957e10 
        M_sol = 1.989e33
        L_sol = 3.839e33
        G = 6.674079999999999e-08 
        sigma_sb = 5.6703669999999995e-05
        
        # data in cgs
        
        m = data[0][1]
        r = data[0][2]
        l = data[0][3]
            
        # calculate Teff [K], logg [dex]
        
        Teff = (l/(sigma_sb * 4*np.pi*r**2))**0.25
        
        g_surf = G*m/(r**2)
        
        logg = np.log10(g_surf)
        
        # convert to solar units
        
        if units=='SOLAR':
            m = m/M_sol
            r = r/R_sol
            l = l/L_sol
            
        return {'M': m,
                'R': r,
                'L': l,
                'Teff': Teff,
                'logg': logg}

    
    
    def read_phot_coeffs_h5(self, filter_x, phot_file):
        """
        This function is old and uses a synspec
        intensity spectrum calculated from iOri's 
        tlusty atmosphere.
        """
        
        l_min = 0
        dl = 1
        l_max = self.resp_coeffs.data['l_max']
        n_l = np.ceil( (l_max - l_min)/dl ) + 1

        l_range = l_min + dl*np.arange(n_l)
        
        # read intensity file
        
        I = ac.atm_coeffs(phot_file)
        I = I.data
        
        # Evaluate photometric data
        I_x_l = {}
        dI_dlnT_x_l = {}
        dI_dlng_x_l = {}
        
        for l in l_range.astype(int):
            I_x_l[l] = I[f'I_{filter_x}_{l}'][0]
            dI_dlnT_x_l[l] = I[f'dlnTeff_{filter_x}_{l}'][0]
            dI_dlng_x_l[l] = I[f'dlng_{filter_x}_{l}'][0]
            
        return {f'I_{filter_x}': I_x_l,
                f'dI_dlnT_{filter_x}': dI_dlnT_x_l,
                f'dI_dlng_{filter_x}': dI_dlng_x_l}

    
    def read_phot_coeffs_msg(self, filter_x): 
        
        # Set atmosphere parameters dict
    
        Teff = self.params['Teff']
        logg = self.params['logg']

        dx = {'logT': np.log10(Teff), 'logg': logg}
        pg = pymsg.PhotGrid(f"{os.environ['GYRELC_DIR']}/grid/{filter_x}.h5")
        
        # Set up intensity moment range
        
        l_min = 0
        dl = 1
        l_max = self.resp_coeffs.data['l_max']
        n_l = np.ceil( (l_max - l_min)/dl ) + 1

        l_range = l_min + dl*np.arange(n_l)
        
        # Evaluate photometric data
        I_x_l = {}
        dI_dlnT_x_l = {}
        dI_dlng_x_l = {}
        
        for l in l_range.astype(int):
            I_x_l[l] = pg.D_moment(dx, l)
            dI_dlnT_x_l[l] = pg.D_moment(dx, l, deriv={'logT':True})/np.log(10.)
            dI_dlng_x_l[l] = pg.D_moment(dx, l, deriv={'logg':True})/np.log(10.)
            
        return {f'I_{filter_x}': I_x_l, 
                f'dI_dlnT_{filter_x}': dI_dlnT_x_l, 
                f'dI_dlng_{filter_x}': dI_dlng_x_l}
    
    
    def make_phot_coeffs_pt_mass(self, filter_x):
        return {f'I_{filter_x}': np.array([0.]), 
                f'dI_dlnT_{filter_x}':np.array([0.]), 
                f'dI_dlng_{filter_x}': np.array(0.)}
        
    
    def read_phot_coeffs(self, filter_x, phot_file=None):
        
        if self.inlist['comp_model_type']=='MESA':
            if phot_file==None:
                self.phot_coeffs.update( self.read_phot_coeffs_msg(filter_x) )
            else:
                self.phot_coeffs.update( self.read_phot_coeffs_h5(filter_x, phot_file) )
                
        elif self.inlist['comp_model_type']=='PT_MASS':
            self.phot_coeffs.update( self.make_phot_coeffs_pt_mass(filter_x) ) 
            
        return
    
    
    def D_moment(self, l, filter_x, deriv=None):
        I = self.phot_coeffs
        if deriv=='lnT':
            return I[f'dI_dlnT_{filter_x}'][l]
        elif deriv=='lng':
            return I[f'dI_dlng_{filter_x}'][l]
        elif deriv=='logT':
            return I[f'dI_dlnT_{filter_x}'][l] * np.log(10.)
        elif deriv=='logT':
            return I[f'dI_dlng_{filter_x}'][l] * np.log(10.)
        else:
            return I[f'I_{filter_x}'][l]
    
    def R_xl(self, filter_x, l):
        I = self.phot_coeffs
        return (2 + l)*(1 - l)*I[f'I_{filter_x}'][l] / I[f'I_{filter_x}'][0]
    
    def T_xl(self, filter_x, l):
        I = self.phot_coeffs
        return I[f'dI_dlnT_{filter_x}'][l] / I[f'I_{filter_x}'][0] 
    
    def G_xl(self, filter_x, l):
        I = self.phot_coeffs
        return I[f'dI_dlng_{filter_x}'][l] / I[f'I_{filter_x}'][0]
    
    def disk_intg_factor (self, filter_x, l):
        I = self.phot_coeffs
        b_l = I[f'I_{filter_x}'][l]/I[f'I_{filter_x}'][0]
        return b_l
    
    
    def eval_fourier_moment (self, filter_x, theta, phi, l,m,k):
        
        i_l = l
        i_m = m + self.resp_coeffs.data['l_max']
        i_k = k

        dR_lmk = self.resp_coeffs.R_lmk[i_l,i_m,i_k]
        dT_lmk = self.resp_coeffs.T_lmk[i_l,i_m,i_k]
        dG_lmk = self.resp_coeffs.G_lmk[i_l,i_m,i_k]
                
        Y_lm = sph_harm(m, l, phi, theta)
                
        R_xl = self.R_xl(filter_x, l)
        T_xl = self.T_xl(filter_x, l) 
        G_xl = self.G_xl(filter_x, l)
        
        if k == 0:
            if m == 0:
                kappa = 0.5
            elif m >= 1:
                kappa = 1.
            else:
                kappa = 0.
        else:
            kappa = 1.
            
        return 2*kappa*(dR_lmk*R_xl + dT_lmk*T_xl + dG_lmk*G_xl)*Y_lm
    

    def eval_fourier (self, filter_x, theta, phi):
        
        if self.inlist['comp_model_type']=='PT_MASS':
            f, A = np.array([0]), np.array([0])
        else:
            resp_coeffs = self.resp_coeffs
            I = self.phot_coeffs
            
            # Initialize the frequencies/amplitudes arrays
            
            f = np.arange(resp_coeffs.data['k_max']+1)*resp_coeffs.data['Omega_orb']
            
            A = np.zeros(resp_coeffs.data['k_max']+1, dtype=complex)
                
            # Loop over l, m and k
        
            for l in np.arange(2, resp_coeffs.data['l_max']+1).astype(int):
                for m in np.arange(-l, l+1).astype(int):
                    for k in np.arange(0, resp_coeffs.data['k_max']+1).astype(int):
        
                        # Add the Fourier contribution * spherical harmonic
            
                        A[k] += self.eval_fourier_moment(filter_x, theta, phi, l,m,k)
                    
        # Return data
        return f, A  

