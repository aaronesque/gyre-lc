import numpy as np
import f90nml as nml

from scipy.integrate import quad as integrate
from scipy.special import lpmv, sph_harm
from scipy.optimize import fsolve

from star_func import star

### Class definitions

class binary:
    
    def __init__ (self, list_path):
        
        # should check if 'inlist', this just checks 'str'
        
        if isinstance(list_path, str):
            self.nml = nml.read(list_path)
        else: raise Exception('Inlist file error')
        
        self.star = {1: star(list_path, 1),\
                     2: star(list_path, 2) }
        
        self.orbit = self.read_orbit(self.nml)
            
        
    def read_orbit(self, inlist):
        
        orb_list = inlist['orbit']
        
        params = {}
        
        #placeholder unit routine
        
        if orb_list['omega_orb_units']=='CYC_PER_DAY':
            omega_orb_units = 1 
        else: omega_orb_units = 1
        
        if orb_list['a_units']=='CM':
            a_units = 6.957e10 #cm per r_sol
        else: a_units = 1
            
        # eventually, we want to get these from resp_data
        # resp_data currently does not contain this output tho
        
        params['a'] = orb_list['a']*a_units
        params['e'] = orb_list['e']
        params['Omega_orb'] = orb_list['omega_orb']*omega_orb_units
        
        return params

#

class irradiation:
    
    def __init__ (self, bin_data, star_number):
        
        # prep for find_bin_sep()
        
        self.e = bin_data.orbit['e']
        self.a = bin_data.orbit['a']
        self.Omega_orb = bin_data.orbit['Omega_orb']
        
        # prep for eval_irrad()
        
        if int(star_number)==1:
            star_neighbor = 2
        elif int(star_number)==2:
            star_neighbor = 1
        else:
            raise Exception('Star unspecified')
        
        self.atm_data = bin_data.star[star_number].atm.data
        self.res_data = bin_data.star[star_number].res.data
        
        self.L1 = bin_data.star[star_number].par['L']
        self.R1 = bin_data.star[star_number].par['R']
        self.L2 = bin_data.star[star_neighbor].par['L']
        
            
    def eval_ramp (self, l, m):
        
        term1 = 2* np.sqrt( ((2*l+1)/(4*np.pi)) \
                           *(np.math.factorial(l-m)/np.math.factorial(l+m)) )
        
        if (np.abs(m)==1):
            term2 = np.pi/2
        else:
            term2 = np.cos(m*np.pi/2)/(1-m**2)
        
        # Can I do this integral using symbolic mu, (as I do below)
        # or should I stick to integrating only those mu which I have data for?
        
        term3, term3_err = integrate(lambda mu: np.sqrt(1-mu**2)*lpmv(m,l, mu), -1, 1)
        
        Z_lm = term1*term2*term3
        
        return Z_lm
    
    
    def find_disk_intg_factor (self, l, x):
        
        # I still don't know if I need a bandpass correction, or if defining b_l
        # using the bandpass corrected intensities (as I do here) is sufficient.
    
        b_l = self.atm_data[f'I_{x}_{l}'][:]/self.atm_data[f'I_{x}_0'][:]
    
        return b_l[0]
        
    
    def find_mean_anom (self, t, t_peri=0):
        
        return self.Omega_orb*(t - t_peri)
        
    
    def find_ecce_anom (self, M):
    
        Keppler = lambda E : E - self.e*np.sin(E) - M
        
        return fsolve(Keppler, 0)[0]
        
    
    def find_true_anom (self, E):
    
        return 2*np.arctan( ((1+self.e)/(1-self.e))*np.tan(E/2) )
    
    
    def convert_t_to_f (self, t, t_peri=0):
    
        M = self.find_mean_anom(t, t_peri)
        
        E = self.find_ecce_anom(M)
        
        f = self.find_true_anom(E)
        
        return f
    
    
    def find_bin_sep (self, t, t_peri=0):
        
        f = self.convert_t_to_f(t, t_peri)#* 2*np.pi
        
        D = self.a*(1-self.e**2)/(1+self.e*np.cos(f))
        
        return D
    
    
    def find_irrad (self, x, t, t_peri=0):
        
        Dt = self.find_bin_sep(t, t_peri)
        
        rel_dJ = np.zeros_like(t)
        
        n_l = self.res_data['l_max']
        
        for l in range(2, n_l+1):
            for m in range(-l, l+1):
                
                Z_lm = self.eval_ramp(l, m)
                b_l = self.find_disk_intg_factor(l, x)
                rel_dJ += b_l*Z_lm*(self.L2/self.L1)*(self.R1/Dt)**2
        
        return rel_dJ
        
            #dJ/unperturbed_observed_flux 
            # = bandpass_correction_coefficient*(L2/L1)*(R1/D(t))**2
            #   * disc_integral_factor_(Burkart_used_Eddington_limb_darkening)
