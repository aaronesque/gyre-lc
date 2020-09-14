import numpy as np
import h5py
from scipy.integrate import quad as integrate
from scipy.special import lpmv, sph_harm

### Class definition

class irradiation:
    
    def __init__ (self, atm_data, resp_data, bin_data, filter_x):
    
        self.atm_data = atm_data
        self.resp_data = resp_data
        self.x = filter_x
        
        self.omega_orb = bin_data.omega_orb
        self.e = bin_data.e
        self.a = bin_data.a
        
        self.L1 = bin_data.L1
        self.R1 = bin_data.R1
        
        self.L2 = bin_data.L2
        self.R2 = bin_data.R2
        
            
            
    def eval_ramp (self, l, m):
        
        term1 = 2* np.sqrt( ((2*l+1)/(4*np.pi)) \
                           *(np.math.factorial(l-m)/np.math.factorial(l+m)) )
        
        if (np.abs(m)==1):
            term2 = np.pi/2
        else:
            term2 = np.cos(m*np.pi/2)/(1-m**2)
        
        # Can I do this integral using symbolic mu,
        
        term3, term3_err = integrate(lambda mu: np.sqrt(1-mu**2)*lpmv(m,l, mu), -1, 1)
        
        # or should I stick to integrating only those mu which I have data for?
        Z_lm = term1*term2*term3
        
        return Z_lm
    
    
    def find_disk_intg_factor (self, l):
        
        atm = self.atm_data
        
        # I still don't know if I need a bandpass correction, or if defining b_l
        # using the bandpass corrected intensities (as I do here) is sufficient.
    
        b_l = atm.data[f'I_{self.x}_{l}'][:]/atm.data[f'I_{self.x}_0'][:]
    
        return b_l[0]
        
    
    def find_mean_anom (self, t, t_peri=0):
        
        return self.omega_orb*(t - t_peri)
        
    
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
    
    
    def eval_irrad (self, t, t_peri=0):
        
        Dt = self.find_bin_sep(t, t_peri)
        
        rel_dJ = np.zeros_like(t)
        
        n_l = self.resp_data.data['l_max']
        
        for l in range(2, n_l+1):
            for m in range(-l, l+1):
                
                Z_lm = self.eval_ramp(l, m)
                b_l = self.find_disk_intg_factor(l)
                #rel_dJ += b_l*Z_lm*(self.L2/self.L1)*(self.R1/Dt)**2
                rel_dJ += b_l*Z_lm*(self.L1/self.L2)*(self.R2/Dt)**2
        
        return rel_dJ
        
            #dJ/unperturbed_observed_flux 
            # = bandpass_correction_coefficient*(L2/L1)*(R1/D(t))**2
            #   * disc_integral_factor_(Burkart_used_Eddington_limb_darkening)
