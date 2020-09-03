import numpy as np
import h5py
from scipy.integrate import quad as integrate
from scipy.special import legendre, lpmv, sph_harm
import f90nml as nml

### Class definition

class irradiation:
    
    def __init__ (self, resp_data, atm_data):
    
        self.resp_data = resp_data
        self.atm_data = atm_data
        #self.binary_data = binary_data
        
    
    def solve_true_anom (t, omega_orb, e):
    
        #def convert_t_to_phase
    
        def eval_true_anom (t, omega_orb, e):
    
            return lambda f : omega_orb*t \
                - 2.*np.arctan((1.-e)*np.tan(f/2)/np.sqrt(1.-e**2)) \
                + e*np.sqrt(1.-e**2)*np.sin(f)/(1.+e*np.cos(f))
    
        soln = np.array([fsolve(eval_true_anom(tau,omega_orb,e), tau%np.pi) for tau in t])
    
        return soln
    
    
    def eval_ramp_func (self, l, m, mu):
        
        term1 = 2 * np.sqrt((2*l+1)/(4*np.pi))
        
        term2 = np.sqrt( np.math.factorial(l-m)/np.math.factorial(l+m) )
        
        term3 = np.cos(m*np.pi/2)/(1-m**2)
        
        # Can I do this integral using symbolic mu,
        
        term4 = integrate(lambda mu: np.sqrt(1-mu**2)*lpmv(m,l, mu), -1, 1)
        
        # or should I stick to integrating only those mu which I have data for?
        
        Z_lm = term1*term2*term3*term4
        
        return Z_lm
    
    
    def eval_di_factor (self, x, l):
        
        atm = self.atm_data
        
        # I still don't know if I need a bandpass correction, or if defining b_l
        # using the bandpass corrected intensities (as I do here) is sufficient.
    
        b_l = atm.data[f'I_{x}_{l}'][:]/atm.data['I_{x}_0'][:]
    
        return b_l
    
    
    #def eval_irrad (self):
    
        #for l in range(2, n_l+1):
            #for m in range(-n_l, n_l+1):
        
            #dJ/unperturbed_observed_flux 
            # = bandpass_correction_coefficient*(L2/L1)*(R1/D(t))**2
            #   * disc_integral_factor_(Burkart_used_Eddington_limb_darkening)
