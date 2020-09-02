import numpy as np
import h5py
from scipy.special import legendre, sph_harm
import f90nml as nml

### Class definition

class observer:
    
    def __init__ (self, resp_data, atm_data):
    
        self.resp_data = resp_data
        self.atm_data = atm_data
        #self.binary_data = binary_data
    
    
    def convert_coords (self, inc, omega):

        theta = inc/180 * np.pi
        phi = (90-omega)/180 * np.pi
        
        return theta, phi

    
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
            
        
    
    def eval_fourier (self, x, inc, omega):
        
        resp = self.resp_data
        atm = self.atm_data
        
        theta, phi = self.convert_coords(inc, omega)
        
        # Initialize the frequencies/amplitudes arrays
        
        f = np.arange(resp.data['k_max']+1)*resp.data['Omega_orb']
        
        A = np.zeros(resp.data['k_max']+1, dtype=np.complex)
        

        # Loop over l, m and k
        
    
        for l in range(2, resp.data['l_max']+1):
            for m in range(-l, l+1):
                for k in range(0, resp.data['k_max']+1):

                    i_l = l
                    i_m = m + resp.data['l_max']
                    i_k = k

                    dR_lmk = resp.R_lmk[i_l,i_m,i_k]
                    dT_lmk = resp.T_lmk[i_l,i_m,i_k]
                    dG_lmk = resp.G_lmk[i_l,i_m,i_k]

                    Y_lm = sph_harm(m, l, phi, theta)
                    
                    R_xlm = atm.R_xl[x][l]
                    T_xlm = atm.T_xl[x][l]
                    G_xlm = atm.T_xl[x][l]
                
                    if k == 0:
                        if m == 0:
                            kappa = 0.5
                        elif m >= 1:
                            kappa = 1.
                        else:
                            kappa = 0.
                    else:
                        kappa = 1.
    
                    # Add the Fourier contribution * spherical harmonic
                    
                    A[i_k] += 2*kappa*(dR_lmk*R_xlm + dT_lmk*T_xlm + dG_lmk*G_xlm)*Y_lm
        
        # Return data
        
        return f, A
    

    def find_flux (self, x, inc, omega, t):
        
        resp = self.resp_data
        Omega_orb = resp.data['Omega_orb']
        
        f, A = self.eval_fourier(x, inc, omega)
        
        # Initialize the frequencies/amplitudes arrays
        
        if isinstance(t, np.ndarray):
            
            diff_flux = np.zeros_like(t)
        
        elif isinstance(t, list):
            
            t = np.array(t)
            diff_flux = np.zeros_like(t)
            
        else: #How to isinstance any mathable t? float, int, etc?
            
            diff_flux = 0.

        # Add contributions from each frequency component

        n = len(A)

        for i in range(n):
            
            diff_flux += np.real(A[i] * np.exp(1j*f[i]*2*np.pi*t))
            
        return diff_flux
    
