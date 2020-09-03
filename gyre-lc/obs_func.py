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
                    
                    R_xl = atm.R_xl[x][i_l]
                    T_xl = atm.T_xl[x][i_l]
                    G_xl = atm.T_xl[x][i_l]
                
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
                    
                    A[i_k] += 2*kappa*(dR_lmk*R_xl + dT_lmk*T_xl + dG_lmk*G_xl)*Y_lm
        
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
    
