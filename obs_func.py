# grid class

import numpy as np
import h5py
from scipy.special import legendre, sph_harm
import f90nml as nml

### Class definition

class observer:
    
    def __init__ (self, resp_data, atm_data):
    
        inc = 62.9
        omega = 122.2

        theta = inc/180 * np.pi
        phi = (90-omega)/180 * np.pi
        
        self.eval_fourier (resp_data, atm_data,'B', theta,phi)
    
    
    def eval_fourier (resp_data, atm_data, x, theta_0, phi_0):
        
        resp = resp_data
        atm = atm_data
        
        # Initialize the frequencies/amplitudes arrays
        f = np.arange(resp.data['k_max']+1)*resp.data['Omega_orb']
        A = np.zeros(resp.data['k_max']+1, dtype=np.complex)

        # Loop over l, m and k
    
        for l in range(2, resp.data['l_max']+1):
            for m in range(-l, l+1):
                for k in range(0, resp.data['k_max']+1):

                    i_l = l
                    i_m = m + l
                    i_k = k

                    dR_lmk = resp.R_lmk[i_l,i_m,i_k]
                    dT_lmk = resp.T_lmk[i_l,i_m,i_k]
                    dG_lmk = resp.G_lmk[i_l,i_m,i_k]
                
                    Y_lm = sph_harm(m, l, phi_0, theta_0)

                    R_xlm = atm.R_xl[x][l] * Y_lm
                    T_xlm = atm.T_xl[x][l] * Y_lm
                    G_xlm = atm.T_xl[x][l] * Y_lm
                
                if k == 0:
                    if m == 0:
                        kappa = 0.5
                    elif m >= 1:
                        kappa = 1.
                    else:
                        kappa = 0.
                else:
                    kappa = 1.

                # Add the Fourier contribution
                
                A[i_k] += 2*kappa*(dR_lmk*R_xlm + dT_lmk*T_xlm)
        
        # Return data
        
        return f, A

