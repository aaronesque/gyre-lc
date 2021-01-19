import numpy as np
import h5py
from scipy.special import sph_harm

from irrad_func import irradiation

### Class definition

class observer:
    
    def __init__ (self, bin_data):
    
        self.bin_data = bin_data
        
        self.irr = {1: irradiation(bin_data, 1),\
                    2: irradiation(bin_data, 2)}
        
    
    def convert_coords (self, inc, omega):

        theta = inc/180 * np.pi
        phi = (90-omega)/180 * np.pi
        
        return theta, phi
    
    
    def calc_bol_I (self, I_0, l):
        
        if l == 0:
            I_l = 1
        elif l == 1:
            I_l = 1/3
        elif l == 2:
            I_l = 1/8
        elif l == 3:
            I_l = 0
        elif l == 4:
            I_l = -1/48
        elif l == 5:
            I_l = 0
        elif l == 6:
            I_l = 1/128
        else:
            raise Exception('No I_l defined')
            
        return I_l*I_0
            
    
    def eval_fourier (self, star_number, inc, omega, x, bol=False):
        
        res = self.bin_data.star[star_number].res
        atm = self.bin_data.star[star_number].atm
        
        theta, phi = self.convert_coords(inc, omega)
        
        # Initialize the frequencies/amplitudes arrays
        
        f = np.arange(res.data['k_max']+1)*res.data['Omega_orb']
        
        A = np.zeros(res.data['k_max']+1, dtype=np.complex)
            

        # Loop over l, m and k
        I = atm.data
    
        for l in range(2, res.data['l_max']+1):
            for m in range(-l, l+1):
                for k in range(0, res.data['k_max']+1):

                    i_l = l
                    i_m = m + res.data['l_max']
                    i_k = k

                    dR_lmk = res.R_lmk[i_l,i_m,i_k]
                    dT_lmk = res.T_lmk[i_l,i_m,i_k]
                    dG_lmk = res.G_lmk[i_l,i_m,i_k]

                    Y_lm = sph_harm(m, l, phi, theta)
                    
                    R_xl = (2 + l)*(1 - l)*I[f'I_{x}_{l}'][:] / I[f'I_{x}_0'][:]
                    T_xl = I[f'dlnTeff_{x}_{l}'][:] / I[f'I_{x}_0'][:] 
                    G_xl = I[f'dlng_{x}_{l}'][:] / I[f'I_{x}_0'][:]
                    
                    if bol==True:
                        I_l = self.calc_bol_I( I[f'I_{x}_0'][:], l)
                        R_xl = (2 + l)*(1 - l)*I_l / I[f'I_{x}_0'][:]
                        T_xl = 4*I_l / I[f'I_{x}_0'][:] 
                        G_xl = 0
                    
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
    

    def find_star_flux (self, star_number, inc, omega, x, t, t_peri=0, bol=False):
        
        res = self.bin_data.star[star_number].res
        Omega_orb = res.data['Omega_orb']
        
        f, A = self.eval_fourier(star_number, inc, omega, x, bol)
        
        # Initialize the frequencies/amplitudes arrays
        
        if isinstance(t, np.ndarray):
            
            star_flux = np.zeros_like(t)
        
        elif isinstance(t, list):
            
            t = np.array(t)
            star_flux = np.zeros_like(t)
            
        else: #How to isinstance any mathable t? float, int, etc?
            
            star_flux = 0.

        # Add contributions from each frequency component

        n_A = len(A)

        for i in range(n_A):
            
            star_flux += np.real(A[i] * np.exp(-1j*f[i]*2*np.pi*(t - t_peri)))
            
        return star_flux
    
    
    def find_flux (self, inc, omega, x, t, t_peri=0, bol=False):
        
        irr_1 = self.irr[1].find_irrad(x, t, t_peri)
        irr_2 = self.irr[2].find_irrad(x, t, t_peri)
        
        star_1 = self.find_star_flux(1, inc, omega, x, t, t_peri, bol)
        star_2 = self.find_star_flux(2, inc, omega, x, t, t_peri, bol)
        
        return irr_1 + star_1 + irr_2 + star_2
