import numpy as np
import h5py

## Class definition

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

