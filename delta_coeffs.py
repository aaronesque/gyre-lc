# grid class

import numpy as np
import h5py

### Class definition

class delta_coeffs:
    
    def __init__ (self, response_file):
        
        data = self.read_response(response_file[0])
        
        # Sanity check
        
        if len(data) == 0:
            raise Exception('Empty response file')
        if not isinstance(data, dict):
            raise Exception('Response file type error')
            
        # Setup axes
        
        n_l = data['l_max'] + 1
        n_m = 2*data['l_max'] + 1
        n_k = data['k_max'] + 1
        
        self.R_coeffs = np.empty([n_l, n_m, n_k], dtype=complex)
        self.T_coeffs = np.empty([n_l, n_m, n_k], dtype=complex)
        self.G_coeffs = np.empty([n_l, n_m, n_k], dtype=complex)

        for i_l in range(2, n_l):
            for i_m in range(0, n_m):
                for i_k in range(0, n_k):

                    self.R_coeffs[i_l,i_m,i_k] = self.Delta_R(data, i_l,i_m,i_k)
                    self.T_coeffs[i_l,i_m,i_k] = self.Delta_T(data, i_l,i_m,i_k)
                    self.G_coeffs[i_l,i_m,i_k] = self.Delta_g(data, i_l,i_m,i_k)
        
        
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
    
    
    def Delta_R(self, data, l,m,k):
        
        return np.sqrt(4.*np.pi) * data['xi_r_ref'][k,m,l]
       
        
    def Delta_T(self, data, l,m,k):
            
        xi_r_ref = data['xi_r_ref'][k,m,l]

        lag_L_ref = data['lag_L_ref'][k,m,l]
        
        return np.sqrt(4.*np.pi)*(lag_L_ref - 2*xi_r_ref)/4.
      
    
    def Delta_g(self, gyre_data, l,m,k):
    
        xi_r_ref = gyre_data['xi_r_ref'][k,m,l]
        
        omega = -k*gyre_data['Omega_orb'] - m*gyre_data['Omega_rot']
        
        return np.sqrt(4*np.pi)*(-omega**2 - 2)*xi_r_ref
