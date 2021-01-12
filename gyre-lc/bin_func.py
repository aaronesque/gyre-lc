import numpy as np
import h5py
import f90nml as nml

### Class definition

class binary:
    
    def __init__ (self, list_path):
        
        # should check if 'inlist', this just checks 'str'
        
        if isinstance(list_path, str):
            bin_list = nml.read(list_path)
            self.nml = bin_list
        else: raise Exception('Inlist file error')
        
        self.star = {}
    
        self.star = self.read_star(bin_list, 1)
        self.star.update( self.read_star(bin_list, 2) )
        
        self.orbit = self.read_orbit(bin_list)
        
            
    def read_star(self, inlist, star_number):
        
        star_list = inlist[f'star_{star_number}']
        params = {}
        
        model_path = star_list['star_model_path']
        
        params.update( self.read_mesa(model_path) )
        
        return {star_number: params}
        
        
    def read_mesa(self, mesa_model_path, units='SOL'):
        
        data = ascii.read(mesa_model_path, data_start=0, data_end=1)
        
        # cgs constants
        
        R_sol = 6.957e10 
        M_sol = 1.989e33
        L_sol = 3.839e33
        G = 6.674079999999999e-08 
        sigma_sb = 5.6703669999999995e-05
        
        # data in cgs
        
        r = data[0][1]
        m = data[0][2]
        l = data[0][3]
            
        # calculate Teff [K], logg [dex]
        
        Teff = (l/(sigma_sb * 4*np.pi*r**2))**0.25
        
        g_surf = (m/r**2)/(M_sol/R_sol**2)
        logg = np.log10(g_surf)
        
        # convert to solar units
        
        if units=='SOL':
            r = r/R_sol
            m = m/M_sol
            l = l/L_sol
            
        return {'R': r,
                'M': m,
                'L': l,
                'Teff': Teff,
                'logg': logg}
    
    
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
        
        #q = f.attrs['q']
        #e = f.attrs['e']
        
        f.close()
        
        return {'xi_r_ref': xi_r_ref_re + 1j*xi_r_ref_im,
                'lag_L_ref': lag_L_ref_re + 1j*lag_L_ref_im,
                'k_max': k_max,
                'l_max': l_max,
                'Omega_rot': Omega_rot,
                'Omega_orb': Omega_orb}#,
                #'q': q,
                #'e': e}
            
            
    def read_orbit(self, inlist):
        
        orb_list = inlist['orbit']
        resp_list = self.read_response(inlist['star_1']['tide_model_path'])
        
        params = {}
        
        #placeholder unit routine
        
        if orb_list['omega_orb_units']=='CYC_PER_DAY':
            omega_orb_units = 1 
        else: omega_orb_units = 1
        
        if orb_list['Omega_orb']=='DEFAULT':
            params['Omega_orb'] = resp_list['Omega_orb']
        else: params['Omega_orb'] = orb_list['omega_orb']
            
        params['Omega_orb'] *= omega_orb_units
        
        if orb_list['a_units']=='CM':
            a_units = 6.957e10 #cm per r_sol
        else: a_units = 1
            
        # eventually, we want to get these from resp_data
        # resp_data currently does not contain this output tho
        
        params['a'] = orb_list['a']*a_units
        params['e'] = orb_list['e']
        
        return params
