import numpy as np
import h5py
import f90nml as nml

### Class definition

class binary:
    
    def __init__ (self, resp_data, bin_data):
        
        self.resp_data = resp_data
        
        # should check if 'inlist', this just checks 'str'
        
        if isinstance(bin_data, str):
            bin_dict = self.read_inlist(bin_data)
        else: bin_dict = bin_data
            
        self.omega_orb = bin_dict['omega_orb']
        self.a = bin_dict['a']
        self.e = bin_dict['e']
        
        self.L1 = bin_dict['L1']
        self.R1 = bin_dict['R1']
        self.M1 = bin_dict['M1']
        
        self.L2 = bin_dict['L2']
        self.R2 = bin_dict['R2']
        self.M2 = bin_dict['M2']
            
            
    def read_inlist(self, filename):
        
        inlist = nml.read(filename)
        
        data = {}
            
        def read_star(inlist, star_data, star_number):
            
            # this needs a better 'defaults' routine
            # this also needs a more general handling of 'units' 
            
            star_list = inlist[f'star{star_number}']
            
            star_data[f'Teff{star_number}'] = star_list['teff'] #K
            star_data[f'logg{star_number}'] = star_list['logg'] #dex
            
            if star_list['lum_units']=='SOL':
                lum_units = 3.839e33 #erg per l_sol
            else: lum_units = 1
                
            if star_list['r_units']=='SOL':
                r_units = 6.957e10 #cm per r_sol
            else: r_units = 1
                    
            if star_list['m_units']=='SOL':
                m_units = 1.989e33 #g per m_sol
            else: m_units = 1
            
            star_data[f'L{star_number}'] = star_list['lum']*lum_units
            star_data[f'R{star_number}'] = star_list['r']*r_units
            star_data[f'M{star_number}'] = star_list['m']*m_units
            
        def read_orbit(inlist, orb_data):
            
            orb_list = inlist['orbit']
            
            #placeholder unit routine
            
            if orb_list['omega_orb_units']=='CYC_PER_DAY':
                omega_orb_units = 1 
            else: omega_orb_units = 1
            
            if orb_list['Omega_orb']=='DEFAULT':
                orb_data['omega_orb'] = self.resp_data.data['Omega_orb']
            else: orb_data['omega_orb'] = orb_list['omega_orb']
                
            orb_data['omega_orb'] *= omega_orb_units
            
            if orb_list['a_units']=='SOL':
                a_units = 6.957e10 #cm per r_sol
            else: a_units = 1
                
            orb_data['a'] = orb_list['a']*a_units
            
            orb_data['e'] = orb_list['e']
            
        read_star(inlist, data, 1)
        read_star(inlist, data, 2)
        read_orbit(inlist, data)
        
        return data
