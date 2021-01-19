import numpy as np
import h5py
import f90nml as nml

### Class definition

class binary:
    
    def __init__ (self, list_path):
        
        # should check if 'inlist', this just checks 'str'
        
        if isinstance(list_path, str):
            self.nml = nml.read(list_path)
        else: raise Exception('Inlist file error')
        
        self.star = {1: star(list_path, 1),\
                     2: star(list_path, 2) }
        
        self.orbit = self.read_orbit(self.nml)
            
        
    def read_orbit(self, inlist):
        
        orb_list = inlist['orbit']
        
        params = {}
        
        #placeholder unit routine
        
        if orb_list['omega_orb_units']=='CYC_PER_DAY':
            omega_orb_units = 1 
        else: omega_orb_units = 1
        
        if orb_list['a_units']=='CM':
            a_units = 6.957e10 #cm per r_sol
        else: a_units = 1
            
        # eventually, we want to get these from resp_data
        # resp_data currently does not contain this output tho
        
        params['a'] = orb_list['a']*a_units
        params['e'] = orb_list['e']
        params['Omega_orb'] = orb_list['omega_orb']*omega_orb_units
        
        return params
