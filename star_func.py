import numpy as np
import h5py
import f90nml as nml
from astropy.io import ascii

from resp_coeffs import resp_coeffs
from atm_coeffs import atm_coeffs

### Class definitions

class star:

    def __init__ (self, list_path, star_number):

        # should check if 'inlist', this just checks 'str'

        if isinstance(list_path, str):
            self.nml = nml.read(list_path)[f'star_{star_number}']
        else: raise Exception('Inlist file error')

        self.par = self.read_mesa_params(self.nml['star_model_path'])

        self.res = resp_coeffs(self.nml['tide_model_path'])

        self.atm = self.read_atm_coeffs(self.nml['star_model_path'])


    def read_mesa_params(self, mesa_model_path, units='SOL'):

        data = ascii.read(mesa_model_path, data_start=0, data_end=1)

        # cgs constants

        R_sol = 6.957e10
        M_sol = 1.989e33
        L_sol = 3.839e33
        G = 6.674079999999999e-08
        sigma_sb = 5.6703669999999995e-05

        # data in cgs

        m = data[0][1]
        r = data[0][2]
        l = data[0][3]

        # calculate Teff [K], logg [dex]

        Teff = (l/(sigma_sb * 4*np.pi*r**2))**0.25

        g_surf = G*m/(r**2)

        logg = np.log10(g_surf)

        # convert to solar units

        if units=='SOL':
            m = m/M_sol
            r = r/R_sol
            l = l/L_sol

        return {'M': m,
                'R': r,
                'L': l,
                'Teff': Teff,
                'logg': logg}


    def read_atm_coeffs(self, inlist):

        Teff = self.par['Teff']
        logg = self.par['logg']

        #atm = atm_coeffs(Teff,logg)

        return atm_coeffs('lc-data/iOri-partials.h5')
