import numpy as np
import h5py
import f90nml as nml
from scipy.special import sph_harm
from astropy.io import ascii

import sys
import os
sys.path.insert(0, os.path.join(os.environ['MSG_DIR'], 'lib'))
import pymsg

import resp_coeffs as rc
import atm_coeffs as ac

### Class definitions

class Star:
    """Takes path to mesa model specifying stellar
    parameters necessary for plotting a light curve.

    :param mesa_model: Path to MESA stellar model, defaults to None
        if blank or if model type is a point mass
    :type mesa_model: str, optional

    :param gyre_model: Path to GYRE tides model, defaults to zero
        tides if blank or if model type is a point mass
    :type gyre_model: str, optional

    :param synspec_model : Path to SYNSPEC atmospheric spectrum model,
        defaults to using MSG for spectra as desired
    :type synspec_model: str, optional

    :param mass: Stellar mass is specified by user if model type is
        a point mass, else is auto-set according to mesa_model
    :type mass: float, optional

    :param radius: Stellar radius can be specified by user if model
        type is a point mass, else is auto-set according to mesa_model
    :type radius: float, optional

    :param luminosity: Bolometric luminosity can be specified by user
        if model type is a point mass, else is auto-set according to
        mesa_model
    :type luminosity: float, optional

    :param Teff: Effective temperature can be specified by user if
        neither mesa_model or synspec_model are specified and MSG
        spectra are desired
    :type Teff: float, optional

    :param logg: Effective surface gravity can be specified by user
        if neither mesa_model or synspec_model are specified and
        MSG spectra are desired
    :type logg: float, optional

    :param units: User may choose between 'CGS' and 'SOLAR'
    :type units: str, default='SOLAR'
    """

    def __init__ (self, **kwargs):
        """Constructor method
        """
        # takes user inputs from limited list
        allowed_keys = {'mass', 'radius', 'luminosity', 'Teff', 'logg', 'units',\
                        'mesa_model', 'gyre_model', 'synspec_model'}
        self.__dict__.update((k, v) for k, v in kwargs.items() if k in allowed_keys)

        # if mesa model is specified, read mesa params
        if kwargs.get('mesa_model'):
            self.read_mesa_params(self.mesa_model)

            # also check for gyre model for response coefficients
            if kwargs.get('gyre_model'):
                self.resp_coeffs = rc.resp_coeffs( self.gyre_model )
            else: self.resp_coeffs = rc.resp_coeffs('')

        # else check for a point mass, read user-specified params
        elif kwargs.get('mass'):

            # make note of model as point mass type
            self.point_mass_model = True
            self.read_pt_mass_params()

            # no gyre_model allowed
            if kwargs.get('gyre_model'):
                raise Exception("A point mass cannot experience tides")
            self.resp_coeffs = rc.resp_coeffs('')


        else: raise Exception("Star() must take valid 'mesa_model' or 'mass' argument")

        self.phot_coeffs = {}


    def read_pt_mass_params(self):
        """Checks for user-specified point mass parameters.
        Does unit conversions as needed.
        """
        # cgs constants
        M_sol = 1.989e33
        L_sol = 3.839e33
        R_sol = 6.957e10
        G = 6.674079999999999e-08
        sigma_sb = 5.6703669999999995e-05
        
        if self.__dict__.get('mass'):
            pass
        else: self.__dict__['mass'] = 0.
        if self.__dict__.get('luminosity'):
            pass
        else: self.__dict__['luminosity'] = 0.
        if self.__dict__.get('radius'):
            pass
        else: self.__dict__['radius'] = 0.

        if self.__dict__.get('units'):
            if upper(self.__dict__['units'])=='SOLAR':
                self.__dict__['__mass_units'] = M_sol
                self.__dict__['__luminosity_units'] = L_sol
                self.__dict__['__radius_units'] = R_sol
            elif upper(self.__dict__['units'])=='CGS':
                self.__dict__['__mass_units'] = 1
                self.__dict__['__luminosity_units'] = 1
                self.__dict__['__radius_units'] = 1
            else:
                raise Exception(f'user-specified `units` unrecognized')
        else: # assume user assumed 'solar'
            self.__dict__['units'] = 'SOLAR'
            self.__dict__['__mass_units'] = M_sol
            self.__dict__['__luminosity_units'] = L_sol
            self.__dict__['__radius_units'] = R_sol

        self.__dict__['__mass'] = self.__dict__['mass']*self.__dict__['__mass_units']
        self.__dict__['__luminosity'] = self.__dict__['luminosity']*self.__dict__['__luminosity_units']
        self.__dict__['__radius'] = self.__dict__['radius']*self.__dict__['__radius_units']
        
        # calculate Teff [K], logg [dex]
        if self.__dict__.get('Teff'): 
            pass
        else:
            try:
                self.Teff = (self.luminosity/(sigma_sb * 4*np.pi*self.radius**2))**0.25
            except ZeroDivisionError:
                self.Teff = 0
                print(f'Invalid radius, r={self.radius} <= 0. Teff set to 0.')
        
        if self.__dict__.get('logg'): 
            pass
        else:
            try:
                g_surf = G*self.mass/(self.radius**2)
            except ZeroDivisionError:
                self.logg = 0
                print(f'Invalid radius, r={self.radius} <= 0. log(g) set to 0.')
            else:
                self.logg = np.log10(g_surf)
        
        return

    def make_pt_mass_response(self):
        """Makes a zero-response for the point mass.
        """
        return {'xi_r_ref': 0+0j,
                'lag_L_ref': 0+0j,
                'k_max': 0.,
                'l_max': 0.,
                'Omega_rot': 0.,
                'Omega_orb': 0.}


    def read_mesa_params(self, mesa_model_path):
        """Reads and builds stellar parameters from
        user-specified MESA model. Converts to solar units
        assuming MESA output is CGS.
        """
        data = ascii.read(mesa_model_path, data_start=0, data_end=1)

        # cgs constants

        M_sol = 1.989e33
        L_sol = 3.839e33
        R_sol = 6.957e10
        G = 6.674079999999999e-08
        sigma_sb = 5.6703669999999995e-05

        m = data[0][1]
        r = data[0][2]
        l = data[0][3]

        # calculate Teff [K], logg [dex]
        self.__dict__['Teff'] = (l/(sigma_sb * 4*np.pi*r**2))**0.25
        g_surf = G*m/(r**2)
        self.__dict__['logg'] = np.log10(g_surf)

        self.__dict__['units'] = 'SOLAR'
        self.__dict__['mass'] = m/M_sol
        self.__dict__['radius'] = r/R_sol
        self.__dict__['luminosity'] = l/L_sol
        return


    def read_phot_coeffs_h5(self, filter_x, synspec_model):
        """
        This function is old and was tested on a
        synspec intensity spectrum calculated from
        iOri's tlusty atmosphere.
        """

        l_min = 0
        dl = 1
        l_max = self.resp_coeffs.data['l_max']
        n_l = np.ceil( (l_max - l_min)/dl ) + 1

        l_range = l_min + dl*np.arange(n_l)

        # read intensity file

        I = ac.atm_coeffs(synspec_model)
        I = I.data

        # Evaluate photometric data
        I_x_l = {}
        dI_dlnT_x_l = {}
        dI_dlng_x_l = {}

        for l in l_range.astype(int):
            I_x_l[l] = I[f'I_{filter_x}_{l}'][0]
            dI_dlnT_x_l[l] = I[f'dlnTeff_{filter_x}_{l}'][0]
            dI_dlng_x_l[l] = I[f'dlng_{filter_x}_{l}'][0]

        return {f'I_{filter_x}': I_x_l,
                f'dI_dlnT_{filter_x}': dI_dlnT_x_l,
                f'dI_dlng_{filter_x}': dI_dlng_x_l}


    def read_phot_coeffs_msg(self, filter_x):

        # Set atmosphere parameters dict

        Teff = self.Teff
        logg = self.logg

        dx = {'logT': np.log10(Teff), 'logg': logg}
        pg = pymsg.PhotGrid(f"{os.environ['GYRELC_DIR']}/grid/{filter_x}.h5")

        # Set up intensity moment range

        l_min = 0
        dl = 1
        l_max = self.resp_coeffs.data['l_max']
        n_l = np.ceil( (l_max - l_min)/dl ) + 1

        l_range = l_min + dl*np.arange(n_l)

        # Evaluate photometric data
        I_x_l = {}
        dI_dlnT_x_l = {}
        dI_dlng_x_l = {}

        for l in l_range.astype(int):
            I_x_l[l] = pg.D_moment(dx, l)
            dI_dlnT_x_l[l] = pg.D_moment(dx, l, deriv={'logT':True})/np.log(10.)
            dI_dlng_x_l[l] = pg.D_moment(dx, l, deriv={'logg':True})/np.log(10.)

        return {f'I_{filter_x}': I_x_l,
                f'dI_dlnT_{filter_x}': dI_dlnT_x_l,
                f'dI_dlng_{filter_x}': dI_dlng_x_l}


    def make_phot_coeffs_pt_mass(self, filter_x):
        return {f'I_{filter_x}': np.array([0.]),
                f'dI_dlnT_{filter_x}':np.array([0.]),
                f'dI_dlng_{filter_x}': np.array(0.)}


    def read_phot_coeffs(self, filter_x, synspec_model=None):
        """
        """
        if self.luminosity==0.:
            self.phot_coeffs.update( self.make_phot_coeffs_pt_mass(filter_x) )

        if self.__dict__.get('mesa_model'):
            if synspec_model==None:
                self.phot_coeffs.update( self.read_phot_coeffs_msg(filter_x) )
            else:
                self.phot_coeffs.update( self.read_phot_coeffs_h5(filter_x, synspec_model) )
        elif self.__dict__.get('point_mass_model'):
            self.phot_coeffs.update( self.make_phot_coeffs_pt_mass(filter_x) )
        else: raise Exception(f"Invalid component model type must be 'mesa_model' or 'point_mass_model'.")
        
        return


    def D_moment(self, l, filter_x, deriv=None):
        I = self.phot_coeffs
        if deriv=='lnT':
            return I[f'dI_dlnT_{filter_x}'][l]
        elif deriv=='lng':
            return I[f'dI_dlng_{filter_x}'][l]
        elif deriv=='logT':
            return I[f'dI_dlnT_{filter_x}'][l] * np.log(10.)
        elif deriv=='logT':
            return I[f'dI_dlng_{filter_x}'][l] * np.log(10.)
        else:
            return I[f'I_{filter_x}'][l]


    def R_xl(self, filter_x, l):
        I = self.phot_coeffs
        return (2 + l)*(1 - l)*I[f'I_{filter_x}'][l] / I[f'I_{filter_x}'][0]

    def T_xl(self, filter_x, l):
        I = self.phot_coeffs
        return I[f'dI_dlnT_{filter_x}'][l] / I[f'I_{filter_x}'][0]

    def G_xl(self, filter_x, l):
        I = self.phot_coeffs
        return I[f'dI_dlng_{filter_x}'][l] / I[f'I_{filter_x}'][0]

    def disk_intg_factor (self, filter_x, l):
        I = self.phot_coeffs
        b_l = I[f'I_{filter_x}'][l]/I[f'I_{filter_x}'][0]
        return b_l


    def eval_fourier_moment (self, filter_x, theta, phi, l,m,k):

        i_l = l
        i_m = m + self.resp_coeffs.data['l_max']
        i_k = k

        dR_lmk = self.resp_coeffs.R_lmk[i_l,i_m,i_k]
        dT_lmk = self.resp_coeffs.T_lmk[i_l,i_m,i_k]
        dG_lmk = self.resp_coeffs.G_lmk[i_l,i_m,i_k]

        Y_lm = sph_harm(m, l, phi, theta)

        R_xl = self.R_xl(filter_x, l)
        T_xl = self.T_xl(filter_x, l)
        G_xl = self.G_xl(filter_x, l)

        if k == 0:
            if m == 0:
                kappa = 0.5
            elif m >= 1:
                kappa = 1.
            else:
                kappa = 0.
        else:
            kappa = 1.

        return 2*kappa*(dR_lmk*R_xl + dT_lmk*T_xl + dG_lmk*G_xl)*Y_lm


    def eval_fourier (self, filter_x, theta, phi):

        if self.__dict__.get('point_mass_model'):
            f, A = np.array([0]), np.array([0])
        else:
            resp_coeffs = self.resp_coeffs
            I = self.phot_coeffs

            # Initialize the frequencies/amplitudes arrays

            f = np.arange(resp_coeffs.data['k_max']+1)*resp_coeffs.data['Omega_orb']

            A = np.zeros(resp_coeffs.data['k_max']+1, dtype=complex)

            # Loop over l, m and k

            for l in np.arange(2, resp_coeffs.data['l_max']+1).astype(int):
                for m in np.arange(-l, l+1).astype(int):
                    for k in np.arange(0, resp_coeffs.data['k_max']+1).astype(int):
                        # Add the Fourier contribution * spherical harmonic

                        A[k] += self.eval_fourier_moment(filter_x, theta, phi, l,m,k)

        # Return data
        return f, A

