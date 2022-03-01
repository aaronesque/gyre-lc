import numpy as np
import h5py
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
    """Takes parameters necessary for plotting a light curve.

    Parameters are taken in as kwargs. To function properly, either 
    `mesa_model` or `mass` (if modeling a point mass) must be taken.

    Attributes:
        mesa_model (str): Path to MESA stellar model, defaults to None
            if blank or if model type is a point mass
        gyre_model (str): Path to GYRE tides model, defaults to zero
            tides if blank or if model type is a point mass
        photgrid (:py:class:`pymsg.PhotGrid`): Photometric grid class
            object, produced using MSG
        mass (float): Stellar mass is specified by user if model type is
            a point mass, else is auto-set according to mesa_model
        radius (float): Stellar radius can be specified by user if model
            type is a point mass, else is auto-set according to mesa_model
        luminosity (float): Bolometric luminosity can be specified by user
            if model type is a point mass, else is auto-set according to
            mesa_model
        Teff (float): Effective temperature can be specified by user if
            neither mesa_model or photgrid are specified and MSG
            spectra are desired
        logg (float): Effective surface gravity can be specified by user
            if neither mesa_model or photgrid are specified and
            MSG spectra are desired
        units (str): User may choose between 'CGS' and 'SOLAR'
        resp_coeffs (dict): A dictionary containing the tidal response
            coefficients from GYRE-tides output
        phot_coeffs (dict): A dictionary containing the photometric 
            coefficients from `gyrelc.Observer`
        point_mass_model (bool): Gets set to True when `mesa_model` is 
            unspecified

    """

    def __init__ (self, mesa_model=None, gyre_model=None, 
            photgrid=None, mass=None, radius=None, 
            luminosity=None, Teff=None, logg=None, units=None):
        """Constructor method
        """
        # if mesa model is specified, read mesa params
        if mesa_model is not None:
            self.model = mesa_model
            self.model_type = 'MESA'
            self.read_mesa_params()

            # also check for gyre model for response coefficients
            if gyre_model is not None:
                self.gyre_model = gyre_model
                self.resp_coeffs = rc.resp_coeffs( self.gyre_model )
            else: self.resp_coeffs = rc.resp_coeffs('')

        # else check for a point mass, read user-specified params
        elif mass is not None:
            # make note of model as point mass type
            self.model = f'{mass} M_sol'
            self.model_type = 'point mass'
            self.read_pt_mass_params()

            # no gyre_model allowed
            if gyre_model is not None:
                raise Exception("A point mass cannot experience tides")
            self.resp_coeffs = rc.resp_coeffs('')

        else: raise Exception("Star() must take valid 'mesa_model' or 'mass' argument")

        self.photgrid = photgrid
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


    def read_mesa_params(self):
        """Reads and builds stellar parameters from
        user-specified MESA model. Converts to solar units
        assuming MESA output is CGS.
        """
        data = ascii.read(self.model, data_start=0, data_end=1)

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

    
    def read_phot_coeffs(self, photgrid):
        """Selects appropriate photometric coefficient routine
        depending on the stellar model
        """
        if self.luminosity==0.:
            self.phot_coeffs.update( self.make_phot_coeffs_pt_mass(photgrid) )

        if self.model_type=='MESA':
            if isinstance(photgrid, pymsg.PhotGrid):
                self.phot_coeffs.update( self.read_phot_coeffs_msg(photgrid) )
            elif isinstance(photgrid, str):
                self.phot_coeffs.update( self.read_phot_coeffs_h5(photgrid) )
            elif isinstance(photgrid, None):
                raise Exception("Must specify photgrid")
            else: raise Exception(f"Invalid photgrid type")
        elif self.model_type=='point mass':
            self.phot_coeffs.update( self.make_phot_coeffs_pt_mass(photgrid) )
        else: raise Exception(f"Invalid component model type must be 'mesa_model' or 'point_mass_model'.")
        
        return


    def read_phot_coeffs_h5(self, synspec_model):
        """This function is deprecated and was tested on a
        synspec intensity spectrum calculated from
        iOri's tlusty atmosphere.
        """
        # Read intensity file
        I = ac.atm_coeffs(synspec_model)
        I = I.data

        # Evaluate photometric data
        I_x_l = {}
        dI_dlnT_x_l = {}
        dI_dlng_x_l = {}
       
        # Use filename as filter name
        filter_x = synspec_model.split(",")[0]

        for l in range(self.resp_coeffs.data['l_max']+1):
            I_x_l[l] = I[f'I_{filter_x}_{l}'][0]
            dI_dlnT_x_l[l] = I[f'dlnTeff_{filter_x}_{l}'][0]
            dI_dlng_x_l[l] = I[f'dlng_{filter_x}_{l}'][0]

        return {f'I_{filter_x}': I_x_l,
                f'dI_dlnT_{filter_x}': dI_dlnT_x_l,
                f'dI_dlng_{filter_x}': dI_dlng_x_l}


    def read_phot_coeffs_msg(self, photgrid):
        """Creates and stores photometric coefficients from
        MSG for the desired filter
        """
        # Set atmosphere parameters dict

        Teff = self.Teff
        logg = self.logg

        dx = {'logT': np.log10(Teff), 'logg': logg}
        pg = photgrid

        # Evaluate photometric data
        I_x_l = {}
        dI_dlnT_x_l = {}
        dI_dlng_x_l = {}

        for l in range(self.resp_coeffs.data['l_max']+1):
            I_x_l[l] = pg.D_moment(dx, l)
            dI_dlnT_x_l[l] = pg.D_moment(dx, l, deriv={'logT':True})/np.log(10.)
            dI_dlng_x_l[l] = pg.D_moment(dx, l, deriv={'logg':True})/np.log(10.)

        return {f'I_x': I_x_l,
                f'dI_dlnT_x': dI_dlnT_x_l,
                f'dI_dlng_x': dI_dlng_x_l}


    def make_phot_coeffs_pt_mass(self, photgrid):
        """Creates and stores zeros as photometric
        coefficients for a point mass model
        """
        return {f'I_x': np.array([0.]),
                f'dI_dlnT_x':np.array([0.]),
                f'dI_dlng_x': np.array(0.)}


    def D_moment(self, l, deriv=None):
        """Makes reading photometric coefficients from
        a gyrelc.Star instance less painful
        """
        I = self.phot_coeffs
        if deriv=='lnT':
            return I[f'dI_dlnT_x'][l]
        elif deriv=='lng':
            return I[f'dI_dlng_x'][l]
        elif deriv=='logT':
            return I[f'dI_dlnT_x'][l] * np.log(10.)
        elif deriv=='logT':
            return I[f'dI_dlng_x'][l] * np.log(10.)
        else:
            return I[f'I_x'][l]


    def R_xl(self, l):
        I = self.phot_coeffs
        return (2 + l)*(1 - l)*I[f'I_x'][l] / I[f'I_x'][0]

    def T_xl(self, l):
        I = self.phot_coeffs
        return I[f'dI_dlnT_x'][l] / I[f'I_x'][0]

    def G_xl(self, l):
        I = self.phot_coeffs
        return I[f'dI_dlng_x'][l] / I[f'I_x'][0]

    def disk_intg_factor (self, l):
        I = self.phot_coeffs
        b_l = I[f'I_x'][l]/I[f'I_x'][0]
        return b_l


    def eval_fourier_moment (self, theta, phi, l,m,k):

        i_l = l
        i_m = m + self.resp_coeffs.data['l_max']
        i_k = k

        dR_lmk = self.resp_coeffs.R_lmk[i_l,i_m,i_k]
        dT_lmk = self.resp_coeffs.T_lmk[i_l,i_m,i_k]
        dG_lmk = self.resp_coeffs.G_lmk[i_l,i_m,i_k]

        Y_lm = sph_harm(m, l, phi, theta)

        R_xl = self.R_xl(l)
        T_xl = self.T_xl(l)
        G_xl = self.G_xl(l)

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


    def eval_fourier (self, theta, phi):

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

                        A[k] += self.eval_fourier_moment(theta, phi, l,m,k)

        # Return data
        return f, A

