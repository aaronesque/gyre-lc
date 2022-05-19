# Import standard modules
import numpy as np
import sys
import os
#import Star

# Import special modules
import h5py
from scipy.special import lpmv, sph_harm
from scipy.integrate import quad as integrate
from scipy.optimize import fsolve
from astropy.io import ascii


### Class definitions

class Star:
    """Takes parameters necessary for plotting a light curve.

    Parameters are taken in as kwargs. To function properly, either 
    `mesa_model` or `mass` (if modeling a point mass) must be taken.

    Attributes:
        resp_coeffs (dict): A dictionary containing the tidal response
            coefficients from GYRE-tides output
        phot_coeffs (dict): A dictionary containing the photometric 
            coefficients from :py:class:`pymsg.PhotGrid`
        params (dict): A dictionary containing the stellar parameters
        model_type (str): Denotes the type of stellar model this is, 
            'MESA' or 'point mass'
    """

    def __init__ (self, mesa_model=None, gyre_model=None, 
            photgrid=None, pt_mass_model=None):
        #mass=None, radius=None, luminosity=None, Teff=None, logg=None, units=None):
        """Constructor method
        """
        self.params = {}
        # if mesa model is specified, read mesa params
        if mesa_model is not None:
            self.model_type = 'MESA'
            self.read_mesa_params(mesa_model)

            # also check for gyre model for response coefficients
            if gyre_model is not None:
                self.resp_coeffs = Response( gyre_model )
            else: self.resp_coeffs = Response(None)

        # else check for a point mass, read user-specified params
        elif pt_mass_model is not None:
            # make note of model as point mass type
            self.model_type = 'point mass'
            self.read_pt_mass_params(pt_mass_model)
            #self.model = f'{mass} M_sol'

            # no gyre_model allowed
            if gyre_model is not None:
                raise Exception("A point mass cannot experience tides")
            self.resp_coeffs = Response(None)

        else: raise Exception("Star() must take valid 'mesa_model' or 'pt_mass_model' argument")

        dx = {'logT': np.log10(self.params['Teff']), 'logg': self.params['logg']}
        self.phot_coeffs = Photosphere(self.resp_coeffs, photgrid, dx)


    def read_mesa_params(self, mesa_model):
        """Reads and builds stellar parameters from
        user-specified MESA model. Converts to solar units
        assuming MESA output is CGS.
        """
        data = ascii.read(mesa_model, data_start=0, data_end=1)

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
        self.params['Teff'] = (l/(sigma_sb * 4*np.pi*r**2))**0.25
        g_surf = G*m/(r**2)
        self.params['logg'] = np.log10(g_surf)

        self.params['units'] = 'SOLAR'
        self.params['mass'] = m/M_sol
        self.params['radius'] = r/R_sol
        self.params['luminosity'] = l/L_sol
        print(self.params)
        return
    

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

###

class Response:
    
    def __init__ (self, response_file):
       
        self.meta = {'model': response_file}
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
        
        if filename == None:
            
            # Fabricate dummy data
            xi_r_ref_re = np.zeros((0,0,0))
            xi_r_ref_im = np.zeros((0,0,0))
            lag_L_ref_re = np.zeros((0,0,0))
            lag_L_ref_im = np.zeros((0,0,0))
            l_max = 0
            k_max = 0
            Omega_rot = 0.
            Omega_orb = 0.
            
        else:
            
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

###

class Photosphere:
    
    def __init__ (self, resp_coeffs, photgrid, dx):
        
        self.resp_coeffs = resp_coeffs
        
        self.coeffs = {}
        self.read_phot_coeffs(photgrid, dx) 
    
    
    def read_phot_coeffs(self, photgrid, dx):
        """Selects appropriate photometric coefficient routine
        depending on the stellar model
        """

        if str(type(photgrid))=="<class 'pymsg.PhotGrid'>":
            self.coeffs.update( self.read_phot_coeffs_msg(photgrid, dx) )
        elif isinstance(photgrid, str):
            self.coeffs.update( self.read_phot_coeffs_h5(photgrid, dx) )
        elif isinstance(photgrid, None):
            raise Exception("Must specify photgrid")
        else: raise Exception(f"Invalid photgrid type")
        
        return
    

    def read_phot_coeffs_msg(self, photgrid, dx):
        """Creates and stores photometric coefficients from
        MSG for the desired filter
        """
        # Set atmosphere parameters dict
        dx = {'logT': np.log10(dx['Teff']), 'logg': dx['logg']}
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


    def R_xl(self, l):
        I = self.coeffs
        return (2 + l)*(1 - l)*I[f'I_x'][l] / I[f'I_x'][0]

    def T_xl(self, l):
        I = self.coeffs
        return I[f'dI_dlnT_x'][l] / I[f'I_x'][0]

    def G_xl(self, l):
        I = self.coeffs
        return I[f'dI_dlng_x'][l] / I[f'I_x'][0]

    def disk_intg_factor (self, l):
        I = self.coeffs
        b_l = I[f'I_x'][l]/I[f'I_x'][0]
        return b_l



    def read_phot_coeffs_h5(self, synspec_model):
        """This function is deprecated and was tested on a
        synspec intensity spectrum calculated from
        iOri's tlusty atmosphere.
        """
        # Read intensity file
        I = self.read_intensity(synspec_model)
        I = I.data

        # Evaluate photometric data
        I_x_l = {}
        dI_dlnT_x_l = {}
        dI_dlng_x_l = {}
       
        # Use filename as filter name
        filter_x = synspec_model.split(",")[0]

        for l in range(Response.resp_coeffs.data['l_max']+1):
            I_x_l[l] = I[f'I_{filter_x}_{l}'][0]
            dI_dlnT_x_l[l] = I[f'dlnTeff_{filter_x}_{l}'][0]
            dI_dlng_x_l[l] = I[f'dlng_{filter_x}_{l}'][0]

        return {f'I_{filter_x}': I_x_l,
                f'dI_dlnT_{filter_x}': dI_dlnT_x_l,
                f'dI_dlng_{filter_x}': dI_dlng_x_l}
    

    
    def read_intensity(self, filename):
       
        self.meta = self.read_info(filename)

        return h5py.File(filename, 'r')
    
    
    def read_info(self, filename):
        
        moments = self.read_intensity(filename)
        colors = []
        ells = []

        for moment in list(moments):
            words = moment.split('_')
            
            if len(words)==2:
                if words[-1] not in colors:
                    colors.append(words[-1])
                    
            if len(words)==3:
                if int(words[-1]) not in ells:
                    ells.append(int(words[-1]))
                    
        l_min, l_max = min(ells), max(ells)
        info = {}
        
        for color in colors:
            info[color] = {'l_min':l_min, 'l_max':l_max}
        
        return info
    
            
    def find_coeffs(self, C, x, l=None):
        
        I = self.data
        
        def find_coeffs_C(I, C, x, l):
            if C=='R': 
                return (2 + l)*(1 - l)*I[f'I_{x}_{l}'][:] / I[f'I_{x}_0'][:]
            if C=='T':
                return I[f'dlnTeff_{x}_{l}'][:] / I[f'I_{x}_0'][:]
            if C=='G':
                return I[f'dlng_{x}_{l}'][:] / I[f'I_{x}_0'][:]
        C_x = {} 
        
        if l==None:
            n_l = self.info[x]['l_max'] # - self.info[x]['l_min'] 
            #this may break if l_min =/= 0
        
            for l in np.arange(n_l+1,dtype=int):
                C_x[l] = find_coeffs_C(I,C,x,l)[0]
            return C_x
        
        elif isinstance(l,int):
            C_x[l] = find_coeffs_C(I,C,x,l)[0]
            return C_x
        
        else:
            raise Exception('TypeError: l must be int or int array')

###
