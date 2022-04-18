# Import standard modules
import numpy as np
import sys
import os

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
            coefficients from :py:class:`pymsg.PhotGrid`
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
                self.resp_coeffs = resp_coeffs( self.gyre_model )
            else: self.resp_coeffs = resp_coeffs('')

        # else check for a point mass, read user-specified params
        elif mass is not None:
            # make note of model as point mass type
            self.model = f'{mass} M_sol'
            self.model_type = 'point mass'
            self.read_pt_mass_params()

            # no gyre_model allowed
            if gyre_model is not None:
                raise Exception("A point mass cannot experience tides")
            self.resp_coeffs = resp_coeffs('')

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
            if str(type(photgrid))=="<class 'pymsg.PhotGrid'>":
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

###

class Irradiation:
    """This class contains the routines necessary for modeling
    irradiation effects.

    This is especially important in close or high eccentricity
    binaries. :py:class:`gyrelc.Irradiation` is turned off by 
    default for :py:func:`gyrelc.Observer.find_flux()` to 
    preserve computational resources.

    Attributes:
        component_1 (:py:class:`gyrelc.Star`): Class representations of
            the binary components whose tidal interactions will be
            simulated and visualized
        component_2 (:py:class:`gyrelc.Star`): Class representations of
            the binary components whose tidal interactions will be
            simulated and visualized
        omega_orb (float): The binary's orbital angular velocity
        omega_orb_units (str): Units of ``omega_orb``. Default is 'CYC_PER_DAY'
        a (float): The binary separation
        a_units (str): Units of ``a``. Default is 'SOLAR' but user may specify
        'CGS' (centimeters) instead
        e (float): The binary eccentricity. :math:`0\leq e <1`
    
    """

    def __init__ (self):
        """Constructor method
        """
        pass
            
    def eval_ramp_coeff (self, l, m):
        """Evaluates the ramp function

        Args:
            l (int): degree of spherical harmonic
            m (int): order of spherical harmonic

        Returns:
            The float resulting from evaluating the ramp function :math:`Z_{lm}`
        """
        term1 = 2* np.sqrt( ((2*l+1)/(4*np.pi)) \
                           *(np.math.factorial(l-m)/np.math.factorial(l+m)) )
        if (np.abs(m)==1):
            term2 = np.pi/2
        else:
            term2 = np.cos(m*np.pi/2)/(1-m**2)
        
        term3, term3_err = integrate(lambda mu: np.sqrt(1-mu**2)*lpmv(m,l, mu), -1, 1)
        Z_lm = term1*term2*term3
        return Z_lm
    
    
    def find_disk_intg_factor (self, star_number, l):
        """Finds the disk integral factor :ads_citep:`Burkart:2012` for
        a given spherical degree l and passband filter x

        Args:
            star_number (int): index denoting primary or secondary star
            l (int): degree of spherical harmonic

        Returns:
            The disk integral factor :math:`b_l`, a float
        """
        if self.component[star_number].luminosity==0.:
            return 0.
        else:
            b_l = self.component[star_number].disk_intg_factor(l)
            return b_l
        
    
    def find_mean_anom (self, t, t_peri=0):
        """Finds the mean anomaly given the observation time,
        time at periastron, and orbital angular velocity

        Args:
            t (float, np.array): observation time(s)
            t_peri (float): the time at periastron

        Returns:
            The mean anomaly
        """
        return self.omega_orb*(t - t_peri)*(2*np.pi)
        
    
    def find_ecce_anom (self, M):
        """Finds the eccentric anomaly from the mean anomaly
        and orbital eccentricity

        Args:
            M (float, np.array): The mean anomaly

        Returns:
            The eccentric anomaly
        """
        K_soln = np.empty_like(M)
        for i, M_val in enumerate(M): 
            Keppler = lambda E : E - self.e*np.sin(E) - M_val
            K_soln[i] = fsolve(Keppler, 0)[0]
        return K_soln
        
    
    def find_true_anom (self, E):
        """Finds the true anomaly from the eccentric 
        anomaly and eccentricity

        Args:
            E (float, np.array): The eccentric anomaly

        Returns:
            The true anomaly
        """
        return 2*np.arctan( ((1+self.e)/(1-self.e))*np.tan(E/2) )
    
    
    def convert_t_to_f (self, t, t_peri=0):
        """Converts from observation time `t` to
        true anomaly `f`

        Args:
            t (float, np.array): The observation time(s)
            t_peri (float): The time of periastron

        Returns:
            The true anomaly
        """
        M = self.find_mean_anom(t, t_peri)
        E = self.find_ecce_anom(M)
        f = self.find_true_anom(E)
        return f
    
    
    def find_bin_sep (self, t, t_peri=0):
        """Finds the binary separation :math:`D(t)` for the
        binary at time `t`, given some `t_peri`

        Args:
            t (float, np.array): The observation time(s)
            t_peri (float): The time of periastron
           
        Returns:
            The binary separation at time(s) `t`
        """
        f = self.convert_t_to_f(t, t_peri)
        D = self.a*(1-self.e**2)/(1+self.e*np.cos(f))
        return D
    
    
    def setup_irrad (self, star_number):
        """Acquires the parameters from either :py:class:`gyrelc.Star`
        component necessary to evaluate the additional flux from
        the binary due to irradiation effects :ads_citep:`Burkart:2012`

        Args:
            star_number (int): The index denoting the primary or secondary component

        Returns:
            The irradiated star's intrinsic luminosity, its radius, its
            photometric coefficients, its tidal response coefficients,
            and its companion's luminosity
        """
        if int(star_number)==1:
            star_neighbor = 2
        elif int(star_number)==2:
            star_neighbor = 1
        else:
            raise Exception('Star unspecified')
        
        phot_coeffs = self.component[star_number].phot_coeffs
        resp_data = self.component[star_number].resp_coeffs.data
        
        L1 = self.component[star_number].luminosity
        R1 = self.component[star_number].radius
        L2 = self.component[star_neighbor].luminosity
        
        return phot_coeffs, resp_data, L1, R1, L2
    
    
    def find_irrad (self, star_number, theta, phi, t, t_peri=0):
        """
        """
        # Set up for sum
        
        phot_coeffs, resp_data, L1, R1, L2 = self.setup_irrad(star_number)
       
        if L1==0.:
            return np.zeros_like(t)
        else:
            pass
        
        Dt = self.find_bin_sep(t, t_peri)
        
        ft = self.convert_t_to_f(t, t_peri)
        
        rel_dJ = np.zeros_like(t)
        
        for l in range(resp_data['l_max']+1):
            for m in np.arange(-l, l+1).astype(int):
                
                Z_lm = self.eval_ramp_coeff(l, m)
                Y_lm = sph_harm(m, l, phi, theta)
                
                Zt_lm = Z_lm*Y_lm*np.exp(-1j*m*ft) # Replace Z_lm*Y_lm w amplitude, ft*m w frequency?
                    
                b_l = self.find_disk_intg_factor(star_number, l)
                
                rel_dJ += b_l*Zt_lm.real*(L2/L1)*(R1/Dt)**2
                
        return rel_dJ

###

class Binary(Irradiation):
    """This is a class representation of a tidally interacting
    binary system. 

    It is a combination of two :class:`gyrelc.Star`
    classes, augmented by several user-input binary parameters
    and the optional irradiation scheme (Burkart 2016?) as
    inherited by :class:`gyrelc.Irradiation`

    Attributes:
        component_1 (:py:class:`gyrelc.Star`): Class representations of
            the binary components whose tidal interactions will be
            simulated and visualized
        component_2 (:py:class:`gyrelc.Star`): Class representations of
            the binary components whose tidal interactions will be
            simulated and visualized
        omega_orb (float): The binary's orbital angular velocity
        omega_orb_units (str): Units of ``omega_orb``. Default is 'CYC_PER_DAY'
        a (float): The binary separation
        a_units (str): Units of ``a``. Default is 'SOLAR' but user may specify
        'CGS' (centimeters) instead
        e (float): The binary eccentricity. :math:`0\leq e <1`
    
    """

    def __init__ (self, component_1, component_2, omega_orb, a, e, 
            omega_orb_units=None, a_units=None):
        """Constructor method
        """
        self.omega_orb = omega_orb
        self.a = a
        self.e = e
        self.component = {1: component_1,\
                          2: component_2 }
        super().__init__()

    def read_units(self):
        """Unit routine
        """
        # assumed conversion factors
        self.__omega_orb_units = 1 # cyc / cyc
        self.__a_units =  6.957e10 # cm / R_sol

        if self.__dict__.get('omega_orb_units'):
            if self.__dict__['omega_orb_units']=='CYC_PER_DAY':
                self.__dict__['__omega_orb_units'] = 1
            else:
                raise Exception(f'user-specified `omega_orb` units unrecognized')
        else:
            self.__dict__['omega_orb_units'] = 'CYC_PER_DAY'

        if self.__dict__.get('a_units'):
            if self.__dict__['a_units']=='SOLAR':
                self.__dict__['__a_units'] = 6.957e10 #cm per r_sol
            else:
                raise Exception(f'user-specified `a` units unrecognized')
        else: self.__dict__['a_units'] = 'SOLAR'

        return


    def eval_fourier (self, inc, omega, t=0, reflection=True):

        omega_orb = self.omega_orb

        f_1, A_1 = self.component[1].eval_fourier(inc, omega)
        f_2, A_2 = self.component[2].eval_fourier(inc, omega+180)

        #if reflection==True:
        #    rf_1, rA_1 = self.system.eval_fourier_irrad(1, self.filter_x, inc, omega, t, t_peri)
        #    rf_2, rA_2 = self.system.eval_fourier_irrad(2, self.filter_x, inc, omega+180, t, t_peri)

        #    return [f_1/omega_orb, rf_1, f_2/omega_orb, rf_2], [A_1, rA_1, A_2, rA_2]
        #else:

        return f_1/omega_orb, np.abs(A_1 + A_2)

###

class Observer:
    """This is a class representation of observer, whose position,
    duration, and instrument determines the light curve observed.

    Attributes:
        system (:py:class:`gyrelc.Binary`): Class representations of
            the binary whose tidal interactions will be
            simulated and visualized
        inc (float): The binary's inclination relative to the observer
        omega (float): The binary's argument of periastron
    
    """
    def __init__ (self, star_system, inc, omega, photgrid=None):
        
        self.system = star_system
        self.inc = inc
        self.omega = omega
        self.photgrid = photgrid

        if isinstance(star_system, Binary):
            
            self.system_type = 'binary'
            # If photgrid not specified here, take from gyrelc.Star
            if self.photgrid is None:
                self.photgrid = self.system.component[1].photgrid

            if self.photgrid is not None:
                self.system.component[1].read_phot_coeffs(self.photgrid)
                self.system.component[2].read_phot_coeffs(self.photgrid)
            # If photgrid not specified at any point, raise Exception
            else: raise Exception('Input error: photgrid not specified')

        elif isinstance(star_system, Star):
            
            self.system_type = 'single'
            # If photgrid not specified here, take from gyrelc.Star
            if self.photgrid is None:
                self.photgrid = self.system.photgrid

            if self.photgrid is not None:
                self.system.read_phot_coeffs(self.photgrid)
            else: raise Exception('Input error: photgrid not specified')

        else: raise Exception(f'Input error: {star_system} must be of class Binary() or Star()')
        
    
    def convert_coords (self, inc, omega):
        theta = inc/180 * np.pi
        phi = (90-omega)/180 * np.pi
        return theta, phi
    
    
    def eval_flux_single (self, star, inc, omega, t, t_peri=0):
        
        theta, phi = self.convert_coords(inc, omega)
        f, A = star.eval_fourier(theta, phi)
        # Initialize the frequencies/amplitudes arrays
        
        if isinstance(t, np.ndarray):
            star_flux = np.zeros_like(t)
        elif isinstance(t, list):
            t = np.array(t)
            star_flux = np.zeros_like(t)
        else: 
            star_flux = 0.

        # Add contributions from each frequency component

        star_dict = {}

        for k in range(len(A)):
            
            star_flux += np.real(A[k] * np.exp(-1j*f[k]*2*np.pi*(t - t_peri)))
            #star_dict[k] = np.real(A[k] * np.exp(-1j*f[k]*2*np.pi*(t - t_peri)))
            
        return star_flux#, star_dict
    
    
    def eval_flux_binary (self, inc, omega, t, t_peri=0, reflection=True):
        
        resp_1 = self.eval_flux_single(self.system.component[1], inc, omega, t, t_peri)
        L1 = self.system.component[1].luminosity
        
        resp_2 = self.eval_flux_single(self.system.component[2], inc, omega+180, t, t_peri)
        L2 = self.system.component[2].luminosity
        
        if reflection==True:
            
            refl_1 = self.system.find_irrad(1, inc, omega, t, t_peri)
            refl_2 = self.system.find_irrad(2, inc, omega+180, t, t_peri)
            
            flux = L1*(refl_1 + resp_1) + L2*(refl_2 + resp_2)
        
        else:
            flux = L1*resp_1 + L2*resp_2

        return flux/(L1+L2)
    
    
    def find_flux (self, t, t_peri=0, reflection=True):
        
        if self.system_type=='binary':
            return self.eval_flux_binary(self.inc, self.omega, t, t_peri, reflection)
        elif self.system_type=='single':
            return self.eval_flux_single(self.system, self.inc, self.omega, t, t_peri)
        

    def find_fourier (self, t=0, reflection=True):
        
        if self.system_type=='binary':
            return self.system.eval_fourier(self.inc, self.omega, t, reflection)
        if self.system_type=='single':
            return self.system.eval_fourier(self.inc, self.omega)
