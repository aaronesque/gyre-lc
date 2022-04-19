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
