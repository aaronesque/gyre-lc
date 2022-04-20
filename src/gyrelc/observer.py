# Import standard modules
import numpy as np
import sys
import os
from gyrelc.binary import Binary
from gyrelc.star import Star

# Import special modules
import h5py
from scipy.special import lpmv, sph_harm
from scipy.integrate import quad as integrate
from scipy.optimize import fsolve
from astropy.io import ascii


### Class definitions

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

        if isinstance(star_system, Binary):
            self.system_type = 'binary'

        elif isinstance(star_system, Star):
            self.system_type = 'single'
        
        else: raise Exception(f'Input error: {star_system} must be of class Binary() or Star()')
        
    
    def convert_coords (self, inc, omega):
        theta = inc/180 * np.pi
        phi = (90-omega)/180 * np.pi
        return theta, phi
    
    
    def eval_flux_single (self, star, inc, omega, t, t_peri=0):
        
        theta, phi = self.convert_coords(inc, omega)
        f, A = self.eval_fourier(star, theta, phi)
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
            
        return star_flux
    
    
    def eval_flux_binary (self, inc, omega, t, t_peri=0, reflection=True):
        
        resp_1 = self.eval_flux_single(self.system.component[1], inc, omega, t, t_peri)
        L1 = self.system.component[1].params['luminosity']
        
        resp_2 = self.eval_flux_single(self.system.component[2], inc, omega+180, t, t_peri)
        L2 = self.system.component[2].params['luminosity']
        
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
        
###
    
    def eval_fourier (self, star, theta, phi):

        if self.__dict__.get('point_mass_model'):
            f, A = np.array([0]), np.array([0])
        else:
            resp_coeffs = star.resp_coeffs
            I = star.phot_coeffs

            # Initialize the frequencies/amplitudes arrays

            f = np.arange(resp_coeffs.data['k_max']+1)*resp_coeffs.data['Omega_orb']

            A = np.zeros(resp_coeffs.data['k_max']+1, dtype=complex)

            # Loop over l, m and k

            for l in np.arange(2, resp_coeffs.data['l_max']+1).astype(int):
                for m in np.arange(-l, l+1).astype(int):
                    for k in np.arange(0, resp_coeffs.data['k_max']+1).astype(int):
                        
                        # Add the Fourier contribution * spherical harmonic
                        A[k] += self.eval_fourier_moment(star, theta, phi, l,m,k)

        # Return data
        return f, A

###

    def eval_fourier_moment (self, star, theta, phi, l,m,k):

        i_l = l
        i_m = m + star.resp_coeffs.data['l_max']
        i_k = k

        dR_lmk = star.resp_coeffs.R_lmk[i_l,i_m,i_k]
        dT_lmk = star.resp_coeffs.T_lmk[i_l,i_m,i_k]
        dG_lmk = star.resp_coeffs.G_lmk[i_l,i_m,i_k]

        Y_lm = sph_harm(m, l, phi, theta)

        R_xl = star.phot_coeffs.R_xl(l)
        T_xl = star.phot_coeffs.T_xl(l)
        G_xl = star.phot_coeffs.G_xl(l)

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

###

    def eval_fourier_binary (self, inc, omega, t=0, reflection=True):

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

    def find_fourier (self, t=0, reflection=True):
        
        if self.system_type=='binary':
            return self.eval_fourier_binary(self.inc, self.omega, t, reflection)
        if self.system_type=='single':
            return self.eval_fourier_single(self.inc, self.omega)
