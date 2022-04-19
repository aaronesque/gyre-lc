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
