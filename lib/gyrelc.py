# Import standard modules

import numpy as np
import sys
import os

# Import pymsg

sys.path.insert(0, os.path.join(os.environ['MSG_DIR'], 'lib'))
import pymsg

# Import gyrelc modules

sys.path.insert(0, os.path.join(os.environ['GYRELC_DIR'], 'src'))
import atm_coeffs as ac
import resp_coeffs as rc
from star_class import Star
from binary_class import Irradiation, Binary


### Class definitions

class Observer:

    def __init__ (self, comp_system, filter_x, use_msg=True, phot_file=None):

        self.filter_x = filter_x
        self.use_msg = use_msg
        self.system = comp_system

        if isinstance(comp_system, Binary):

            self.system_type = 'binary'
            if use_msg==True:
                self.system.component[1].read_phot_coeffs(self.filter_x)
                self.system.component[2].read_phot_coeffs(self.filter_x)
            else:
                self.system.component[1].read_phot_coeffs(self.filter_x, 1)
                self.system.component[2].read_phot_coeffs(self.filter_x, 2)

        elif isinstance(comp_system, Star):

            self.system_type = 'single'
            if use_msg==True:
                self.system.read_phot_coeffs(self.filter_x)

        else: raise Exception(f'Input error: {comp_system} must be of class Binary() or Star()')


    def convert_coords (self, inc, omega):

        theta = inc/180 * np.pi
        phi = (90-omega)/180 * np.pi

        return theta, phi


    def eval_flux_single (self, star, inc, omega, t, t_peri=0):

        theta, phi = self.convert_coords(inc, omega)
        f, A = star.eval_fourier(self.filter_x, theta, phi)

        # Initialize the frequencies/amplitudes arrays

        if isinstance(t, np.ndarray):
            comp_flux = np.zeros_like(t)
        elif isinstance(t, list):
            t = np.array(t)
            comp_flux = np.zeros_like(t)
        else: 
            comp_flux = 0.

        # Add contributions from each frequency component

        for k in range(len(A)):
            comp_flux += np.real(A[k] * np.exp(-1j*f[k]*2*np.pi*(t - t_peri)))

        return comp_flux


    def eval_flux_binary (self, inc, omega, t, t_peri=0, reflection=True):

        resp_1 = self.eval_flux_single(self.system.component[1], inc, omega, t, t_peri)
        L1 = self.system.component[1].params['L']

        resp_2 = self.eval_flux_single(self.system.component[2], inc, omega+180, t, t_peri)
        L2 = self.system.component[2].params['L']

        if reflection==True:

            refl_1 = self.system.find_irrad(1, self.filter_x, inc, omega, t, t_peri)
            refl_2 = self.system.find_irrad(2, self.filter_x, inc, omega+180, t, t_peri)

            flux = L1*(refl_1 + resp_1) + L2*(refl_2 + resp_2)

        else:
            flux = L1*resp_1 + L2*resp_2

        return flux/(L1+L2)


    def find_flux (self, inc, omega, t, t_peri=0, reflection=True):

        if self.system_type=='binary':
            return self.eval_flux_binary(inc, omega, t, t_peri, reflection)
        elif self.system_type=='single':
            return self.eval_flux_single(self.system, inc, omega, t, t_peri)


    def find_fourier (self, inc, omega, t=0, reflection=True):

        if self.system_type=='binary':
            return self.system.eval_fourier(self.filter_x, inc, omega, t, reflection)
        if self.system_type=='single':
            return self.system.eval_fourier(self.filter_x, inc, omega)
