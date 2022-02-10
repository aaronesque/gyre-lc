import numpy as np
from scipy.integrate import quad as integrate
from scipy.special import lpmv, sph_harm
from scipy.optimize import fsolve
from star_class import Star

### Class definitions

class Irradiation:
    
    def __init__ (self):
        """Constructor method
        """
        pass
            
    def eval_ramp_coeff (self, l, m):
        """
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
    
    
    def find_disk_intg_factor (self, star_number, filter_x, l):
        """
        """
        if self.component[star_number].luminosity==0.:
            return 0.
        else:
            b_l = self.component[star_number].disk_intg_factor(filter_x, l)
            return b_l
        
    
    def find_mean_anom (self, t, t_peri=0):
        """
        """
        return self.omega_orb*(t - t_peri)*(2*np.pi)
        
    
    def find_ecce_anom (self, M):
        """
        """
        K_soln = np.empty_like(M)
        for i, M_val in enumerate(M): 
            Keppler = lambda E : E - self.e*np.sin(E) - M_val
            K_soln[i] = fsolve(Keppler, 0)[0]
        return K_soln
        
    
    def find_true_anom (self, E):
        """
        """
        return 2*np.arctan( ((1+self.e)/(1-self.e))*np.tan(E/2) )
    
    
    def convert_t_to_f (self, t, t_peri=0):
        """
        """
        M = self.find_mean_anom(t, t_peri)
        E = self.find_ecce_anom(M)
        f = self.find_true_anom(E)
        return f
    
    
    def find_bin_sep (self, t, t_peri=0):
        """
        """
        f = self.convert_t_to_f(t, t_peri)
        D = self.a*(1-self.e**2)/(1+self.e*np.cos(f))
        return D
    
    
    def setup_irrad (self, star_number):
        """
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
    
    
    def find_irrad (self, star_number, filter_x, theta, phi, t, t_peri=0):
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
        
        # Set up intensity moment range
        
        l_min = 0
        dl = 1
        l_max = resp_data['l_max']
        n_l = np.ceil( (l_max - l_min)/dl ) + 1

        l_range = l_min + dl*np.arange(n_l)
        
        for l in l_range.astype(int):
            for m in np.arange(-l, l+1).astype(int):
                
                Z_lm = self.eval_ramp_coeff(l, m)
                Y_lm = sph_harm(m, l, phi, theta)
                
                Zt_lm = Z_lm*Y_lm*np.exp(-1j*m*ft) # Replace Z_lm*Y_lm w amplitude, ft*m w frequency?
                    
                b_l = self.find_disk_intg_factor(star_number, filter_x, l)
                
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
            omega_orb_units=False, a_units=False):
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


    def eval_fourier (self, filter_x, inc, omega, t=0, reflection=True):

        omega_orb = self.omega_orb

        f_1, A_1 = self.component[1].eval_fourier(filter_x, inc, omega)
        f_2, A_2 = self.component[2].eval_fourier(filter_x, inc, omega+180)

        #if reflection==True:
        #    rf_1, rA_1 = self.system.eval_fourier_irrad(1, self.filter_x, inc, omega, t, t_peri)
        #    rf_2, rA_2 = self.system.eval_fourier_irrad(2, self.filter_x, inc, omega+180, t, t_peri)

        #    return [f_1/omega_orb, rf_1, f_2/omega_orb, rf_2], [A_1, rA_1, A_2, rA_2]
        #else:

        return f_1/omega_orb, np.abs(A_1 + A_2)
