import numpy as np
import f90nml as nml
from scipy.integrate import quad as integrate
from scipy.special import lpmv, sph_harm
from scipy.optimize import fsolve
from star_class import Star

### Class definitions

class Irradiation:

    def __init__ (self):
        pass

    def eval_ramp_coeff (self, l, m):

        term1 = 2* np.sqrt( ((2*l+1)/(4*np.pi)) \
                           *(np.math.factorial(l-m)/np.math.factorial(l+m)) )

        if (np.abs(m)==1):
            term2 = np.pi/2
        else:
            term2 = np.cos(m*np.pi/2)/(1-m**2)

        term3, term3_err = integrate(lambda mu: np.sqrt(1-mu**2)*lpmv(m,l, mu), -1, 1)

        Z_lm = term1*term2*term3

        return Z_lm


    def find_disk_intg_factor (self, comp_number, filter_x, l):

        b_l = self.component[comp_number].disk_intg_factor(filter_x, l)

        return b_l


    def find_mean_anom (self, t, t_peri=0):

        return self.orbit_params['Omega_orb']*(t - t_peri)*(2*np.pi)


    def find_ecce_anom (self, M):

        K_soln = np.empty_like(M)

        for i, M_val in enumerate(M):

            Keppler = lambda E : E - self.orbit_params['e']*np.sin(E) - M_val
            K_soln[i] = fsolve(Keppler, 0)[0]

        return K_soln


    def find_true_anom (self, E):

        return 2*np.arctan( ((1+self.orbit_params['e'])/(1-self.orbit_params['e']))*np.tan(E/2) )


    def convert_t_to_f (self, t, t_peri=0):

        M = self.find_mean_anom(t, t_peri)

        E = self.find_ecce_anom(M)

        f = self.find_true_anom(E)

        return f


    def find_bin_sep (self, t, t_peri=0):

        f = self.convert_t_to_f(t, t_peri)

        D = self.orbit_params['a']*(1-self.orbit_params['e']**2)/(1+self.orbit_params['e']*np.cos(f))

        return D


    def setup_irrad (self, comp_number):

        if int(comp_number)==1:
            comp_neighbor = 2
        elif int(comp_number)==2:
            comp_neighbor = 1
        else:
            raise Exception('Component unspecified')

        phot_coeffs = self.component[comp_number].phot_coeffs
        resp_data = self.component[comp_number].resp_coeffs.data

        L1 = self.component[comp_number].params['L']
        R1 = self.component[comp_number].params['R']
        L2 = self.component[comp_neighbor].params['L']

        return phot_coeffs, resp_data, L1, R1, L2


    def find_irrad (self, comp_number, filter_x, theta, phi, t, t_peri=0):

        if self.component[comp_number].inlist['comp_model_type']=='PT_MASS':
            return t*0
        else: None

        # Set up for sum

        phot_coeffs, resp_data, L1, R1, L2 = self.setup_irrad(comp_number)

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

                b_l = self.find_disk_intg_factor(comp_number, filter_x, l)

                rel_dJ += b_l*Zt_lm.real*(L2/L1)*(R1/Dt)**2

        return rel_dJ

###

class Binary(Irradiation):
    
    def __init__ (self, inlist_path):
        
        # should check if 'inlist', this just checks 'str'
        
        if isinstance(inlist_path, str):
            self.nml = nml.read(inlist_path)
        else: raise Exception('Inlist file error')
        
        self.orbit_params = self.read_orbit()
        
        self.component = {1: Star(inlist_path, 1),\
                          2: Star(inlist_path, 2) }
        
        super().__init__()
        
    def read_orbit(self):
        
        orbit_inlist = self.nml['orbit']
        
        params = {}
        
        #placeholder unit routine
        
        if orbit_inlist['omega_orb_units']=='CYC_PER_DAY':
            omega_orb_units = 1 
        else: omega_orb_units = 1
        
        if orbit_inlist['a_units']=='CM':
            a_units = 6.957e10 #cm per r_sol
        else: a_units = 1
            
        # eventually, we want to get these from resp_data
        # resp_data currently does not contain this output tho
        
        params['a'] = orbit_inlist['a']*a_units
        params['e'] = orbit_inlist['e']
        params['Omega_orb'] = orbit_inlist['omega_orb']*omega_orb_units
        
        return params
    
    
    def eval_fourier (self, filter_x, inc, omega, t=0):
            
        omega_orb = self.orbit_params['Omega_orb']
        
        f_1, A_1 = self.component[1].eval_fourier(filter_x, inc, omega)
        f_2, A_2 = self.component[2].eval_fourier(filter_x, inc, omega+180)
    
        return f_1/omega_orb, np.abs(A_1) + np.abs(A_2)
