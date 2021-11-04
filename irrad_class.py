import numpy as np
from scipy.integrate import quad as integrate
from scipy.special import lpmv, sph_harm
from scipy.optimize import fsolve

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


    def find_disk_intg_factor (self, star_number, filter_x, l):

        b_l = self.component[star_number].disk_intg_factor(filter_x, l)

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


    def setup_irrad (self, star_number):

        if int(star_number)==1:
            star_neighbor = 2
        elif int(star_number)==2:
            star_neighbor = 1
        else:
            raise Exception('Star unspecified')

        phot_coeffs = self.component[star_number].phot_coeffs
        resp_data = self.component[star_number].resp_coeffs.data

        L1 = self.component[star_number].params['L']
        R1 = self.component[star_number].params['R']
        L2 = self.component[star_neighbor].params['L']

        return phot_coeffs, resp_data, L1, R1, L2


    def find_irrad (self, star_number, filter_x, theta, phi, t, t_peri=0):

        # Set up for sum

        phot_coeffs, resp_data, L1, R1, L2 = self.setup_irrad(star_number)

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
