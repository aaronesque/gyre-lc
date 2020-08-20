#!/usr/bin/env python

from resp_delta import resp_coeffs
from atm_coeffs import atm_coeffs
from obs_func import observer
import sys

if __name__ == '__main__':
    #filename = sys.argv[1:]
    #if len(filename) == 0:
    #    raise Exception("Syntax: main.py [filename]")

    resp = resp_coeffs('../gyre_run/response.001.h5')
    atm = atm_coeffs('iOri-partials-2.h5', 'partials.in')

    inc = 62.9
    omega = 122.2
    x = 'B'

    obs = observer(resp, atm)

    obs.find_flux('B', inc,omega, -2.2)
