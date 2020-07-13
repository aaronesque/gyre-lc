#!/usr/bin/env python

import grid
import sys
            
if __name__ == '__main__':

    filenames_list = sys.argv[1:]

    if len(filenames_list) == 0:
        raise Exception("Syntax: main.py [filename list]")

    grid = grid.from_file(filenames_list, debug=True)

    grid.show_topology()

    try:
        
        user_i = 3 #int(input('Enter i=7 for Teff=30000K: '))
        user_j = 3 #int(input('Enter j=3 for logg=4.0: '))
        print('Entering i={:d} for Teff={:.0f}K \nEntering j={:d} for logg={:3.1f}'.format(user_i, grid.Teff_axis[user_i], user_j, grid.logg_axis[user_j]))

        derivs = grid.find_derivs((user_i,user_j), 'cross', True)
        
        print(derivs)

        # Do data stenciling

        data, Teff_axis, logg_axis = grid.recon_stencil([user_i,user_j], show=True)

        print('Data from stenciling:\n', data)

    except IndexError:
        print("\nOpe, that right there's gonna be a problem. How's 'bout we try a different (Teff,logg)?\n")
        raise
