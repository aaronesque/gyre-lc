#!/usr/bin/env python
#
# Unit testing for grid.py

import sys
import numpy as np
import grid as gr
import node as nd

def data_func (Teff, logg):

    return Teff*logg**3 + logg**2


def deriv_func (Teff, logg):

    return logg**3, 3*Teff*logg**2 + 2*logg, 3*logg**2


def bound_func (Teff, logg):

    return Teff**4/10**logg/1E15 < 1.


def build_grid (type=None):

    # Set up axes

    if type == 'fine':
        Teff_axis = np.linspace(2500., 50000., 100)
        logg_axis = np.linspace(2.5, 4.5, 100)
    else:
        Teff_axis = [2500., 5000., 7500., 10000., 15000., 20000., 25000., 30000., 40000., 50000.]
        logg_axis = [2.5, 3.0, 3.5, 4.0, 4.5]

    # Build the grid

    grid = gr.from_func(Teff_axis, logg_axis, data_func, bound_func=bound_func)

    return grid
    

def test_locate_nodes ():

    print('Testing locate() at nodes')

    # Build the grid

    grid = build_grid()

    # Check locate at nodes

    for i in range(grid.n_Teff):

        for j in range(grid.n_logg):

            i_loc, j_loc = grid.locate(grid.Teff_axis[i], grid.logg_axis[j])

            i_chk = i-1 if i == grid.n_Teff-1 else i
            j_chk = j-1 if j == grid.n_logg-1 else j
            
            if i_loc != i_chk or j_loc != j_chk:
                raise Exception(f'locate() index mismatch: ({i_loc},{j_loc}) != ({i_chk},{j_chk})')

            
def test_locate_centers ():

    print('Testing locate() at centers')

    # Build the grid

    grid = build_grid()

    # Check locate at centers

    for i in range(grid.n_Teff-1):

        Teff = 0.5*(grid.Teff_axis[i] + grid.Teff_axis[i+1])

        for j in range(grid.n_logg-1):

            logg = 0.5*(grid.logg_axis[j] + grid.logg_axis[j+1])

            i_loc, j_loc = grid.locate(Teff, logg)

            i_chk = i
            j_chk = j

            if i_loc != i_chk or j_loc != j_chk:
                raise Exception(f'locate() index mismatch: ({i_loc},{j_loc}) != ({i_chk},{j_chk})')


def test_find_neighbors (verbose=False):

    print('Testing find_neighbors()')

    # Build the grid

    grid = build_grid()

    # Check find_neighbors
    
    for i in range(grid.n_Teff):

        for j in range(grid.n_logg):
            
            nbrs = grid.find_neighbors((i,j))

            if verbose:  print(f'Node({i},{j})')

            for ni in range(-1,2):

                for nj in range(-1,2):
                    
                    try:
                        
                        i_chk = i+ni
                        j_chk = j+nj

                        # First find the correct answer to: does a node exist at this index?

                        if verbose:  print(f'checking neighbor({ni+1},{nj+1}):', end=' ')

                        if (i_chk >= 0) and (i_chk < grid.n_Teff) and (j_chk >= 0) and (j_chk < grid.n_logg)\
                                and isinstance(grid.nodes[i_chk,j_chk], nd.Node):  correct = True
                        
                        else:  correct = False

                        # Next, check nbrs[index] == correct at index

                        if nbrs[ni+1,nj+1] == correct:  
                            if verbose:  print(f'{nbrs[ni+1,nj+1]}=={correct}')
                        
                        else:  raise Exception(f'Invalid result: {nbrs[ni+1,nj+1]} != {correct}')

                    except Exception as err_i:
                        print(err_i)
                        pass
                            
            if verbose:  print(' ')


def test_recon_stencil(verbose=False):

    print('Testing recon_stencil()')

    # Build the grid

    grid = build_grid()

    # Check recon_stencil
    
    for i in range(grid.n_Teff):

        for j in range(grid.n_logg):
            
            if verbose:  print(f'Node({i},{j})')

            data, Teff_axis, logg_axis = grid.recon_stencil((i,j))


def test_find_derivs (verbose=False):

    print('Testing find_derivs()')

    # Build the grid

    grid = build_grid(type='fine')

    # Set the error tolerances

    tol_dTeff = 1E-2
    tol_dlogg = 1E-2
    tol_cross = 1E-2

    # Check find_derivs

    for i in range(grid.n_Teff):

        for j in range(grid.n_logg):
            
            if verbose:  print(f'Node({i},{j})')

            if isinstance(grid.nodes[i,j], nd.Node):

                ddata_dTeff = grid.find_derivs((i,j), 'dTeff')
                ddata_dlogg = grid.find_derivs((i,j), 'dlogg')
                ddata_cross = grid.find_derivs((i,j), 'cross')

                Teff = grid.Teff_axis[i]
                logg = grid.logg_axis[j]

                ddata_dTeff_chk, ddata_dlogg_chk, ddata_cross_chk = deriv_func(Teff, logg)

                err_dTeff = (ddata_dTeff - ddata_dTeff_chk)/ddata_dTeff_chk
                err_dlogg = (ddata_dlogg - ddata_dlogg_chk)/ddata_dlogg_chk
                err_cross = (ddata_cross - ddata_cross_chk)/ddata_cross_chk

                if err_dTeff > tol_dTeff:
                    raise Exception(f'dTeff derivative error {err_dTeff} exceeds tolerance at node ({i},{j})')
                
                if err_dlogg > tol_dlogg:
                    raise Exception(f'dlogg derivative error {err_dlogg} exceeds tolerance at node ({i},{j})')
                
                if err_cross > tol_cross:
                    raise Exception(f'cross derivative error {err_cross} exceeds tolerance at node ({i},{j})')
                
                if verbose:
                    print('  errs:', err_dTeff, err_dlogg, err_cross)


if __name__ == '__main__':

    # Run tests

    test_locate_nodes()
    test_locate_centers()
    test_find_neighbors()
    test_recon_stencil()
    test_find_derivs()

