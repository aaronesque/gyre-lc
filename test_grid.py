#!/usr/bin/env python
#
# Unit testing for grid.py

import sys
import numpy as np
import grid as gr
import node as nd
import copy as cp

def data_func (Teff, logg):

    t = Teff/1e4
    g = logg

    return np.sin(t)*np.cos(g)


def derivs_func (Teff, logg):

    t = Teff/1e4
    g = logg

    dt = 1./1e4
    dg = 1.

    return {
        'T': np.cos(t)*np.cos(g)*dt,
        'g': -np.sin(t)*np.sin(g)*dg,
        'Tg': -np.cos(t)*np.sin(g)*dt*dg,
        'T_scale': dt,
        'g_scale': dg }


def bound_func (Teff, logg):

    return Teff**4/10**logg/1E15 < 1.


def build_grid (type=None, ragged=True):

    # Set up axes

    if type == 'fine':
        Teff_axis = np.linspace(2500., 50000., 1000)
        logg_axis = np.linspace(2.5, 4.5, 1000)
    elif type == '3x3':
        Teff_axis = [9000., 10000., 11000.]
        logg_axis = [3.9, 4.0, 4.1]
    else:
        Teff_axis = [2500., 5000., 7500., 10000., 15000., 20000., 25000., 30000., 40000., 50000.]
        logg_axis = [2.5, 3.0, 3.5, 4.0, 4.5]

    # Build the grid

    if ragged:
        grid = gr.from_func(Teff_axis, logg_axis, data_func, bound_func=bound_func)
    else:
        grid = gr.from_func(Teff_axis, logg_axis, data_func)

    return grid
    

def test_locate_nodes ():

    print('Testing locate() at nodes')

    # Build the grid

    grid = build_grid(ragged=True)

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

    grid = build_grid(ragged=True)

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

    grid = build_grid(ragged=True)

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

    grid = build_grid('3x3')

    # Check recon_stencil

    for k in range(512):

        # Copy the grid

        grid_copy = cp.deepcopy(grid)

        # Delete selected nodes

        for i in range(3):
            for j in range(3):
                if 2**(i+3*j) & k:
                    grid_copy.nodes[i,j] = None

        # Do the stencil reconstruction

        data, Teff_axis, logg_axis = grid_copy.recon_stencil((1,1))

        if data is not None:

            # Check the reconstructed values

            for i in range(3):
                for j in range(3):

                    # Calculate what the data should be

                    data_chk = data_func(Teff_axis[i], logg_axis[j])

                    # Evaluate the error

                    err = np.abs(data[i,j] - data_chk)

                    # Evaluate the expected error (using 2nd-order
                    # Taylor series estimates)

                    derivs = derivs_func(Teff_axis[1], logg_axis[1])

                    dTeff = np.abs(Teff_axis[i] - Teff_axis[1])
                    dlogg = np.abs(logg_axis[j] - logg_axis[1])

                    if i == 1 and j == 1:

                        # Center point; always should have zero error

                        err_chk = 0.

                    elif i == 1 and (j == 0 or j == 2):

                        # Face

                        err_chk = derivs['g_scale']**2*dlogg**2

                    elif (i == 0 or i == 2) and j == 1:

                        # Face

                        err_chk = derivs['T_scale']**2*dTeff**2

                    else:

                        # Corner

                        err_chk = 2*derivs['T_scale']*derivs['g_scale']*dTeff*dlogg + derivs['T_scale']**2*dTeff**2 + derivs['g_scale']**2*dlogg**2

                    # Check it is within tolerances

                    if err > err_chk:
                        raise Exception(f'Error {err} outside toerance {err_chk} at k={k},  i={i}, j={j}')



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
                    print(ddata_cross, ddata_cross_chk)
                    raise Exception(f'cross derivative error {err_cross} exceeds tolerance at node ({i},{j})')
                
                if verbose:
                    print('  errs:', err_dTeff, err_dlogg, err_cross)


if __name__ == '__main__':

    # Run tests

    test_locate_nodes()
    test_locate_centers()
    test_find_neighbors()
    test_recon_stencil()
#    test_find_derivs()

