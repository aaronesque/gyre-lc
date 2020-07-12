#!/usr/bin/env python
#
# Unit testing for grid.py

import sys
import numpy as np
import grid as gr
import node as nd

def build_grid ():

    # Set up axes

    Teff_axis = [2500., 5000., 7500., 10000., 15000., 20000., 25000., 30000., 40000., 50000.]
    logg_axis = [2.5, 3.0, 3.5, 4.0, 4.5]

    # Define data funcs

    def func (Teff, logg):
        return Teff*logg**3 + logg**2

    def bound_func (Teff, logg):
        return Teff**4/10**logg/1E15 < 1.

    # Build the grid

    grid = gr.from_func(Teff_axis, logg_axis, func, bound_func=bound_func)

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


def test_find_neighbors (verbose=True):

    print('Testing find_neighbors()')

    # Build the grid

    grid = build_grid()

    # Check locate at neighbors
    
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



def take_first_deriv (nbrs_data, nbrs_axis):

    x_0 = nbrs_axis[0]
    x_1 = nbrs_axis[-1]

    I_0 = nbrs_data[0]
    I_1 = nbrs_data[-1]

    return (I_1 - I_0)/(x_1 - x_0)



def test_recon_stencil(ij):

    print('Testing recon_stencil()')

    nbrs_stencil, nbrs_Teff, nbrs_logg = grid.recon_stencil(ij)

    return


def test_find_derivs (ij):

    print('Testing find_derivs()')

    # Build the grid

    grid = build_grid()

    # Testing against true deriv of data_function() in create_grid.py

    nbrs_stencil, nbrs_Teff, nbrs_logg = grid.recon_stencil(ij)

    return


if __name__ == '__main__':

    # Run tests

    test_locate_nodes()
    test_locate_centers()
    test_find_neighbors(verbose=False)
    test_find_derivs((7,3))

