#!/usr/bin/env python
#
# Unit testing for grid.py

import sys
import numpy as np
import grid as gr

def test_locate_nodes (grid):

    print('Testing locate() at nodes')

    # Check locate at nodes

    for i in range(grid.n_Teff):

        for j in range(grid.n_logg):

            i_loc, j_loc = grid.locate(grid.Teff_axis[i], grid.logg_axis[j])

            i_chk = i-1 if i == grid.n_Teff-1 else i
            j_chk = j-1 if j == grid.n_logg-1 else j
            
            if i_loc != i_chk or j_loc != j_chk:
                raise Exception(f'locate() index mismatch: ({i_loc},{j_loc}) != ({i_chk},{j_chk})')

def test_locate_centers (grid):

    print('Testing locate() at centers')

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


def test_find_neighbors (grid, verbose=True):

    print('Testing find_neighbors()')

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
                                and isinstance(grid.nodes[i_chk,j_chk], gr.Node):  correct = True
                        
                        else:  correct = False

                        # Next, check nbrs[index] == correct at index

                        if nbrs[ni+1,nj+1] == correct:  
                            if verbose:  print(f'{nbrs[ni+1,nj+1]}=={correct}')
                        
                        else:  raise Exception(f'Invalid result: {nbrs[ni+1,nj+1]} != {correct}')

                    except Exception as err_i:
                        print(err_i)
                        pass
                            
            if verbose:  print(' ')


if __name__ == '__main__':

    filenames = sys.argv[1:]

    if len(filenames) == 0:
        raise Exception("Syntax: grid.py [filename list]")

    # Construct grid

    grid = gr.Grid(filenames)

    # Run tests

    test_locate_nodes(grid)
    test_locate_centers(grid)
    test_find_neighbors(grid, verbose=False)
