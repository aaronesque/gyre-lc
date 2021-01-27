#!/home/alopez/anaconda3/bin/python

import h5py
import numpy as np
import f90nml as nml
from shutil import copyfile
import sys
import os
import tempfile as tmp
import atlas as at



#**** Functions ****#

##### Read and list params from model grid file

def leaf_param_list(model_file):

    atm = at.read_atmos(model_file)
    
    dirs, Teffs, loggs = [],[],[]
    
    for leaf,params in atm.items():
        
        dirs.append(leaf)
        
        for param,val in params.items():
            
            if param=='Teff':
                Teffs.append(float(val))
                
            if param=='logg':
                loggs.append(float(val))
                
    return dirs, Teffs, loggs


## Make axes on the fly

def make_axis(param_list, reverse=False):

    # Identify unique values for axes ticks, sort descending order

    param_axis = np.array(list(set(param_list)))

    param_axis.sort()

    if reverse==True:

        param_axis = param_axis[::-1]

    return param_axis


## Indexing params

def index_param_lists (list1, list2):

    #* Identifying data abscissa + range

    axis1 = make_axis(list1)
    axis2 = make_axis(list2)

    #* Mapping data coords to index coords

    # first list value pairs as coordinates

    coord_list = [(i,j) for i in list1 for j in list2]

    # Identify indices

    i_list = [i for coord in coord_list for i,t_ in enumerate(axis1) if t_==coord[0] ]
    j_list = [j for coord in coord_list for j,g_ in enumerate(axis2) if g_==coord[1] ]

    return i_list, j_list


# Algorithm that returns an index pair's coordinate in the param grid

def return_params( index_pair, x_axis,y_axis, Teff_axis,logg_axis):
    
    x,y = index_pair
    
    ### NOTE: we should be able to remove x_axis,y_axis dependence here

    index_grid = [(i,j) for i in x_axis for j in y_axis]
    
    # Check index pair exists in grid
    
    if (x,y) in index_grid:
        
        Teff = [Teff_ for i, Teff_ in enumerate(Teff_axis) if i==x][0]
        
        logg = [logg_ for i, logg_ in enumerate(logg_axis) if i==y][0]
        
        return (Teff,logg)
        
    else: 
        
        print(f'Atmosphere (i_Teff={x}, i_logg={y}) missing in grid.')
        
        return None


# Algorithm that verifies point's existence in the index list

def check_index( index_pair, Teff_axis, logg_axis, dirs, verbose=False ):
    
    x,y = index_pair

    Teff, logg = Teff_axis[x], logg_axis[y]

    leaf = make_leaf(Teff, logg)
    
    if leaf in dirs:
        
        if (verbose==True):
            print(f'Atmosphere EXISTS at index ({x},{y})')
            
        return True
    
    else:
        
        if (verbose==True):
            print(f'Atmosphere MISSING at index ({x},{y})')
        
        return False

# Algorithm that returns an index pair's coordinate in the param grid

def return_index( Teff, logg, Teff_axis, logg_axis ):
    
    #* Creating index grid for neighbor identification algorithm

    coord_grid = [(i,j) for i in Teff_axis for j in logg_axis]
    
    # Check index pair exists in index list
    if (Teff,logg) in coord_grid:
        
        x = [i for i, Teff_ in enumerate(Teff_axis) if Teff==Teff_][0]
        
        y = [i for i, logg_ in enumerate(logg_axis) if logg==logg_][0]
        
        return (x,y)
        
    else: 
        
        print(f'Atmosphere (Teff={Teff}, logg={logg}) missing in grid.')
        
        return None

        
# Algorithm that verifies point's existence in the index list

def check_params( coord_pair, Teff_list, logg_list, verbose=False ):
    
    Teff,logg = coord_pair
    
    if (Teff,logg) in list(zip(Teff_list, logg_list)):
        
        if (verbose==True):
            print(f'Atmosphere EXISTS at param coordinate ({Teff},{logg})')
            
        return True
    
    else:
        
        if (verbose==True):
            print(f'Atmosphere MISSING at param coordinate ({Teff},{logg})')
        
        return False


## Quickly make filename "leaf"

def make_leaf(Teff, logg):
    return f't{int(Teff):05d}g{int(logg*100)}'


## Identify and return nearest neighbors

def find_neighbors ( Teff, logg, Teff_list, logg_list, dirs, return_as='params'):
    
    # Identifying data abscissa + range
    
    Teff_axis = make_axis(Teff_list, reverse=True)
    logg_axis = make_axis(logg_list, reverse=True)

    # Mapping data coords to index coords

    x_list, y_list = index_param_lists(Teff_axis, logg_axis)

    # Identify unique indices, sort ascending order

    x_axis = make_axis(x_list)
    y_axis = make_axis(y_list)
                          
    # Index param coordinate
                          
    x,y = return_index( Teff, logg, Teff_axis, logg_axis)
    
    
    neighbors = {}
                          
    # Check index pair exists in index list
    if check_index( (x,y), Teff_axis, logg_axis, dirs):
        None        
    else:
        print('Index pair does not match any Teff, logg in the atmosphere list.')
        return
        
    for j in [-1,0,1]:
        for i in [-1,0,1]:
            
            # Check neighbor point exists
            if check_index( (x+i,y+j), Teff_axis, logg_axis, dirs):
                neighbors[f'{1+i}{1+j}'] = (x+i,y+j)
            else: # we're at a boundary. use self as point, or omit altogether?
                neighbors[f'{1+i}{1+j}'] = (x,y)
            
    if return_as=='params':
        for k in neighbors.keys():
            index_pair = neighbors[k]
            neighbors[k] = return_params( index_pair , x_axis,y_axis, Teff_axis,logg_axis )
            
    return neighbors


## Takes first derivative of I wrt Teff and logg

def first_derivs(I, neighbors, color):

    Teff = I.attrs['Teff']
    logg = I.attrs['logg']

    a_ = {'dx':{}, 'dy':{}, 'dxdy':{}}

    diff_ = {'dx':{}, 'dy':{}, 'dxdy':{}}

    derivs_ = {'dx':{}, 'dy':{}, 'dxdy':{}}

    # d/dx

    Teff_ = { k : neighbors[k][0] for k in neighbors.keys() if k in ['01','21','10','12']}
    
    logg_ = { k : neighbors[k][1] for k in neighbors.keys() if k in ['01','21','10','12']}

   # d/dx
    
    if abs( Teff_['21'] - Teff_['01'] ) > 0:
        
        for i in [-1,1]:
            
            a_['dx'][i] = abs(Teff_[f'{1+i}1'] - Teff)/abs(Teff_['21'] - Teff_['01'])
            
            a = a_['dx'][i]
                
            diff_['dx'][i] = Teff_[f'{1+i}1'] - Teff
                
            with tmp.TemporaryDirectory() as tdir:
                    
                leaf = make_leaf(Teff_[f'{1+i}1'],logg_[f'{1+i}1'])
                    
                orig_file = f'../{leaf}/color_moments.h5'
                tmp_file = f'{tdir}/{leaf}.h5'
                copyfile( orig_file , tmp_file )
                    
                I_i = h5py.File(tmp_file, 'r')
                
                derivs_['dx'][i] = np.sign(a) * a * (I_i[f'I_{color}'][:] - I[f'I_{color}'][:])/diff_['dx'][i]
                
                I_i.close()
                
        derivs_['dx'] = np.sum( [derivs_['dx'][i] for i in [-1,1]], axis=0)
        
    else:
        
        derivs_['dx'] = np.zeros_like( I_i[f'I_{color}'][:] )
    
    # d/dy
    
    if abs( logg_['12'] - logg_['10'] ) > 0:
        
        for i in [-1,1]:
            
            a_['dy'][i] = abs(logg_[f'1{1+i}'] - logg)/abs(logg_['12'] - logg_['10'])
            
            a = a_['dy'][i]
                
            diff_['dy'][i] = logg_[f'1{1+i}'] - logg
                    
            with tmp.TemporaryDirectory() as tdir:
                    
                leaf = make_leaf(Teff_[f'1{1+i}'],logg_[f'1{1+i}'])
                    
                orig_file = f'../{leaf}/color_moments.h5'
                tmp_file = f'{tdir}/{leaf}.h5'
                copyfile( orig_file , tmp_file )
                    
                I_i = h5py.File(tmp_file, 'r')
                
                derivs_['dy'][i] = np.sign(a) * a * (I_i[f'I_{color}'][:] - I[f'I_{color}'][:])/diff_['dy'][i]
                
                I_i.close()
                    
        derivs_['dy'] = np.sum( [derivs_['dy'][i] for i in [-1,1]], axis=0)
                
    else:
        
        derivs_['dy'] = np.zeros_like( I_i[f'I_{color}'][:] )
    
    return derivs_, a_, diff_



## Takes cross derivative of I wrt Teff and logg

def cross_derivs(I_, neighbors, weights, diffs, color):

    Teff = I_.attrs['Teff']
    logg = I_.attrs['logg']
    I = I_[f'I_{color}'][:]

    a_ = weights

    derivs_ = {'dx':{}, 'dy':{}, 'dxdy':{}}

    Teff_ = { k : neighbors[k][0] for k in neighbors.keys() if k in ['00','02','20','22']}

    logg_ = { k : neighbors[k][1] for k in neighbors.keys() if k in ['00','02','20','22']}

    for i in [-1,1]:
        for j in [-1,1]:

            if (a_['dx'][i] and a_['dy'][j]) > 0:

                dTeff = np.abs(Teff_[f'{1+i}{1+j}'] - Teff)
                dlogg = np.abs(logg_[f'{1+i}{1+j}'] - logg)

                if (dTeff and dlogg) != 0:

                    ai, aj = a_['dx'][i], a_['dy'][j]

                    with tmp.TemporaryDirectory() as tdir:

                        leaf = make_leaf(Teff_[f'{1+i}{1+j}'],logg_[f'{1+i}{1+j}'])
                        orig_file = f'../{leaf}/color_moments.h5'
                        tmp_file = f'{tdir}/{leaf}.h5'
                        copyfile( orig_file , tmp_file )

                        I_ij_ = h5py.File(tmp_file, 'r')
                        I_ij = I_ij_[f'I_{color}'][:]

                        leaf = make_leaf(Teff_[f'{1+i}0'],logg_[f'{1+i}0'])
                        orig_file = f'../{leaf}/color_moments.h5'
                        tmp_file = f'{tdir}/{leaf}.h5'
                        copyfile( orig_file , tmp_file )

                        I_i0_ = h5py.File(tmp_file, 'r')
                        I_i0 = I_i0_[f'I_{color}'][:]

                        leaf = make_leaf(Teff_[f'0{1+j}'],logg_[f'0{1+j}'])
                        orig_file = f'../{leaf}/color_moments.h5'
                        tmp_file = f'{tdir}/{leaf}.h5'
                        copyfile( orig_file , tmp_file )

                        I_0j_ = h5py.File(tmp_file, 'r')
                        I_0j = I_0j_[f'I_{color}'][:]

                        derivs_['dxdy'][i+j*1j] = ai*aj*i*(I_ij - I_i0 - I_0j + I)/(dTeff*dlogg)

                        I_ij_.close()
                        I_i0_.close()
                        I_0j_.close()

                else:

                    derivs_['dxdy'][i+j*1j] = np.zeros_like(I)
            else:

                derivs_['dxdy'][i+j*1j] = np.zeros_like(I)

    derivs_['dxdy'] = np.sum( [derivs_['dxdy'][i] for i in derivs_['dxdy'].keys()], axis=0)

    return derivs_['dxdy']



## Takes a log derivative

def ln_derivs(I, neighbors, color, l):

    Teff = I.attrs['Teff']
    logg = I.attrs['logg']

    a_ = {'dlnx':{}, 'dlny':{}}

    diff_ = {'dlnx':{}, 'dlny':{}}

    derivs_ = {'dlnx':{}, 'dlny':{}}

    Teff_ = { k : neighbors[k][0] for k in neighbors.keys() if k in ['01','21','10','12']}
    
    logg_ = { k : neighbors[k][1] for k in neighbors.keys() if k in ['01','21','10','12']}
    
    # d/dx
    
    if abs( Teff_['21'] - Teff_['01'] ) > 0:
        
        for i in [-1,1]:
            
            a_['dlnx'][i] = abs(Teff_[f'{1+i}1'] - Teff)/abs(Teff_['21'] - Teff_['01'])
            
            a = a_['dlnx'][i]
                
            diff_['dlnx'][i] = np.log(Teff_[f'{1+i}1']) - np.log(Teff)
            
            with tmp.TemporaryDirectory() as tdir:
                    
                leaf = make_leaf(Teff_[f'{1+i}1'],logg_[f'{1+i}1'])

                I_i = h5py.File('../'+leaf+'/color_moments.h5', 'r')
                
                derivs_['dlnx'][i] = np.sign(a) * a *\
                (I_i[f'I_{color}_{l}'][:] - I[f'I_{color}_{l}'][:])/diff_['dlnx'][i]
                
                I_i.close()
                
        derivs_['dlnx'] = np.sum( [derivs_['dlnx'][i] for i in [-1,1]], axis=0)
        
    else:
        
        derivs_['dlnx'] = np.zeros_like( I_i[f'I_{color}_{l}'][:] )
    
    # d/dy
    
    if abs( logg_['12'] - logg_['10'] ) > 0:
        
        for i in [-1,1]:
            
            a_['dlny'][i] = abs(logg_[f'1{1+i}'] - logg)/abs(logg_['12'] - logg_['10'])
            
            a = a_['dlny'][i]
                
            diff_['dlny'][i] = np.log(10.0)*(logg_[f'1{1+i}'] - logg)
                
            with tmp.TemporaryDirectory() as tdir:
                    
                leaf = make_leaf(Teff_[f'1{1+i}'],logg_[f'1{1+i}'])
                I_i = h5py.File('../'+leaf+'/color_moments.h5', 'r')
                
                derivs_['dlny'][i] = np.sign(a) * a *\
                (I_i[f'I_{color}_{l}'][:] - I[f'I_{color}_{l}'][:])/diff_['dlny'][i]
                
                I_i.close()
                
        derivs_['dlny'] = np.sum( [derivs_['dlny'][i] for i in [-1,1]], axis=0)
                
    else:
        
        derivs_['dlny'] = np.zeros_like( I_i[f'I_{color}'][:] )
    
    return derivs_


##### Main function

def main(inlist_file, model_file):

    ## Inlist

    if os.path.exists(inlist_file):
        #print(f'opening {inlist_file}')
        None
    else:
        print(f'error: inlist {inlist_file} not found')
        return

    inlist = nml.read(inlist_file)

    input_h5 = inlist['params']['output_filename']
    output_h5 = inlist['params']['partials_filename']

    copyfile(input_h5,output_h5)

    # Read input intensity file

    I_ = h5py.File(input_h5, 'r')

    # Create and read ouptut intensity file
    
    dI_ = h5py.File(output_h5, 'r+')

    # Read stellar parameters

    Teff = inlist['star']['teff']
    logg = inlist['star']['logg']

    # Identify leafs + values

    dirs, Teff_list, logg_list = leaf_param_list(model_file)

    # Find neighbors

    neighbors = find_neighbors(Teff, logg, Teff_list, logg_list, dirs)

    # Loop through colors, loop through modes,
    # find partials of I wrt Teff, logg, and
    # for each mode

    for color in inlist['color']:

        x = color['filter']

        derivs, weights, diffs = first_derivs(I_, neighbors, x)

        derivs['dxdy'] = cross_derivs(I_, neighbors, weights, diffs, x)

        # Create storage for partial colors

        if f'dTeff_{x}' not in list(dI_):
            dI_.create_dataset(f'dTeff_{x}',
                                   (I_[f'I_{x}'].shape[0],) )

        if f'dlogg_{x}' not in list(dI_):
            dI_.create_dataset(f'dlogg_{x}',
                                   (I_[f'I_{x}'].shape[0],) )

        if f'dTeff-dlogg_{x}' not in list(dI_):
            dI_.create_dataset(f'dTeff-dlogg_{x}',
                                   (I_[f'I_{x}'].shape[0],) )

        dI_[f'dTeff_{x}'][:] = derivs['dx']
        dI_[f'dlogg_{x}'][:] = derivs['dy']
        dI_[f'dTeff-dlogg_{x}'][:] = derivs['dxdy']

        ## Next, find partial color moments

        # Check if multiple or single l given
        
        if 'l_min' and 'l_max' in color:
            lrange = [*range(color['l_min'],color['l_max']+1)]
            
        elif 'l' in color:
            lrange = color['l']
            
        else:
            print('Please specify Legendre moments l')
            return
            
        # For all moments l,
        
        for l in lrange:
            
            # Use log identities to facilitate finite diff calculations

            derivs = ln_derivs(I_, neighbors, x, l)

            # ... create storage for partial quantities
            if f'dlng_{x}_{l}' not in list(dI_):
                dI_.create_dataset(f'dlng_{x}_{l}',
                                   (I_[f'I_{x}_{l}'].shape[0],) )

            if f'dlnTeff_{x}_{l}' not in list(dI_):
                dI_.create_dataset(f'dlnTeff_{x}_{l}',
                                   (I_[f'I_{x}_{l}'].shape[0],) )

            # ... and store intensity moments for posterity

            if f'I_{x}_{l}' in list(dI_):
                #del dI_[f'I_{x}_{l}']
                None

            # Finally, save results

            dI_[f'dlnTeff_{x}_{l}'][:] = derivs['dlnx']
            dI_[f'dlng_{x}_{l}'][:] = derivs['dlny']
           

    # Close h5 files

    dI_.close()
    I_.close()

    # Completed without errors for log

    leaf = make_leaf(Teff, logg)
    with open(f'{leaf}.sterr', 'w') as err:
        err.write('done')

    print('done')

    return


########
#**** Run ****#

inlist_file = sys.argv[1]
model_file = sys.argv[2]

try:
    main(inlist_file, model_file)
except Exception as e:
    f = nml.read(inlist_file)
    leaf = make_leaf(f['star']['Teff'],f['star']['logg'])
    with open(f'{leaf}.sterr','w') as err:
        err.write(f'Error: \n {e}')
