# grid class

import numpy as np

import node as nd

### Class definition

class Grid:

    def __init__ (self, nodes_list, debug=False):

        # Sanity check

        if len(nodes_list) == 0:
            raise Exception('Empty nodes list')

        if not all(isinstance(node, nd.Node) for node in nodes_list):
            raise Exception('Nodes list contains invalid instance')

        # Extract axes

        self.Teff_axis = np.unique([node.Teff for node in nodes_list])
        self.logg_axis = np.unique([node.logg for node in nodes_list])

        self.n_Teff = len(self.Teff_axis)
        self.n_logg = len(self.logg_axis)

        # Store the nodes

        self.nodes = np.empty([self.n_Teff,self.n_logg], dtype=object)

        for node in nodes_list:

            # Find where the node belongs in the list

            i_Teff = np.where(node.Teff == self.Teff_axis)[0][0]
            i_logg = np.where(node.logg == self.logg_axis)[0][0]

            self.nodes[i_Teff,i_logg] = node

        # Set the debug flag

        self.debug = debug

        
    def lookup (self, Teff, logg):

        if Teff in self.Teff_axis and logg in self.logg_axis:

            # Find location in array

            i_Teff = np.where(Teff == self.Teff_axis)[0][0]
            i_logg = np.where(logg == self.logg_axis)[0][0]

            # Return the node if it exists

            node = self.nodes[i_Teff,i_logg]
            
            if isinstance(node, nd.Node):
                return node
            else:
                return None

        else:

            return None
        

    def show_topology (self):

        for i_logg in reversed(range(self.n_logg)):

            for i_Teff in range(self.n_Teff):

                if isinstance(self.nodes[i_Teff,i_logg], nd.Node):
                    print("X", end=' ')
                else:
                    print(".", end=' ')

            print()
            print()
            

    def locate (self, Teff, logg):
        
        # First check Teff, logg within grid

        if (Teff < self.Teff_axis[0]) or (Teff > self.Teff_axis[-1]):
            raise IndexError(f'Teff={Teff} outside axis range')
        if (logg < self.logg_axis[0]) or (logg > self.logg_axis[-1]):
            raise IndexError(f'logg={logg} outside axis range')

        # Find indices

        def locate_x (x, x_axis):

            # Use bisection to find i such that
            # x_axis[k] <= x < x_axis[k+1].

            # Special case if x == x_axis[-1]

            if x == x_axis[0]:

                k = 0

            elif x == x_axis[-1]:

                k = len(x_axis) - 2

            else:

                k_0 = -1
                k_1 = len(x_axis)

                while (k_1 - k_0) > 1:
    
                    k = (k_0 + k_1)//2

                    if x >= x_axis[k]:
                        k_0 = k
                    else:
                        k_1 = k

                k = k_0

            return k

        i = locate_x(Teff, self.Teff_axis)
        j = locate_x(logg, self.logg_axis)

        # Return tuple (i,j) if all is well

        if self.debug:
            print(f'\nLocated grid point (i,j)=({i},{j}) such that \n')
            print(f'Teff_axis[i]={self.Teff_axis[i]}  <=  Teff={Teff}  <  Teff_axis[i+1]={self.Teff_axis[i+1]}')
            print(f'logg_axis[j]={self.logg_axis[j]}  <=  logg={logg}  <  logg_axis[j+1]={self.logg_axis[j+1]}\n')

        return (i,j)


    def find_neighbors (self, ij, show=False):
        
        stencil = np.zeros([3,3], dtype=bool)
        
        # Find where node has neighbors
        
        for k in range(-1,2):
            
            for h in range(-1,2):
                
                j_Teff = ij[0] + h
                j_logg = ij[1] + k
                
                # if neighbor exists, store True in bool stencil
                
                if ((0 <= j_Teff) and (j_Teff < self.n_Teff)) and \
                ((0 <= j_logg) and (j_logg < self.n_logg)):
                
                    if isinstance(self.nodes[j_Teff,j_logg], nd.Node):
                        if show==True: print("X", end=' ')
                        stencil[h+1, k+1] = True
                    else:
                        if show==True: print(".", end=' ')
                    
            if show==True:    
                print()
                print()
        
        return stencil

    

    def find_derivs (self, ij, which_deriv, show=False):
  
        nbrs_stencil, nbrs_Teff, nbrs_logg  = self.recon_stencil(ij, show)

        def take_first_deriv (nbrs_data, nbrs_axis):

            x_0 = nbrs_axis[0]
            x_1 = nbrs_axis[-1]

            I_0 = nbrs_data[0]
            I_1 = nbrs_data[-1]

            return (I_1 - I_0)/(x_1 - x_0)

        def take_cross_deriv (nbrs_data, nbrs_Teff, nbrs_logg):

            dx = nbrs_Teff[-1] - nbrs_Teff[0]
            dy = nbrs_logg[-1] - nbrs_logg[0]

            I_a = nbrs_data[2][2]
            I_b = nbrs_data[0][2]
            I_c = nbrs_data[0][0]
            I_d = nbrs_data[2][0]

            return (I_a - I_b + I_c - I_d)/(dx*dy)
        
        if which_deriv=='dTeff':
            return take_first_deriv(nbrs_stencil[:,1], nbrs_Teff)        
        
        elif which_deriv=='dlogg': 
            return take_first_deriv(nbrs_stencil[1,:], nbrs_logg)    
        
        elif which_deriv=='cross':
            return take_cross_deriv(nbrs_stencil, nbrs_Teff, nbrs_logg)

        elif which_deriv=='dlnTeff':
            return take_first_deriv(nbrs_stencil[:,1], np.log(nbrs_Teff))
        
        elif which_deriv=='dlng':
            return take_first_deriv(nbrs_stencil[1,:], nbrs_logg)/np.log(10)

        else:
            raise Exception("Specify which deriv for find_deriv(): 'dTeff','dlogg','cross','dlnTeff','dlng'.")


    def recon_stencil (self, ij, show=False):

        # Grab the neighbor data

        nbrs = self.find_neighbors(ij)

        # Check whether we have enough info to reconstruct a stencil
        # covering the neighbors: at least a 2x2 square of neighbors
        # (including the point itself)

        if not (np.all(nbrs[0:2,0:2]) or np.all(nbrs[1:3,0:2]) or
                np.all(nbrs[0:2,1:3]) or np.all(nbrs[1:3,1:3])):
            return None, None, None

        # Perform the reconstruction

        # First, create the axes -- copied across or extrapolated from
        # the grid

        def create_axis (k, axis):
            if k == 0:
                return np.array([2*axis[0]-axis[1], axis[0], axis[1]])
            elif k == len(axis)-1:
                return np.array([axis[-2], axis[-1], 2*axis[-1]-axis[-2]])
            else:
                return axis[k-1:k+2]

        Teff_axis = create_axis(ij[0], self.Teff_axis)
        logg_axis = create_axis(ij[1], self.logg_axis)

        # Now reconstruct the data. First, copy over existing data

        data = np.empty([3,3], dtype=object)

        for di in range(-1,2):

            i = ij[0] + di

            for dj in range(-1,2):

                j = ij[1] + dj

                if nbrs[di+1,dj+1]:
                    data[di+1,dj+1] = self.nodes[i,j].data

        # Now fill in missing data using extrapolation

        def extrap_in_Teff(i, j, i_ex):
            assert nbrs_ex[i[0],j] and nbrs_ex[i[1],j]
            w = (Teff_axis[i_ex] - Teff_axis[i[0]])/(Teff_axis[i[1]] - Teff_axis[i[0]])
            return (1-w)*data_ex[i[0],j] + w*data_ex[i[1],j]

        def extrap_in_logg(i, j, j_ex):
            assert nbrs_ex[i,j[0]] and nbrs_ex[i,j[1]]
            w = (logg_axis[j_ex] - logg_axis[j[0]])/(logg_axis[j[1]] - logg_axis[j[0]])
            return (1-w)*data_ex[i,j[0]] + w*data_ex[i,j[1]]

        def extrap_in_both(i, j, i_ex, j_ex):
            return 0.5*(extrap_in_Teff(i, j_ex, i_ex) + extrap_in_logg(i_ex, j, j_ex))

        data_ex = np.copy(data)
        nbrs_ex = np.copy(nbrs)

        # Extrapolate face values

        # Bottom

        if not nbrs[1,0]:
            data_ex[1,0] = extrap_in_logg(1, [1,2], 0)
            nbrs_ex[1,0] = True
            if show:
                print('Extrapolated [1,0] face')

        # Right

        if not nbrs[2,1]:
            data_ex[2,1] = extrap_in_Teff([0,1], 1, 2)
            nbrs_ex[2,1] = True
            if show:
                print('Extrapolated [2,1] face')

        # Top

        if not nbrs[1,2]:
            data_ex[1,2] = extrap_in_logg(1, [0,1], 2)
            nbrs_ex[1,2] = True
            if show:
                print('Extrapolated [1,2] face')

        # Left

        if not nbrs[0,1]:
            data_ex[0,1] = extrap_in_Teff([1,2], 1, 0)
            nbrs_ex[0,1] = True
            if show:
                print('Extrapolated [0,1] face')

        # Extrapolate corner values (first using original corners)

        # Bottom-left

        if not nbrs[0,0]:
            if nbrs[2,0] and not nbrs[0,2]:
                data_ex[0,0] = extrap_in_Teff([1,2], 0, 0)
                if show:
                    print('Extrapolated [0,0] corner along Teff')
            elif not nbrs[2,0] and nbrs[0,2]:
                data_ex[0,0] = extrap_in_logg(0, [1,2], 0)
                if show:
                    print('Extrapolated [0,0] corner along logg')
            elif nbrs[2,0] and nbrs[0,2]:
                data_ex[0,0] = extrap_in_both([1,2], [1,2], 0, 0)
                if show:
                    print('Extrapolated [0,0] corner along both')
            nbrs_ex[0,0] = nbrs[2,0] or nbrs[0,2]

        # Bottom-right

        if not nbrs[2,0]:
            if nbrs[0,0] and not nbrs[2,2]:
                data_ex[2,0] = extrap_in_Teff([0,1], 0, 2)
                if show:
                    print('Extrapolated [2,0] corner along Teff')
            elif not nbrs[0,0] and nbrs[2,2]:
                data_ex[2,0] = extrap_in_logg(2, [1,2], 0)
                if show:
                    print('Extrapolated [2,0] corner along logg')
            elif nbrs[0,0] and nbrs[2,2]:
                data_ex[2,0] = extrap_in_both([0,1], [1,2], 2,  0)
                if show:
                    print('Extrapolated [2,0] corner along both')
            nbrs_ex[2,0] = nbrs[0,0] or nbrs[2,2]

        # Top-right

        if not nbrs[2,2]:
            if nbrs[0,2] and not nbrs[2,0]:
                data_ex[2,2] = extrap_in_Teff([0,1], 2, 2)
                if show:
                    print('Extrapolated [2,2] corner along Teff')
            elif not nbrs[0,2] and nbrs[2,0]:
                data_ex[2,2] = extrap_in_logg(2, [0,1], 2)
                if show:
                    print('Extrapolated [2,2] corner along logg')
            elif nbrs[0,2] and nbrs[2,0]:
                data_ex[2,2] = extrap_in_both([0,1], [0,1], 2,  2)
                if show:
                    print('Extrapolated [2,2] corner along both')
            nbrs_ex[2,2] = nbrs[0,2] or nbrs[2,0]

        # Top-left

        if not nbrs[0,2]:
            if nbrs[2,2] and not nbrs[0,0]:
                data_ex[0,2] = extrap_in_Teff([1,2], 2, 0)
                if show:
                    print('Extrapolated [0,2] corner along Teff')
            elif not nbrs[2,2] and nbrs[0,0]:
                data_ex[0,2] = extrap_in_logg(0, [0,1], 2)
                if show:
                    print('Extrapolated [0,2] corner along logg')
            elif nbrs[2,2] and nbrs[0,0]:
                data_ex[0,2] = extrap_in_both([1,2], [0,1], 0, 2)
                if show:
                    print('Extrapolated [0,2] corner along both')
            nbrs_ex[0,2] = nbrs[2,2] or nbrs[0,0]

        # Extrapolate the one remaining corner value

        if not nbrs_ex[0,0]:
            data_ex[0,0] = extrap_in_both([1,2], [1,2], 0, 0)
            nbrs_ex[0,0] = True 
            if show:
                print('Extrapolated [0,0] corner along both')
        elif not nbrs_ex[2,0]:
            data_ex[2,0] = extrap_in_both([0,1], [1,2], 2,  0)
            nbrs_ex[2,0] = True
            if show:
                print('Extrapolated [2,0] corner along both')
        elif not nbrs_ex[2,2]:
            data_ex[2,2] = extrap_in_both([0,1], [0,1], 2,  2)
            nbrs_ex[2,2] = True
            if show:
                print('Extrapolated [2,2] corner along both')
        elif not nbrs_ex[0,2]:
            data_ex[0,2] = extrap_in_both([0,1], [1,2], 2,  0)
            nbrs_ex[0,2] = True
            if show:
                print('Extrapolated [0,2] corner along both')

        # Sanity check
        
        assert np.all(nbrs_ex)

        # Return the data and axes

        return data_ex, Teff_axis, logg_axis

### Factory methods

def from_file (filenames_list, debug=False):

    # Read the nodes

    nodes_list = []

    for filename in filenames_list:
        nodes_list += [nd.from_file(filename)]

    # Return a new Grid

    return Grid(nodes_list, debug)


def from_func (Teff_axis, logg_axis, func, debug=False):

    # Create the nodes

    nodes_list = []

    for Teff in Teff_axis:
        for logg in logg_axis:
            nodes_list += [nd.Node(Teff, logg, func(Teff, logg))]

    # Return a new Grid

    return Grid(nodes_list, debug)
