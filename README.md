# atmosphere-grid

### I_from_grid.py

- finds specific intensity from pre-computed (TLUSTY, SYNSPEC) grid of specific intensities

### I_from_func.py

- finds specific intensity from analytical mass, radius, $T_eff$ relations

### grid.py
- creates `Node` class
- creates `Grid` class

### test_grid.py
For unit testing.

# gyre-lc

### gyre-lc.py

- takes inlist 'bin_params.py' and returns timeseries flux data to screen.
- to-do: hook to color grid maker via 'colormoment()' function
- to-do: print data to file
- to-do: consolidate objects into fewer (or single) libraries
- to-do: get 'run_plotter()' to work
