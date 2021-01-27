# gyre-lc

### gyre-lc.py

Run with
> ./gyre-lc.py bin_params.in

- takes inlist 'bin_params.in' and returns timeseries flux data to screen.

To-do:
- [ ] hook to color grid maker via 'colormoment()' function
- [ ] print data to file
- [ ] consolidate objects into fewer (or single) libraries
- [ ] get 'run_plotter()' to work


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
