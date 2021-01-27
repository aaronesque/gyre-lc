# gyre-lc


Run with
> ./gyre-lc.py bin_params.in

A description of the various modules follows.

### gyre-lc.py

- takes inlist 'bin_params.in' and returns timeseries flux data to screen.

To-do:
- [ ] hook to color grid maker via 'colormoment()' function
- [ ] print data to file
- [ ] get 'run_plotter()' to work

### obs_func.py

- creates `observer` class with functions
- - `find_star_flux()`
- - `find_flux()`

### orbit_func.py

- creates `binary` class
- creates `irradiation` class

### star_func.py

- creates `star` class

### atm_coeffs.py

- creates `atm_coeffs` class

### resp_coeffs.py

- creates `resp_coeffs` class


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
