# gyre-lc

GYRE-lc is a Python library for the production of synthetic light curves for tidally distorted binary systems. It relies on GYRE-tides (Sun et al., 2021), an as-yet unreleased extension of GYRE to model stellar tides. It requires at least one MESA stellar model and its corresponding GYRE tide model as inputs. A model spectrum is also required– GYRE-lc works best with MSG interpolated spectra for speed, ease of use, accuracy, and reliability.

A description of the various modules follows. See https://gyre-lc.readthedocs.io/en/latest/ for the latest version of documentation.

### star.py

- creates `Star` class
- creates `Response` class
- creates `Photosphere` class

### binary.py

- creates `Binary` class
- creates `Irradiation` class

### observer.py

- creates `Observer` class with functions
  - `find_flux()`
