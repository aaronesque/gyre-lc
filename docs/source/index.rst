.. gyre-lc documentation master file, created by
   sphinx-quickstart on Tue Dec 14 13:12:55 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

###################
GYRE-lc
###################

**GYRE-lc** is a Python library for the production of synthetic light curves for tidally distorted binary systems. It relies on GYRE-tides :ads_citep:`Sun:2021`, an as-yet unreleased extension of GYRE to model stellar tides. It requires at least one `MESA <mesa.sourceforge.net>`_ stellar model and its corresponding `GYRE <https://gyre.readthedocs.io/en/stable/>`_ tide model as inputs. A model spectrum is also required-- GYRE-lc works best with `MSG <http://www.astro.wisc.edu/~townsend/resource/docs/msg/>`_ interpolated spectra for speed, ease of use, accuracy, and reliability.

.. make sure you include the build_spectrum script in the bundle

.. note:: This project is under active development.

.. toctree::
   :caption: User Guide
   :name: user-guide
   :maxdepth: 2

   user-guide/preliminaries.rst
   user-guide/quick-start.rst
   user-guide/python-walkthrough.rst
   
.. toctree::
   :caption: Reference Guide
   :name: ref-guide
   :maxdepth: 2

   ref-guide/installation.rst
   ref-guide/how-gyre-lc-works.rst
   ref-guide/python-interface.rst

