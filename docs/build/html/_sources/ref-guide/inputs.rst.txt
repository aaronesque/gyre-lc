.. _inputs:

.. gyre-lc documentation master file, created by
   sphinx-quickstart on Tue Dec 14 13:12:55 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

#############################
Inputs
#############################

**********************
Namelist Input Files
**********************

GYRE-lc reads parameters from an *inlist*, which is an input file that defines a number of "namelist" groups. Inlists are an input format designed for Fortran, and GYRE-lc is presently written entirely in Python. However, the Fortran-heavy workflow of GYRE and MESA makes inlists a naturally convenient way for users to categorize groups of inputs for GYRE-lc.

========
&comp_1
========

The ``&comp_1`` namelist includes all parameters specifically relevant to the primary component of the binary. GYRE-lc does not currently support rapid rotation. 

``comp_model_type``
  Type of binary component model to use:
  
  - ``'MESA'`` : Evolutionary model read from a MESA GYRE-format text file
  - ``'PT_MASS'`` : Point mass model

``comp_model_path``
  Path to binary component model file (when ``comp_model_type`` = ``'MESA'``)

``tide_model_path``
  Path to tide model file (when ``comp_model_type`` = ``'MESA'``). Note that no "tide_model_type" parameter exists currently-- GYRE-lc remains exclusively compatible with GYRE-format pulsation models

``mass``
  This component's mass (when ``comp_model_type`` = ``'PT_MASS'``)

``mass_units``
  The units of the component's mass:
  
  - ``'SOLAR'`` : solar mass units, :math:`{M_\odot}` (default)
  - ``'CGS'`` : grams
 
``lum``
  This component's luminosity (when ``comp_model_type`` = ``'PT_MASS'``). If left blank, a component of zero luminosity and zero "reflection" is assumed

``lum_units``
  The units of the component's luminosity:

  - ``'SOLAR'`` : solar luminosity units, :math:`{L_\odot}` (default)
  - ``'CGS'`` : ergs per second

========
&comp_2
========

Same parameters available as ``&comp_1``, but for the second component of the binary.

=======
&orbit
=======

The ``&orbit`` namelist includes the parameters that specify the orbital configuration of the binary. 

``a``
  The binary separation distance; the length of the axis connecting the centers of mass of both components

``a_units``
  The units of the binary separation:
  
  - ``'SOLAR'`` : solar radii, :math:`{R_\odot}` (default)
  - ``'CGS'`` : centimeters
	
``e``
  The binary's orbital eccentricity (:math:`0 < e \leq 1`)
 
``omega_orb``
  The binary's orbital velocity

``omega_orb_units``
  The units of the orbital velocity

  - ``'CYC_PER_DAY'`` : cycles per day (default)

============
&observer
============

.. note:: This project is under active development.


