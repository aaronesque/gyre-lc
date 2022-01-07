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

===================
&star_1 and &star_2
===================
	
``bc_model_type``
  Type of binary component model to use:
  
  - ``'MESA'`` : Evolutionary model read from a MESA GYRE-format text file
  - ``'PT_MASS'`` : Point mass model

``bc_model_path``
  Path to binary component model file (when ``bc_model_type`` = ``'MESA'``)

``tide_model_path``
  Path to tide model file (when ``bc_model_type`` = ``'MESA'``). Note that no "tide_model_type" parameter exists currently-- GYRE-lc remains exclusively compatible with GYRE-format pulsation models

``mass``
  This component's mass (when ``bc_model_type`` = ``'PT_MASS'``).

``mass_units``
  The 
  

=======
&orbit
=======



============
&observer
============

.. note:: This project is under active development.




