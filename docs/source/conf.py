# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
import re

sys.path.insert(0, os.path.abspath('exts'))
sys.path.insert(0, os.path.abspath('../../src/'))

# -- Project information -----------------------------------------------------

project = 'GYRE-lc'
author = 'Aaron Lopez'
version = 'master'
release = 'master'
copyright = '2022, Aaron Lopez'
author = 'Aaron Lopez'



# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
#extensions = ['sphinxcontrib.bibtex']
#bibtex_bibfiles = ['refs.bib']
#bibtex_reference_style = 'author_year'

# Numbered figures
numfig = True

extensions = [
        'sphinx_rtd_theme',
        'sphinx.ext.extlinks',
        'sphinx.ext.intersphinx',
        'sphinx.ext.autodoc',
        'sphinx.ext.mathjax',
        'sphinx.ext.napoleon',
        'ads_cite'
        ]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme' #'alabaster'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# -- Options for EPUB output -------------------------------------------------
epub_show_urls = 'footnote'

# -- Additional configuration ------------------------------------------------

# Set master doc
master_doc = 'index'

# Mathjax + LaTeX config
macros = {}

with open('macros.def', encoding='utf-8') as f:
    line = f.readline()
    while line:
        macro, defn = line.rstrip().split('\t')
        macros[macro] = defn
        line = f.readline()

mathjax_macros = {}

for macro, defn in macros.items():
    argnums = re.findall('#(\d)', defn)
    if argnums:
        mathjax_macros[macro] = [defn, int(max(argnums))]
    else:
        mathjax_macros[macro] = defn

mathjax3_config = { 
    'loader': {
        'load': ['[tex]/mathtools']
    },
    'tex': {
        'packages': {
            '[+]': ['mathtools']
        },
        'macros': mathjax_macros
    }
}

latex_macros = ''

for macro, defn in macros.items():
    argnums = re.findall('#(\d)', defn)
    if argnums:
        latex_macros += f'\\def\\{macro}[#{int(max(argnums))}]{{{defn}}}\n'
    else:
        latex_macros += f'\\def\\{macro}{{{defn}}}\n'

latex_elements = {
    'preamble': '\\usepackage{mathtools}\n'+latex_macros
}

# Set logo
html_logo = 'gyre-lc-logo-1.png'

# Set up intersphinx
intersphinx_mapping = {
    'numpy': ('http://docs.scipy.org/doc/numpy/', None),
    'matplotlib': ('https://matplotlib.org/stable/', None)
}

# Set up extlinks
extlinks = {'ads': ('https://ui.adsabs.harvard.edu/abs/%s/abstract', '')}

# Set up autodoc
autoclass_content = 'class'
autodoc_member_order = 'bysource'

# Set up napoleon
napoleon_google_docstring = True
napoleon_include_init_with_doc = True

