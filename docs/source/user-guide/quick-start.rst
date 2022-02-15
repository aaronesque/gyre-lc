.. _quick-start:

==============================
Quick Start
==============================

GYRE-lc presumes a basic familiarity with introductory python, which includes the ability to call functions and make simple 1-dimensional plots.

To get started with GYRE-lc, follow these steps:

- download & install the `MESA Software Development Kit (SDK) <http://www.astro.wisc.edu/~townsend/static.php?ref=mesasdk>`_;
- download, install, & setup `MSG <http://www.astro.wisc.edu/~townsend/resource/docs/msg/>`_;
- download & unpack the `GYRE-lc <https://github.com/aaronesque/gyre-lc>`_ source code;
- set the ``GYRELC_DIR`` environment variable to point to the newly created source directory;
- in a Python environment, import standard packages ``sys`` and ``os`` 
- finally, implement GYRE-lc with ``sys.path.insert(0, os.path.join(os.environ['GYRELC_DIR'], 'lib'))`` and ``import gyrelc``

For a more in-depth installation guide, refer to the Installation chapter. Otherwise, proceed to the next chapter where you’ll learn to run your first GYRE-lc calculation.


.. note:: This project is under active development.

