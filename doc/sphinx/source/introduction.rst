.. libecpint intro file

.. _`sec:introduction`:

=============
Introduction
=============  

Overview
========

libecpint is a C++ library for the efficient evaluation of integrals over ab initio effective core potentials, using a mixture of generated, recursive code and Gauss-Chebyshev quadrature. It is designed to be standalone and generic, but is currently in development and may not be completely stable. If you experience any problems please raise an issue here; contributions and suggestions are also welcome.

This assumes ECPs and basis sets of the form usually seen in electronic structure calculations, namely those expanded in terms of Gaussian functions. The angular momentum of function that can be treated is in theory arbitrary, but is limited by your choice of maximum when the library is built.

Citing libecpint
================

If you use this library in your program and find it helpful, that's great! Any feedback would be much appreciated. If you publish results using this library, please consider citing the following paper detailing the implementation:

R. A. Shaw, J. G. Hill, J. Chem. Phys. 147, 074108 (2017); doi: [10.1063/1.4986887](http://dx.doi.org/10.1063/1.4986887)

and the JOSS paper that presents the library:

R. A. Shaw, J. G. Hill, J. Open Source Softw. 6, 3039 (2021); doi: [10.21105/joss.03039](https://doi.org/10.21105/joss.03039)

Full bibtex citations can be found in CITATION in the main directory.

Please also cite the ECPs and basis sets you use. 

Requirements
============

For the library
^^^^^^^^^^^^^^^

- A modern C++ compiler, at least C++11 standard library is required. This has currently only been tested with GCC (6.3.0 and above, but will in theory work with any > 4.9) and clang (9.0.0 and above). Intel compilers have been known to cause issues.
- CMake/CTest build tools v. >= 3.12 
- Python (2.7 or above, including 3 and higher)

For the docs
^^^^^^^^^^^^

- Doxygen
- Sphinx
- Breathe
- Exhale

For radial code regeneration
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- Python 3.6 or above
- numpy
- sympy

License
=======

libecpint is available under an MIT License, allowing for free and open use, reproduction, and modification of the library, so long as the copyright and license notices are preserved. The authors hold no liability for, and give no warranty against, results of the use of this software.

Support
=======

If you have any problems or would like to make suggestions for improvements, please raise an issue on the github repo. We will endeavour to get back to you as soon as possible, but as "we" is predominantly just "me" (Robert), it may take a while. 

Help is always welcome, and if you wish to make contributions to the code yourself, please take a look at the library API docs and have a go. Send a pull request with any enhancements!

.. toctree::
   :hidden:
