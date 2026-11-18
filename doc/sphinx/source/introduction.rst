.. libecpint intro file

.. _`sec:introduction`:

=============
Introduction
=============

Overview
========

libecpint is a C++ library for the efficient evaluation of integrals over ab initio effective
core potentials (ECPs), and their first and second derivatives, over Gaussian basis functions. It
combines generated and recursive code for the angular parts with Gauss-Legendre quadrature for
the radial parts, and is designed to be standalone and generic. If you experience any problems
please raise an issue on GitHub; contributions and suggestions are also welcome.

It treats ECPs and basis sets of the form usual in electronic structure calculations, expanded in
Gaussian functions, with angular momentum up to 8 (l functions) for both, including the order of
derivative for the basis.

What's new in version 1.1
=========================

- **Faster**: typically 3--15x faster than v1.0.7 for whole molecules (7.5x on average over the
  benchmark), and 45x for the gradient of a 386-atom molecule, through a shared radial grid,
  screening at every level (including derivatives), and leaner angular integrals. See
  :ref:`sec:usage` (Performance) for the figures.
- **More accurate and stable**: fixed (non-adaptive) quadrature, more accurate Bessel functions
  and fixed numerical failures (e.g. type 1 integrals for steep ECPs, and off-centre type 2
  integrals for high angular momenta and diffuse functions) give integrals converged to ~1e-13 and
  smooth in the geometry. Every result is checked against an independent reference
  implementation.
- **New features**: gradients and Hessians contracted directly with a density, without storing the
  derivative matrices; derivatives with respect to a magnetic field (for gauge-including atomic
  orbitals); angular momentum up to 8, fixed rather than a build option; and a built-in ECP
  library expanded from 6 to 91 families, read from MOLPRO-format files.
- **Build changes**: some build options and dependencies have been removed; see
  :ref:`sec:upgrading`.

The public API (``ECPIntegrator`` and ``ECPIntegral``) is unchanged, apart from additions.

Python bindings
===============

`pylibecpint <https://github.com/robashaw/pylibecpint>`_ provides Python bindings for the main
functionality of this library, with all inputs and outputs as NumPy arrays, and is installable
via pip:

.. code-block:: bash

   pip install pylibecpint

It wraps the high-level ``ECPIntegrator`` API and the lower-level ``ECPIntegral`` shell-pair
interface. See its repository for documentation and examples.

Citing libecpint
================

If you use this library in your program and find it helpful, that's great! Any feedback would be
much appreciated. If you publish results using this library, please consider citing the paper
detailing the implementation:

R. A. Shaw, J. G. Hill, J. Chem. Phys. 147, 074108 (2017);
doi: `10.1063/1.4986887 <https://doi.org/10.1063/1.4986887>`_

and the paper presenting the library:

R. A. Shaw, J. G. Hill, J. Open Source Softw. 6, 3039 (2021);
doi: `10.21105/joss.03039 <https://doi.org/10.21105/joss.03039>`_

Full BibTeX citations can be found in ``CITATION`` in the main directory. Please also cite the
ECPs and basis sets you use: the files in the ECP library record the references for each atom.

Requirements
============

For the library
^^^^^^^^^^^^^^^

- A C++17 compiler, e.g. GCC (8 and above) or clang (7 and above). Intel compilers have been
  known to cause issues.
- CMake (3.12 and above)
- Python (any version), used to configure the code generator

For the docs
^^^^^^^^^^^^

- Doxygen
- Sphinx, with the packages in ``doc/sphinx/requirements.txt``

License
=======

libecpint is available under an MIT License, allowing for free and open use, reproduction, and
modification of the library, so long as the copyright and license notices are preserved. The
authors hold no liability for, and give no warranty against, results of the use of this software.

Support
=======

If you have any problems or would like to make suggestions for improvements, please raise an
issue on the GitHub repository. We will endeavour to get back to you as soon as possible, but as
"we" is predominantly just "me" (Robert), it may take a while.

Help is always welcome: if you wish to contribute to the code yourself, see
:ref:`sec:developer`, and send a pull request with any enhancements!

.. toctree::
   :hidden:
