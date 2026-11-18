.. libecpint install file

.. _`sec:installation`:

=====================================
Installation
=====================================   

Obtaining libecpint
===================

The latest stable release of libecpint can always be found at the Github Repo_. 

.. _Repo: https://www.github.com/robashaw/libecpint.git

It can be downloaded directly from there, or you can clone it locally using git with the command

.. code-block:: bash

	git clone https://github.com/robashaw/libecpint.git 
	
If you are a developer looking to make changes to the code, please fork the repo into your own version, and make a pull request when you think your changes are production ready. We will not accept any attempts to push directly into main.

Building
========

To build the library, do the following in the top of the source tree (out-of-source build is required!):

.. code-block:: bash

	mkdir build
	cd build
	cmake [options] ..
	make [-jn]

The ``-jn`` flag tells make to use ``n`` threads while compiling (e.g. ``-j4`` would use four threads), and is highly recommended if your computer can cope, as the generated code files can all be compiled independently of one another.

CMake Options
^^^^^^^^^^^^^

The options after cmake above can be included using the syntax `-DOPTION=value`. The pertinent options are as follows:

+----------------------+---------------+------------------------------------------------------+
| Option       	       | Default       |  Description                                         |    
+======================+===============+======================================================+
| CMAKE_CXX_FLAGS      | N/A           | Flags to pass to the C++ compiler. These will depend |
|                      |               | on which compiler you're using, but in general we    |
|                      |               | *strongly* recommend passing optimisation flags,     |
|                      |               | specifically ``-O2`` or ``-O3``.                     |
+----------------------+---------------+------------------------------------------------------+
| CMAKE_INSTALL_PREFIX | /usr/local    | The directory where the library will be installed.   |
|                      |               | You must have permissions to edit this directory.    |
+----------------------+---------------+------------------------------------------------------+
| LIBECPINT_MAX_UNROL  | 1             | The max. angular momentum for which the type 2       |
|                      |               | angular integration is generated as unrolled code;   |
|                      |               | all other cases use general code. 1 is fastest: 0    |
|                      |               | is a few percent slower for s and p shells, while 2  |
|                      |               | makes the library slower (up to 2x for d shells,     |
|                      |               | 1.5x for some gradients) and the build ~20x longer;  |
|                      |               | 3 takes hours to build. Leave this at its default.   |
+----------------------+---------------+------------------------------------------------------+

The maximum angular momentum is fixed, and is no longer a build option (``LIBECPINT_MAX_L`` in
``config.hpp``): 8 for the orbital basis plus the order of derivative, so second derivatives up to
i functions, and 8 for the ECPs, including the local channel. The library throws
``std::invalid_argument`` beyond it.

	  
Documentation
^^^^^^^^^^^^^

This documentation can be generated locally via CMake by running

.. code-block:: bash
	
	make docs
	
This requires the following to be available:

- Doxygen
- Sphinx
- Breathe
- Exhale

Closed-form radial integrals
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The type 2 radial integrals are calculated by quadrature. The generator for the closed forms of
[Shaw2017], which earlier versions used for low angular momenta, is kept, unused, in
``src/generated/radial``; its ``README.md`` explains why it is not used and how to check it.


Testing
=======

To run the tests, in the build directory run

.. code-block:: bash

	ctest -j8 -LE reference

This runs unit tests of the library's components, tests of the radial integrals and of cases with closed forms, and checks of the integrals and their first and second derivatives for a set of test systems against values from an independent reference engine, stored in ``tests/data/reference``. If any test fails, the reason can be found in ``build/Testing/Temporary/LastTest.log``, or by running ``ctest --output-on-failure``; a failing comparison names the worst element by atom, shell and Cartesian component. If you have followed these instructions and the tests are still failing, please raise an issue on Github, giving details of your environment, any options you gave, and the relevant contents of LastTest.log.

Configuring with ``-DLIBECPINT_TEST_REFERENCE=ON`` adds slower tests that compare against the reference engine computed at test time (``ctest -L reference``). ``tests/README.md`` describes how the tests are organised, how to diagnose failures, and how to regenerate the stored values or add tests.

Benchmark
^^^^^^^^^

To see how efficient your build is, you can run the timing benchmark (not run by ``ctest``):

.. code-block:: bash

	make Benchmark

This times the shell-pair integrals for each combination of angular momenta, and the integrals and their first and second derivatives for several molecules, reporting the median and spread of repeated runs. ``benchmarks/ecpint_bench scaling`` times iodoalkane chains of increasing length, with estimated scaling exponents, and ``benchmarks/compare.py`` compares results between builds. ``benchmarks/README.md`` describes them.

Installation
============

To install the library and share directory, run

.. code-block:: bash

	make install
	
which will create the following files/directories:

.. code-block:: bash

	${CMAKE_INSTALL_PREFIX}/lib/libecpint.a
	${CMAKE_INSTALL_PREFIX}/include/libecpint.hpp
	${CMAKE_INSTALL_PREFIX}/include/libecpint
	${CMAKE_INSTALL_PREFIX}/share/libecpint


.. toctree::
   :hidden:
