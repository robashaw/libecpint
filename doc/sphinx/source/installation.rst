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
	
If you are a developer looking to make changes to the code, please fork the repo into your own version, and make a pull request when you think your changes are production ready. We will not accept any attempts to push directly into master.

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
| LIBECPINT_MAX_L      | 5             | The maximum angular momentum (in the orbital and ECP |
|                      |               | basis) that the code will be able to handle. The     |
|                      |               | higher this is, the longer the build stage will take,|
|                      |               | although not significantly so. NOTE: If you want     |
|                      |               | derivatives of up to ``L``, this must be ``L+n``     |
|                      |               | where ``n`` is the order of derivative (1 or 2).     |
+----------------------+---------------+------------------------------------------------------+
| LIBECPINT_MAX_UNROL  | 1             | NOW REDUNDANT. The max. angular momentum the code is |
|                      |               | unrolled up to. Increasing this will make the build  |
|                      |               | much slower but no longer gives any noticeable       |
|                      |               | advantage. We *strongly* recommend leaving this at   |
|                      |               | its default value.                                   |
+----------------------+---------------+------------------------------------------------------+

	  
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

Regenerating the radial code
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The recursive radial integral code has been pre-generated, as the current setting has been calibrated to balance accuracy and efficiency. If you would like to experiment (warning: after reading the paper), go into the directory ``src/generated/radial``. Edit the top line of  ``unrol_radial.py`` to change ``MAX_UNROL_AM``, the maximum angular momentum to be unrolled. Then do the following:

.. code-block:: bash
	
	python3 unrol_radial.py
	./generate.sh

This will generate the simplified recursive integrals and then piece together the ``radial_gen.cpp`` file and place it in the correct location. It should be very safe (but not very efficient) to decrease ``MAX_UNROL_AM``, but be prepared for things to break if you increase it too much. 


Testing
=======

To run all the tests, in the build directory run

.. code-block:: bash
	
	make test

This will give results that look as follows:

.. code-block:: bash

	 Running tests...
	 Test project [build-dir]
	       Start  1: MathUtil
	  1/16 Test  #1: MathUtil .........................   Passed    0.01 sec
	       Start  2: MultiArray
	  2/16 Test  #2: MultiArray .......................   Passed    0.00 sec
	       Start  3: Bessel
	  3/16 Test  #3: Bessel ...........................   Passed    0.02 sec
	       Start  4: GaussianShell
	  4/16 Test  #4: GaussianShell ....................   Passed    0.01 sec
	       Start  5: GaussianECP
	  5/16 Test  #5: GaussianECP ......................   Passed    0.02 sec
	       Start  6: GaussQuad
	  6/16 Test  #6: GaussQuad ........................   Passed    0.01 sec
	       Start  7: Generator
	  7/16 Test  #7: Generator ........................   Passed    0.00 sec
	       Start  8: IntTest1
	  8/16 Test  #8: IntTest1 .........................   Passed    0.04 sec
	       Start  9: IntTest2
	  9/16 Test  #9: IntTest2 .........................   Passed    0.02 sec
	       Start 10: DerivTest1
	 10/16 Test #10: DerivTest1 .......................   Passed    0.04 sec
	       Start 11: DerivTest2
	 11/16 Test #11: DerivTest2 .......................   Passed    0.12 sec
	       Start 12: HessTest1
	 12/16 Test #12: HessTest1 ........................   Passed    0.12 sec
	       Start 13: HessTest2
	 13/16 Test #13: HessTest2 ........................   Passed    0.11 sec
	       Start 14: APITest1
	 14/16 Test #14: APITest1 .........................   Passed    0.03 sec
	       Start 15: APITest2
	 15/16 Test #15: APITest2 .........................   Passed    0.13 sec
	       Start 16: Type1Test
	 16/16 Test #16: Type1Test ........................   Passed    0.03 sec
	 100% tests passed, 0 tests failed out of 16
	 

If any of these tests fails, the reasons for the failure can then be found in ``build/Testing/Temporary/LastTest.log``. This should help you troubleshoot and problems. If you have followed these instructions and the tests are still failing, please raise an issue on Github, giving details of your environment, any options you gave, and the relevant contents of LastTest.log. 

Stress test
^^^^^^^^^^^

To see how efficient your build is, you can make an additional test as follows:

.. code-block:: bash

	make StressTest
	
This will compute integrals and derivatives for increasing clusters of silver atoms, giving timings for each. At the end it will give estimated scaling exponents. Our latest build (clang 12.0.0, -O3 flag, single core on MacBook Pro 2017) gave these results:

.. code-block:: bash

	N: 2
	  Initialisation...   done. TIME TAKEN:        0.158172 seconds
	       Integrals...   done. TIME TAKEN:        0.116014 seconds
	      1st derivs...   done. TIME TAKEN:        0.321756 seconds
	      2nd derivs...   done. TIME TAKEN:        0.785262 seconds
	N: 4
	  Initialisation...   done. TIME TAKEN:        0.138962 seconds
	       Integrals...   done. TIME TAKEN:         1.23218 seconds
	      1st derivs...   done. TIME TAKEN:         4.16765 seconds
	N: 6
	  Initialisation...   done. TIME TAKEN:        0.128387 seconds
	       Integrals...   done. TIME TAKEN:         4.34048 seconds
	      1st derivs...   done. TIME TAKEN:         15.7498 seconds
	Scaling of integrals: N**3.31
	Scaling of 1st derivs: N**3.56

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
