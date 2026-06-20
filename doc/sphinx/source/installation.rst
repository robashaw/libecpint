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

This will run unit tests and integration tests covering integrals, derivatives, Hessians, the high-level API, type-1 integrals, reference atom values, and shell screening. If any of these tests fails, the reasons for the failure can then be found in ``build/Testing/Temporary/LastTest.log``. This should help you troubleshoot any problems. If you have followed these instructions and the tests are still failing, please raise an issue on Github, giving details of your environment, any options you gave, and the relevant contents of LastTest.log.

Stress test
^^^^^^^^^^^

To see how efficient your build is, you can make an additional test as follows:

.. code-block:: bash

	make StressTest

This will compute integrals and derivatives for increasing clusters of silver atoms, giving timings for each. At the end it will give estimated scaling exponents. As of v1.1, shell screening reduces the effective scaling from roughly cubic to quadratic for larger molecules.

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
