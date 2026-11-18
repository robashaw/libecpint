.. libecpint install file

.. _`sec:installation`:

=====================================
Installation
=====================================

Obtaining libecpint
===================

The latest stable release of libecpint can always be found on GitHub_. It can be downloaded
directly from there, or cloned with

.. _GitHub: https://github.com/robashaw/libecpint

.. code-block:: bash

   git clone https://github.com/robashaw/libecpint.git

If you are a developer looking to make changes to the code, please fork the repository, and make
a pull request when you think your changes are ready (see :ref:`sec:developer`).

Building
========

From the top of the source tree (an out-of-source build is required):

.. code-block:: bash

   mkdir build
   cd build
   cmake [options] ..
   make -j8
   ctest -j8 -LE reference   # optional, ~10 s
   make install

The ``-j`` flag sets the number of parallel build jobs.

CMake options
^^^^^^^^^^^^^

Options are passed to ``cmake`` as ``-DOPTION=value``:

+-------------------------------+-------------+--------------------------------------------------+
| Option                        | Default     | Description                                      |
+===============================+=============+==================================================+
| ``CMAKE_BUILD_TYPE``          | Release     | ``Release`` builds with optimisation (``-O3``),  |
|                               |             | which you almost certainly want.                 |
+-------------------------------+-------------+--------------------------------------------------+
| ``CMAKE_INSTALL_PREFIX``      | /usr/local  | Where ``make install`` puts the library.         |
+-------------------------------+-------------+--------------------------------------------------+
| ``LIBECPINT_BUILD_TESTS``     | ON          | Build the tests, the benchmark and the example.  |
+-------------------------------+-------------+--------------------------------------------------+
| ``LIBECPINT_TEST_REFERENCE``  | OFF         | Add the slower tests against the reference       |
|                               |             | engine computed at test time (see                |
|                               |             | :ref:`sec:developer`).                           |
+-------------------------------+-------------+--------------------------------------------------+
| ``LIBECPINT_BUILD_DOCS``      | ON          | Add the ``docs`` target (never built by default; |
|                               |             | see below).                                      |
+-------------------------------+-------------+--------------------------------------------------+
| ``LIBECPINT_RADIAL_CACHE_MB`` | 4           | Per-thread cache for the radial factors of       |
|                               |             | off-centre type 2 integrals, in MB: one of 1, 2, |
|                               |             | 4, 8, 16, 32, 64. Larger values can help large   |
|                               |             | basis sets; results do not depend on it.         |
+-------------------------------+-------------+--------------------------------------------------+
| ``LIBECPINT_MAX_UNROL``       | 1           | The largest angular momentum for which the type  |
|                               |             | 2 angular sums are generated as unrolled code.   |
|                               |             | Leave this at its default: 0 is a few percent    |
|                               |             | slower for s and p shells, while 2 makes the     |
|                               |             | library slower (up to 2x for d shells) and the   |
|                               |             | build ~20x longer; 3 takes hours to build.       |
+-------------------------------+-------------+--------------------------------------------------+

The maximum angular momentum is fixed at 8, for the ECPs and for the orbital basis plus the order
of derivative; the library throws ``std::invalid_argument`` beyond it.

Installation
============

``make install`` creates

.. code-block:: text

   ${CMAKE_INSTALL_PREFIX}/lib/libecpint.a
   ${CMAKE_INSTALL_PREFIX}/lib/cmake/ecpint/        # CMake package files
   ${CMAKE_INSTALL_PREFIX}/include/libecpint.hpp
   ${CMAKE_INSTALL_PREFIX}/include/libecpint/
   ${CMAKE_INSTALL_PREFIX}/share/libecpint/raw/     # the ECP library

The path to ``share/libecpint`` is what ``set_ecp_basis_from_library`` needs.

Using libecpint from CMake
^^^^^^^^^^^^^^^^^^^^^^^^^^

After installation:

.. code-block:: cmake

   find_package(ecpint CONFIG REQUIRED)
   target_link_libraries(myprogram PRIVATE ECPINT::ecpint)

adding the install prefix to ``CMAKE_PREFIX_PATH`` if it is not a standard location. To build
libecpint as part of your project instead, ``add_subdirectory`` its source directory (perhaps with
``-DLIBECPINT_BUILD_TESTS=OFF``) and link to the same ``ECPINT::ecpint`` target. The target brings
its include directories and the C++17 requirement with it.

.. _`sec:upgrading`:

Upgrading from version 1.0
==========================

The public API is unchanged, apart from additions, but the build has changed:

- **C++17** is required (previously C++11).
- ``LIBECPINT_MAX_L`` is no longer an option: angular momentum up to 8 is always supported.
  Remove it from your configuration; code that read it from ``config.hpp`` should use
  ``LIBECPINT_MAX_L`` there, which is now always 8.
- The **pugixml, Faddeeva and libcerf dependencies have been removed**, with the options
  ``LIBECPINT_USE_PUGIXML`` and ``LIBECPINT_USE_CERF`` and the macro ``HAS_PUGIXML``. Only
  ``libecpint.a`` needs linking (there is no longer a ``libFaddeeva.a``), and
  ``set_ecp_basis_from_library`` is always available.
- The **ECP library** is now ``share/libecpint/raw``, in MOLPRO format, replacing the XML files
  in ``share/libecpint/xml``; ``share_dir`` is still the ``share/libecpint`` directory.
  ``ECPBasis::addECP_from_file`` now reads MOLPRO-format files (e.g. an export from the Basis Set
  Exchange), and throws ``std::runtime_error`` if the file or atom is missing.
- The generated header ``qgen.hpp`` is internal: it is written to the build tree and not
  installed.
- ``init``'s ``smallGrid`` and ``bigGrid`` arguments are accepted but no longer have any effect.

Some public but internal classes and functions (the radial, angular and quadrature classes, for
example) have changed or been removed. These are not part of the documented API.

Documentation
=============

This documentation can be built from the build directory with

.. code-block:: bash

   pip install -r ../doc/sphinx/requirements.txt
   make docs

which writes it to ``doc/html/index.html`` in the build directory. This needs Doxygen, and Sphinx
with the packages in ``doc/sphinx/requirements.txt``. If ``doxygen`` or ``sphinx-build`` is not on
your path, pass its location to CMake with ``-DDOXYGEN_EXECUTABLE=...`` or
``-DSPHINX_EXECUTABLE=...``.

Testing
=======

To check your build, run in the build directory

.. code-block:: bash

   ctest -j8 -LE reference

This checks the library's components, its radial integrals and closed-form cases, and the
integrals and their derivatives for a set of molecules against stored values from an independent
reference implementation. If a test fails, ``ctest --output-on-failure`` (or
``Testing/Temporary/LastTest.log``) gives the reason, naming the worst element. If tests still
fail after following these instructions, please raise an issue on GitHub, giving details of your
environment, the options you used, and the relevant part of ``LastTest.log``.

:ref:`sec:developer` describes the tests in more detail, and the timing benchmark
(``make Benchmark``).

.. toctree::
   :hidden:
