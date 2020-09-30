
.. _program_listing_file__Users_robertshaw_devfiles_libecpint_new_README.md:

Program Listing for File README.md
==================================

|exhale_lsh| :ref:`Return to documentation for file <file__Users_robertshaw_devfiles_libecpint_new_README.md>` (``/Users/robertshaw/devfiles/libecpint_new/README.md``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: markdown

   ## Libecpint
   
   [![Build Status](https://dev.azure.com/robertshaw383/libecpint/_apis/build/status/robashaw.libecpint?branchName=readecp)](https://dev.azure.com/robertshaw383/libecpint/_build/latest?definitionId=2&branchName=readecp)
   [![codecov](https://codecov.io/gh/robashaw/libecpint/branch/readecp/graph/badge.svg)](https://codecov.io/gh/robashaw/libecpint)
   
   Libecpint is a C++ library for the efficient evaluation of integrals over ab initio effective core potentials, using a mixture of generated, recursive code and Gauss-Chebyshev quadrature. It is designed to be standalone and generic, but is currently in development and may not be completely stable. If you experience any problems please raise an issue here; contributions and suggestions are also welcome.
   
   ### Applications
   
   This assumes ECPs and basis sets of the form usually seen in electronic structure calculations, namely those expanded in terms of Gaussian functions. The angular momentum of function that can be treated is in theory arbitrary, but is limited by your choice of maximum when the library is built.
   
   
   ## Building and testing
   
   ### Dependencies
   
   - A modern C++ compiler, at least C++11 standard library is required. This has currently only been tested with GCC (6.3.0 and above, but will in theory work with any > 4.9) and clang (9.0.0 and above). Intel compilers have been known to cause issues.
   - CMake/CTest build tools
   - Python (2.7 or above, including 3 and higher)
   
   Additionally, if you wish to regenerate the radial code (see below),  Python3 is required with numpy and sympy.
   
   ### Build instructions
   
   To build the library, fork the repo locally and do the following in the top of the source tree (out-of-source build is required!):
   
   ```
   mkdir build
   cd build
   cmake [options] ..
   make [-jn]
   ```
   
   The `-jn` flag tells make to use `n` threads while compiling (e.g. `-j4` would use four threads), and is highly recommended if your computer can cope, as the generated code files can all be compiled independently of one another.
   
   The options after cmake above can be included using the syntax `-DOPTION=value`. The pertinent options are as follows:
   - `LIBECPINT_MAX_L` = the maximum angular momentum (in both orbital and ECP basis) required. The default is 5 (i.e. h-type functions), but this can easily be increased. Note that the higher this value, the longer the code generation will take (especially if optimization flags have not been added - see below), but it will not greatly affect compilation time.
   - `CMAKE_CXX_FLAGS` = any additional flags to be passed to the compiler. It is _strongly recommended_ that you provide optimization flags, e.g. at least `-O2` if not `-O3` for gcc/clang.
   - `LIBECPINT_MAX_UNROL` = the maximum angular momentum for which the whole integral is unrolled. The default is 2. It is _strongly recommended_ that you do not increase this past 4, as the compilation time and file sizes increase significantly. For reference, the following table gives compilation times and max. file sizes with `-O3` optimization flags and GCC 6.3.0:
   
   |  LIBECPINT_MAX_UNROL |  Compilation time (minutes)  | Max. file size (MB)  |
   |---|---|---|
   |  1  |  2  |  0.02  |
   |  2  |  10  | 1.0  |
   |  3  |  54  |  24.1  |
   |  4  |  191  |  556.2  |
   
   
   ### Testing and installation
   
   To test and install the build, do
   ```
   make test
   make install
   ```
   
   ## Documentation
   
   Code documentation can be generated using doxygen in the folder `doc/doxygen`. Examples of how to use the library can be found in the `tests` directory.
   
   The usage documentation is currently spare and under development - apologies, we hope to rectify this soon, along with providing a better API!
   
   ## Performance
   
   The angular momenta that have been fully unrolled will evaluate very rapidly, but the higher angular momenta will be noticeably slower. We are currently looking at ways to reduce this cost, but as these only make up a very small amount of the total number of integrals (which in turn are only a fraction of the total computation time in an actual calculation), it is very unlikely the ECP integrals will ever become a bottleneck.
   
   
   ## Acknowledging usage
   
   If you use this library in your program and find it helpful, that's great! Any feedback would be much appreciated. If you publish results using this library, please consider citing the following paper detailing the implementation:
   
   R. A. Shaw, J. G. Hill, J. Chem. Phys. 147, 074108 (2017); doi: [10.1063/1.4986887](http://dx.doi.org/10.1063/1.4986887)
   
   A full bibtex citation can be found in CITATION in the main directory.
   
   ## Work in progress
   
   ### Currently in the development version but not yet stable
   - First and second derivatives
   - GoogleTest unit testing suite, to work with continuous integration
   
   ### In the near future
   - An improved API
   - A built in ECP library
   
   ## Regenerating the radial code
   
   The recursive radial integral code has been pre-generated, as the current setting has been calibrated to balance accuracy and efficiency. If you would like to experiment (warning: after reading the paper cited above), go into the directory `src/generated/radial`. Edit the top line of  `unrol_radial.py` to change `MAX_UNROL_AM`, the maximum angular momentum to be unrolled. Then do the following:
   
   ```
   python3 unrol_radial.py
   ./generate.sh
   ```
   
   This will generate the simplified recursive integrals and then piece together the `radial_gen.cpp` file and place it in the correct location. It should be very safe (but not very efficient) to decrease `MAX_UNROL_AM`, but be prepared for things to break if you increase it too much. 
