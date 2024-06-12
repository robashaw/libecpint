# Libecpint 1.0.7

[![Build Status](https://dev.azure.com/robertshaw383/libecpint/_apis/build/status/robashaw.libecpint?branchName=master)](https://dev.azure.com/robertshaw383/libecpint/_build/latest?definitionId=2&branchName=master)
[![codecov](https://codecov.io/gh/robashaw/libecpint/branch/master/graph/badge.svg)](https://codecov.io/gh/robashaw/libecpint)
[![Documentation Status](https://readthedocs.org/projects/libecpint/badge/?version=latest)](https://libecpint.readthedocs.io/en/latest/index.html)
[![Code Quality](https://www.code-inspector.com/project/15206/status/svg)]()

[![Zenodo](https://zenodo.org/badge/DOI/10.5281/zenodo.4694353.svg)](https://doi.org/10.5281/zenodo.4694353)
[![JOSS](https://joss.theoj.org/papers/10.21105/joss.03039/status.svg)](https://doi.org/10.21105/joss.03039)

Libecpint is a C++ library for the efficient evaluation of integrals over ab initio effective core potentials, using a mixture of generated, recursive code and Gauss-Chebyshev quadrature. It is designed to be standalone and generic, and is now in its first stable release. If you experience any problems please raise an issue here; contributions and suggestions are also welcome.

## Contributing

Contributions are welcomed, either in the form of raising issues or pull requests on this repo. Please take a look at the Code of Conduct before interacting, which includes instructions for reporting any violations.

## New in first full release

- Analytical 1st and 2nd derivatives;
- Integration now >10x faster;
- New, high level API, with ECP library;
- Automated testing suite.

### Patch 1

- Bug fix in screening of on-ECP type 2 integrals
- Improvements in CMake build steps, thanks to nabbelbabbel/moritzBens

### Patch 2

- Fix for memory leaks in derivative routines
- Minor changes to CMake files

### Patch 3

- Fix bug in radial type 1 integrals where quadrature could fail to converge
- Const correctness throughout, should allow for parallelisation
- Minor updates to docs

### Patch 4

- Code generation now takes considerably less time and memory; MAX_L=8 takes ~35 seconds, peaking at 1.5GB of memory (joint effort with Thomas Dresselhaus and Peter Bygrave)
- This will be the final patch before v1.1

## Dependencies

- A modern C++ compiler, at least C++11 standard library is required. This has been tested with:
  * gcc (v6.3.0 and above)
  * clang (v10.0.0 and above), you may need the CXX flag "-std=c++14"
  * icpc (v20.2.1), may also need the CXX flag "-std=c++14"
- CMake/CTest build tools (v3.12 and higher)
- Python (2.7 or above, including 3 and higher)

Additionally, if you wish to regenerate the radial code (see below),  Python >=3.6 is required with numpy and sympy.

## Documentation

Please refer to the main documentation [here](https://libecpint.readthedocs.io/en/latest/index.html).

## Examples

There is also a working example in the example folder, with instructions of how to build and link against the library. Please also the API tests in tests/lib/

## Acknowledging usage

If you use this library in your program and find it helpful, that's great! Any feedback would be much appreciated. If you publish results using this library, please consider citing the following papers detailing the implementation and the library, respectively:

R. A. Shaw, J. G. Hill, J. Chem. Phys. 147, 074108 (2017); doi: [10.1063/1.4986887](http://dx.doi.org/10.1063/1.4986887)

R. A. Shaw, J. G. Hill, J. Open Source Softw. 6, 3039 (2021); doi: [10.21105/joss.03039](https://doi.org/10.21105/joss.03039)

Full bibtex citations can be found in CITATION in the main directory.
