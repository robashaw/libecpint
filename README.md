# Libecpint 1.1.0

[![Build Status](https://dev.azure.com/robertshaw383/libecpint/_apis/build/status/robashaw.libecpint?branchName=main)](https://dev.azure.com/robertshaw383/libecpint/_build/latest?definitionId=2&branchName=main)
[![codecov](https://codecov.io/gh/robashaw/libecpint/branch/main/graph/badge.svg)](https://codecov.io/gh/robashaw/libecpint)
[![Documentation Status](https://readthedocs.org/projects/libecpint/badge/?version=latest)](https://libecpint.readthedocs.io/en/latest/index.html)
[![Code Quality](https://www.code-inspector.com/project/15206/status/svg)]()

[![Zenodo](https://zenodo.org/badge/DOI/10.5281/zenodo.4694353.svg)](https://doi.org/10.5281/zenodo.4694353)
[![JOSS](https://joss.theoj.org/papers/10.21105/joss.03039/status.svg)](https://doi.org/10.21105/joss.03039)

Libecpint is a C++ library for the efficient evaluation of integrals over ab initio effective core potentials, using a mixture of generated, recursive code and Gauss-Chebyshev quadrature. It is designed to be standalone and generic, and is now in its first stable release. If you experience any problems please raise an issue here; contributions and suggestions are also welcome.

## Python bindings

**pylibecpint** provides Python bindings for this library, installable via pip:

```bash
pip install pylibecpint
```

See the [pylibecpint repository](https://github.com/robashaw/pylibecpint) for documentation and usage examples.

## Contributing

Contributions are welcomed, either in the form of raising issues or pull requests on this repo. Please take a look at the Code of Conduct before interacting, which includes instructions for reporting any violations.

## New in v1.1

- **Bug fix**: type-1 ECP integrals no longer collapse to zero for steep ECPs
- **Performance**: shell screening in `compute_integrals` reduces scaling from roughly cubic to quadratic for larger molecules
- **Performance**: Bessel function evaluation optimised with scratch arrays and improved memory access (~5% speedup)
- **Performance**: radial quadrature allocation and indirection cleanups
- **Performance**: contraction weight hoisting out of inner loops in qgen (~3% integrals, ~7% derivatives)
- **API**: screening tolerances exposed in `ECPIntegralFactory` and propagated to 1st/2nd derivative integrals
- **Build**: restored CI (Azure image + GoogleTest), fixed GoogleTest lib64 path, fixed clang build on Azure
- **Build**: C++ standard migrated to targets; minimal Windows install and macOS arm64 cross-compile support
- **Build**: removed the unused Faddeeva/libcerf dependency and the `LIBECPINT_USE_CERF` option
- **Build/API**: removed the pugixml dependency; the ECP library is read natively from the MOLPRO-format files in `share/libecpint/raw`, so `LIBECPINT_USE_PUGIXML`, `HAS_PUGIXML`, the XML files and `parseecp.py` are gone. `set_ecp_basis_from_library` is always available, and `addECP_from_file` takes a MOLPRO-format file (e.g. a Basis Set Exchange MOLPRO export) and throws `std::runtime_error` if the file or atom is missing
- **Build**: example built as part of tests
- **Code style**: codebase reformatted with clang-format (Google style)
- **Correctness**: Bessel small-z guard, FAC guard, FAST_POW assert, VLA removal
- **Build/API**: `LIBECPINT_MAX_L` is no longer a build option. Angular momentum is supported up to 8, for the orbital basis plus the order of derivative (second derivatives up to i functions) and for the ECPs, and the library throws `std::invalid_argument` beyond it (release builds used to return wrong derivatives silently). Only the small unrolled type 2 kernels are generated, and the generated `qgen.hpp` is written to the build tree rather than the source tree
- **Performance**: the angular integral tables store only the entries that are used, cutting the memory taken at initialisation by about 50x (e.g. 550 MB to 10 MB for second derivatives with g functions)

## Dependencies

- A C++17 compiler, e.g. gcc (v8 and above) or clang (v7 and above)
- CMake/CTest build tools (v3.12 and higher)
- Python (2.7 or above, including 3 and higher)

Additionally, if you wish to regenerate the radial code (see below),  Python >=3.6 is required with numpy and sympy.

## Documentation

Please refer to the main documentation [here](https://libecpint.readthedocs.io/en/latest/index.html).

## Examples

There is also a working example in the example folder, with instructions of how to build and link against the library.

## Testing

After building, run the tests from the build directory with

```bash
ctest -j8 -LE reference
```

This checks the library's components, its radial integrals and closed-form cases, and the integrals and their first and second derivatives for a set of test systems, against values from an independent reference engine stored in `tests/data/reference`. Configuring with `-DLIBECPINT_TEST_REFERENCE=ON` adds slower tests against the reference engine computed at test time (`ctest -L reference`). After changing a test system or the reference engine, regenerate the stored values with `cmake --build . --target regenerate_references`.

See [tests/README.md](tests/README.md) for how the tests are organised, how to diagnose a failure, and how to add tests.

## Acknowledging usage

If you use this library in your program and find it helpful, that's great! Any feedback would be much appreciated. If you publish results using this library, please consider citing the following papers detailing the implementation and the library, respectively:

R. A. Shaw, J. G. Hill, J. Chem. Phys. 147, 074108 (2017); doi: [10.1063/1.4986887](http://dx.doi.org/10.1063/1.4986887)

R. A. Shaw, J. G. Hill, J. Open Source Softw. 6, 3039 (2021); doi: [10.21105/joss.03039](https://doi.org/10.21105/joss.03039)

Full bibtex citations can be found in CITATION in the main directory.
