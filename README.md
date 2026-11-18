# Libecpint 1.1.0

[![Build Status](https://dev.azure.com/robertshaw383/libecpint/_apis/build/status/robashaw.libecpint?branchName=main)](https://dev.azure.com/robertshaw383/libecpint/_build/latest?definitionId=2&branchName=main)
[![codecov](https://codecov.io/gh/robashaw/libecpint/branch/main/graph/badge.svg)](https://codecov.io/gh/robashaw/libecpint)
[![Documentation Status](https://readthedocs.org/projects/libecpint/badge/?version=latest)](https://libecpint.readthedocs.io/en/latest/index.html)
[![Zenodo](https://zenodo.org/badge/DOI/10.5281/zenodo.4694352.svg)](https://doi.org/10.5281/zenodo.4694352)
[![JOSS](https://joss.theoj.org/papers/10.21105/joss.03039/status.svg)](https://doi.org/10.21105/joss.03039)

Libecpint is a C++ library for the efficient evaluation of integrals over ab initio effective core potentials (ECPs), and their first and second derivatives, using a mixture of generated, recursive code and Gauss-Legendre quadrature. It is designed to be standalone and generic. If you experience any problems please raise an issue here; contributions and suggestions are also welcome.

## Python bindings: pylibecpint

**[pylibecpint](https://github.com/robashaw/pylibecpint)** provides Python bindings for this library, with NumPy inputs and outputs:

```bash
pip install pylibecpint
```

See the [pylibecpint repository](https://github.com/robashaw/pylibecpint) for documentation and examples.

## What's new in v1.1

Version 1.1 is a substantial rework of the internals since v1.0.7. The public API (`ECPIntegrator` and `ECPIntegral`) is unchanged, apart from additions.

**Performance.** Typically 3-15x faster than v1.0.7 for whole molecules (7.5x on average over the benchmark systems), with single shell pairs about 15x faster. The radial integrals share one quadrature grid per ECP across shell pairs, with per-thread caches; work is screened at every level, from whole shells down to primitives, now including derivatives; and the angular integral tables are about 50x smaller. Gradients now scale almost linearly with molecular size: for a 386-atom iodoalkane, the gradient takes 2.4 s against 107 s with v1.0.7 (one core, Apple M3 Pro).

**Accuracy and stability.** The radial integrals now use fixed, non-adaptive quadrature throughout, in place of adaptive quadrature and closed forms that could converge spuriously or lose precision for diffuse functions, so integrals are converged to about 1e-13 and vary smoothly with the geometry. Several numerical failures of v1.0.7 are fixed, including type 1 integrals collapsing to zero for steep ECPs and inaccurate off-centre type 2 integrals for high angular momenta and diffuse functions: even for ordinary basis sets, v1.0.7's second derivatives could be off by up to ~4e-4 relative to their largest element, and by far more for difficult cases. Every integral and derivative is now checked against an independent reference implementation, over a range of molecules and ECPs, with tests that isolate the screening, the quadrature and the bookkeeping.

**New features.**

- `compute_gradient` and `compute_hessian` contract the derivatives with a density matrix directly, without storing the derivative matrices (which could take many GB).
- Derivatives with respect to a magnetic field, for gauge-including atomic orbitals (`compute_bfield_derivs`).
- Angular momentum up to 8 (l functions) for the ECPs, and for the basis including the order of derivative (so integrals up to l, first derivatives up to k, second derivatives up to i), always: it is no longer a build option.
- The built-in ECP library has grown from 6 to 91 families (Stuttgart/Cologne, Karlsruhe, LANL, CRENBL/CRENBS, SBKJC and more), read from MOLPRO-format files, which can also be exported directly from the Basis Set Exchange.

**Build changes that may affect downstream users.**

- C++17 is now required.
- The pugixml, Faddeeva and libcerf dependencies have been removed, along with the `LIBECPINT_USE_PUGIXML`, `LIBECPINT_USE_CERF` and `LIBECPINT_MAX_L` options. Only `libecpint.a` needs linking.
- The ECP library is in `share/libecpint/raw` (MOLPRO format) instead of `share/libecpint/xml`, and `addECP_from_file` reads MOLPRO-format files.
- Some public but internal classes and function signatures (e.g. the radial, angular and quadrature classes) have changed or been removed; these are not part of the documented API.

See [Upgrading from version 1.0](https://libecpint.readthedocs.io/en/latest/installation.html#upgrading-from-version-1-0) in the documentation for details.

## Installation

You need a C++17 compiler (e.g. GCC 8+ or clang 7+), CMake 3.12+ and Python (any version). Then

```bash
mkdir build && cd build
cmake ..
make -j8
ctest -j8 -LE reference   # optional, ~10 s
make install
```

and in your own CMake project, `find_package(ecpint CONFIG REQUIRED)` and link to `ECPINT::ecpint`. The documentation lists the build options.

## Documentation

The full documentation, including usage of the API, is on [Read the Docs](https://libecpint.readthedocs.io/en/latest/index.html). There is also a working example in the `example` folder, with instructions of how to build and link against the library.

## Testing

After building, run the tests from the build directory with

```bash
ctest -j8 -LE reference
```

This checks the library's components, its radial integrals and closed-form cases, and the integrals and their first, second and magnetic field derivatives for a set of molecules, against values from an independent reference engine stored in `tests/data/reference`. See [tests/README.md](tests/README.md) for how the tests are organised, how to diagnose a failure, and how to add tests, and [benchmarks/README.md](benchmarks/README.md) for the timing benchmark.

## Contributing

Contributions are welcomed, either in the form of raising issues or pull requests on this repo. Please take a look at the Code of Conduct before interacting, which includes instructions for reporting any violations, and at [ARCHITECTURE.md](ARCHITECTURE.md) and the developer guide in the documentation for how the code is organised.

## Acknowledging usage

If you use this library in your program and find it helpful, that's great! Any feedback would be much appreciated. If you publish results using this library, please consider citing the following papers detailing the implementation and the library, respectively:

R. A. Shaw, J. G. Hill, J. Chem. Phys. 147, 074108 (2017); doi: [10.1063/1.4986887](https://doi.org/10.1063/1.4986887)

R. A. Shaw, J. G. Hill, J. Open Source Softw. 6, 3039 (2021); doi: [10.21105/joss.03039](https://doi.org/10.21105/joss.03039)

Full bibtex citations can be found in CITATION in the main directory.
