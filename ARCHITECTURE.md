# Architecture

This document describes the high-level architecture of libecpint, to help any new contributors find their way around the codebase.

## High-level overview

The code is roughly divided into four chunks:

- the interface,
- the fixed backend,
- the generated backend,
- the code generator

As a general rule, the API will not change. If functionality is added to the library, it should be exposed through ```api.hpp```, following the style of other functions in that interface. Changing any existing parts of the API will break compatibility, so we will only do this if it is *absolutely* necessary.

The other three parts are discussed below.

## Code Map

### Fixed backend

This is found in ```src/lib``` and can be further divided into the following functionalities:

- Quadrature: the fixed Gauss-Legendre rules (```gaussquad.cpp```), the type 1 and type 2 radial integrals (```radial.cpp```), and the radial grid shared between shell pairs, with its per-thread caches (```shared_grid.cpp```)
- Screening: every bound and estimate used to skip negligible work, from whole shells down to primitive triples (```screening.cpp```)
- Assembly: packing per-ECP shell-pair derivatives into derivatives with respect to the centres, and the loops that place shell-pair blocks in the molecular matrices (```assembly.cpp```)
- Analytic angular integrals (```angular.cpp```)
- Containers (```ecp.cpp, gshell.cpp```)
- Mathematical utilities (```bessel.cpp, mathutil.cpp```)

In brief, this is where you should put changes of any of the above kinds. These are all things that are used by the library as part of core functionality, and that do not change depending on how the generated part of the code is initialised.

Each file generally describes one class (for example, the ```ECP``` object in ```ecp.cpp``` or the ```AngularIntegral``` object in ```angular.cpp```) _OR_ one set of functionalities (e.g. the screening bounds, in ```screening.cpp```). Please stick to this convention.

As a general rule, any changes in the fixed backend will only affect other parts of the fixed backend, and will not affect that API or generated code. The exception is ```ecpint.cpp``` which describes an ```ECPIntegral``` object, but connects the fixed and generated code together. The screening it applies is defined in ```screening.cpp```.

### Generated backend

This is found in ```src/generated``` and ```src/generated/radial```. It is where any generated code, or templates required for code generation, are located. Currently there are two "part" files, templates for the generated header ```qgen.hpp``` and the kernel lookup ```ecpint_gen.cpp```, which are written to the build tree (```<build>/src/generated```) with the generated kernels.

The radial folder holds the generator for closed forms of the primitive radial integrals, unrolled from their recurrence relations. It is not used at present (the radial integrals are calculated by quadrature); see the README.md there.

### Code generator

The generated code is handled by ```src/generate.cpp``` and the utility functions in ```include/generate.hpp```. The former works by finding all the non-zero integral terms for a particular radial integral class using the ```SumTerm``` objects in the header. These are then sorted and made unique, before the unrolled code is written into the generated backend described above. Only the small cases are generated, Q(LA, LB, lam) with LA <= LB <= ```LIBECPINT_MAX_UNROL``` and LA + LB + lam <= 3 ```LIBECPINT_MAX_UNROL```, looked up by ```qgen::unrolled```. All other type 2 integrals are computed by the general ```qgen::rolled_up```, from the radial integrals that the angular integrals select, which ```ECPIntegral``` finds once when it is constructed (```buildType2Radials```, the same selection as the generator's). So the maximum angular momentum (```LIBECPINT_MAX_L```, fixed at 8) does not depend on code generation; it is set by the power and Gamma function tables in ```mathutil.hpp```, which ```ecpint.cpp``` checks with ```static_assert```.

This part is compiled and run as a separate project before the main library. Therefore if any changes you make depend on the generated code, and in particular ```qgen::unrolled```, you should think very carefully about whether it will be affected by any changes in the code generation. If so, functionality should be added to the generated backend, in the same way that ```ecpint_gen.cpp``` is handled. 

### Tests and benchmarks

The tests are in ```tests```: unit tests of single components (```tests/units```), the type 2 radial integrals (```tests/radial```), closed-form cases (```tests/analytic```), and whole test systems against stored values from an independent reference engine (```tests/suite```, with the support code in ```tests/support``` and the data in ```tests/data```). ```tests/README.md``` describes the design, how to diagnose failures, and how to add tests. A change to the fixed backend should keep them passing, and comes with a test that fails without it. The timing benchmark, outside the tests, is in ```benchmarks``` (```make Benchmark```; see ```benchmarks/README.md```).

## Known limitations

- **Type 2 precision far from the ECP.** The type 2 integrals expand each shell about the ECP
  centre (the binomial coefficients `CA`, `CB` used by `qgen::rolled_up` and the unrolled kernels). For a primitive
  with exponent a at distance d from the ECP, the terms of that expansion exceed the result by up
  to ~(1 + d sqrt(2a))^l, l including the angular momentum added by derivatives, so the relative
  error of a block grows as ~1e-16 (1 + d sqrt(2a))^(2l): about 1e-15 for p shells, but 6e-8 for
  an h shell (exponent 8.7) 2 bohr away, and more at the limit of l = 8 (i shells with second
  derivatives). The affected blocks are small in absolute terms; the
  test suite allows for this (`CONDITIONING_TOL` in `tests/support/compare.hpp`). Removing the
  limit needs a formulation without the expansion, e.g. writing Cartesian Gaussians as
  derivatives of s-type ones with respect to their centre. To be revisited.
