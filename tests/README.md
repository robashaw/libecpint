# libecpint tests

The tests are built in layers, so that a failure points at one part of the library: the
bookkeeping (packing and assembly of derivatives), the screening, the radial quadrature, or the
integrals as a whole. Every test uses GoogleTest and runs under `ctest`.

## Running

From the build directory:

```bash
ctest -j8 -LE reference     # everything except the live reference tests: ~10 s
ctest -j8                   # everything (live tests only if configured, see below)
ctest -R Suite.Stored       # one group, by regular expression
ctest -R Radial --output-on-failure
```

The slower tests compare the library with the reference engine as it runs. They exist only when
the build is configured with `-DLIBECPINT_TEST_REFERENCE=ON`, and carry the label `reference`:

```bash
cmake -DLIBECPINT_TEST_REFERENCE=ON .
ctest -j8 -L reference      # Suite.Live.*: ~20 s each
```

Run them when the reference engine changes, or to find where a stored comparison fails. Day to
day, `ctest -LE reference` is enough: the stored comparisons check the same systems in seconds.

A single test of a GoogleTest binary can also be run directly, e.g.
`tests/suite/suite_test --gtest_filter=Systems/Stored.*/la2o3` or
`tests/radial/radial_test --gtest_filter=Type2.TightPrimitives`. Most accuracy tests print
`[ info ]` lines with their worst errors, so that a drift shows before a test fails.

## Layout

| Directory | ctest names | What it checks |
|---|---|---|
| `units/` | `MathUtil`, `Bessel`, `GaussQuad`, `Assembly`, `DerivativeRecurrences`, `Screening`, `SharedGrid`, `Api`, `BField`, `TestSupport`, ... | single components, mostly against exact or independent values |
| `radial/` | `Radial.Type2.*`, `Radial.Smoothness.*` | type 2 radial integrals against independent quadrature, case by case |
| `analytic/` | `Analytic` | integrals with all centres coincident, against closed forms, for both the library and the reference engine |
| `suite/` | `Suite.<check>.<system>` | whole systems: integrals, first and second derivatives, magnetic field derivatives |
| `support/` | | the test support library `ecpint_testsupport` (below) |
| `data/` | | the systems, basis sets, ECPs and stored reference values |

The timing benchmark is in `../benchmarks`, outside ctest (`make Benchmark`; see its README).

## Design

### The reference engine

`support/integral_reference.hpp` (`intref::Engine`) computes ECP integrals over contracted
Cartesian shells, and their first and second derivatives with respect to the shell centres, by a
different method from the library's. With `bfield`, it also computes the derivatives with respect
to a magnetic field over gauge-including orbitals, by multiplying each basis function pointwise by
its factor `[(A - C) x r]_k`, where the library raises the shells by one. It integrates the semi-local projections numerically on the
sphere, with grids aligned to each Gaussian, and the radial parts on composite Gauss-Legendre
panels. It shares no code with the library, expands nothing about the ECP centre (the source of
the library's main precision limit, below), and screens nothing beyond terms below ~1e-40. It is
slow, so its results are stored.

It computes only the blocks for each (ECP, shell pair): the integrals and the derivatives with
respect to the two shell centres. Derivatives with respect to the ECP centre (translational
invariance) and the placement of blocks in molecular matrices are left to the library's own code
(`assembly.hpp`: `pack_first_derivatives`, `pack_second_derivatives` and the assembly loops), which
has its own unit tests (`units/assembly`, `units/derivatives`). So a difference between the
library and the reference in the suite is a difference in the integrals, not in bookkeeping.

`support/radial_reference.hpp` gives independent radial integrals (`radref::G`) and scaled Bessel
functions (`radref::K`, from their power series), for `radial/` and `units/`.

### Systems and stored references

`data/systems.txt` lists the systems, one per line:

```
name=la2o3  xyz=systems/la2o3.xyz  basis=basis/la2o3-regression.gbs  ecp=La:lib/def2-ecp  deriv=2 live=yes
```

- `xyz`: geometry (angstrom, or bohr if the comment line says so)
- `basis`: a Gaussian94 basis file, trimmed to the elements used; blocks keyed by element or by
  atom number (1-based, for that atom only)
- `ecp`: `ELEMENT:lib/NAME` for `share/libecpint/raw/NAME.ecp`, or `ELEMENT:data/NAME` for
  `data/ecp/NAME.ecp` (Molpro format)
- `deriv`: the highest derivative order checked
- `live`: whether to add a `Suite.Live` test (ag3_lanl2dz and auh2_highl are `no`: their
  references take ~145 s and ~100 s)

`data/reference/<system>.ref` holds the reference engine's values for each system: the integrals
in full, and the first, second and magnetic field derivatives (the last computed for every system
with `deriv` >= 1) contracted with two fixed test densities (dense
derivatives would make the files too large; the live tests compare them in full). Each value has
a conditioning scale (below). The header records a hash of the system's inputs; if the geometry,
basis or ECP changes, the stored tests fail with "references stale" rather than with a misleading
difference.

### The suite checks

For each system:

- **Stored**: the library, with default settings, through `compute_integrals`,
  `compute_gradient`, `compute_hessian` and `compute_bfield_derivs`, against the stored values.
- **Live** (label `reference`): the library against the reference engine computed now, every
  element of the integrals and the dense derivatives, with the worst element named by atom, shell
  and Cartesian component.
- **Screening**: the library with its default screening against the library with none. A
  difference is a screening error, whatever the accuracy of the integrals.
- **Threads**: every (ECP, shell pair) block of every per-block kernel (integrals, first and
  second derivatives up to the system's `deriv`, magnetic field derivatives) computed by several
  threads sharing one `ECPIntegral` is identical to the serial result; a failure names the kernel.

### Tolerances

Each value must satisfy

    |value - reference| <= 1e-13 * max |reference| + 1e-14 * scale

with the maximum over the set (integrals, first, second or magnetic field derivatives). The first term is the
expected agreement between two accurate methods (`DEFAULT_REL_TOL` in `support/compare.hpp`).

The second allows for the library's precision limit for type 2 integrals: it expands each shell
about the ECP centre, and for a primitive with exponent a at distance d from it the terms of that
expansion exceed the result by up to ~(1 + d sqrt(2a))^l (l including the angular momentum that
derivatives add), so small relative errors in the radial integrals are amplified by that much.
The scale of a value estimates the size of those terms: for each block contributing to it, the
largest element of the integral over the shells raised to l + o (for derivatives of order o) with
coefficients |c| (2a)^o (1 + d sqrt(2a))^(l + o). The magnetic field derivatives are built from
the shells raised by one with no factor 2a, weighted by the distances in the formula (see
`reference_integrals` in `support/reference.hpp`). `CONDITIONING_TOL` (1e-14) is calibrated on
every element of the live systems; see its comment, and "Known limitations" in `ARCHITECTURE.md`.

The unit, radial and analytic tests set their own tolerances, each with the measured error in a
comment, typically ~50x below the tolerance.

## Diagnosing a failure

| Fails | Passes | Likely cause |
|---|---|---|
| `Assembly` | | packing of derivatives or the assembly loops |
| `DerivativeRecurrences` | `Assembly` | the derivative recurrences of `ECPIntegral` |
| `BField`, or only the magnetic field derivatives in `Suite.*` | `Assembly`, the other derivatives | `compute_shell_pair_derivative_Bfield` (its special cases are named by the `BField` tests) |
| `Screening` (unit) or `Suite.Screening` | radial, analytic | a screening bound or estimate (`screening.hpp`) |
| `Radial.*`, `SharedGrid`, `Bessel`, `GaussQuad` | | radial quadrature, the shared grid, Bessel functions |
| `Analytic` for the library only | | on-centre integrals (type 1, or the both-on-ECP closed form) |
| `Analytic` for the reference engine | | the reference engine itself: do not trust stored values |
| `Suite.Live` | all of the above | the angular part or assembly of the integrals: a true integration bug |
| `Suite.Stored` only | `Suite.Live` | stale reference values: regenerate (below) |
| `Suite.Threads` | | a data race (run under ThreadSanitizer) |

A failure report names the worst element, its error and its tolerance, e.g.

    second derivatives: 20 values beyond tolerance (1.65e-12 = rel * max |reference| 1.651e+01,
    plus any conditioning allowance); worst error 6.014e-12 (tolerance 3.02e-12) at
    d2/dO3_z dO3_z V [O3 shell 4 s, O4 shell 5 d(zz)]

To narrow a failure down, compare the same element with screening off (`System::setup` with
tolerances `-1.0, 0.0`), per ECP and per primitive, and against the reference engine with tighter
grids. Check the reference itself on a case with a known answer before blaming the library: its
convergence matters as much as the library's.

## Regenerating the stored references

```bash
cmake --build <build> --target regenerate_references
```

This runs the reference engine on all hardware threads for every system in the manifest and
writes `data/reference/*.ref` into the source tree (~3 minutes, most of it ag3_lanl2dz and auh2_highl).
`<build>/tests/suite/regenerate_refs NAME...` regenerates only the named systems. Regenerate
after changing a system's inputs or the reference engine, never to make a library change pass.

## Adding tests

**A system.** Add the geometry to `data/systems/`, a Gaussian94 basis trimmed to its elements to
`data/basis/` (e.g. from the Basis Set Exchange), and any ECP not in the library to `data/ecp/`.
Add a line to `data/systems.txt`, keeping it small enough for the reference engine (a few hundred
basis functions at most), re-run CMake, then regenerate its reference with
`regenerate_refs NAME`. Commit the `.ref` file with the system.

**A radial regression case.** Add a `TEST(Type2, ...)` to `radial/radial.cpp` with a comment
recording the failure it reproduces, comparing against `radref::G`. Measure the error first, and
set the tolerance from it. Check that the test fails with the bug reintroduced.

**A screening bound.** Document its formula in `screening.hpp`, and add a property test to
`units/screening/screening.cpp` comparing it with the exact quantity it bounds, over a sweep of
parameters and over the manifest's systems.

**A unit test.** Add a directory under `units/` with a `CMakeLists.txt` like the others
(link `gtest` and `ecpint`, or `ecpint_testsupport` to use the support library), and add it to
`units/CMakeLists.txt`.
