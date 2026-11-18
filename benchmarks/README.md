# libecpint benchmark

`ecpint_bench` times the library, and `compare.py` compares its results between builds, e.g. before
and after a change, or with different settings such as `LIBECPINT_MAX_UNROL`. Neither runs under
ctest. The benchmark is built with the tests (`LIBECPINT_BUILD_TESTS`, the default), and uses the
test systems in `tests/data`.

## Running

From the build directory:

```bash
make Benchmark                           # kernels and systems, ~20 s
benchmarks/ecpint_bench kernels          # one group
benchmarks/ecpint_bench scaling --max-carbons 128 --wasted
benchmarks/ecpint_bench --reps 21 --csv out.csv --tag mybuild
```

Use a release build, on an otherwise idle machine.

## What it measures

- **kernels**: `ECPIntegral::compute_shell_pair` per call, for each pair of angular momenta up to
  g, against the def2 ECP of Au (s, p and d projectors, f local channel). `off` has both shells
  away from the ECP, about 3.8 bohr from it: the type 2 code that unrolling changes. `A on ECP` and
  `both on ECP` use other code paths that unrolling does not change, so they show how much the
  timings vary between builds that should be equally fast. The same pair is called repeatedly, so
  the library's per-thread radial caches are warm, and the angular part weighs more than in a
  whole calculation.
- **systems**: `ECPIntegrator` initialisation, `compute_integrals`, `compute_gradient` and
  `compute_hessian` (contracted with a test density), for test systems from
  `tests/data/systems.txt` and an iodoalkane chain C16H34 with four iodines (def2-SVP, ecp28mdf).
  Every calculation runs on a fresh thread, so the per-thread caches start empty, as for the first
  calculation at a geometry. The setup is not included in the times, except in `init`.
- **scaling**: `compute_integrals` and `compute_gradient` for the iodoalkane chains C8, C16, ...,
  with local exponents of the growth with the number of atoms. `--wasted` also counts the
  shell-pair calls of `compute_integrals` that pass the screening but return only values below the
  tolerance: the slack in the screening.

## Reading the results

Each measurement is run once and discarded, then repeated (`--reps`; 11 for kernels, 5 for systems,
3 for scaling). A system-level sample is the mean of as many calculations as take at least 50 ms.
The report gives:

- `median` and `min` of the repetitions. Compare medians. The minimum is the least disturbed run,
  but is more sensitive to luck.
- `spread`: interquartile range / median. A difference smaller than about twice the spread is not
  significant.
- `checksum`: a weighted sum of every value computed. Builds that compute the same integrals agree
  to ~1e-15 relative; a larger difference means a change in the values, not just the speed.

`--csv` writes one row per measurement (tag, group, name, unit, reps, median, min, spread, checksum).

## Comparing builds

`compare.py` compares results written with `--csv`, grouped by `--tag`. Several files with the same
tag are repeated runs of one build. Alternate the builds' runs rather than running each build in
turn, so that drift (e.g. thermal throttling) affects all of them alike:

```bash
for round in 1 2 3; do
  for b in old new; do
    build-$b/benchmarks/ecpint_bench --csv $b-$round.csv --tag $b
  done
done
benchmarks/compare.py old-*.csv new-*.csv --base old
```

For each measurement, it prints the base build's median over its runs and the ratio of each other
build to it (below 1 is faster). Each ratio comes with the noise of the two builds: the larger of
the spread of the per-run medians and the typical spread within a run. Geometric means summarise
the `off` kernels, the controls (which should be ~1 unless the change affects them), the systems
and the scaling. A `!` marks a checksum that differs from the base's. `--out` writes the
comparison as CSV.
