.. libecpint developer guide

.. _`sec:developer`:

===============
Developer guide
===============

This page is for those working on libecpint itself: how it computes the integrals, how the work
is screened, and how to test and benchmark changes. ``ARCHITECTURE.md`` in the repository maps
the source files; ``tests/README.md`` and ``benchmarks/README.md`` give more detail on testing
and benchmarking.

Methodology
===========

The method is that of [Shaw2017]_, with the radial integrals now computed entirely by
quadrature. In outline:

The integrals
^^^^^^^^^^^^^

A semi-local ECP centred on C is

.. math::

   U(\mathbf{r}) = U_L(r) + \sum_{l=0}^{L-1} \sum_{m=-l}^{l} |lm\rangle U_l(r) \langle lm|,
   \qquad U_l(r) = \sum_k d_k r^{n_k - 2} e^{-\eta_k r^2},

giving **type 1** integrals over the local part :math:`U_L` and **type 2** integrals over the
projectors. Each Cartesian shell is expanded about C (binomial coefficients in the components of
:math:`\mathbf{A} - \mathbf{C}`), so that each integral is a sum of products of

- **angular integrals** over real spherical harmonics, which depend only on angular momenta.
  ``ECPIntegral`` tabulates them once, when it is constructed, storing only the entries that are
  used; and
- **radial integrals** of the form

  .. math::

     G(N, l_1, l_2) = \int_0^\infty r^{N} U_l(r)\, e^{-a(r-A)^2} K_{l_1}(2aAr)\,
     e^{-b(r-B)^2} K_{l_2}(2bBr)\, dr

  with :math:`K_l(z) = e^{-z} i_l(z)` the scaled modified spherical Bessel functions, and A, B
  the distances of the shells from C.

Shells on the ECP centre take simpler paths (only :math:`K_0` survives), and when both are on it
the type 2 integrals have a closed form.

The angular sums for type 2 integrals are done by general code (``qgen::rolled_up``), except for
the smallest classes, which are unrolled by a code generator run at build time
(``LIBECPINT_MAX_UNROL``, default 1; higher values make the library *slower*, see
:ref:`sec:installation`).

Radial quadrature
^^^^^^^^^^^^^^^^^

All radial integrals are computed by fixed (non-adaptive) composite Gauss-Legendre rules: 48
points per panel, with panels no wider than :math:`11/\sqrt{\zeta}` for a Gaussian of exponent
:math:`\zeta`. The rules are converged to ~1e-14 relative for any power of r up to 40, and, being
fixed, give integrals that are smooth in the geometry (adaptive quadrature, as in earlier
versions, could stop on a spurious convergence, and the closed forms they used for low angular
momenta lost precision for diffuse primitives).

For the type 2 integrals, the work is shared between shell pairs:

- each ECP has one **shared radial grid**, on which every primitive with exponent up to 30 is
  integrated;
- on it, the contracted factors :math:`F_{A,l}(r) = \sum_a c_a e^{-a(r-A)^2} K_l(2aAr)` depend
  only on the shell and its distance from the ECP, so they are computed once and reused for every
  partner shell and every ECP channel, instead of evaluating Bessel functions for every
  primitive triple;
- tighter primitives (core-like functions) are integrated per primitive triple, each on its own
  window of r outside which the integrand is negligible.

Grids and factors are kept in small per-thread least-recently-used caches, keyed by content
(``LIBECPINT_RADIAL_CACHE_MB`` per thread, default 4), so the results depend neither on the
cache state nor on the number of threads.

The Bessel functions are tabulated (with a Taylor series) for arguments below 16, and computed by
recurrence above, upwards while stable and downwards (Miller's method) otherwise, to ~1e-16
(absolute) for every order the library uses.

Derivatives
^^^^^^^^^^^

Derivatives with respect to the shell centres are integrals over shells of raised and lowered
angular momentum,

.. math::

   \partial_{A_x} \phi_{l_x} = l_x \phi_{l_x - 1} - 2a\, \phi_{l_x + 1},

applied once or twice; derivatives with respect to the ECP centre follow from translational
invariance. ``assembly.hpp`` packs the per-ECP shell-pair derivatives into derivatives with
respect to atoms (including the cases where centres coincide), and holds the loops that place
blocks into molecular matrices, or contract them directly with a density
(``compute_gradient``, ``compute_hessian``). The magnetic field derivatives are integrals over
shells raised by one, weighted by the components of :math:`(\mathbf{A} - \mathbf{C}) \times
\mathbf{r}`.

Screening
^^^^^^^^^

Work is skipped at four levels, from the coarsest to the finest. All the bounds are in
``screening.hpp``, each documented with its formula and derivation, and kept apart from the
integration code so that each can be tested on its own (a bound must bound).

+----------------------+-----------------------------+------------------------------------------+
| Level                | Quantity                    | Skipped when                             |
+======================+=============================+==========================================+
| shell against ECP    | ``shell_ecp_bound`` times   | below ``tol``: every pair it forms would |
|                      | the largest partner's       | be skipped too                           |
+----------------------+-----------------------------+------------------------------------------+
| shell pair against   | product of the two shells'  | below ``tol`` (default 1e-12)            |
| ECP                  | ``shell_ecp_bound``         |                                          |
+----------------------+-----------------------------+------------------------------------------+
| ECP channel, within  | ``channel_bounds``          | below 1e-2 x ``radial_tol``              |
| a shell pair         | (rigorous bound)            | (default 1e-15)                          |
+----------------------+-----------------------------+------------------------------------------+
| primitive triple, in | ``type2_bound``             | below its share of ``radial_tol``        |
| a radial integral    | (rigorous bound)            |                                          |
+----------------------+-----------------------------+------------------------------------------+

**Shells and shell pairs.** For a shell of angular momentum l, smallest exponent
:math:`\alpha`, coefficients :math:`c_i` and exponents :math:`a_i`, at distance A from an ECP
with most diffuse exponent :math:`\eta` and coefficients :math:`d_k`,

.. math::

   \mathrm{shell\_ecp\_bound} = \mathrm{shell\_bound}(l + o, \alpha, A^2, \eta/2)
   \sum_i |c_i| (2a_i + l + 1)^{o} \sqrt{\textstyle\sum_k |d_k|}

for derivatives of order o, where ``shell_bound`` bounds
:math:`\max_r |r - A|^l e^{-\alpha |r - A|^2} e^{-\eta r^2 / 2}` (a Gaussian product with an
optimised split of the exponent). Each shell takes half of the ECP's Gaussian, so the product of
two shells' values estimates the integral over the pair; it is computed once per (shell, ECP)
when the integrals are requested. Pairs whose product is below ``tol`` are skipped. Before that,
a whole shell is dropped against an ECP if its value times the largest value of any shell against
that ECP is below ``tol``: then every pair it forms would be skipped, so this test only saves
work. The estimates are not proven bounds, so they are tested against the exact integrals they
skip, over a sweep of parameters and over every test system.

These bounds build on the shell-bounding Gaussians of McKenzie et al. [McKenzie2018]_:
``shell_bound`` is their shell factor, and ``channel_bounds`` (below) their three-centre projected
bound, each extended as described in ``screening.hpp``. The pair estimate is our own, chosen so
that it can be tabulated once per (shell, ECP). Their shell test, which compares a shell's factor
with a bound on any partner shell over all space, is replaced by the largest actual partner value,
because that bound does not cover the coefficients and derivative factors in
``shell_ecp_bound``.

**Channels.** Within a shell pair, ``channel_bounds`` bounds the full contribution of each ECP
channel (polynomial factors and angular sums included), extending the bound of [McKenzie2018]_
to ECPs with powers of r other than :math:`r^0`; channels below 1e-2 times
``radial_tol`` are skipped. The factor 1e-2 allows for the tens of blocks, each possibly missing
skipped channels, that are summed into one element of the second derivatives.

**Primitive triples.** In each radial integral, ``type2_bound`` bounds :math:`|G|` for every
power N and Bessel order at once, after a cheap check of its logarithm (one square root, no
exponentials) that discards most negligible triples. A triple is skipped below

.. math::

   \frac{\mathrm{radial\_tol}}{(2 + 2\max(A, B))^{L_A + L_B}\, n_\mathrm{triples}},

which allows for the binomial factors applied downstream and keeps the skipped triples together
below the tolerance. The integration window (``integrand_window``) is derived in the same way.

Precision limits
^^^^^^^^^^^^^^^^

Expanding each shell about the ECP centre makes the terms of the expansion larger than the result
by up to :math:`(1 + d\sqrt{2a})^{l}` for a primitive with exponent a at distance d (l
including the angular momentum added by derivatives), so the relative error of a block grows as
:math:`\sim 10^{-16} (1 + d\sqrt{2a})^{2l}`: about 1e-15 for p shells, but 6e-8 for an h shell with
exponent 8.7 two bohr away, and more at the limit of l = 8. Removing it needs a formulation
without the expansion; this is noted in ``ARCHITECTURE.md`` as future work. The test suite allows
for it explicitly (below).

Testing
=======

All tests use GoogleTest and run under ``ctest``. From the build directory:

.. code-block:: bash

   ctest -j8 -LE reference            # everything except the live reference tests: ~10 s
   ctest -R Suite.Stored              # one group, by regular expression
   ctest -R Radial --output-on-failure

Configuring with ``-DLIBECPINT_TEST_REFERENCE=ON`` adds slower tests (label ``reference``) that
compare against the reference engine as it runs; run them when the reference engine changes, or
to locate a failure.

Layout
^^^^^^

The tests are built in layers, so that a failure points at one part of the library:

+--------------------+------------------------------------+---------------------------------------+
| Directory          | ctest names                        | What it checks                        |
+====================+====================================+=======================================+
| ``tests/units``    | ``MathUtil``, ``Bessel``,          | single components, against exact or   |
|                    | ``GaussQuad``, ``Assembly``,       | independent values                    |
|                    | ``DerivativeRecurrences``,         |                                       |
|                    | ``Screening``, ``SharedGrid``,     |                                       |
|                    | ``Api``, ``BField``, ...           |                                       |
+--------------------+------------------------------------+---------------------------------------+
| ``tests/radial``   | ``Radial.*``                       | type 2 radial integrals against       |
|                    |                                    | independent quadrature                |
+--------------------+------------------------------------+---------------------------------------+
| ``tests/analytic`` | ``Analytic``                       | all centres coincident, against       |
|                    |                                    | closed forms                          |
+--------------------+------------------------------------+---------------------------------------+
| ``tests/suite``    | ``Suite.<check>.<system>``         | whole molecules: integrals, first and |
|                    |                                    | second derivatives, magnetic field    |
|                    |                                    | derivatives                           |
+--------------------+------------------------------------+---------------------------------------+

The **reference engine** (``tests/support/integral_reference.hpp``) computes the same integrals
by a different method: it integrates the projections numerically on the sphere, with grids
aligned to each Gaussian, and the radial parts on composite panels. It shares no code with the
library, expands nothing about the ECP centre, and screens nothing beyond ~1e-40. It is slow, so
its values for each test system (listed in ``tests/data/systems.txt``) are stored in
``tests/data/reference``.

For each system, the suite runs four checks:

- **Stored**: the library, with default settings, against the stored reference values;
- **Live**: the same against the reference engine computed now, element by element;
- **Screening**: the library with its default screening against the library with none, which
  isolates screening errors from integration errors;
- **Threads**: every kernel run from several threads sharing one ``ECPIntegral`` gives results
  identical to the serial ones.

**Tolerances.** Each value must satisfy
:math:`|v - v_\mathrm{ref}| \le 10^{-13} \max|v_\mathrm{ref}| + 10^{-14}\, s`, where the first
term is the agreement expected between two accurate methods, and the second allows for the
precision limit above: the scale s of each value estimates the size of the terms in its expansion
about the ECP. The unit, radial and analytic tests set their own tolerances, typically ~50x the
measured error.

Diagnosing a failure
^^^^^^^^^^^^^^^^^^^^

+--------------------------------+---------------------------+--------------------------------------+
| Fails                          | Passes                    | Likely cause                         |
+================================+===========================+======================================+
| ``Assembly``                   |                           | packing of derivatives or the        |
|                                |                           | assembly loops                       |
+--------------------------------+---------------------------+--------------------------------------+
| ``DerivativeRecurrences``      | ``Assembly``              | the derivative recurrences           |
+--------------------------------+---------------------------+--------------------------------------+
| ``BField``                     | ``Assembly``              | the magnetic field derivatives       |
+--------------------------------+---------------------------+--------------------------------------+
| ``Screening`` or               | radial, analytic          | a screening bound or estimate        |
| ``Suite.Screening``            |                           |                                      |
+--------------------------------+---------------------------+--------------------------------------+
| ``Radial.*``, ``SharedGrid``,  |                           | radial quadrature, the shared grid,  |
| ``Bessel``, ``GaussQuad``      |                           | Bessel functions                     |
+--------------------------------+---------------------------+--------------------------------------+
| ``Analytic``, library only     |                           | on-centre integrals                  |
+--------------------------------+---------------------------+--------------------------------------+
| ``Analytic``, reference engine |                           | the reference engine itself: do not  |
|                                |                           | trust the stored values              |
+--------------------------------+---------------------------+--------------------------------------+
| ``Suite.Live``                 | all of the above          | the angular part or assembly of the  |
|                                |                           | integrals                            |
+--------------------------------+---------------------------+--------------------------------------+
| ``Suite.Stored`` only          | ``Suite.Live``            | stale stored values (regenerate)     |
+--------------------------------+---------------------------+--------------------------------------+
| ``Suite.Threads``              |                           | a data race (run under               |
|                                |                           | ThreadSanitizer)                     |
+--------------------------------+---------------------------+--------------------------------------+

A failing comparison names the worst element by atom, shell and Cartesian component, with its
error and tolerance. Most accuracy tests also print ``[ info ]`` lines with their worst errors,
so that a drift shows before a test fails.

Stored references and new tests
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

After changing a test system or the reference engine (never to make a library change pass),
regenerate the stored values with

.. code-block:: bash

   cmake --build . --target regenerate_references   # ~3 minutes

To add a system, add its geometry, basis (Gaussian94 format) and any ECP not in the library to
``tests/data``, add a line to ``tests/data/systems.txt``, re-run CMake, and generate its reference
with ``tests/suite/regenerate_refs NAME``. ``tests/README.md`` describes how to add radial,
screening and unit tests.

Benchmarking
============

The timing benchmark, ``ecpint_bench``, is built with the tests but does not run under ``ctest``.
Use a release build, on an otherwise idle machine. From the build directory:

.. code-block:: bash

   make Benchmark                                     # kernels and systems, ~20 s
   benchmarks/ecpint_bench scaling --max-carbons 128 --wasted
   benchmarks/ecpint_bench --csv out.csv --tag mybuild

It has three groups:

- **kernels**: ``compute_shell_pair`` for each pair of angular momenta up to g, against the def2
  ECP of Au, with the shells off the ECP centre, and (as controls) on it;
- **systems**: ``init``, ``compute_integrals``, ``compute_gradient`` and ``compute_hessian`` for
  test systems and an iodoalkane chain, each from cold caches;
- **scaling**: integrals and gradients for iodoalkane chains of doubling length, with the local
  exponent of the growth (``--wasted`` also counts the shell-pair calls that pass the screening
  but return only negligible values).

Each measurement reports the median, minimum and spread of repeated runs, and a checksum of the
values computed, so that a change in the results (rather than the speed) is caught. To compare
builds, run them interleaved, then use ``benchmarks/compare.py``:

.. code-block:: bash

   for round in 1 2 3; do
     for b in old new; do build-$b/benchmarks/ecpint_bench --csv $b-$round.csv --tag $b; done
   done
   benchmarks/compare.py old-*.csv new-*.csv --base old

Latest results
^^^^^^^^^^^^^^

v1.1 against v1.0.7, measured in November 2026 on one core of an Apple M3 Pro, both built with
GCC 16 (``-O3``) and their default settings, over three interleaved rounds (noise 1--9%). Ratios
below 1 are faster. v1.0.7 has no ``compute_gradient`` or ``compute_hessian``, so for it the
benchmark ran the same contraction loops over v1.0.7's own kernels, unscreened as in its
``compute_first_derivs``. This reproduces v1.0.7's own derivatives to 1e-14, and, by avoiding
storing the dense matrices, if anything flatters v1.0.7.

+-------------------------------------+-----------+-----------+---------------+
| Measurement                         | v1.0.7    | v1.1      | Ratio         |
+=====================================+===========+===========+===============+
| kernels, off the ECP (15 classes,   | 67 µs --  | 6.7 µs -- | 0.056 -- 0.10 |
| ss to gg)                           | 6.5 ms    | 0.44 ms   | (mean 0.066)  |
+-------------------------------------+-----------+-----------+---------------+
| kernels, on the ECP (controls)      |           |           | mean 0.053    |
+-------------------------------------+-----------+-----------+---------------+
| Ag\ :sub:`3` (valence basis):       | 54 / 186  | 5.4 / 19  | 0.10 / 0.10 / |
| integrals / gradient / Hessian (ms) | / 551     | / 51      | 0.091         |
+-------------------------------------+-----------+-----------+---------------+
| La\ :sub:`2`\ O\ :sub:`3`:          | 13 / 38   | 1.7 / 5.3 | 0.13 / 0.14 / |
| integrals / gradient / Hessian (ms) | / 93      | / 14      | 0.14          |
+-------------------------------------+-----------+-----------+---------------+
| CH\ :sub:`3`\ I, def2-SVP:          | 71 / 128  | 4.7 / 10  | 0.065 / 0.080 |
| integrals / gradient / Hessian (ms) | / 285     | / 23      | / 0.080       |
+-------------------------------------+-----------+-----------+---------------+
| Au\ :sub:`2`, def2-TZVPP:           | 221 / 571 | 17 / 52   | 0.078 / 0.091 |
| integrals / gradient (ms)           |           |           |               |
+-------------------------------------+-----------+-----------+---------------+
| C\ :sub:`16` iodoalkane (50 atoms): | 0.30 /    | 0.060 /   | 0.19 / 0.26 / |
| integrals / gradient / Hessian (s)  | 0.95 /    | 0.25 /    | 0.29          |
|                                     | 2.7       | 0.79      |               |
+-------------------------------------+-----------+-----------+---------------+
| C\ :sub:`128` iodoalkane (386       | 4.1 / 107 | 0.55 /    | 0.13 / 0.022  |
| atoms): integrals / gradient (s)    |           | 2.4       |               |
+-------------------------------------+-----------+-----------+---------------+
| all systems (geometric mean)        |           |           | 0.134         |
+-------------------------------------+-----------+-----------+---------------+

Local scaling exponents for the iodoalkane chains (time ~ N\ :sup:`p` in the number of atoms):

+--------------+-----------------------+-----------------------+
| Atoms        | Integrals: v1.0.7 /   | Gradient: v1.0.7 /    |
|              | v1.1                  | v1.1                  |
+==============+=======================+=======================+
| 26 -- 50     | 1.37 / 1.38           | 1.57 / 1.45           |
+--------------+-----------------------+-----------------------+
| 50 -- 98     | 1.24 / 1.12           | 1.80 / 1.17           |
+--------------+-----------------------+-----------------------+
| 98 -- 194    | 1.24 / 1.09           | 2.35 / 1.06           |
+--------------+-----------------------+-----------------------+
| 194 -- 386   | 1.33 / 1.04           | 2.76 / 1.05           |
+--------------+-----------------------+-----------------------+

Every checksum differs between the versions, because v1.0.7's results are less accurate. Against
v1.1 (which agrees with the reference engine to ~1e-13), v1.0.7's largest errors, relative to
the largest value, are:

+----------------------------------------------+-------------+-------------+-------------+
| System                                       | Integrals   | Gradient    | Hessian     |
+==============================================+=============+=============+=============+
| Ag\ :sub:`3`, ecp28mdf                       | 1e-9        | 6e-7        | 1e-5        |
+----------------------------------------------+-------------+-------------+-------------+
| Ag\ :sub:`3`, LANL2DZ                        | 3e-8        | 9e-5        | 4e-4        |
+----------------------------------------------+-------------+-------------+-------------+
| CH\ :sub:`3`\ I, def2-SVP                    | 7e-9        | 1e-7        | 1e-5        |
+----------------------------------------------+-------------+-------------+-------------+
| Au\ :sub:`2`, def2-TZVPP                     | 4e-9        | 5e-6        |             |
+----------------------------------------------+-------------+-------------+-------------+
| La\ :sub:`2`\ O\ :sub:`3` (regression basis) | 2e-5        | 4e-3        | 0.19        |
+----------------------------------------------+-------------+-------------+-------------+
| testpair (diffuse f shells)                  | 2e-4        | 0.06        | 24          |
+----------------------------------------------+-------------+-------------+-------------+

The last two systems are designed to be hard (tight and diffuse functions of high angular
momentum away from the ECP). For a single atom (Pt, def2-TZVP), where only the on-centre code is
used, the two versions agree to 1e-14.

Contributing
============

Contributions are welcome, as issues or pull requests on GitHub. A change to the library should
keep the tests passing, and come with a test that fails without it. Code is formatted with
``clang-format`` (the configuration is in the repository, and CI checks it), so run it before
committing. Functionality for users should be exposed through ``api.hpp``, without changing the
existing interface.

The documentation is built with ``make docs`` (see :ref:`sec:installation`), from the
reStructuredText files in ``doc/sphinx/source`` and Doxygen comments in the headers.

The generator for closed forms of the radial integrals, used by earlier versions, is kept, unused,
in ``src/generated/radial``; its ``README.md`` explains why.

.. [Shaw2017] R. A. Shaw, J. G. Hill, J. Chem. Phys. 147, 074108 (2017);
   doi: `10.1063/1.4986887 <https://doi.org/10.1063/1.4986887>`_

.. [McKenzie2018] S. C. McKenzie, E. Epifanovsky, G. M. J. Barca, A. T. B. Gilbert,
   P. M. W. Gill, J. Phys. Chem. A 122 (2018);
   doi: `10.1021/acs.jpca.7b12679 <https://doi.org/10.1021/acs.jpca.7b12679>`_

.. toctree::
   :hidden:
