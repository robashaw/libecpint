.. libecpint usage file

.. _`sec:usage`:

=====
Usage
=====

libecpint can be used in two ways:

- the **high-level API**, ``ECPIntegrator``: you give it the basis set, the ECPs and the
  geometry, and it returns whole matrices of integrals, or their derivatives contracted with a
  density. This is what most programs need, and the place to start.
- the **low-level API**, ``ECPIntegral``: you compute the integrals over one shell pair and one
  ECP at a time, and assemble them yourself. This suits programs with their own integral loops,
  e.g. to parallelise over shell pairs.

Both are included with ``#include <libecpint.hpp>``, and live in the ``libecpint`` namespace.

Quick start
===========

A complete program, for HI with the iodine ECP from the built-in library, is in
``example/example.cpp`` (build instructions in ``example/README.md``). In outline:

.. code-block:: c++

   #include <libecpint.hpp>
   using namespace libecpint;

   ECPIntegrator factory;

   // 1. The orbital basis: one entry per shell (coordinates in bohr), primitives concatenated
   factory.set_gaussian_basis(nshells, g_coords, g_exps, g_coefs, g_ams, g_lengths);

   // 2. The ECPs, here from the built-in library
   std::vector<std::string> names = {"ecp28mdf"};
   int charges[1] = {53};
   factory.set_ecp_basis_from_library(1, I_coords, charges, names, "/path/to/share/libecpint");

   // 3. Initialise, for derivatives up to second order
   factory.init(2);

   // 4. Compute
   factory.compute_integrals();
   std::shared_ptr<std::vector<double>> V = factory.get_integrals();  // ncart x ncart
   std::vector<double> gradient = factory.compute_gradient(D);        // 3 * natoms
   std::vector<double> hessian = factory.compute_hessian(D);          // 3 natoms x 3 natoms

where ``D`` is a density matrix (``ncart x ncart``, row-major). The rest of this page explains
each step.

Conventions
===========

Read this section before passing in your own basis sets: most problems come from one of these.

- **Units.** All coordinates are in bohr.
- **Cartesian functions.** All integrals and derivatives are over Cartesian Gaussians, in
  canonical order within a shell (for d: ``xx, xy, xz, yy, yz, zz``), with shells in the order
  given to ``set_gaussian_basis``. Converting to spherical harmonics is left to the calling
  program. ``factory.ncart`` is the total number of Cartesian functions.
- **Matrices** are stored row-major: element ``(i, j)`` is ``[i * ncart + j]``.
- **Normalisation.** Contraction coefficients are used exactly as given, so they must include
  the normalisation of the primitives. One coefficient is used for every Cartesian component
  of a shell, so (as is usual) only the ``x^l`` component is normalised by it.
- **Atoms** are identified by position: shells and ECPs within 1e-4 bohr of each other are on
  the same atom. Atoms are numbered in the order they first appear in ``set_gaussian_basis``,
  then any ECP-only centres; this is the atom order of every derivative.
- **ECP powers.** Each ECP primitive is ``d r^(n-2) exp(-a r^2)``: ``n`` *includes* the
  ``r^2`` of the spherical Jacobian, as in the Stuttgart/Cologne ECPs and the MOLPRO format, so
  ``n`` is normally 0, 1 or 2 (never more than 4). If your ECP's powers are mostly -2, -1 or 0,
  it follows the other convention: add 2. If in doubt, compare with the built-in library.
- **Angular momentum** is supported up to 8 (l functions) for the ECPs, and for the orbital
  basis *including* the order of derivative: integrals up to l (L = 8), first derivatives up to
  k (L = 7), and second derivatives up to i (L = 6). Beyond that the library throws
  ``std::invalid_argument``.

High-level API
==============

Setting up
^^^^^^^^^^

An ``ECPIntegrator`` is set up in three steps, in this order: the orbital basis, the ECPs, then
``init``.

.. code-block:: c++

   factory.set_gaussian_basis(nshells, g_coords, g_exps, g_coefs, g_ams, g_lengths);
   factory.set_ecp_basis(necps, u_coords, u_exps, u_coefs, u_ams, u_ns, u_lengths);
   factory.init(deriv_order);

+-------------------------+---------------------------------------------------------------------+
| Argument                | Meaning                                                             |
+=========================+=====================================================================+
| ``nshells``, ``necps``  | number of shells in the orbital basis; number of ECP centres        |
+-------------------------+---------------------------------------------------------------------+
| ``g_coords``,           | coordinates of each shell / ECP centre, ``x, y, z`` in turn         |
| ``u_coords``            | (length 3 x number of shells / centres)                             |
+-------------------------+---------------------------------------------------------------------+
| ``g_exps``, ``g_coefs`` | exponents and coefficients of all primitives, shell after shell     |
+-------------------------+---------------------------------------------------------------------+
| ``g_ams``               | angular momentum of each shell                                      |
+-------------------------+---------------------------------------------------------------------+
| ``g_lengths``           | number of primitives in each shell                                  |
+-------------------------+---------------------------------------------------------------------+
| ``u_exps``, ``u_coefs``,| exponent, coefficient, angular momentum and power ``n`` of every    |
| ``u_ams``, ``u_ns``     | ECP primitive, centre after centre (the local channel has the       |
|                         | highest ``l`` of its ECP)                                           |
+-------------------------+---------------------------------------------------------------------+
| ``u_lengths``           | number of primitives on each ECP centre                             |
+-------------------------+---------------------------------------------------------------------+
| ``deriv_order``         | highest derivative order needed: 0, 1 or 2                          |
+-------------------------+---------------------------------------------------------------------+

``init`` also takes optional screening settings (see Settings_). ECP input is validated: an
angular momentum or power out of range throws ``std::invalid_argument``.

Instead of ``set_ecp_basis``, ECPs can be taken from the library supplied with libecpint (see
`The built-in ECP library`_):

.. code-block:: c++

   factory.set_ecp_basis_from_library(necps, u_coords, u_charges, u_names, share_dir);

where ``u_charges`` are the atomic numbers of the ECP atoms, ``u_names`` the name of each one's
ECP (e.g. ``"ecp28mdf"``), and ``share_dir`` the path to the installed ``share/libecpint``
directory. It throws ``std::runtime_error`` if an ECP is not found.

Integrals
^^^^^^^^^

.. code-block:: c++

   factory.compute_integrals();
   std::shared_ptr<std::vector<double>> V = factory.get_integrals();
   double Vij = (*V)[i * factory.ncart + j];

Gradients and Hessians
^^^^^^^^^^^^^^^^^^^^^^

Most uses of the derivatives (energy gradients, the integral part of the Hessian) only need
them contracted with a density matrix. This is the recommended way to get them: the derivative
matrices are never stored, so the memory needed is about that of the integrals.

.. code-block:: c++

   // D: ncart x ncart, row-major, in the order of the integrals; it need not be symmetric
   std::vector<double> gradient = factory.compute_gradient(D);  // needs init(1) or init(2)
   std::vector<double> hessian = factory.compute_hessian(D);    // needs init(2)

- ``gradient[3 * n + i]`` is the sum over ``k, l`` of ``D[k * ncart + l]`` times the derivative
  of integral ``(k, l)`` with respect to coordinate ``i`` (0, 1, 2 for x, y, z) of atom ``n``.
- ``hessian`` is the full, symmetric ``3 natoms x 3 natoms`` matrix of the corresponding second
  derivatives, row-major: ``hessian[(3 * m + i) * 3 * natoms + 3 * n + j]``.

Derivative matrices
^^^^^^^^^^^^^^^^^^^

If the derivative matrices themselves are needed (e.g. for coupled-perturbed equations), they
can be stored in full. Beware of the memory: ``3 natoms`` matrices for the first derivatives
(about 7 GB for 77 atoms and 1900 basis functions), and ``3 natoms (3 natoms + 1) / 2`` for the
second (about 800 GB for the same molecule).

.. code-block:: c++

   factory.compute_first_derivs();   // needs init(1) or init(2)
   auto d1 = factory.get_first_derivs();
   double dVij_dy_n = (*d1[3 * n + 1])[i * factory.ncart + j];

   factory.compute_second_derivs();  // needs init(2)
   auto d2 = factory.get_second_derivs();

The first derivatives are ordered ``x, y, z`` for each atom in turn. The second derivatives are
stored once per pair of atoms, in blocks ``AA, AB, AC, ..., BB, BC, ..., CC, ...``; within a
block the components are ``xx, xy, xz, yy, yz, zz`` for a diagonal block (e.g. ``AA``), and
``xx, xy, xz, yx, yy, yz, zx, zy, zz`` for an off-diagonal block (e.g. ``AB``). The macro
``H_START(a, b, natoms)`` (``assembly.hpp``) gives the start of the block for atoms ``a < b``;
add 3 for a diagonal block (``a == b``). For example, the derivative with respect to ``y`` of atom
0 and ``x`` of atom 2 is

.. code-block:: c++

   auto d2_Ay_Cx = d2[H_START(0, 2, factory.natoms) + 3];  // yx is component 3

Magnetic field derivatives
^^^^^^^^^^^^^^^^^^^^^^^^^^

For a basis of gauge-including atomic orbitals (London orbitals), the derivatives of the ECP
integrals with respect to a uniform magnetic field (e.g. for NMR shieldings) are

.. code-block:: c++

   factory.compute_bfield_derivs();  // needs init(1) or init(2)
   auto dB = factory.get_bfield_derivs();

giving three matrices, for ``x``, ``y`` and ``z``. Each is real and antisymmetric, the derivative
being ``i/2`` times it, and the ECP's projectors are gauge-transformed to the ECP centre. See
``ECPIntegral::compute_shell_pair_derivative_Bfield`` in ``ecpint.hpp`` for the precise
definition.

New geometries
^^^^^^^^^^^^^^

To move the atoms (e.g. in a geometry optimisation), pass the new coordinates, in the same order
as before, then compute again:

.. code-block:: c++

   factory.update_gaussian_basis_coords(nshells, g_coords);
   factory.update_ecp_basis_coords(necps, u_coords);

The pointers returned by the ``get_*`` functions refer to copies, so call them again after each
computation.

.. _Settings:

Settings
^^^^^^^^

.. code-block:: c++

   factory.init(deriv_order, tol, radial_tol);

+-----------------+-----------+------------------------------------------------------------------+
| Parameter       | Default   | Meaning                                                          |
+=================+===========+==================================================================+
| ``deriv_order`` | 0         | highest derivative order needed (0, 1 or 2)                      |
+-----------------+-----------+------------------------------------------------------------------+
| ``tol``         | 1e-12     | screening tolerance for shells and shell pairs: their            |
|                 |           | contributions are skipped when estimated to be below it          |
+-----------------+-----------+------------------------------------------------------------------+
| ``radial_tol``  | 1e-15     | screening tolerance within each shell pair (ECP channels and     |
|                 |           | primitives); 0 skips nothing                                     |
+-----------------+-----------+------------------------------------------------------------------+

Two further arguments, ``smallGrid`` and ``bigGrid``, are accepted for compatibility with earlier
versions but no longer have any effect: the radial quadrature is fixed, and converged for every
case (see :ref:`sec:developer`).

The defaults give integrals and derivatives accurate to about 1e-13 relative to their largest
value, for typical basis sets. Loosening ``tol`` (e.g. to 1e-10) speeds up large molecules at the
cost of accuracy; there is no need to tighten it. Integrals over very high angular momentum shells
far from an ECP are less precise than this (see `Accuracy`_).

Performance
^^^^^^^^^^^

Version 1.1 is much faster than v1.0.7. Timings on one core of an Apple M3 Pro (GCC 16,
``-O3``), each calculation starting from cold caches, with gradients and Hessians contracted with
a density:

+-----------------------------------------+-------------+-----------+---------+
| Calculation                             | v1.1        | v1.0.7    | Speedup |
+=========================================+=============+===========+=========+
| CH\ :sub:`3`\ I, def2-SVP: integrals    | 4.7 ms      | 71 ms     | 15x     |
+-----------------------------------------+-------------+-----------+---------+
| CH\ :sub:`3`\ I, def2-SVP: Hessian      | 23 ms       | 285 ms    | 12x     |
+-----------------------------------------+-------------+-----------+---------+
| Au\ :sub:`2`, def2-TZVPP: integrals     | 17 ms       | 221 ms    | 13x     |
+-----------------------------------------+-------------+-----------+---------+
| Au\ :sub:`2`, def2-TZVPP: gradient      | 52 ms       | 571 ms    | 11x     |
+-----------------------------------------+-------------+-----------+---------+
| Iodoalkane, 50 atoms: Hessian           | 0.79 s      | 2.7 s     | 3.4x    |
+-----------------------------------------+-------------+-----------+---------+
| Iodoalkane, 386 atoms: integrals        | 0.55 s      | 4.1 s     | 7.5x    |
+-----------------------------------------+-------------+-----------+---------+
| Iodoalkane, 386 atoms: gradient         | 2.4 s       | 107 s     | 45x     |
+-----------------------------------------+-------------+-----------+---------+

Over all the benchmark calculations, v1.1 is 7.5x faster (geometric mean), and single shell pairs
are about 15x faster. Gradients now grow almost linearly with the size of the molecule (as
N\ :sup:`1.05` for the largest chains, against N\ :sup:`2.8` for v1.0.7, which did not screen
derivatives). The results are also more accurate: for ordinary basis sets, v1.0.7's Hessians were
off by 1e-5 to 4e-4 relative to their largest element, against ~1e-13 now. The
:ref:`sec:developer` gives the full results.

``ECPIntegrator`` runs on one thread. To use several, call the low-level API from your own
threads: one ``ECPIntegral`` may be shared, and its ``compute_*`` functions are safe to call
concurrently (each thread keeps its own caches, and the results do not depend on the number of
threads).

Accuracy
^^^^^^^^

Every integral and derivative is checked against an independent numerical reference engine,
over a range of molecules and ECPs (see :ref:`sec:developer`). The agreement is about 1e-13
relative to the largest value, with one known limitation: the method expands each shell about the
ECP centre, which loses precision for high angular momentum shells with tight exponents far from
the ECP: an h shell with exponent 8.7, 2 bohr from the ECP, gives blocks of integrals accurate
to ~1e-7 relative, against ~1e-15 for p shells. These blocks are small in absolute terms, and
for the angular momenta and exponents of typical valence basis sets the effect is negligible.

The built-in ECP library
^^^^^^^^^^^^^^^^^^^^^^^^

The ECPs supplied with libecpint are in ``share/libecpint/raw``, one MOLPRO-format file per
family. The Stuttgart/Cologne energy-consistent pseudopotentials are named ``ecpnxy`` after
their ECPnXY label: n is the number of core electrons, X is S (single-electron fit) or M
(multi-electron fit), and Y is HF (non-relativistic), WB (quasi-relativistic) or DF (fully
relativistic). They come from the Stuttgart/Cologne library (http://www.tc.uni-koeln.de/PP),
except where the ECPs of the cc-pVnZ-PP basis sets on the Basis Set Exchange differ (currently
only Y in ``ecp28mdf``), in which case the Basis Set Exchange version is used. The other families
come from the Basis Set Exchange. The Karlsruhe ``lcecp-n`` sets are the Stuttgart 4f-in-core
ECPnMWB potentials grouped by valence, except that Nd in ``lcecp-1`` has different parameters
from ``ecp49mwb``. Each file records the source and references for every atom.

.. list-table::
   :widths: 25 75
   :header-rows: 1

   * - Name
     - Atoms
   * - ``ecp2mhf``
     - Ne
   * - ``ecp2mwb``
     - B -- Ne
   * - ``ecp2sdf``
     - Li -- F
   * - ``ecp10mdf``
     - K -- Kr
   * - ``ecp10mhf``
     - Ar, Ca -- Ga
   * - ``ecp10mwb``
     - Al -- Ca, Cu -- Ga
   * - ``ecp10sdf``
     - Na -- Cl
   * - ``ecp18mhf``
     - Ti
   * - ``ecp18mwb``
     - Ga
   * - ``ecp18sdf``
     - K, Ca
   * - ``ecp18shf``
     - K, Ca
   * - ``ecp28mdf``
     - Ge -- Br, Rb -- Xe
   * - ``ecp28mhf``
     - Kr, Sr -- In, Ce -- Yb
   * - ``ecp28mwb``
     - Zn -- In, La -- Lu
   * - ``ecp28oldmdf``
     - I
   * - ``ecp28sdf``
     - Cu -- Br
   * - ``ecp28shf``
     - Cu, Ge
   * - ``ecp36mwb``
     - In
   * - ``ecp36sdf``
     - Rb, Sr
   * - ``ecp36shf``
     - Rb, Sr
   * - ``ecp46mdf``
     - Sn -- I, Cs, Ba
   * - ``ecp46mhf``
     - Xe, Ba, La
   * - ``ecp46mwb``
     - Cd -- Ce
   * - ``ecp46sdf``
     - Ag -- I
   * - ``ecp46shf``
     - Ag
   * - ``ecp47mhf``
     - La, Ce
   * - ``ecp47mwb``
     - La -- Pr
   * - ``ecp48mhf``
     - Ce, Pr
   * - ``ecp48mwb``
     - Ce -- Nd
   * - ``ecp49mhf``
     - Pr, Nd
   * - ``ecp49mwb``
     - Pr, Nd
   * - ``ecp50mhf``
     - Nd, Pm
   * - ``ecp50mwb``
     - Nd, Pm
   * - ``ecp51mhf``
     - Pm, Sm
   * - ``ecp51mwb``
     - Pm, Sm
   * - ``ecp52mhf``
     - Sm, Eu
   * - ``ecp52mwb``
     - Sm, Eu
   * - ``ecp53mhf``
     - Eu, Gd
   * - ``ecp53mwb``
     - Eu -- Tb
   * - ``ecp54mhf``
     - Gd, Tb
   * - ``ecp54mwb``
     - Gd -- Dy
   * - ``ecp54sdf``
     - Cs, Ba
   * - ``ecp54shf``
     - Cs
   * - ``ecp55mhf``
     - Tb, Dy
   * - ``ecp55mwb``
     - Tb, Dy
   * - ``ecp56mhf``
     - Dy, Ho
   * - ``ecp56mwb``
     - Dy, Ho
   * - ``ecp57mhf``
     - Ho, Er
   * - ``ecp57mwb``
     - Ho, Er
   * - ``ecp58mhf``
     - Er, Tm
   * - ``ecp58mwb``
     - Er, Tm
   * - ``ecp59mhf``
     - Tm, Yb
   * - ``ecp59mwb``
     - Tm, Yb
   * - ``ecp60mdf``
     - Yb, Hf -- Rn, Ac -- U
   * - ``ecp60mhf``
     - Yb -- Tl, Ac -- Lr
   * - ``ecp60mwb``
     - Yb -- Tl, Ac -- Lr
   * - ``ecp60oldmdf``
     - Au, Hg
   * - ``ecp68mdf``
     - Yb
   * - ``ecp68mwb``
     - Tl
   * - ``ecp78mdf``
     - Pb -- At, Fr, Ra
   * - ``ecp78mhf``
     - Hg -- Rn
   * - ``ecp78mwb``
     - Hg -- Rn, Ac -- U
   * - ``ecp78sdf``
     - Au
   * - ``ecp78shf``
     - Au
   * - ``ecp79mwb``
     - Th -- U
   * - ``ecp80mwb``
     - Pa -- Np
   * - ``ecp81mwb``
     - U -- Pu
   * - ``ecp82mwb``
     - Np -- Am
   * - ``ecp83mwb``
     - Pu, Am
   * - ``ecp84mwb``
     - Pu -- Cm
   * - ``ecp85mwb``
     - Am -- Bk
   * - ``ecp86mwb``
     - Cm -- Cf
   * - ``ecp87mwb``
     - Bk, Cf
   * - ``ecp88mwb``
     - Cf, Es
   * - ``ecp89mwb``
     - Es, Fm
   * - ``ecp90mwb``
     - Fm, Md
   * - ``ecp91mwb``
     - Md, No
   * - ``ecp92mdf``
     - Db
   * - ``ecp92mdfb``
     - Cn -- Og
   * - ``ecp92mdfq``
     - Rg -- Og
   * - ``ecp92mhf``
     - Db
   * - ``ecp92mwb``
     - No, Lr, Db
   * - ``def2-ecp``
     - def2-ECP (Karlsruhe; def2-SD in Turbomole): Rb -- Rn
   * - ``dhf-ecp``
     - dhf-ECP (Karlsruhe): Rb -- Ba, Hf -- Rn
   * - ``lanl2dz``
     - LANL2DZ (Los Alamos): Na -- La, Hf -- Bi, U -- Pu
   * - ``crenbl``
     - CRENBL (large core): Li -- Og
   * - ``crenbs``
     - CRENBS (small core): Sc -- Kr, Y -- Xe, La, Hf -- Rn, Rf -- Og
   * - ``sbkjc``
     - SBKJC: Li -- Ce, Hf -- Rn
   * - ``lcecp-0``
     - lcecp-0 (Karlsruhe; 4f-in-core, divalent): La -- Yb
   * - ``lcecp-1``
     - lcecp-1 (Karlsruhe; 4f-in-core, trivalent): La -- Lu
   * - ``lcecp-2``
     - lcecp-2 (Karlsruhe; 4f-in-core, tetravalent): Ce -- Nd, Tb, Dy

**Adding an ECP.** Save it in MOLPRO format as ``share/libecpint/raw/NAME.ecp``, with ``NAME``
on the first line, and exponents written with ``E`` rather than ``D``. A MOLPRO export from the
`Basis Set Exchange <https://www.basissetexchange.org>`_ can be used as it is: any orbital basis
and spin-orbit terms in the file are ignored. Please consider a pull request, so that everyone
can use it. A file outside the library can be read with ``ECPBasis::addECP_from_file``.

Low-level API
=============

Shells and ECPs
^^^^^^^^^^^^^^^

The low-level functions take ``GaussianShell`` and ``ECP`` objects:

.. code-block:: c++

   double A[3] = {Ax, Ay, Az};       // bohr
   GaussianShell shellA(A, l);
   shellA.addPrim(exponent, coef);   // for each primitive

   double C[3] = {Cx, Cy, Cz};
   ECP U(C);
   U.addPrimitive(n, l, exponent, coef);  // for each primitive (n as in Conventions)

The same conventions apply as for the high-level API. ``ECPBasis::addECP_from_file`` reads an
ECP from a MOLPRO-format file.

The ECPIntegral object
^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: c++

   ECPIntegral ecpint(maxLB, maxLU, deriv_order);

where ``maxLB`` and ``maxLU`` are the highest angular momenta in the orbital basis and in the ECPs,
and ``deriv_order`` the highest derivative order needed (default 0). An optional fourth argument
is the radial screening tolerance (``radial_tol`` above, default 1e-15). Construct it once and
reuse it: the constructor builds the angular integral tables. Its compute functions are
``const``, and can be called concurrently from several threads.

Integrals
^^^^^^^^^

.. code-block:: c++

   TwoIndex<double> values;
   ecpint.compute_shell_pair(U, shellA, shellB, values);
   double v = values(i, j);  // or values.data[i * ncartB + j]

``values`` has dimensions ``(ncartA, ncartB)``, in the order described in Conventions.

First derivatives
^^^^^^^^^^^^^^^^^

.. code-block:: c++

   std::array<TwoIndex<double>, 9> d1;
   ecpint.compute_shell_pair_derivative(U, shellA, shellB, d1);

gives the derivatives with respect to ``Ax, Ay, Az, Bx, By, Bz, Cx, Cy, Cz``, where ``A``, ``B``
and ``C`` are the centres of ``shellA``, ``shellB`` and ``U``. They are designed to be added: if
``A == B``, the derivative with respect to that atom's x coordinate is ``Ax + Bx``.

Second derivatives
^^^^^^^^^^^^^^^^^^

.. code-block:: c++

   std::array<TwoIndex<double>, 45> d2;
   ecpint.compute_shell_pair_second_derivative(U, shellA, shellB, d2);

gives 45 matrices, in the order

.. code-block:: text

   AxAx, AxAy, AxAz, AyAy, AyAz, AzAz,
   AxBx, AxBy, AxBz, AyBx, AyBy, AyBz, AzBx, AzBy, AzBz,
   AxCx, AxCy, AxCz, AyCx, AyCy, AyCz, AzCx, AzCy, AzCz,
   BxBx, BxBy, BxBz, ByBy, ByBz, BzBz,
   BxCx, BxCy, BxCz, ByCx, ByCy, ByCz, BzCx, BzCy, BzCz,
   CxCx, CxCy, CxCz, CyCy, CyCz, CzCz

These too are additive, but when two or three of the centres coincide, the symmetric blocks need
care. ``pack_second_derivatives`` in ``assembly.hpp`` does this, and the assembly loops there
show how the high-level API places blocks in molecular matrices.

The magnetic field derivatives over one shell pair are given by
``compute_shell_pair_derivative_Bfield``.

Examples of the low-level API are in the tests: ``tests/analytic/analytic.cpp``
(``compute_shell_pair``), ``tests/units/derivatives/derivatives.cpp`` (the derivatives), and
``tests/suite/suite.cpp`` (all kernels, including from several threads).

.. toctree::
   :hidden:
