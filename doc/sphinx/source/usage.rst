.. libecpint API file

.. _`sec:usage`:

=====
Usage
=====  

There are two main ways to use the libecpint library. In the high-level API, you pass details of the system (basis set, coordinates, ECPs) to libecpint, and it handles the computation of all the integrals and/or derivatives automatically, returning arrays of the (Cartesian) integrals. In the low-level API, you control the calculation of the integrals yourself, calling the relevant routines as and when you need them. We envisage that for most purposes, the high-level API will be more appropriate, and is easier to use. 

Important note about ECP definitions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

ECP powers (``u_ns``) are assumed to be including the ``r^2`` from the spherical Jacobian, as is the convention followed with the Stuttgart-Dresden ECPs. In the code, we subtract two from this (see the ECP constructor) as integration is done in Cartesian coordinates. In general, the inputted ``n`` should only ever be 0, 1, or (most commonly) 2. If your ECP definitions are instead -2 or -1, or predominantly 0, that suggests you are following the other convention of not including the Jacobian, so you should add 2. Yes, this is annoying; if in doubt check the built-in ECP library for examples.

High Level API
==============

A complete example of using this API is ``example/example.cpp`` (HI, with the iodine ECP from the built-in library: integrals and first derivatives), with build instructions in ``example/README.md``. ``tests/units/api/api.cpp`` also shows ``update_gaussian_basis_coords``, ``update_ecp_basis_coords`` and the ``get_*`` functions. 

The ECPIntegrator Object
^^^^^^^^^^^^^^^^^^^^^^^^

The first step in using this API is to include the libecpint header and create an ``ECPIntegrator`` object:

.. code-block:: c++

	include <libecpint.hpp>
	ECPIntegrator factory;
	
This object will form the main interface to all of the subsequent routines, and is described in detail in the Library API.  

Initialisation
^^^^^^^^^^^^^^

There are three steps to initialising the ``ECPIntegrator`` before it can be used to calculate integrals. These are:

	1) specifying the Gaussian basis set
	2) specifying the ECP basis
	3) calling the ``init`` routine
	
These steps are performed as follows:

.. code-block:: c++
	
	factory.set_gaussian_basis(N_shells, g_coords, g_exps, g_coefs, g_ams, g_lengths);
	factory.set_ecp_basis(N_ecps, u_coords, u_exps, u_coefs, u_ams, u_ns, u_lengths);
	factory.init(deriv_order, tol, radial_tol, smallGrid, bigGrid);

where ``N_shells`` and ``N_ecps`` are the numbers of shells in the Gaussian basis, and the number of ECP centers, respectively, while ``deriv_order`` is the maximum derivatives needed (0, 1, or 2). The ``init`` parameters after ``deriv_order`` are optional and control screening and quadrature (see Settings_ below). The rest of the parameters are:

 - the Cartesian coordinates *in Bohr* (``g_coords, u_coords``);
 - the exponents, coefficients, and angular momenta (and powers, ``u_ns``, for the ECPs) comprising the basis sets
 - the number of exponents per shell (``g_lengths, u_lengths``) 
 
 ``example/example.cpp`` shows how to specify these for HI. 
 
**NOTE**: The atom order given in ``set_gaussian_basis`` fixes the atom order for all the derivatives, as will be described later. 

The ECP library
^^^^^^^^^^^^^^^

Alternatively, step 2 can be replaced by reading the ECPs from the built-in library provided with libecpint. This can be found in ``share/libecpint``. To do this, you call:

.. code-block:: c++
	
	factory.set_ecp_basis_from_library(N_ecps, u_coords, u_charges, u_names, share_dir);
	
The new parameters are:
	 
	 - ``u_charges`` a list of atomic numbers for the ECPs, corresponding to the centers in u_coords;
	 - ``u_names`` the ECP names for each ECP, e.g. ``ecp10mdf``;
	 - ``share_dir`` the absolute path to the share/libecpint directory, which must be passed by you.
	 
The ECPs in the library, and the atoms they are available for, are given below. The Stuttgart/Cologne energy-consistent pseudopotentials are named ``ecpnxy`` after their ECPnXY label: n is the number of core electrons, X is S (single-electron fit) or M (multi-electron fit), and Y is HF (non-relativistic), WB (quasi-relativistic) or DF (fully relativistic). They come from the Stuttgart/Cologne library (http://www.tc.uni-koeln.de/PP), except where the ECPs of the cc-pVnZ-PP basis sets on the Basis Set Exchange differ (currently only Y in ``ecp28mdf``), in which case the Basis Set Exchange version is used. The remaining families come from the Basis Set Exchange. The Karlsruhe ``lcecp-n`` sets are the Stuttgart 4f-in-core ECPnMWB potentials grouped by valence, except that Nd in ``lcecp-1`` has different parameters from ``ecp49mwb``. Each file records the source and references for every atom.

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
	  
** TO ADD AN ECP TO THE LIBRARY ** 

	1) Put the ECP in MOLPRO format in ``share/libecpint/raw`` as ``NAME.ecp``, where ``NAME`` is the name of the ECP; make sure that any exponents are with ``E`` (C-convention) not ``D`` (Fortran-convention). The MOLPRO export from the `Basis Set Exchange <https://www.basissetexchange.org>`_ can be used as it is; any orbital basis and spin-orbit terms in the file are ignored.
	2) Make the top line of ``NAME.ecp`` be the ``NAME``.

	The ECP will now be available for use by libecpint. Please consider creating a pull request so that everyone can benefit from the addition!
	
*Note* that the ``n`` value for each ECP primitive should typically be 0, 1, or 2 (for the Stuttgart-Dresden ECPs, for example, it is `always` 2). Some input formats follow a convention of subtracting or adding 2 to this.  
	

Computing integrals
^^^^^^^^^^^^^^^^^^^

Computing integrals over all shell pairs is then very simple:

.. code-block:: c++

	factory.compute_integrals()
	
To retrieve these you then create a shared pointer to a vector:
	
.. code-block:: c++
	
	std::shared_ptr<std::vector<double>> integrals = factory.get_integrals();
	double I00 = (*ints)[0]; // example for accessing element (0, 0)

These are stored in row-order, and are in a *Cartesian Gaussian basis*. Typically these would be converted to a spherical harmonic Gaussian basis (we might add the ability to do this later). We follow canonical Cartesian order, so for a d-type function this would be ``xx, xy, xz, yy, yz, zz``, and the order of the shells is the same as when you called ``set_gaussian_basis``. The total number of Cartesian gaussians is stored in ``factory.ncart``; you can access the ij-th integral as

.. code-block:: c++

	double Iij = (*ints)[i*factory.ncart + j]
	

Derivatives contracted with a density matrix
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Most uses of the derivatives, such as energy gradients and the integral part of the Hessian, only need them contracted with a density matrix. This is the recommended way to get them, as the derivative matrices are never stored:

.. code-block:: c++

	// D: ncart x ncart, row-major, in the same order as the integrals
	std::vector<double> gradient = factory.compute_gradient(D.data());  // needs deriv_order >= 1
	std::vector<double> hessian = factory.compute_hessian(D.data());    // needs deriv_order >= 2

``gradient[3*n + i]`` is the sum over k, l of ``D[k*ncart + l]`` times the derivative of integral (k, l) with respect to coordinate i (x, y, z) of atom n, with atoms in the order given in ``set_gaussian_basis``. ``hessian`` is the full, symmetric ``3*natoms`` by ``3*natoms`` matrix of the corresponding second derivatives, in row-major order. ``D`` need not be symmetric. The memory needed is about that of the integrals themselves, whereas the full derivative matrices below need ``3*natoms`` (first) or ``3*natoms*(3*natoms+1)/2`` (second) times as much.

First derivatives
^^^^^^^^^^^^^^^^^

If the derivative matrices themselves are needed (e.g. for coupled-perturbed equations), they can be stored in full. *Note* that the memory needed grows as ``natoms * ncart^2``: about 7 GB for a 77-atom molecule with 1900 Cartesian functions.

First derivatives are similarly calculated by calling

.. code-block:: c++

	factory.compute_first_derivs()
	
*Note* that this will only work if ``init`` was called with ``deriv_order > 0``. This will return an array of ``3*factory.natoms`` shared pointers to the integral derivatives with respect to each coordinate. The order is x, y, z, and the order of atoms matches that specified in ``set_gaussian_basis``. For example, to get the array of integral derivatives with respect to the y-coordinate of the n-th atom, you would do:

.. code-block:: c++

	std::vector<std::shared_ptr<std::vector<double>>> first_derivs = factory.get_first_derivs();
	I_dy_atom_n_00 = (*first_derivs[3*n+1])[0];

The order of the elements in each array is identical to that from ``compute_integrals``.

Second derivatives
^^^^^^^^^^^^^^^^^^

*Note* that storing all second derivative matrices needs memory growing as ``natoms^2 * ncart^2``, which is prohibitive beyond small systems (about 800 GB for the example above); use ``compute_hessian`` where possible.

As for first derivatives, second derivatives are computed as

.. code-block:: c++

	factory.compute_second_derivs()
	
and are provided as a vector of shared pointers to arrays. The order of these derivatives is somewhat more complicated though, and takes full advantage of symmetry. If the atoms are A, B, C, ... as specified in ``set_gaussian_basis``, then they are blocked as follows: 

.. code-block:: bash
	
	AA, AB, AC, ..., BB, BC, ..., CC, ...

Within each block the order is ``xx, xy, xz, yy, yz, z`` on the diagonal (e.g. AA, BB, CC, ...) and ``xx, xy, xz, yx, yy, yz, zx, zy, zz`` on the off-diagonal (e.g. AB, AC, BC, ...). There is a helper macro for this, ``H_MACRO``, defined in ``api.hpp``. So for example to get the derivative with respect to Ay (atom index 0) and Cx (atom index = 2) in a system with ``N`` atoms, you would do

.. code-block:: c++

	int deriv_index = H_START(0, 2, N) + 3; // yx is index 3 for off-diagonal blocks
	std::shared_ptr<std::vector<double>> h_Ay_Cx = factory.get_second_derivs()[deriv_index];
	
Hopefully this doesn't give you too much of a headache working out.

Magnetic field derivatives
^^^^^^^^^^^^^^^^^^^^^^^^^^

For a basis of gauge-including atomic orbitals (London orbitals), the derivatives of the ECP integrals with respect to a uniform magnetic field (e.g. for NMR shieldings) are computed by

.. code-block:: c++

	factory.compute_bfield_derivs();  // needs deriv_order >= 1
	std::vector<std::shared_ptr<std::vector<double>>> bfield_derivs = factory.get_bfield_derivs();

giving three matrices, for x, y and z, in the order of the integrals. Each is real and antisymmetric, the derivative being ``i/2`` times it; the ECP's projectors are gauge-transformed to the ECP centre. See ``ECPIntegral::compute_shell_pair_derivative_Bfield`` in ``ecpint.hpp`` for the precise definition.

Updating coordinates
^^^^^^^^^^^^^^^^^^^^

To update the coordinates for the basis and ECPs (for example, after a step in a geometry optimisation), simply pass the new coordinates *in the same order they were given when initialised*. This is done as:

.. code-block:: c++

	 factory.update_gaussian_basis_coords(N_shells, g_coords);
	 factory.update_ecp_basis_coords(N_ecps, u_coords);
	 
You then call the compute routines when you need the new integrals and/or derivatives. 

**NOTE** you will need to re-get the pointers using the get routines every time you recompute integrals/derivatives.

.. _Settings:

Settings
^^^^^^^^

The ``init`` method accepts optional parameters that control screening and quadrature accuracy:

.. code-block:: c++

	factory.init(deriv_order, tol, radial_tol, smallGrid, bigGrid);

+-----------------+-----------+----------------------------------------------------------+
| Parameter       | Default   |  Description                                             |
+=================+===========+==========================================================+
| deriv_order     | 0         | Maximum derivative order to compute (0, 1, or 2).       |
+-----------------+-----------+----------------------------------------------------------+
| tol             | 1e-12     | Screening tolerance for skipping shell-pair              |
|                 |           | contributions. Also determines ``shell_tolerance``       |
|                 |           | for skipping entire shell contributions based on         |
|                 |           | the maximum angular momentum in the ECP basis.           |
+-----------------+-----------+----------------------------------------------------------+
| radial_tol      | 1e-15     | Convergence tolerance for the radial quadrature.         |
+-----------------+-----------+----------------------------------------------------------+
| smallGrid       | 256       | Maximum number of quadrature points for the small        |
|                 |           | radial integration grid.                                 |
+-----------------+-----------+----------------------------------------------------------+
| bigGrid         | 1024      | Maximum number of quadrature points for the large        |
|                 |           | radial integration grid.                                 |
+-----------------+-----------+----------------------------------------------------------+

The ``tolerance`` and ``shell_tolerance`` fields on ``ECPIntegrator`` are set by ``init`` and are used during integral, first derivative, and second derivative evaluation to skip negligible contributions. For most applications the defaults are appropriate; loosening ``tol`` (e.g. to 1e-10) can speed up calculations on larger molecules at the cost of some accuracy.


Low Level API
=============

Examples of using this API can be found in the tests: ``compute_shell_pair`` in ``tests/analytic/analytic.cpp``, ``compute_shell_pair_derivative`` in ``tests/radial/radial.cpp`` (``Type2.HighOrderDerivatives``), and the derivative recurrences in ``tests/units/derivatives/derivatives.cpp``. 

The ECPIntegral Object
^^^^^^^^^^^^^^^^^^^^^^

To be able to calculate integrals and derivatives at the shell-pair level, you need to create an ``ECPIntegral`` object instead. This is done as follows:

.. code-block:: c++

	#include <libecpint/ecpint.hpp>
	ECPIntegral ecpint(maxLB, maxLU, deriv_order);

where ``maxLB`` is the maximum angular momentum in the Gaussian basis, ``maxLU`` is the maximum angular momentum in the ECP basis, and ``deriv_order`` is the maximum order of derivative needed (defaults to 0). 

Making shells and ECPs
^^^^^^^^^^^^^^^^^^^^^^

The compute functions in this API require ``ECP`` and ``GaussianShell`` objects representing the ECP and Gaussian basis functions, respectively. These are populated for the ECP as follows:

.. code-block:: c++

	#include <libecpint/ecp.hpp>
	double C[3] = {Cx, Cy, Cz};
	ECP newU(C);
	newU.addPrimitive(n, l, x, c);
	// addPrimitive for each primitive in ECP

where ``C`` is the coordinates of the center of the ECP (in Bohr), and ``n, l, x, c`` are the power, angular momentum, exponent, and coefficient of each primitive in that ECP. For the Gaussian basis instead:

.. code-block:: c++	
	
	#include <libecpint/gshell.hpp>
	double A[3] = {Ax, Ay, Az};
	GaussianShell shellA(A, l); 
	shellA.addPrim(x, c);
	// addPrim for each primitive in this shell

where ``A`` is the coordinates of the center of the shell (in Bohr), and ``l, x, c`` are as above. You need to either create these objects for every shell and ECP as they are needed, or store them, so that they can be passed to the compute routines below. 
 
Computing integrals
^^^^^^^^^^^^^^^^^^^

Computing integrals over a shell pair is then fairly simple. You call

.. code-block:: c++

	#include <libecpint/mathutil.hpp>
	TwoIndex<double> results;
	ecpint.compute_shell_pair(newU, shellA, shellB, results);
	
This will place the integrals in a matrix, ``results``, with dimensions ``(ncartA, ncartB)``. The order in each is the canonical Cartesian order as described earlier for the high-level API. Elements of the results matrix can be accessed in two ways:

.. code-block:: c++
	 
	 double Iij = results(i, j);
	 // or
	 double Iij = results.data[i*ncartB + j];

First derivatives
^^^^^^^^^^^^^^^^^

Similarly, first derivatives are calculated using

.. code-block:: c++
	
	std::vector<TwoIndex<double>> results;
	ecpint.compute_shell_pair_derivative(newU, shellA, shellB, results);
	
Now ``results`` is an array of 9 matrices (the matrices are ordered the same as for ``compute_shell_pair`` above). These are derivative matrices with respect to ``Ax, Ay, Az, Bx, By, Bz, Cx, Cy, Cz``, where ``A, B, C`` are the centers of ``shellA``, ``shellB``, and ``newU`` respectively.

**NOTE** These derivatives are designed to be additive. Thus, if ``A==B``, then the derivative for the Ax coordinate will be the sum ``Ax+Bx``, etc. 

Second derivatives
^^^^^^^^^^^^^^^^^^

The second derivatives are calculated in the same simple manner, but their ordering is, as with the high-level API, much more complicated. 

.. code-block:: c++

	std::vector<TwoIndex<double>> results;
	ecpint.compute_shell_pair_second_derivative(newU, shellA, shellB, results);
	
Now ``results`` contains **45** different derivative matrices. These are, in order:

.. code-block:: bash

	 AxAx, AxAy, AxAz, AyAy, AyAz, AzAz,
	 AxBx, AxBy, AxBz, AyBx, AyBy, AyBz, AzBx, AzBy, AzBz,
	 AxCx, AxCy, AxCz, AyCx, AyCy, AyCz, AzCx, AzCy, AzCz,
	 BxBx, BxBy, BxBz, ByBy, ByBz, BzBz,
	 BxCx, BxCy, BxCz, ByCx, ByCy, ByCz, BzCx, BzCy, BzCz,
	 CxCx, CxCy, CxCz, CyCy, CyCz, CzCz

where ``A, B, C`` are as described earlier. These are again additive but THERE IS NOW SYMMETRY TO CONSIDER. I *strongly* recommend that you look at the code in ``api.cpp`` for the construction of the full Hessian, to see how this symmetry has to be taken account of when any of the three centers are the same, as it is too complicated to describe here. 
 
Settings
^^^^^^^^

The ``ECPIntegral`` constructor accepts the following parameters:

.. code-block:: c++

	ECPIntegral ecpint(maxLB, maxLU, deriv_order, thresh, smallGrid, bigGrid);

+-----------------+-----------+----------------------------------------------------------+
| Parameter       | Default   |  Description                                             |
+=================+===========+==========================================================+
| maxLB           | (required)| Maximum angular momentum in the Gaussian basis.          |
+-----------------+-----------+----------------------------------------------------------+
| maxLU           | (required)| Maximum angular momentum in the ECP basis.               |
+-----------------+-----------+----------------------------------------------------------+
| deriv_order     | 0         | Maximum derivative order to compute (0, 1, or 2).       |
+-----------------+-----------+----------------------------------------------------------+
| thresh          | 1e-15     | Convergence tolerance for the radial quadrature.         |
+-----------------+-----------+----------------------------------------------------------+
| smallGrid       | 256       | Maximum number of quadrature points for the small        |
|                 |           | radial integration grid.                                 |
+-----------------+-----------+----------------------------------------------------------+
| bigGrid         | 1024      | Maximum number of quadrature points for the large        |
|                 |           | radial integration grid.                                 |
+-----------------+-----------+----------------------------------------------------------+

.. toctree::
   :hidden:
