.. libecpint API file

.. _`sec:usage`:

=====
Usage
=====  

There are two main ways to use the libecpint library. In the high-level API, you pass details of the system (basis set, coordinates, ECPs) to libecpint, and it handles the computation of all the integrals and/or derivatives automatically, returning arrays of the (Cartesian) integrals. In the low-level API, you control the calculation of the integrals yourself, calling the relevant routines as and when you need them. We envisage that for most purposes, the high-level API will be more appropriate, and is easier to use. 

High Level API
==============

Examples of using this API can be found in ``tests/lib/api_test1`` and ``tests/lib/api_test2``. 

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
	factory.init(deriv_order);
	
where ``N_shells`` and ``N_ecps`` are the numbers of shells in the Gaussian basis, and the number of ECP centers, respectively, while ``deriv_order`` is the maximum derivatives needed (0, 1, or 2). The rest of the parameters are:

 - the Cartesian coordinates *in Bohr* (``g_coords, u_coords``);
 - the exponents, coefficients, and angular momenta (and powers, ``u_ns``, for the ECPs) comprising the basis sets
 - the number of exponents per shell (``g_lengths, u_lengths``) 
 
 ``tests/lib/api_test1`` shows how to specify these for HBr in the AVDZ/AVDZ-PP basis. 
 
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
	 
The currently available ECPs in the library (more being added soon), and the atoms they are available for, are given below: 

.. list-table::
	:widths: 25 75
	:header-rows: 1
	
	* - Name
	  - Atoms
	* - ECP10MDF
	  - K -- Kr (Z = 19 -- 36)
	* - ECP28MDF
	  - Rb -- Xe (Z = 37 -- 54)
	* - ECP46MDF
	  - Cs, Ba (Z = 55, 56)
	* - ECP60MDF
	  - Hf -- Rn (Z = 72 -- 86), Ac -- U (Z = 89 -- 92)
	* - ECP78MDF
	  - Fr, Ra (Z = 87, 88)
	* - LANL2DZ
	  - Na -- La (Z = 11 -- 57), Hf -- Bi (Z = 72 -- 83), U -- Pu (Z = 92 -- 94)
	  
** TO ADD AN ECP TO THE LIBRARY ** 

	1) Put the ECP in MOLPRO format in ``share/libecpint/raw`` as ``NAME.ecp``, where ``NAME`` is the name of the ECP; make sure that any exponents are with ``E`` (C-convention) not ``D`` (Fortran-convention).  
	2) Make the top line of ``NAME.ecp`` be the ``NAME``.
	3) In ``share/libecpint`` run ``python3 parseecp.py NAME`` (Python >=3.6 required, with lxml module). This will create ``NAME.xml`` in the ``share/libecpint/xml`` folder, and this ECP will now be available for use by libecpint. Please consider creating a pull request so that everyone can benefit from the addition!
	
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
	

First derivatives
^^^^^^^^^^^^^^^^^

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
 
Updating coordinates
^^^^^^^^^^^^^^^^^^^^

To update the coordinates for the basis and ECPs (for example, after a step in a geometry optimisation), simply pass the new coordinates *in the same order they were given when initialised*. This is done as:

.. code-block:: c++

	 factory.update_gaussian_basis_coords(N_shells, g_coords);
	 factory.update_ecp_basis_coords(N_ecps, u_coords);
	 
You then call the compute routines when you need the new integrals and/or derivatives. 

**NOTE** you will need to re-get the pointers using the get routines every time you recompute integrals/derivatives.

Settings
^^^^^^^^

**TODO** detail optional settings that can be passed to ECPIntegrator


Low Level API
=============

Examples of using this API can be found in ``tests/lib/[name]_test[number]`` where [name] is int (integrals), deriv (first derivatives), or hess (second derivatives), and [number] is 1 or 2. 

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

**TODO** detail optional settings that can be passed to ECPIntegral

.. toctree::
   :hidden:
