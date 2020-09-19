/* 
 *      Copyright (c) 2020 Robert Shaw
 *		This file is a part of Libecpint.
 *
 *      Permission is hereby granted, free of charge, to any person obtaining
 *      a copy of this software and associated documentation files (the
 *      "Software"), to deal in the Software without restriction, including
 *      without limitation the rights to use, copy, modify, merge, publish,
 *      distribute, sublicense, and/or sell copies of the Software, and to
 *      permit persons to whom the Software is furnished to do so, subject to
 *      the following conditions:
 *
 *      The above copyright notice and this permission notice shall be
 *      included in all copies or substantial portions of the Software.
 *
 *      THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 *      EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 *      MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 *      NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
 *      LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 *      OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 *      WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#ifndef ECPINT_HEAD
#define ECPINT_HEAD

#include <vector>
#include <array>
#include "multiarr.hpp"
#include "gaussquad.hpp"
#include "ecp.hpp"
#include "bessel.hpp"
#include "radial.hpp"
#include "angular.hpp"
#include "gshell.hpp"

#include "config.hpp"

namespace libecpint {

#define N_INDEX(l, m) ((l+m)*(l+m+1)/2 + m)

	/** 
	* \class ECPIntegral
	* \brief Calculates ECP integrals. 
	* 
	* Given an ECP basis, and orbital bases, this will calculate the ECP integrals over all ECP centers. 
	* 
	* REFERENCES: 
	* (Shaw2017) R. A. Shaw, J. G. Hill, J. Chem. Phys. 147 (2017), 074108
	* (Flores06) R. Flores-Moreno et al., J. Comput. Chem. 27 (2006), 1009
	* (MM81) L. E. McMurchie and E. R. Davidson, J. Comp. Phys. 44 (1981), 289 - 301
	*/
	class ECPIntegral
	{
	private:
		RadialIntegral radInts; ///< The interface to the radial integral calculation
		AngularIntegral angInts; ///< The angular integrals, which can be reused over all ECP centers
	
		/// Worker functions for calculating binomial expansion coefficients
		double calcC(int a, int m, double A) const;
		
		/// Array of function pointers to generated integral evaluators in qgen
		static void(*QGEN[LIBECPINT_MAX_L+1][LIBECPINT_MAX_L+1][LIBECPINT_MAX_L+1])(ECP&, GaussianShell&, GaussianShell&,
		 FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double,
		 RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);

	public:
		
		/** 
		  * Constructs the coefficients in the binomial expansion (see REF. Shaw2017)
		  * @param C - reference to a FiveIndex array to store the results in
		  * @param L - maximum angular momentum to go up to in expansion
		  * @param A - xyz coordinates for the center to calculate over
		  */ 
		void makeC(FiveIndex<double> &C, int L, double *A);
		
		/**
		  * Creates an ECP integrator, initialising the radial and angular parts for subsequent calculations.
		  * @param maxLB - the maximum angular momentum in the orbital basis
		  * @param maxLU - the maximum angular momentum in the ECP basis
		  * @param deriv - the maximum order of derivative to be calculated (TODO: derivs currently being implemented)
		  */ 
		ECPIntegral(int maxLB, int maxLU, int deriv=0);
	
		/**
		  * Calculates the type 1 integrals for the given ECP center over the given shell pair, using quadrature
		  * @param U - reference to ECP
		  * @param shellA - the first basis shell (rows in values)
		  * @param shellB - the second basis shell (cols in values)
		  * @param data - wrapper for data about shell pair
		  * @param CA - binomial expansion coefficients for shellA, made with makeC
		  * @param CB - binomial expansion coefficients for shellB, made with makeC
		  * @param values - array in which results are returned
		  */
		void type1(ECP& U, GaussianShell &shellA, GaussianShell &shellB, ShellPairData &data, FiveIndex<double> &CA, FiveIndex<double> &CB, TwoIndex<double> &values);
		
		/**
		  * Calculates the type 2 integrals for the given ECP center over the given shell pair
		  * @param l - angular momentum shell of ECP to calculate over
		  * @param U - reference to ECP
		  * @param shellA - the first basis shell (rows in values)
		  * @param shellB - the second basis shell (cols in values)
		  * @param data - wrapper for data about shell pair
		  * @param CA - binomial expansion coefficients for shellA, made with makeC
		  * @param CB - binomial expansion coefficients for shellB, made with makeC
		  * @param values - array in which results are returned
		  */
		void type2(int l, ECP& U, GaussianShell &shellA, GaussianShell &shellB, ShellPairData &data, FiveIndex<double> &CA, FiveIndex<double> &CB, ThreeIndex<double> &values);
	
		/**
		  * Computes the overall ECP integrals over the given ECP center and shell pair. This is currently the main interface.
		  * Results are returned with rows corresponding to shellA and cols to shellB, with the Cartesian functions in alpha order e.g.
		  * {xxx, xxy, xxz, xyy, xyz, xzz, yyy, yyz, yzz, zzz} l = 3
		  * 
		  * @param U - reference to the ECP to calculate the integral over
		  * @param shellA - the first basis shell (rows in values) 
		  * @param shellB - the second basis shell (cols in values)
		  * @param values - reference to TwoIndex array where the results will be stored
		  */ 
		void compute_shell_pair(ECP &U, GaussianShell &shellA, GaussianShell &shellB, TwoIndex<double> &values, int shiftA = 0, int shiftB = 0);
		
		void left_shell_derivative(ECP &U, GaussianShell &shellA, GaussianShell &shellB, std::array<TwoIndex<double>, 3> &results); 
		void compute_shell_pair_derivative(ECP &U, GaussianShell &shellA, GaussianShell &shellB, std::array<TwoIndex<double>, 9> &results);
		void left_shell_second_derivative(ECP &U, GaussianShell &shellA, GaussianShell &shellB, std::array<TwoIndex<double>, 6> &results); 
		void mixed_second_derivative(ECP &U, GaussianShell &shellA, GaussianShell &shellB, std::array<TwoIndex<double>, 9> &results); 
		void compute_shell_pair_second_derivative(ECP &U, GaussianShell &shellA, GaussianShell &shellB, std::array<TwoIndex<double>, 45> &results);
	};

}
#endif
