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

#ifndef MATHUTILHEADERDEF
#define MATHUTILHEADERDEF

/**
  * \file  mathutil.hpp
  * \brief Mathematical constants, special function tables, and utility functions
  */

#include <vector>
#include "multiarr.hpp"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define MAX_FAC 100
#define MAX_DFAC 200

namespace libecpint {
	
	const double ROOT_PI = 1.772453850905516;

	extern double FAC[MAX_FAC];		///< Array of factorials
	extern double DFAC[MAX_DFAC]; 	///< Array of double factorials
	
	/**
	  * Gamma function tabulation, where GAMMA[i] = Gamma((i+1)/2)
	  * e.g. GAMMA[0] = Gamma(1/2) = sqrt(Pi), GAMMA[1] = 0! = 1, GAMMA[2] = Gamma(3/2) = sqrt(Pi) / 2, etc.
	  */ 
	const double GAMMA[30] = {
		1.7724538509055,
		1.0,
		0.88622692545275,
		1.0,
		1.3293403881791,
		2.0,
		3.3233509704478,
		6.0,
		11.631728396567,
		24.0,
		52.342777784553,
		120.0,
		287.88527781504,
		720.0,
		1871.2543057978,
		5040.0,
		14034.407293483,
		40320.0,
		1.1929246199461e5,
		3.62880e5,
		1.1332783889488e6,
		3.628800e6,
		1.1899423083962e7,
		3.9916800e7,
		1.3684336546556e8, 
		4.79001600e8, 
		1.7105420683196e9, 
		6.227020800e9, 
		2.3092317922314e10,
		8.7178291200e10
	};
	
	/** 
	* Calculates real spherical harmonics S_lm(theta, phi) for all l, m up to lmax
	* @param lmax - the maximum angular momentum needed
	* @param x - cos(theta), where theta is the polar angle in spherical coordinates
	* @param phi - the azimuth angle in spherical coordinates
	* @return a matrix S(l, l+m) of the spherical harmonic values
	*/
	TwoIndex<double> realSphericalHarmonics(int lmax, double x, double phi);  
	
	/**
	 * Initialises the global factorial and double factorial arrays
	 */
	void initFactorials(); 
}

#endif
