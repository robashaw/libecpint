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

#ifndef GENERATE_HEAD
#define GENERATE_HEAD

#include <iostream>
#include <fstream>
#include <string>
#include "angular.hpp"
#include "mathutil.hpp"
#include "multiarr.hpp"
#include "config.hpp"
#include <cmath>

namespace libecpint {
	
	const int maxL = LIBECPINT_MAX_L;
	const int maxUnrol = LIBECPINT_MAX_UNROL; 
	const int maxN = 15; 
	const int CAX1 = maxN + 1;
	const int CAX2 = CAX1 * (maxL + 1);
	const int CAX3 = CAX2 * (maxL + 1);

	/** 
	  * \struct SumTerm
	  * \brief  Stores details of terms in ECP integral expansion
	  */
	struct SumTerm {

		Pair SA; ///< (l, m) for spherical harmonic on shellA
		Pair SB; ///< (l, m) for spherical harmonic on shellB
		Triple radial; ///< (N, l1, l2) radial integral required
		Quintuple CA; ///< (0, cartesian index, x, y, z) for binomial coefficient on shellA
		Quintuple CB; ///< (0, cartesian index, x, y, z) for binomial coefficient on shellB
		double ang;   ///< Value of product of angular integrals
		int mu; ///< Value of mu, where the ECP ang. momentum is lam, and mu can range from -lam .. lam
		int na; ///< Index of cartesian function on shellA, in alpha order
		int nb; ///< Index of cartesian function on shellB, in alpha order
	
		/**
		  * orders SumTerm's by mu, then by radial indices, then by angular integral value
		  * @param s - term to compare with
		  * @return true if this term is less than the given term
		  */ 
		bool operator<(const SumTerm& s) const {
			bool result = (mu < s.mu);
			if (mu == s.mu) {
				result = radial < s.radial; 
				if (radial == s.radial)
					result = ang < s.ang; 
			}
			return result; 
		}
	   
	    /// @return true if this term is less than or equal to the given term
		bool operator<=(const SumTerm& s) const {
			return (*this == s) || (*this < s); 
		}
	
		/// @return true if the mu and radial indices are equal
		bool operator==(const SumTerm& s) const {
			return (mu == s.mu) && (radial == s.radial); 
		}
	
		/// @return the compressed index of CA
		int ca_index() const {
			return std::get<1>(CA) + std::get<2>(CA)*CAX1 
				+ std::get<3>(CA)*CAX2 + std::get<4>(CA)*CAX3; 
		}
	
		/// @return the compressed index of CB
		int cb_index() const {
			return std::get<1>(CB) + std::get<2>(CB)*CAX1 
				+ std::get<3>(CB)*CAX2 + std::get<4>(CB)*CAX3; 
		}
	
		/// Converts term to string with compressed indices
		std::string to_string(bool full = true) {
			std::stringstream ss;
		
			if (full) 
				ss << "\tvalue" << mu << " += " << ang << " * CA[" << ca_index()
					<< "] * CB[" << cb_index() << "] * radials" << radial << " * SA"
						<< SA << " * SB" << SB << ";"; 
			else
				ss << "\ttmp += " << "CA[" << ca_index() << "] * CB[" << cb_index() << "] * SA"
					<< SA << " * SB" << SB << ";"; 
		 
			return ss.str(); 
		}
	
		/**
		  * Compares two SumTerm objects
		  * @param s - the SumTerm to compare with
		  * @return a tuple of equalities (0 = false, 1 = true) in the order {mu, radial, SA, SB, ang, CA, CB}
		  */ 
		Heptuple compare(const SumTerm& s) const {
			int f1 = mu == s.mu ? 1 : 0;
			int f2 = radial == s.radial ? 1 : 0;
			int f3 = SA == s.SA ? 1 : 0;
			int f4 = SB == s.SB ? 1 : 0; 
			int f5 = fabs(ang - s.ang) < 1e-10 ? 1 : 0;
			int f6 = CA == s.CA ? 1 : 0;
			int f7 = CB == s.CB ? 1 : 0; 
		
			return {f1, f2, f3, f4, f5, f6, f7};  
		}
	
		/// Prints out a SumTerm without compressing the indices - currently preferred 
		friend std::ostream& operator<<(std::ostream& os, const SumTerm& s); 

	};

	std::ostream& operator<<(std::ostream& os, const SumTerm& s) {
		os << "values(" << s.na << ", " << s.nb << ", " << s.mu << ") += "
			<< s.ang << " * CA" << s.CA << " * CB" << s.CB
				<< " * radials" << s.radial << " * SA" << s.SA << " * SB" << s.SB << ";";
		return os;
	}

} // end namespace

#endif
