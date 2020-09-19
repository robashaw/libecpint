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

#include "api.hpp"
#include <iostream>
#include <cmath>

namespace libecpint {
	
	void ECPIntegrator::set_gaussian_basis (int nshells, double* coords, double* exponents, double* coefs, int* ams, int* shell_lengths) {
		int ctr = 0;
		for (int i = 0; i < nshells; i++) {
			ncart += (ams[i]+1)*(ams[i]+2)/2;
			std::array<double, 3> cvec = {coords[3*i], coords[3*i+1], coords[3*i+2]};
			GaussianShell newShell(cvec, ams[i]);
			if (ams[i] > maxLB) maxLB = ams[i];
			for (int n = 0; n < shell_lengths[i]; n++) {
				newShell.addPrim(exponents[ctr], coefs[ctr]);
				ctr++;
			}
			shells.push_back(newShell);
		} 
		basis_is_set = true;
	}
	
	void ECPIntegrator::set_ecp_basis(int necps, double* coords, double* exponents, double* coefs, int* ams, int* ns, int* shell_lengths) {
		int ctr = 0;
		for (int i = 0; i < necps; i++) {
			ECP newU(&coords[3*i]);
			for (int n = 0; n < shell_lengths[i]; n++) {
				newU.addPrimitive(ns[ctr], ams[ctr], exponents[ctr], coefs[ctr]); 
				ctr++;
			}
			newU.sort();
			ecps.addECP(newU, 0);
		}
		ecp_is_set = true;
	}
	
	void ECPIntegrator::update_gaussian_basis_coords(int nshells, double* coords) {
		assert(nshells = shells.size());
		
		for (int i = 0; i < nshells; i++){
			shells[i].localCenter[0] = coords[3*i];
			shells[i].localCenter[1] = coords[3*i+1];
			shells[i].localCenter[2] = coords[3*i+2];
		}
	}
	
	void ECPIntegrator::update_ecp_basis_coords(int necps, double* coords) {
		assert(necps = ecps.getN());
		
		for (int i = 0; i < necps; i++) 
			ecps.getECP(i).setPos(coords[3*i], coords[3*i+1], coords[3*i+2]);
	}
	
	void ECPIntegrator::init(int deriv_) {
		assert(ecp_is_set);
		assert(basis_is_set);
		deriv = std::max(0, std::min(2, deriv_));
		ecpint = std::make_shared<ECPIntegral>(maxLB, ecps.getMaxL(), deriv);
	}
	
	void ECPIntegrator::compute() {
		// initialise all to zero
		integrals.assign(ncart, ncart, 0.0);
		if (deriv > 0) {
			 for (auto& r : first_derivs) r.assign(ncart, ncart, 0.0);
			 if (deriv > 1) 
				 for (auto& r : second_derivs) r.assign(ncart, ncart, 0.0);
		}
		
		// loop over shells
		TwoIndex<double> tempValues;
		std::array<TwoIndex<double>, 9> tempGrads;
		std::array<TwoIndex<double>, 45> tempHess;
		int nshells = shells.size();
		
		int n1 = 0;
		for(auto s1=0; s1<nshells; ++s1) {
			GaussianShell& shellA = shells[s1];
			int ncartA = shellA.ncartesian();
			
			int n2 = 0;
			for(auto s2=0; s2<=s1; ++s2) {
				GaussianShell& shellB = shells[s2];
				int ncartB = shellB.ncartesian();
				
				TwoIndex<double> shellPairInts(ncartA, ncartB, 0.0);
				std::array<TwoIndex<double>, 9> shellPairGrads;
				std::array<TwoIndex<double>, 45> shellPairHess;
				if (deriv > 0) {
					for (auto& r : shellPairGrads) r.assign(ncartA, ncartB, 0.0);
					if (deriv > 1)
						for (auto& r : shellPairHess) r.assign(ncartA, ncartB, 0.0);
				}
				
				
				for (int i = 0; i < ecps.getN(); i++) {
					ECP& U = ecps.getECP(i);
					
					ecpint->compute_shell_pair(U, shellA, shellB, tempValues);
					shellPairInts.add(tempValues);
					
					if (deriv > 0) {
						ecpint->compute_shell_pair_derivative(U, shellA, shellB, tempGrads);
						for (int i = 0; i < 9; i++) shellPairGrads[i].add(tempGrads[i]);
						if (deriv > 1) {
							ecpint->compute_shell_pair_second_derivative(U, shellA, shellB, tempHess);
							for (int i = 0; i < 45; i++) shellPairHess[i].add(tempHess[i]);
						}
					}
				}
			
				for (int i = n1; i < n1 + ncartA; i++) {
					for (int j = n2; j < n2 + ncartB; j++) {
						integrals(i, j) = shellPairInts(i-n1, j-n2);
						integrals(j, i) = integrals(i, j);
					}
				}
			
				if (deriv > 0) {
					for (int i = n1; i < n1 + ncartA; i++) {
						for (int j = n2; j < n2 + ncartB; j++) {
							for (int k = 0; k < 9; k++) {
								first_derivs[k](i, j) = shellPairGrads[k](i-n1, j-n2);
								first_derivs[k](j, i) = first_derivs[k](i, j);
							}
						}
					}
					
					if (deriv > 1) {
						for (int i = n1; i < n1 + ncartA; i++) {
							for (int j = n2; j < n2 + ncartB; j++) {
								for (int k = 0; k < 45; k++) {
									second_derivs[k](i, j) = shellPairHess[k](i-n1, j-n2);
									second_derivs[k](j, i) = second_derivs[k](i, j);
								}
							}
						}
					}
				}
			
				n2 += ncartB;
			}
		
			n1 += ncartA; 
		}
	}
}
