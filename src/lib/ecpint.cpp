/* 
 *      Copyright (c) 2017 Robert Shaw
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

#include "ecpint.hpp"
#include <iostream>
#include <cmath>
#include "Faddeeva.hpp"
#include "mathutil.hpp"
#include "qgen.hpp"
#include <cassert>

namespace libecpint {

	ECPIntegral::ECPIntegral(int maxLB, int maxLU, int deriv) { 
		// Make sure library can perform requested integrals
		assert(maxLB <= LIBECPINT_MAX_L); 
		assert(maxLU <= LIBECPINT_MAX_L);
		
		// Initialise singletons
		initFactorials();
		
		// Initialise angular and radial integrators
		angInts.init(maxLB + deriv, maxLU);
		angInts.compute();
		radInts.init(2*(maxLB+deriv) + maxLU, 1e-15, 256, 512);
	};

	double ECPIntegral::calcC(int a, int m, double A) const {
		double value = 1.0 - 2*((a-m) % 2);
		value *= pow(A, a-m);
		value *= FAC[a]/(FAC[m] * FAC[a-m]);
		return value;
	}

	void ECPIntegral::makeC(FiveIndex<double> &C, int L, double *A) {
		int z; double Ck, Cl;
		int na = 0;
		for (int x = L; x >= 0; x--) {
			for (int y = L-x; y >= 0; y--) {
				z = L - x - y;
			
				for (int k = 0; k<= x; k++) {
					Ck = calcC(x, k, A[0]);
					for (int l = 0; l <= y; l++) {
						Cl = calcC(y, l, A[1]);
						for (int m = 0; m <= z; m++) C(0, na, k, l, m) = Ck * Cl * calcC(z, m, A[2]);
					}
				}
				na++;
			}
		}
	}

	void ECPIntegral::type1(ECP &U, GaussianShell &shellA, GaussianShell &shellB, ShellPairData &data, FiveIndex<double> &CA, FiveIndex<double> &CB, TwoIndex<double> &values) { 

		int LA = data.LA; int LB = data.LB;
		int maxLBasis = data.maxLBasis;
	
		// Build radial integrals
		int L = LA + LB;
		TwoIndex<double> temp;
		ThreeIndex<double> radials(L+1, L+1, 2*L+1);
		for (int ix = 0; ix <= L; ix++) {
			radInts.type1(ix, ix, ix % 2, U, shellA, shellB, data, temp);
			for(int l = 0; l <= ix; l++) {
				for (int m = -l; m <= l; m++) radials(ix, l, l+m) = temp(l, l+m);
			}
		}
	
		// Unpack positions
		double Ax = data.A[0]; double Ay = data.A[1]; double Az = data.A[2];
		double Bx = data.B[0]; double By = data.B[1]; double Bz = data.B[2];
	
		// Calculate chi_ab for all ab in shells
		int z1, z2, lparity, mparity, msign, ix, k, l, m;
		double C;
		int na = 0, nb = 0;
		for (int x1 = LA; x1 >= 0; x1--) {
			for (int y1 = LA-x1; y1 >= 0; y1--) {
				z1 = LA - x1 - y1;
				nb = 0;
			
				for (int x2 = LB; x2 >= 0; x2--) {
					for (int y2 = LB-x2; y2 >= 0; y2--) {
						z2 = LB - x2 - y2;
					
						for (int k1 = 0; k1 <= x1; k1++) {
							for (int k2 = 0; k2 <= x2; k2++) {
								k = k1 + k2;
							
								for (int l1 = 0; l1 <= y1; l1++) {
									for (int l2 = 0; l2 <= y2; l2++) {
										l = l1 + l2;
									
										for (int m1 = 0; m1 <= z1; m1++) {
											for (int m2 = 0; m2 <= z2; m2++){
												m = m1 + m2;
												C = CA(0, na, k1, l1, m1) * CB(0, nb, k2, l2, m2);
												if ( fabs(C) > 1e-14 ) {
													// Build radial integrals
													ix = k + l + m;
													
													// Certain terms can be neglected as the angular integrals will always be zero
													// See Flores06 appendix for details. 
													lparity = ix % 2; 
													msign = 1 - 2*(l%2);
													mparity = (lparity + m) % 2;
												
													for (int lam = lparity; lam <= ix; lam+=2) {
														for (int mu = mparity; mu <= lam; mu+=2)
															values(na, nb) += C * angInts.getIntegral(k, l, m, lam, msign*mu) * radials(ix, lam, lam+msign*mu);
													}
								
												}
											}
										}
									}
								}
							}
						}
					
						values(na, nb) *= 4.0 * M_PI;
						nb++;
					}
				}
			
				na++;
			}
		}
	
	}

	void ECPIntegral::type2(int lam, ECP& U, GaussianShell &shellA, GaussianShell &shellB, ShellPairData &data, FiveIndex<double> &CA, FiveIndex<double> &CB, ThreeIndex<double> &values) {
	
		// Unpack some data for convenience
		int LA = data.LA;
		int LB = data.LB;
		int L = LA + LB;	
		int maxLBasis = data.maxLBasis;
	
		double Am = data.Am; double Bm = data.Bm;
		
		// If shellA or shellB are on the same centre as the ECP, simpler integrals can be performed
		bool A_on_ecp = Am < 1e-7;
		bool B_on_ecp = Bm < 1e-7;

		if (A_on_ecp && B_on_ecp) {
			
			// Both on ECP, simplest case - see Shaw2017 supplementary material
			double prefactor = 4.0 * M_PI; 
			int npA = shellA.nprimitive();
			int npB = shellB.nprimitive();
			int npC = U.getN(); 
		
			double zA, zB, zC, dA, dB, dC, p; 
			int nC, z1, z2;
			
			int na = 0;
			for (int x1 = LA; x1 >= 0; x1--) {
				for (int r1 = LA-x1; r1 >= 0; r1--) {
					z1 = LA - x1 - r1; 
			
					int nb = 0;
					for (int x2 = LB; x2 >= 0; x2--) {
						for (int y2 = LB - x2; y2 >= 0; y2--) {
							z2 = LB - x2 - y2; 
						
							double value = 0.0;
							for (int c = 0; c < npC; c++) {
								GaussianECP& g = U.getGaussian(c);
								if (g.l == lam) {
									zC = g.a;
									dC = g.d;
									nC = g.n; 
 
									for (int a = 0; a < npA; a++) {
										zA = shellA.exp(a);
										dA = shellA.coef(a);
									
										for (int b = 0; b < npB; b++) {
											zB = shellB.exp(b);
											dB = shellB.coef(b); 
										
											p = zA + zB + zC;
										
											double o_root_p = 1.0 / sqrt(p);
											int N = 2 + LA + LB + nC;
											value += 0.5*dA*dB*dC*GAMMA[N]*pow(o_root_p, N+1); 
										}
									}
								}
							}
							
							for (int mu = -lam; mu <= lam; mu++) {
						
								double angular = prefactor * angInts.getIntegral(x1, r1, z1, lam, mu, 0, 0) * angInts.getIntegral(x2, y2, z2, lam, mu, 0, 0); 
								values(na, nb, lam+mu) = angular * value; 
							}
							nb++; 
						}
					}
				
					na++;
				}
			}
			
		} else {
			
			// At least one of the shells is not on the ECP, so spherical harmonics will be required
			
			double xA = Am > 0 ? data.A[2] / Am : 0.0;
			double xB = Bm > 0 ? data.B[2] / Bm : 0.0;
			double phiA = atan2(data.A[1], data.A[0]);
			double phiB = atan2(data.B[1], data.B[0]);
			TwoIndex<double> SA = realSphericalHarmonics(lam+LA, xA, phiA);
			TwoIndex<double> SB = realSphericalHarmonics(lam+LB, xB, phiB);
		
			if (A_on_ecp || B_on_ecp) {
				// Radial integrals need to be calculated by a different recursive scheme, or by quadrature
				ThreeIndex<double> radials(L+1, lam + LA + 1, lam + LB + 1); 
				TwoIndex<double> temp;

				for (int N = 0; N < L+1; N++) {
					radInts.type2(lam, 0, lam + LA, 0, lam + LB, N, U, shellA, shellB, data, temp); 
					for (int l1 = 0; l1 < lam + LA + 1; l1++)
						for (int l2 = 0; l2 < lam + LB + 1; l2++)
							radials(N, l1, l2) = temp(l1, l2);
				}
				
				// TODO: Write a version of rolled_up specifically for this case, as a significant number of terms
				// can be neglected a priori - see Shaw2017 supplementary material. 
				qgen::rolled_up(lam, LA, LB, radials, CA, CB, SA, SB, angInts, values);
				
			} else {
				
				// Neither is on the ECP, the full recursive scheme with generated integrals can be used
				// Need LA <= LB, but symmetry means we can just swap the arguments if LB > LA. 
				if (LA <= LB) 
					QGEN[LA][LB][lam](U, shellA, shellB, CA, CB, SA, SB, Am, Bm, radInts, angInts, values);
				else {
					ThreeIndex<double> temp_values(data.ncartB, data.ncartA, 2*U.getL() + 1); 
					QGEN[LB][LA][lam](U, shellB, shellA, CB, CA, SB, SA, Bm, Am, radInts, angInts, temp_values);
					for (int na = 0; na < data.ncartA; na++)
						for (int nb = 0; nb < data.ncartB; nb++)
							for (int nu = 0; nu < 2*U.getL() + 1; nu++)
								values(na, nb, nu) = temp_values(nb, na, nu); 
				}
					
			}			
		}
	}

	void ECPIntegral::compute_shell_pair(ECP &U, GaussianShell &shellA, GaussianShell &shellB, TwoIndex<double> &values, int shiftA, int shiftB) {
	
		ShellPairData data;
		
		// Shift A and B to be relative to U
		const double* C = U.center();
		data.A[0] = shellA.center()[0] - C[0]; 
		data.A[1] = shellA.center()[1] - C[1];
		data.A[2] = shellA.center()[2] - C[2]; 
		data.B[0] = shellB.center()[0] - C[0]; 
		data.B[1] = shellB.center()[1] - C[1];
		data.B[2] = shellB.center()[2] - C[2]; 
	
	  	// Construct data that will be reused everywhere, and takes account of derivative shifts
		data.LA = shellA.am() + shiftA; 
		data.LB = shellB.am() + shiftB;
		data.maxLBasis = data.LA > data.LB ? data.LA : data.LB;
		data.ncartA = (data.LA+1)*(data.LA+2)/2;
		data.ncartB = (data.LB+1)*(data.LB+2)/2;
	
		data.A2 = data.A[0]*data.A[0] + data.A[1]*data.A[1] + data.A[2]*data.A[2];
		data.Am = sqrt(data.A2);
		data.B2 = data.B[0]*data.B[0] + data.B[1]*data.B[1] + data.B[2]*data.B[2];
		data.Bm = sqrt(data.B2);
		double RAB[3] = {data.A[0] - data.B[0], data.A[1] - data.B[1], data.A[2] - data.B[2]};
		data.RAB2 = RAB[0]*RAB[0] + RAB[1]*RAB[1] + RAB[2]*RAB[2];
		data.RABm = sqrt(data.RAB2);
		
		// Prepare the radial integrator
		radInts.buildParameters(shellA, shellB, data);
	
		// Construct coefficients 
		FiveIndex<double> CA(1, data.ncartA, data.LA+1, data.LA+1, data.LA+1);
		FiveIndex<double> CB(1, data.ncartB, data.LB+1, data.LB+1, data.LB+1);
		makeC(CA, data.LA, data.A);
		makeC(CB, data.LB, data.B);
	
		// Calculate type1 integrals, if necessary
		values.assign(data.ncartA, data.ncartB, 0.0);
		if (!U.noType1())
			type1(U, shellA, shellB, data, CA, CB, values);
		
		// Now all the type2 integrals
		ThreeIndex<double> t2vals(data.ncartA, data.ncartB, 2*U.getL() + 1);
		for (int l = 0; l < U.getL(); l++) {
			t2vals.fill(0.0);
			type2(l, U, shellA, shellB, data, CA, CB, t2vals);
		
			for (int m = -l; m <= l; m++) {
				for(int na = 0; na < data.ncartA; na++) {
					for (int nb = 0; nb < data.ncartB; nb++) {
						values(na, nb) += t2vals(na, nb, l+m);
					}
				}
			}
		}
	}
	
	void ECPIntegral::left_shell_derivative(ECP &U, GaussianShell &shellA, GaussianShell &shellB, std::array<TwoIndex<double>, 3> &results) {
		int LA = shellA.am();
		int LB = shellB.am();
		
		int ncartB = (LB+1) * (LB+2) / 2;
		int ncartA = (LA+1) * (LA+2) / 2;
		int ncartA_plus = (LA+2) * (LA+3) / 2;
		int ncartA_minus = LA * (LA+1) / 2;
		TwoIndex<double> Q_minus, Q_plus; 
		
		for (int i = 0; i < 3; i++) results[i].assign(ncartA, ncartB, 0.0); 
		
		if (LA != 0)
			compute_shell_pair(U, shellA, shellB, Q_minus, -1, 0); 
		
		// hack in the exponents to the coefficients
		GaussianShell tempA = shellA.copy();
		for (int i = 0; i < tempA.nprimitive(); i++) 
			tempA.coeffs[i] *= tempA.exps[i];
		compute_shell_pair(U, tempA, shellB, Q_plus, 1, 0); 
		
		// Now compile the derivatives
		if (LA != 0) {
			int nA = 0;
			int nA_minus, nA_plus;
			for (int k=LA; k >= 0; k--) {
				for (int l=LA-k; l>=0; l--) {
					int m = LA - k - l;
						
					for (int nB = 0; nB < ncartB; nB++) {
						nA_minus = nA_plus = N_INDEX(l, m);
						results[0](nA, nB) = -k*Q_minus(nA_minus, nB) + 2.0*Q_plus(nA_plus, nB);
						
						nA_minus = std::max(0, N_INDEX(l-1, m));
						nA_plus  = N_INDEX(l+1, m);
						results[1](nA, nB) = -l*Q_minus(nA_minus, nB) + 2.0*Q_plus(nA_plus, nB);
						
						nA_minus = std::max(0, N_INDEX(l, m-1));
						nA_plus  = N_INDEX(l, m+1);
						results[2](nA, nB) = -m*Q_minus(nA_minus, nB) + 2.0*Q_plus(nA_plus, nB);
					}
					nA += 1;
				}
			}
		} else {
			for (int nB = 0; nB < ncartB; nB++) {
				results[0](0, nB) = 2.0*Q_plus(0, nB);
				results[1](0, nB) = 2.0*Q_plus(1, nB);
				results[2](0, nB) = 2.0*Q_plus(2, nB);
			}
		}
	}
	
	void ECPIntegral::compute_shell_pair_derivative(ECP &U, GaussianShell &shellA, GaussianShell &shellB, std::array<TwoIndex<double>, 9> &results) {		
		// First we check centres
		double A[3], B[3], C[3];
		*A = *shellA.center();
		*B = *shellB.center();
		*C = *U.center();
		
		double dAC = std::abs(A[0] - C[0]) + std::abs(A[1] - C[1]) + std::abs(A[2] - C[2]);
		double dBC = std::abs(B[0] - C[0]) + std::abs(B[1] - C[1]) + std::abs(B[2] - C[2]);
		
		// Calculate shell derivatives
		std::array<TwoIndex<double>, 3> QA, QB;
		if (dAC > 1e-6) 
			left_shell_derivative(U, shellA, shellB, QA);
		if (dBC > 1e-6)
			left_shell_derivative(U, shellB, shellA, QB);
		
		// initialise results matrices
		int ncartA = (shellA.am()+1) * (shellA.am()+2) / 2;
		int ncartB = (shellB.am()+1) * (shellB.am()+2) / 2;
		
		// Now construct the nuclear derivs
		if (dAC > 1e-6) {
			results[0] = QA[0];
			results[1] = QA[1];
			results[2] = QA[2];
			if (dBC > 1e-6) {
				results[3] = QB[0].transpose();
				results[4] = QB[1].transpose();
				results[5] = QB[2].transpose();
				for (int i = 6; i < 9; i++) results[i].assign(ncartA, ncartB, 0.0);
				for (int nA = 0; nA < ncartA; nA++) {
					for (int nB = 0; nB < ncartB; nB++){
						results[6](nA, nB) = -1.0 * (results[0](nA, nB) + results[3](nA, nB));
						results[7](nA, nB) = -1.0 * (results[1](nA, nB) + results[4](nA, nB));
						results[8](nA, nB) = -1.0 * (results[2](nA, nB) + results[5](nA, nB));
					}
				}
			} else {
			   results[3] = results[0]; results[3].multiply(-1.0);
			   results[4] = results[1]; results[4].multiply(-1.0);
			   results[5] = results[2]; results[5].multiply(-1.0);
			}
		} else if (dBC > 1e-6) {
			results[3] = QB[0].transpose();
			results[4] = QB[1].transpose();
			results[5] = QB[2].transpose();
		   	results[0] = results[3]; results[0].multiply(-1.0);
		   	results[1] = results[4]; results[1].multiply(-1.0);
		   	results[2] = results[5]; results[2].multiply(-1.0);
		} else {
			// else everything is zero
			for (int i = 0; i < 9; i++) results[i].assign(ncartA, ncartB, 0.0);
		}
	}

}
