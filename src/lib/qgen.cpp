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

#include "qgen.hpp"
#include <cmath> 
#include "mathutil.hpp"

namespace libecpint {
  namespace qgen {	
	void rolled_up(int lam, int LA, int LB, ThreeIndex<double> &radials, FiveIndex<double> &CA, FiveIndex<double> &CB, TwoIndex<double> &SA, TwoIndex<double> &SB, AngularIntegral &angint, ThreeIndex<double> &values)
	{
		double prefac = 16.0 * M_PI * M_PI;
		int L = LA + LB;	
		
		int z1, z2; 
		double C, val1, val2; 

		// Loop over cartesian shell functions in alpha order, e.g. {xx xy, xz, yy, yz, zz} for l=2
		int na = 0; // Rows are shellA
		for (int x1 = LA; x1 >= 0; x1--) {
			for (int r1 = LA-x1; r1 >= 0; r1--) {
				z1 = LA - x1 - r1; 
			
				int nb = 0; // Cols are shellB
				for (int x2 = LB; x2 >= 0; x2--) {
					for (int y2 = LB - x2; y2 >= 0; y2--) {
						z2 = LB - x2 - y2; 
					
						// Begin full ECP integral expansion
						for (int alpha_x = 0; alpha_x <= x1; alpha_x++) {
							for (int alpha_y = 0; alpha_y <= r1; alpha_y++) {
								for (int alpha_z = 0; alpha_z <= z1; alpha_z++) {
									int alpha = alpha_x + alpha_y + alpha_z; 
								
									for (int beta_x = 0; beta_x <= x2; beta_x++) {
										for (int beta_y = 0; beta_y <= y2; beta_y++) {
											for (int beta_z = 0; beta_z <= z2; beta_z++) {
												int beta = beta_x + beta_y + beta_z; 
												int N = alpha + beta; 
												C = CA(0, na, alpha_x, alpha_y, alpha_z) * CB(0, nb, beta_x, beta_y, beta_z); 
												
												if (std::abs(C) > 1e-15) {
													for (int lam1 = 0; lam1 <= lam + alpha; lam1++) {
														
														int lam2start = (lam1 + N) % 2; 
														for (int lam2 = lam2start; lam2 <= lam + beta; lam2+=2) {
															val1 = prefac * C * radials(N, lam1, lam2);
													
															for (int mu1 = -lam1; mu1 <= lam1; mu1++) {
																for (int mu2 = -lam2; mu2 <= lam2; mu2++) {
															
																	val2 = val1 * SA(lam1, lam1+mu1) * SB(lam2, lam2+mu2);
																																
																	for (int mu = -lam; mu <= lam; mu++)
																		values(na, nb, lam+mu) += val2 * angint.getIntegral(alpha_x, alpha_y, alpha_z, lam, mu, lam1, mu1) * angint.getIntegral(beta_x, beta_y, beta_z, lam, mu, lam2, mu2); 
																}
															}
														}
													}
												}
											
											}
										}
									}
								}
							}
						}
					
						nb++;
					}
				}
			
				na++; 
			}
		}
	}
  }
}
