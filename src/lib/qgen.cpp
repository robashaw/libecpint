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
#include <iostream>

namespace libecpint {
	namespace qgen {	
		void rolled_up(int lam, int LA, int LB, ThreeIndex<double> &radials, FiveIndex<double> &CA, FiveIndex<double> &CB, TwoIndex<double> &SA, TwoIndex<double> &SB, AngularIntegral &angint, ThreeIndex<double> &values)
		{
			double prefac = 16.0 * M_PI * M_PI;
			int L = LA + LB;	
		
			int z1, z2, w_m;
			int w_ax, w_ay, w_az, w_l1, w_m1; 
			int w_bx, w_by, w_bz, w_l2, w_m2; 
			double C, val1, val2; 
			int* mults = angint.getOmegaMults();
			std::vector<double>& omega = angint.getOmegaData();
			int w_lam = lam * mults[3];
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
								w_ax = w_lam + alpha_x*mults[0];
								for (int alpha_y = 0; alpha_y <= r1; alpha_y++) {
									w_ay = w_ax + alpha_y*mults[1];
									for (int alpha_z = 0; alpha_z <= z1; alpha_z++) {
										w_az = w_ay + alpha_z*mults[2];
										int alpha = alpha_x + alpha_y + alpha_z; 
								
										for (int beta_x = 0; beta_x <= x2; beta_x++) {
											w_bx = w_lam + beta_x*mults[0];
											for (int beta_y = 0; beta_y <= y2; beta_y++) {
												w_by = w_bx + beta_y*mults[1];
												for (int beta_z = 0; beta_z <= z2; beta_z++) {
													w_bz = w_by + beta_z*mults[2]; 
													int beta = beta_x + beta_y + beta_z; 
													int N = alpha + beta; 
													C = CA(0, na, alpha_x, alpha_y, alpha_z) * CB(0, nb, beta_x, beta_y, beta_z); 
												
													if (std::abs(C) > 1e-15) {
														for (int lam1 = 0; lam1 <= lam + alpha; lam1++) {
															w_l1 = w_az + lam1*(1+mults[5]);
															int lam2start = (lam1 + N) % 2; 
															for (int lam2 = lam2start; lam2 <= lam + beta; lam2+=2) {
																w_l2 = w_bz + lam2*(1+mults[5]);
																val1 = prefac * C * radials(N, lam1, lam2);
													
																for (int mu1 = -lam1; mu1 <= lam1; mu1++) {
																	w_m1 = w_l1 + mu1;
																	for (int mu2 = -lam2; mu2 <= lam2; mu2++) {
																		w_m2 = w_l2 + mu2;
																		val2 = val1 * SA(lam1, lam1+mu1) * SB(lam2, lam2+mu2);
																		w_m = -mults[4];
																		for (int mu = -lam; mu <= lam; mu++) {
																			w_m += mults[4];
																			values(na, nb, lam+mu) += val2 * omega[w_m+w_m1] * omega[w_m+w_m2];
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
							}
					
							nb++;
						}
					}
			
					na++; 
				}
			}
		}
  
		void rolled_up_special(int lam, int LA, int LB, ThreeIndex<double>& radials, FiveIndex<double>& CB, TwoIndex<double>& SB, AngularIntegral& angint, ThreeIndex<double>& values) {
			double prefac = 8.0 * M_PI * std::sqrt(M_PI);
			int L = LA + LB;	
	
			int z1, z2; 
			double C, val1, val2; 

			int w_bx, w_by, w_bz, w_l2, w_m2, w1, w_m; 
			int* mults = angint.getOmegaMults();
			std::vector<double>& omega = angint.getOmegaData();
			int w_lam = lam * mults[3];
			
			// Loop over cartesian shell functions in alpha order, e.g. {xx xy, xz, yy, yz, zz} for l=2
			int na = 0; // Rows are shellA
			for (int x1 = LA; x1 >= 0; x1--) {
				for (int r1 = LA-x1; r1 >= 0; r1--) {
					z1 = LA - x1 - r1; 
					w1 = w_lam + x1*mults[0] + r1*mults[1] + z1*mults[2];
		
					int nb = 0; // Cols are shellB
					for (int x2 = LB; x2 >= 0; x2--) {
						for (int y2 = LB - x2; y2 >= 0; y2--) {
							z2 = LB - x2 - y2; 
				
							// Begin full ECP integral expansion
							int alpha = x1 + r1 + z1; 
							
							for (int beta_x = 0; beta_x <= x2; beta_x++) {
								w_bx = w_lam + beta_x*mults[0];
								for (int beta_y = 0; beta_y <= y2; beta_y++) {
									w_by = w_bx + beta_y*mults[1];
									for (int beta_z = 0; beta_z <= z2; beta_z++) {
										w_bz = w_by + beta_z*mults[2];
										int beta = beta_x + beta_y + beta_z; 
										int N = alpha + beta; 
										C = CB(0, nb, beta_x, beta_y, beta_z); 
											
										if (std::abs(C) > 1e-15) {
													
											int lam2start =  N % 2; 
											for (int lam2 = lam2start; lam2 <= lam + beta; lam2+=2) {
												w_l2 = w_bz + lam2*(1+mults[5]);
												val1 = prefac * C * radials(N, 0, lam2);
												
														
												for (int mu2 = -lam2; mu2 <= lam2; mu2++) {
													w_m2 = w_l2 + mu2;
													val2 = val1 *  SB(lam2, lam2+mu2);
													w_m = -mults[4];																		
													for (int mu = -lam; mu <= lam; mu++) {
														w_m += mults[4];
														values(na, nb, lam+mu) += val2 * omega[w1+w_m] * omega[w_m2+w_m];
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
