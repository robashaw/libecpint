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

#include <algorithm>
#include <set>

#include "generate.hpp"

/**
* \file generate.cpp
* \brief ECP integral code generator
*/ 

/** 
* For a given ECP integral, Q(LA, LB, lam), generates the integral code.
* This determines which radial integrals are necessary, based on angular integral screening.
* If LA, LB <= maxUnrol and (LA + LB + lam) <= 3*maxUnrol, then it unrols the angular integration as well.
* 
* @param LA - the shellA angular momentum
* @param LB - the shellB angular momentum
* @param lam - the ECP angular momentum
* @param angInts - the angular integrals
*/
void generate_lists(int LA, int LB, int lam, libecpint::AngularIntegral& angInts) { 
	using namespace libecpint;
	
	// Create the code file
	std::string ofname = "generated/Q" + std::to_string(LA) + std::to_string(LB) + std::to_string(lam) + ".cpp"; 
	std::ofstream outfile(ofname); 
	
	if (!outfile.is_open())
		std::cerr << "Problems writing to file!" << std::endl; 
	else {
		
		std::cout << "Generating Q(" << LA << ", " << LB << ", " << lam << ")... " << std::flush; 
		
		// Top matter
		outfile << "// Generated as part of Libecpint, Copyright 2017 Robert A Shaw" << std::endl; 
		outfile << "#include \"qgen.hpp\"" << std::endl; 
		outfile << "namespace libecpint {" << std::endl << "namespace qgen {" << std::endl;
		outfile << "void Q" << LA << "_" << LB << "_" << lam << "(const ECP& U, const GaussianShell& shellA, const GaussianShell& shellB, "
			<< "const FiveIndex<double> &CA, const FiveIndex<double> &CB, const TwoIndex<double> &SA, const TwoIndex<double> &SB, const double Am, const double Bm, "
				<< "const RadialIntegral &radint, const AngularIntegral& angint, const RadialIntegral::Parameters& parameters, ThreeIndex<double> &values) {" << std::endl << std::endl;
		
		// Do we need to unrol the angular integrals too? 
		bool unrolling = LA <= maxUnrol && LB <= maxUnrol && (LA + LB + lam) <= 3*maxUnrol;
		
		// Store the terms and radials if unrolling, just radial indices if not
		std::vector<SumTerm> terms; 
		std::set<Triple> radial_triples; 
		
		int z1, z2, w_m, w_l;
		int w_ax, w_ay, w_az, w_l1; 
		int w_bx, w_by, w_bz, w_l2; 
		double ang;
		double prefac = 16.0 * M_PI * M_PI;
		const int* mults = angInts.getOmegaMults();
		const std::vector<double>& omega = angInts.getOmegaData();
		int w_size = 2*lam+1; 
		double w1_contr[w_size*(lam+LA+1)];
		double w2_contr[w_size*(lam+LB+1)];
		
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
													
												for (int lam1 = 0; lam1 <= lam + alpha; lam1++) {
													w_l = lam1*w_size+lam; 
													w_l1 = w_az + lam1*(1+mults[5]);
													w_m = w_l1-mults[4];
													for (int mu = -lam; mu <= lam; mu++) {
														w_m += mults[4];
														w1_contr[w_l+mu] = 0.0;
														for (int mu1 = -lam1; mu1 <= lam1; mu1++)
															w1_contr[w_l+mu] += omega[w_m+mu1];
													}
												}
									
												for (int lam2 = 0; lam2 <= lam+beta; lam2++) {
													w_l  = lam2*w_size+lam;
													w_l2 = w_bz + lam2*(1+mults[5]);
													w_m = w_l2-mults[4];
													for (int mu = -lam; mu <= lam; mu++) {
														w_m += mults[4];
														w2_contr[w_l+mu] = 0.0;
														for (int mu2 = -lam2; mu2 <= lam2; mu2++) 
															w2_contr[w_l+mu] += omega[w_m+mu2];
													}
												}
														
												for (int lam1=0; lam1 <= lam+alpha; lam1++) {
													w_l1 = lam1*w_size+lam;
													int lam2start = (lam1 + N) % 2; 
													for (int lam2 = lam2start; lam2 <= lam + beta; lam2+=2) {
														w_l2 = lam2*w_size+lam;
														
														ang = 0.0;
														for (int mu = -lam; mu <= lam; mu++) 
															ang += prefac * w1_contr[w_l1+mu] * w2_contr[w_l2+mu];
														
														if (fabs(ang) > 1e-15) {
															radial_triples.insert(Triple{N, lam1, lam2}); 
															if (unrolling) {
																for (int mu = -lam; mu <= lam; mu++) {
																	for (int mu1 = -lam1; mu1 <= lam1; mu1++) {
																		for (int mu2 = -lam2; mu2 <= lam2; mu2++) {
																			SumTerm newTerm; 
																			newTerm.SA = Pair(lam1, lam1+mu1); 
																			newTerm.SB = Pair(lam2, lam2+mu2);
																			newTerm.radial = Triple(N, lam1, lam2);
																			newTerm.CA = Quintuple(0, na, alpha_x, alpha_y, alpha_z); 
																			newTerm.CB = Quintuple(0, nb, beta_x, beta_y, beta_z); 
																			newTerm.ang = prefac * angInts.getIntegral(alpha_x, alpha_y, alpha_z, lam, mu, lam1, mu1);
																			newTerm.ang *= angInts.getIntegral(beta_x, beta_y, beta_z, lam, mu, lam2, mu2); 
																			newTerm.mu = lam+mu; 
																			newTerm.na = na;
																			newTerm.nb = nb;
													
																			terms.push_back(newTerm); 
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
						}
				
						nb++;
					}
				}
		
				na++; 
			}
		}
		
		// Determine the maximum number of base integrals needed across the set of all radial integrals
		int nbase = 0; 
		if (!radial_triples.empty()) {
			const Triple& tmax = *radial_triples.crbegin();
			nbase = std::get<0>(tmax) + std::get<1>(tmax) - 1; 
			nbase = nbase < 0 ? 0 : nbase; 
		}
		
		// Sort the radials into two lists, depending on whether l1 <= l2 (radial_A), or l2 > l1 (radial_B)
		// swapping the order of l1/l2 in the latter case
		std::vector<Triple> radial_A, radial_B; 
		for (const Triple& t : radial_triples) {
			if (std::get<1>(t) <= std::get<2>(t)) radial_A.push_back(t);  
			else radial_B.push_back({std::get<0>(t), std::get<2>(t), std::get<1>(t)});
		}
		
		// Compute the correctly ordered radials first
		outfile << "\tstd::vector<Triple> radial_triples_A = {" << std::endl; 
		bool first = true; 
		for (Triple& t : radial_A) {
			if (!first) outfile << "," << std::endl; 
			else first = false;
			outfile << "\t\t{" + std::to_string(std::get<0>(t)) + ", "
				+ std::to_string(std::get<1>(t)) + ", " 
					+ std::to_string(std::get<2>(t)) + "}"; 
		}
		outfile << "\t};" << std::endl << std::endl;  
		
		outfile << "\tThreeIndex<double> radials(" << lam+LA+LB+1 << ", " << lam+LA+1 << ", " << lam+LB+1 << ");" << std::endl; 
		outfile << "\tradint.type2(radial_triples_A, " << nbase << ", " << lam << ", U, shellA, shellB, Am, Bm, radials);" << std::endl << std::endl;
		
		// Now compute the reverse-ordered radials
		outfile << "\tstd::vector<Triple> radial_triples_B = {" << std::endl; 
		first = true;
		for (Triple& t : radial_B) {
			if (!first) outfile << "," << std::endl; 
			else first = false; 
			outfile << "\t\t{" + std::to_string(std::get<0>(t)) + ", "
				+ std::to_string(std::get<1>(t)) + ", " 
					+ std::to_string(std::get<2>(t)) + "}"; 
		}
		outfile << "\t};" << std::endl << std::endl;  
		
		outfile << "\tThreeIndex<double> radials_B(" << lam+LA+LB+1 << ", " << lam+LB+1 << ", " << lam+LA+1 << ");" << std::endl; 
		outfile << "\tradint.type2(radial_triples_B, " << nbase << ", " << lam << ", U, shellB, shellA, Bm, Am, radials_B);" << std::endl;
		// These need to be compressed into the radials array, with l1/l2 reversed back
		outfile << "\tfor (Triple& t : radial_triples_B) radials(std::get<0>(t), std::get<2>(t), std::get<1>(t)) = radials_B(std::get<0>(t), std::get<1>(t), std::get<2>(t));" << std::endl << std::endl; 
		
		if (unrolling) {
			// Print out the unrolled angular integral code if needed
			std::cout << "unrolling... " << std::flush; 
			for (auto& term : terms) outfile << "\t" << term << std::endl; 
		} else {
			// Just use the generic rolled-up angular integral code
			outfile << "\trolled_up(" << lam << ", " << LA << ", " << LB << ", radials, CA, CB, SA, SB, angint, values);" << std::endl; 
		}
		outfile << "}" << std::endl << "}" << std::endl << "}" << std::endl; 
		
		std::cout << "done." << std::endl; 
		outfile.close();
	}
}


int main(int argc, char* argv[]) {
	
	// Factorial singletons will not have been initialised
	libecpint::initFactorials();
	int maxL = libecpint::maxL;
	
	libecpint::AngularIntegral angInts(maxL, maxL); 
	angInts.compute(); 
	
	// Generate the qgen.hpp header file
	std::string header_name; 
	if (argc > 1) {
		header_name = argv[1]; 
		header_name += "qgen.hpp"; 
	} else {
		header_name = "generated/qgen"; 
	}

	std::ifstream qgen_part("generated/qgen.part");
	std::ofstream qgen_head(header_name); 
	if (!qgen_part.is_open() || !qgen_head.is_open()) 
		std::cerr << "Problem creating qgen header file!" << std::endl; 
	else {
		std::string line; 
		while(!qgen_part.eof()) {
			std::getline(qgen_part, line); 
			qgen_head << line << std::endl;  
		}
		qgen_part.close();
		
		// Loop over all possible (l1, l2, lam) integrals up to l1 = l2 = lam = maxL
		// with l1 <= l2, generating the code and adding the function to the header file.
		for (int j = 0; j <= maxL; j++) {
			for (int i = 0; i <= j; i++) {
				for (int k = 0; k <= maxL; k++) {
					generate_lists(i, j, k, angInts); 
					qgen_head << "\tvoid Q" << i << "_" << j << "_" << k << "("
						<< "const ECP&, const GaussianShell&, const GaussianShell&, const FiveIndex<double>&, const FiveIndex<double>&, "
							<< "const TwoIndex<double>&, const TwoIndex<double>&, double, double, const RadialIntegral&, "
								<< "const AngularIntegral&, const RadialIntegral::Parameters&, ThreeIndex<double>&);" << std::endl;
				}
			}
		}
		qgen_head << std::endl << "}" << std::endl << "}" << std::endl; 
		qgen_head << "#endif" << std::endl; 
		qgen_head.close(); 
		
		// Now generate the function pointer array in ecpint_gen.cpp
		std::ifstream ecpgen_part("generated/ecpint_gen.part"); 
		std::ofstream ecpgen_head("generated/ecpint_gen.cpp"); 
		
		if (!ecpgen_part.is_open() || !ecpgen_head.is_open())
			std::cerr << "Problem reading/writing ecpgen file!" << std::endl;
		else {
			while(!ecpgen_part.eof()) {
				std::getline(ecpgen_part, line); 
				ecpgen_head << line << std::endl;  
			}
			ecpgen_part.close();
			
			for (int i =0; i <= maxL; i++) {
				ecpgen_head << "\t\t{ "; 
				
				for (int j = 0; j<= maxL; j++) {
					ecpgen_head << "\t\t\t{"; 
					
					int I = std::min(i, j);
					int J = std::max(i, j); 
					
					for (int k = 0; k< maxL; k++) 
						ecpgen_head << "qgen::Q" << I << "_" << J << "_" << k << ", ";
					
					ecpgen_head << "qgen::Q" << I << "_" << J << "_" << maxL << "}";
					if (j != maxL) ecpgen_head << ","; 
					ecpgen_head << std::endl;
				}
				
				ecpgen_head << "\t\t}";
				if (i != maxL) ecpgen_head << ","; 
				ecpgen_head << std::endl; 
			}
			
			ecpgen_head << "\t};" << std::endl << "}" << std::endl;
			ecpgen_head.close();  
		}
	}
	return 0; 
}
	