#include "ecp.hpp"
#include "gshell.hpp"
#include "ecpint.hpp"
#include "multiarr.hpp"
#include "testutil.hpp"
#include <iostream>
#include <array>


using namespace libecpint; 

void finite_diff(ECPIntegral& ecpint,  double h, GaussianShell& shellA, GaussianShell& shellB, GaussianShell& s1, GaussianShell& s2,
			ECP& U, int axis1, int axis2, bool shift_s1, bool shift_s2, bool shift_ecp_1, bool shift_ecp_2, std::vector<double> &results) {
	TwoIndex<double> Ipp, Ipm, Imm, Imp;
	double rv;
	double o4h2 = 0.25 / (h * h);

	if (shift_s1) s1.centerVec[axis1] += h;
	if (shift_s2) s2.centerVec[axis2] += h;
	if (shift_ecp_1) U.center_[axis1] += h;
	if (shift_ecp_2) U.center_[axis2] += h;
	ecpint.compute_shell_pair(U, shellA, shellB, Ipp);
	
	if (shift_s1) s1.centerVec[axis1] -= 2*h;
	if (shift_ecp_1) U.center_[axis1] -= 2*h;
	ecpint.compute_shell_pair(U, shellA, shellB, Imp);
	
	if (shift_s2) s2.centerVec[axis2] -= 2*h;
	if (shift_ecp_2) U.center_[axis2] -= 2*h;
	ecpint.compute_shell_pair(U, shellA, shellB, Imm);
	
	if (shift_s1) s1.centerVec[axis1] += 2*h;
	if (shift_ecp_1) U.center_[axis1] += 2*h;
	ecpint.compute_shell_pair(U, shellA, shellB, Ipm);
	
	if (shift_s1) s1.centerVec[axis1] -= h;
	if (shift_s2) s2.centerVec[axis2] += h;
	if (shift_ecp_1) U.center_[axis1] -= h;
	if (shift_ecp_2) U.center_[axis2] += h;
	
	for (int i = 0; i < Ipp.data.size(); i++) {
		rv = Ipp.data[i] - Imp.data[i] - Ipm.data[i] + Imm.data[i];
		results.push_back(rv * o4h2);
	}
}


int main(int argc, char* argv[]) {
	
	double A[3] = {0.5, 1.5, 0.9};
	double B[3] = {-0.5, -0.5, -0.5}; 

	ECP newU(A);
	newU.addPrimitive(2, 4, 1.0, 0.0);
	newU.addPrimitive(2, 0, 2.015185434, 0.437903605);
	newU.addPrimitive(2, 1, 1.225772391, 0.418499875);
	newU.addPrimitive(2, 2, 3.415601371, 10.482661239);
	newU.addPrimitive(2, 3, 3.574403408, -4.478284427, true);
	
	GaussianShell shellA(A, 1); 
	GaussianShell shellB(B, 2); 
	
	shellA.addPrim(0.15, 0.3);
	shellA.addPrim(2.49, 0.9);
	shellA.addPrim(10.11, 1.2); 
	
	shellB.addPrim(0.02, 1.1);
	shellB.addPrim(1.41, 0.9);
	shellB.addPrim(8.71, 0.4); 
	
	ECPIntegral ecpint(2, 4, 2); 
	std::array<TwoIndex<double>, 45> results; 
	std::vector<double> flat_result, flat_finite;
	
	ecpint.compute_shell_pair_second_derivative(newU, shellA, shellA, results);
	for (int i = 0; i < 45; i++) {
		for (auto& v : results[i].data) 
			flat_result.push_back(v); 
	}

	ecpint.compute_shell_pair_second_derivative(newU, shellA, shellB, results);
	for (int i = 0; i < 45; i++) {
		for (auto& v : results[i].data) 
			flat_result.push_back(v); 
	}
	
	ecpint.compute_shell_pair_second_derivative(newU, shellB, shellB, results);
	for (int i = 0; i < 45; i++) {
		for (auto& v : results[i].data) 
			flat_result.push_back(v); 
	}
	
	/*double h = 0.005;
	finite_diff(ecpint, h, shellB, shellB, shellB, shellB, newU, 0, 0, true, true, false, false, flat_finite);
	finite_diff(ecpint, h, shellB, shellB, shellB, shellB, newU, 0, 1, true, true, false, false, flat_finite);
	finite_diff(ecpint, h, shellB, shellB, shellB, shellB, newU, 0, 2, true, true, false, false, flat_finite);
	finite_diff(ecpint, h, shellB, shellB, shellB, shellB, newU, 1, 1, true, true, false, false, flat_finite);
	finite_diff(ecpint, h, shellB, shellB, shellB, shellB, newU, 1, 2, true, true, false, false, flat_finite);
	finite_diff(ecpint, h, shellB, shellB, shellB, shellB, newU, 2, 2, true, true, false, false, flat_finite);
	                                                                 
	finite_diff(ecpint, h, shellB, shellB, shellB, shellA, newU, 0, 0, true, true, false, true, flat_finite);
	finite_diff(ecpint, h, shellB, shellB, shellB, shellA, newU, 0, 1, true, true, false, true, flat_finite);
	finite_diff(ecpint, h, shellB, shellB, shellB, shellA, newU, 0, 2, true, true, false, true, flat_finite);
	finite_diff(ecpint, h, shellB, shellB, shellB, shellA, newU, 1, 0, true, true, false, true, flat_finite);
	finite_diff(ecpint, h, shellB, shellB, shellB, shellA, newU, 1, 1, true, true, false, true, flat_finite);
	finite_diff(ecpint, h, shellB, shellB, shellB, shellA, newU, 1, 2, true, true, false, true, flat_finite);
	finite_diff(ecpint, h, shellB, shellB, shellB, shellA, newU, 2, 0, true, true, false, true, flat_finite);
	finite_diff(ecpint, h, shellB, shellB, shellB, shellA, newU, 2, 1, true, true, false, true, flat_finite);
	finite_diff(ecpint, h, shellB, shellB, shellB, shellA, newU, 2, 2, true, true, false, true, flat_finite);
	                                                                 
	finite_diff(ecpint, h, shellB, shellB, shellA, shellA, newU, 0, 0, true, true, true, true, flat_finite);
	finite_diff(ecpint, h, shellB, shellB, shellA, shellA, newU, 0, 1, true, true, true, true, flat_finite);
	finite_diff(ecpint, h, shellB, shellB, shellA, shellA, newU, 0, 2, true, true, true, true, flat_finite);
	finite_diff(ecpint, h, shellB, shellB, shellA, shellA, newU, 1, 1, true, true, true, true, flat_finite);
	finite_diff(ecpint, h, shellB, shellB, shellA, shellA, newU, 1, 2, true, true, true, true, flat_finite);
	finite_diff(ecpint, h, shellB, shellB, shellA, shellA, newU, 2, 2, true, true, true, true, flat_finite);	                          
	
	std::string names[21] = {"AxAx", "AxAy", "AxAz", "AyAy", "AyAz", "AzAz",
							 "AxCx", "AxCy", "AxCz", "AyCx", "AyCy", "AyCz", "AzCx", "AzCy", "AzCz",
						 	 "CxCx", "CxCy", "CxCz", "CyCy", "CyCz", "CzCz"};
    int next = 0; 
	int skip = 0;
	double v;
	for (int i = 0; i < flat_finite.size(); i++) {
		if (i % 36 == 0)  {
			std::cout << std::endl << names[next++] << std::endl;	
			if (names[next-1] == "AxCx") skip = 324;
			else if (names[next-1] == "CxCx") skip = 864;
		}
		
		if (next == 1) 
			v = flat_result[i] + 2*flat_result[i+216] + flat_result[i+864];
		else if (next == 2)
			v = flat_result[i] + flat_result[i+216] + flat_result[i+288]+ flat_result[i+864];
		else if (next == 3)
			v = flat_result[i] + flat_result[i+216] + flat_result[i+360]+ flat_result[i+864];
		else if (next == 4)
			v = flat_result[i] + 2*flat_result[i+252] + flat_result[i+864];	
		else if (next == 5)
			v = flat_result[i] + flat_result[i+252] + flat_result[i+324]+ flat_result[i+864];
		else if (next == 6)
			v = flat_result[i] + 2*flat_result[i+324]+ flat_result[i+864];		
		else if (next < 16)
			v = flat_result[i+skip] + flat_result[i+skip+540];
		else
			v = flat_result[i+skip];
		
		std::cout << v << " " << flat_finite[i] << " " << v - flat_finite[i] << std::endl;
	}*/
	
	return check_file<double>("hess_test1.output", flat_result);
}
