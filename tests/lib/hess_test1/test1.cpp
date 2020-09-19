#include "ecp.hpp"
#include "gshell.hpp"
#include "ecpint.hpp"
#include "multiarr.hpp"
#include "testutil.hpp"
#include <iostream>
#include <array>


using namespace libecpint; 

void finite_diff(ECPIntegral& ecpint, ECP& U, GaussianShell& shellA, GaussianShell& shellB, double h,
				int axis1, int axis2, GaussianShell& s1, GaussianShell& s2, std::vector<double> &results) {
	TwoIndex<double> Ipp, Imm, Imp, Ipm; 
	double rv;
	double o4h2 = 0.25 / (h * h);

	s1.centerVec[axis1] += h;
	s2.centerVec[axis2] += h;
	ecpint.compute_shell_pair(U, shellA, shellB, Ipp);
	
	s2.centerVec[axis2] -= 2.0*h;
	ecpint.compute_shell_pair(U, shellA, shellB, Ipm);
	
	s1.centerVec[axis1] -= 2.0*h;
	ecpint.compute_shell_pair(U, shellA, shellB, Imm);
	
	s2.centerVec[axis2] += 2.0*h;
	ecpint.compute_shell_pair(U, shellA, shellB, Imp);
	
	s1.centerVec[axis1] += h;
	s2.centerVec[axis2] -= h;

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
		for (auto v : results[i].data) 
			flat_result.push_back(v); 
	}
	

	ecpint.compute_shell_pair_second_derivative(newU, shellA, shellB, results);
	for (int i = 0; i < 45; i++) {
		for (auto v : results[i].data) 
			flat_result.push_back(v); 
	}
	
	ecpint.compute_shell_pair_second_derivative(newU, shellB, shellB, results);
	for (int i = 0; i < 45; i++) {
		for (auto v : results[i].data) 
			flat_result.push_back(v); 
	}
	
	return check_file<double>("hess_test1.output", flat_result);
}
