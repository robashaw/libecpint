#include "ecp.hpp"
#include "gshell.hpp"
#include "ecpint.hpp"
#include "multiarr.hpp"
#include "testutil.hpp"
#include <iostream>
#include <array>

int main(int argc, char* argv[]) {
	using namespace libecpint; 
	
	double A[3] = {0.5, 1.5, 0.9};
	double B[3] = {-0.5, -0.5, -0.5}; 
	double C[3] = {0.0, 0.0, 0.0};
	ECP newU(C);
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
	
	ecpint.compute_shell_pair_second_derivative(newU, shellB, shellB, results);
	/*for (int i = 0; i < 45; i++) {
		for (auto v : results[i].data) 
			flat_result.push_back(v); 
	}*/
	for (int i = 0; i < results[0].data.size(); i++)
		flat_result.push_back(results[0].data[i] + 2.0*results[6].data[i] + results[24].data[i]);

	/*ecpint.compute_shell_pair_second_derivative(newU, shellA, shellB, results);
	for (int i = 0; i < 45; i++) {
		for (auto v : results[i].data) 
			flat_result.push_back(v); 
	}*/
	
	/*ecpint.compute_shell_pair_second_derivative(newU, shellB, shellB, results);
	for (int i = 0; i < 45; i++) {
		for (auto v : results[i].data) 
			flat_result.push_back(v); 
	}*/
	
	
	TwoIndex<double> Iplus, Iminus, I0;
	double rv;
	//xx
	ecpint.compute_shell_pair(newU, shellB, shellB, I0);
	shellB.centerVec[0] += 0.01; 
	ecpint.compute_shell_pair(newU, shellB, shellB, Iplus);
	shellB.centerVec[0] -= 0.02;
	ecpint.compute_shell_pair(newU, shellB, shellB, Iminus);

	for (int i = 0; i < I0.data.size(); i++) {
		rv = 10000.0 * (Iplus.data[i] - 2.0*I0.data[i] + Iminus.data[i]);
		flat_finite.push_back(rv);
	}
	
	for (int i = 0; i < flat_finite.size(); i++) {
		 std::cout << flat_result[i] << "  " << flat_finite[i] << std::endl;
	}
	//return check_file<double>("deriv_test1.output", flat_result);
}
