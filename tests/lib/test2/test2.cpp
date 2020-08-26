#include "ecp.hpp"
#include "gshell.hpp"
#include "ecpint.hpp"
#include "multiarr.hpp"
#include "testutil.hpp"
#include <iostream>

int main(int argc, char* argv[]) {
	using namespace libecpint; 
	
	double A[3] = {0.5, 1.5, 0.9};
	double B[3] = {-0.5, -0.5, -0.5}; 
	double C[3] = {0.0, 0.0, 1.4};
	ECP newU(C);
	newU.addPrimitive(2, 4, 1.0, 0.0);
	newU.addPrimitive(2, 0, 2.015185434, 0.437903605);
	newU.addPrimitive(2, 1, 1.225772391, 0.418499875);
	newU.addPrimitive(2, 2, 3.415601371, 10.482661239);
	newU.addPrimitive(2, 3, 3.574403408, -4.478284427, true);
	
	GaussianShell shellA(A, 0); 
	GaussianShell shellB(B, 2); 
	
	shellA.addPrim(0.005, 0.3);
	shellA.addPrim(2.49, 0.9);
	shellA.addPrim(10.11, 1.2); 
	
	shellB.addPrim(0.02, 1.1);
	shellB.addPrim(1.41, 0.9);
	shellB.addPrim(8.71, 0.4); 
	
	ECPIntegral ecpint(2, 4); 
	TwoIndex<double> result; 
	std::vector<double> flat_result; 
	
	ecpint.compute_shell_pair(newU, shellA, shellA, result);
	for (auto v : result.data) flat_result.push_back(v); 
	
	ecpint.compute_shell_pair(newU, shellA, shellB, result);
	for (auto v : result.data) flat_result.push_back(v); 
	
	ecpint.compute_shell_pair(newU, shellB, shellB, result);
	for (auto v : result.data) flat_result.push_back(v); 

	return check_file<double>("test2.output", flat_result);
}
