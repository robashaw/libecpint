#include "ecp.hpp"
#include "gshell.hpp"
#include "ecpint.hpp"
#include "multiarr.hpp"
#include <iostream>
#include <fstream>
#include <cmath>

int main(int argc, char* argv[]) {
	using namespace libecpint; 
	
	double A[3] = {0.5, 1.5, 0.9};
	double B[3] = {-0.5, -0.5, -0.5}; 
	ECP newU(A);
	newU.addPrimitive(2, 4, 1.0, 0.0);
	newU.addPrimitive(2, 0, 2.015185434, 0.437903605);
	newU.addPrimitive(2, 1, 1.225772391, 0.418499875);
	newU.addPrimitive(2, 2, 3.415601371, 10.482661239);
	newU.addPrimitive(2, 3, 3.574403408, -4.478284427, true);
	
	GaussianShell shellA(A, 1); 
	GaussianShell shellB(B, 3); 
	
	shellA.addPrim(0.15, 0.3);
	shellA.addPrim(2.49, 0.9);
	shellA.addPrim(10.11, 1.2); 
	
	shellB.addPrim(0.02, 1.1);
	shellB.addPrim(1.41, 0.9);
	shellB.addPrim(8.71, 0.4); 
	
	ECPIntegral ecpint(3, 4); 
	TwoIndex<double> result; 
	std::vector<double> flat_result; 
	
	ecpint.compute_shell_pair(newU, shellA, shellA, result);
	for (auto v : result.data) flat_result.push_back(v); 
	
	ecpint.compute_shell_pair(newU, shellA, shellB, result);
	for (auto v : result.data) flat_result.push_back(v); 
	
	ecpint.compute_shell_pair(newU, shellB, shellB, result);
	for (auto v : result.data) flat_result.push_back(v); 
	
	std::ifstream input_file("test1.output"); 
	if (input_file.is_open()) {
		
		std::vector<double> benchmark; 
		std::string line;
		while(!input_file.eof()) {
			std::getline(input_file, line); 
			benchmark.push_back(std::stod(line)); 
		}
		
		if (benchmark.size() != flat_result.size()) {
			std::cerr << "Size of output is incorrect!" << std::endl;
			return 1;
		} else {
			double error = 0.0;
			for (int i = 0; i < benchmark.size(); i++) {
				double abserror = std::abs(benchmark[i] - flat_result[i]);
				if (std::abs(benchmark[i])>1e-10) error += abserror / std::abs(benchmark[i]);
			}
			error /= double(benchmark.size());

			if (error > 5e-6) {
				std::cerr << "Average error in output is " << error << " percent!" << std::endl;
				return 1;
			} else {
				std::cout << "Test1 passed!" << std::endl; 
			}
		}
		
	} else {
		std::cerr << "Problem opening results file for test1!" << std::endl; 
		return 1; 
	}
	
	return 0;
}
