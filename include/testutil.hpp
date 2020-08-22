#ifndef TESTING_HEAD
#define TESTING_HEAD

#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <sstream>
#include <vector>

namespace libecpint {

template <typename T>
int check_file(std::string filename, std::vector<T>& results) {
	std::ifstream input_file(filename); 
	if (input_file.is_open()) {
		
		std::vector<T> benchmark; 
		std::string line;
		while(!input_file.eof()) {
			std::getline(input_file, line); 
			T converted_value; 
			std::stringstream ss(line);
			ss >> converted_value;
			benchmark.push_back(converted_value); 
		}
		
		if (benchmark.size() != results.size()) {
			std::cerr << "Size of output is incorrect!" << std::endl;
			std::cerr << "File has " << benchmark.size() << " records" << std::endl;
			std::cerr << "Results has " << results.size()  << " records" << std::endl;
			return 1;
		} else {
			double error = 0.0;
			for (int i = 0; i < benchmark.size(); i++) {
				double abserror = std::abs(benchmark[i] - results[i]);
				if (std::abs(benchmark[i])>1e-10) error += abserror / std::abs(benchmark[i]);
			}
			error /= double(benchmark.size());
    
			if (error > 5e-6) {
				std::cerr << "Average error in output is " << error << " percent!" << std::endl;
				return 1;
			} else {
				std::cout << "Test passed!" << std::endl; 
			}
		}
		
	} else {
		std::cerr << "Problem opening results file " << filename << std::endl; 
		return 1; 
	}
	
	return 0;
}

} // end namespace
#endif