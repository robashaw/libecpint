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

#ifndef TESTING_HEAD
#define TESTING_HEAD

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <string>
#include <sstream>
#include <vector>

namespace libecpint {

/*! A helper function for tests that reads in a file of values and compares them
	to those provided by the test, returning 0 on success or 1 on failure.
	
	@tparam T - the type of the value; must be pipeable from a stringstream.
	@param filename - the file to read data from
	@param results - reference to the vector of calculated results
	@return 0 if the results agree with the file within 0.00005%, 1 otherwise
 */
template <typename T>
int check_file(std::string filename, std::vector<T>& results, double thresh=1e-5, double precision=1e-10) {
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
				if (abserror > precision) {
					std::cout << std::setw(10) << "Line " << std::setw(5) << i 
						<< std::setw(5) << ":" << std::setw(15) <<  benchmark[i] 
						<< " / " << std::setw(15) << results[i] << std::endl;
				}
				if (std::abs(benchmark[i])>precision) error += abserror / std::abs(benchmark[i]);
			}
			error /= double(benchmark.size());
    
			if (error > thresh) {
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
