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

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
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
int check_file(std::string filename, std::vector<T>& results, double thresh = 1e-5,
               double precision = 1e-10, double maxrel = 1e-3) {
  std::ifstream input_file(filename);
  if (input_file.is_open()) {
    std::vector<T> benchmark;
    std::string line;
    while (!input_file.eof()) {
      std::getline(input_file, line);
      T converted_value;
      std::stringstream ss(line);
      ss >> converted_value;
      benchmark.push_back(converted_value);
    }

    if (benchmark.size() != results.size()) {
      std::cerr << "Size of output is incorrect!" << std::endl;
      std::cerr << "File has " << benchmark.size() << " records" << std::endl;
      std::cerr << "Results has " << results.size() << " records" << std::endl;
      return 1;
    } else {
      // Scale for deciding which elements are physically significant: elements far below
      // the largest entry are dominated by quadrature noise (integrals that are zero by
      // symmetry), so a *relative*-error test on them is meaningless and not reproducible
      // across platforms. Such elements are instead checked by absolute tolerance; the
      // relative (average and worst-case) gates are applied only to significant elements.
      double maxval = 0.0;
      for (int i = 0; i < benchmark.size(); i++)
        maxval = std::max(maxval, std::abs((double)benchmark[i]));
      const double sig_floor = std::max(precision, 1e-6 * maxval);
      const double abs_tol = std::max(precision, 1e-6 * maxval);  // tolerance for ~zero elements

      double error = 0.0;
      double worst = 0.0;      // largest single-element relative error among significant elements
      double worst_abs = 0.0;  // largest absolute error among insignificant (~zero) elements
      int nsig = 0;
      for (int i = 0; i < benchmark.size(); i++) {
        double abserror = std::abs(benchmark[i] - results[i]);
        if (abserror > precision) {
          std::cout << std::setw(10) << "Line " << std::setw(5) << i << std::setw(5) << ":"
                    << std::setw(15) << benchmark[i] << " / " << std::setw(15) << results[i]
                    << std::endl;
        }
        if (std::abs((double)benchmark[i]) > sig_floor) {
          double relerror = abserror / std::abs(benchmark[i]);
          error += relerror;
          worst = std::max(worst, relerror);
          nsig++;
        } else {
          worst_abs = std::max(worst_abs, abserror);
        }
      }
      error /= double(std::max(1, nsig));

      if (error > thresh) {
        std::cerr << "Average relative error over significant elements is " << error << std::endl;
        return 1;
      } else if (worst > maxrel) {
        // Guard against a few badly-wrong significant elements being masked by the average
        std::cerr << "Largest single-element relative error is " << worst << ", exceeding "
                  << maxrel << std::endl;
        return 1;
      } else if (worst_abs > abs_tol) {
        // An element that should be ~zero came out non-negligible
        std::cerr << "Largest absolute error on a near-zero element is " << worst_abs
                  << ", exceeding " << abs_tol << std::endl;
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

}  // namespace libecpint
#endif
