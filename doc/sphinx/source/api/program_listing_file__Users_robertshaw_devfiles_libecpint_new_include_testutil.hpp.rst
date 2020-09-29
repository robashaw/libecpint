
.. _program_listing_file__Users_robertshaw_devfiles_libecpint_new_include_testutil.hpp:

Program Listing for File testutil.hpp
=====================================

|exhale_lsh| :ref:`Return to documentation for file <file__Users_robertshaw_devfiles_libecpint_new_include_testutil.hpp>` (``/Users/robertshaw/devfiles/libecpint_new/include/testutil.hpp``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   /* 
    *      Copyright (c) 2020 Robert Shaw
    *      This file is a part of Libecpint.
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
