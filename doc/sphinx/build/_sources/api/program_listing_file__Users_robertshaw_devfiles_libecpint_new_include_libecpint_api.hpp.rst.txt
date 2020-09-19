
.. _program_listing_file__Users_robertshaw_devfiles_libecpint_new_include_libecpint_api.hpp:

Program Listing for File api.hpp
================================

|exhale_lsh| :ref:`Return to documentation for file <file__Users_robertshaw_devfiles_libecpint_new_include_libecpint_api.hpp>` (``/Users/robertshaw/devfiles/libecpint_new/include/libecpint/api.hpp``)

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
   
   #ifndef API_HEAD
   #define API_HEAD
   
   #include <vector>
   #include <array>
   #include <string>
   #include <memory>
   #include "ecpint.hpp"
   #include "multiarr.hpp"
   
   namespace libecpint {
       
       struct ECPIntegrator {
           std::vector<GaussianShell> shells;
           std::vector<std::array<double, 3>> centers;
           ECPBasis ecps; 
           std::shared_ptr<ECPIntegral> ecpint; 
           int maxLB, deriv, ncart;
           
           bool ecp_is_set, basis_is_set;
           
           TwoIndex<double> integrals; 
           std::array<TwoIndex<double>, 9> first_derivs;
           std::array<TwoIndex<double>, 45> second_derivs;
           
           ECPIntegrator() { ecp_is_set = basis_is_set = false; maxLB = ncart = 0; }
           
           void set_gaussian_basis(int nshells, double* coords, double* exponents, double* coefs, int* ams, int* shell_lengths);
           void set_ecp_basis(int necps, double* coords, double* exponents, double* coefs, int* ams, int* ns, int* shell_lengths);
           //void set_ecp_basis_from_library(int necps, double* coords, int* charges, std::vector<std::string> names);
           void update_gaussian_basis_coords(int nshells, double* coords);
           void update_ecp_basis_coords(int necps, double* coords);
           
           void init(int deriv_ = 0);
           void compute();
           
           std::shared_ptr<std::vector<double>> get_integrals() { return std::make_shared<std::vector<double>>(integrals.data); }
           
           std::array<std::shared_ptr<std::vector<double>>, 9> get_first_derivs() {
               std::array<std::shared_ptr<std::vector<double>>, 9> results;
               for (int i = 0; i < 9; i++) results[i] = std::make_shared<std::vector<double>>(first_derivs[i].data);
               return results;
           }
           
           std::array<std::shared_ptr<std::vector<double>>, 45> get_second_derivs() {
               std::array<std::shared_ptr<std::vector<double>>, 45> results;
               for (int i = 0; i < 45; i++) results[i] = std::make_shared<std::vector<double>>(second_derivs[i].data);
               return results;
           }
       };
       
   }
   
   #endif
