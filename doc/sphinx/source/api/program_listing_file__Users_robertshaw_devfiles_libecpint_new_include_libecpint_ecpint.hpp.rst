
.. _program_listing_file__Users_robertshaw_devfiles_libecpint_new_include_libecpint_ecpint.hpp:

Program Listing for File ecpint.hpp
===================================

|exhale_lsh| :ref:`Return to documentation for file <file__Users_robertshaw_devfiles_libecpint_new_include_libecpint_ecpint.hpp>` (``/Users/robertshaw/devfiles/libecpint_new/include/libecpint/ecpint.hpp``)

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
   
   #ifndef ECPINT_HEAD
   #define ECPINT_HEAD
   
   #include <vector>
   #include <array>
   #include "multiarr.hpp"
   #include "gaussquad.hpp"
   #include "ecp.hpp"
   #include "bessel.hpp"
   #include "radial.hpp"
   #include "angular.hpp"
   #include "gshell.hpp"
   
   #include "config.hpp"
   
   namespace libecpint {
   
   #define N_INDEX(l, m) (((l+m)*(l+m+1))/2 + m)
   
       class ECPIntegral
       {
       private:
           RadialIntegral radInts; 
           AngularIntegral angInts; 
       
           double calcC(int a, int m, double A) const;
           
           static void(*QGEN[LIBECPINT_MAX_L+1][LIBECPINT_MAX_L+1][LIBECPINT_MAX_L+1])(ECP&, GaussianShell&, GaussianShell&,
            FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double,
            RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
   
       public:
           
           void makeC(FiveIndex<double> &C, int L, double *A);
           
           ECPIntegral(int maxLB, int maxLU, int deriv=0);
       
           void type1(ECP& U, GaussianShell &shellA, GaussianShell &shellB, ShellPairData &data, FiveIndex<double> &CA, FiveIndex<double> &CB, TwoIndex<double> &values);
           
           void type2(int l, ECP& U, GaussianShell &shellA, GaussianShell &shellB, ShellPairData &data, FiveIndex<double> &CA, FiveIndex<double> &CB, ThreeIndex<double> &values);
       
           void compute_shell_pair(ECP &U, GaussianShell &shellA, GaussianShell &shellB, TwoIndex<double> &values, int shiftA = 0, int shiftB = 0);
           
           void compute_shell_pair_derivative(ECP &U, GaussianShell &shellA, GaussianShell &shellB, std::array<TwoIndex<double>, 9> &results);
           
           void compute_shell_pair_second_derivative(ECP &U, GaussianShell &shellA, GaussianShell &shellB, std::array<TwoIndex<double>, 45> &results);
           
           void left_shell_derivative(ECP &U, GaussianShell &shellA, GaussianShell &shellB, std::array<TwoIndex<double>, 3> &results); 
           
           void left_shell_second_derivative(ECP &U, GaussianShell &shellA, GaussianShell &shellB, std::array<TwoIndex<double>, 6> &results); 
           
           void mixed_second_derivative(ECP &U, GaussianShell &shellA, GaussianShell &shellB, std::array<TwoIndex<double>, 9> &results); 
           
       };
   
   }
   #endif
