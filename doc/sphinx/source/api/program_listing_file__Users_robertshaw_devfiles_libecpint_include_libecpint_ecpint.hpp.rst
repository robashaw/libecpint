
.. _program_listing_file__Users_robertshaw_devfiles_libecpint_include_libecpint_ecpint.hpp:

Program Listing for File ecpint.hpp
===================================

|exhale_lsh| :ref:`Return to documentation for file <file__Users_robertshaw_devfiles_libecpint_include_libecpint_ecpint.hpp>` (``/Users/robertshaw/devfiles/libecpint/include/libecpint/ecpint.hpp``)

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
   
       static constexpr double tolerance = 1e-12;
       
           double calcC(int a, int m, double A) const;
           
           static void(*QGEN[LIBECPINT_MAX_L+1][LIBECPINT_MAX_L+1][LIBECPINT_MAX_L+1])(
               const ECP&, const GaussianShell&, const GaussianShell&,
           const FiveIndex<double>&, const FiveIndex<double>&,
           const TwoIndex<double>&, const TwoIndex<double>&,
           double, double,
           const RadialIntegral&, const AngularIntegral&, const RadialIntegral::Parameters&,
           ThreeIndex<double>&);
   
       public:
           int skipped, zero, nonzero;
           
           void makeC(FiveIndex<double> &C, int L, const double *A) const;
           
           ECPIntegral(int maxLB, int maxLU, int deriv=0);
       
           void type1(const ECP& U, const GaussianShell &shellA, const GaussianShell &shellB,
                  const ShellPairData &data, const FiveIndex<double> &CA, const FiveIndex<double> &CB,
                  const RadialIntegral::Parameters & parameters, TwoIndex<double> &values) const;
           
           void type2(int l,
                  const ECP& U, const GaussianShell &shellA, const GaussianShell &shellB,
                  const ShellPairData &data, const FiveIndex<double> &CA, const FiveIndex<double> &CB,
                  const RadialIntegral::Parameters & parameters, ThreeIndex<double> &values) const;
           
           void estimate_type2(
           const ECP& U, const GaussianShell &shellA, const GaussianShell &shellB,
           const ShellPairData &data, double* results) const;
       
           void compute_shell_pair(
           const ECP &U, const GaussianShell &shellA, const GaussianShell &shellB,
           TwoIndex<double> &values, int shiftA = 0, int shiftB = 0) const;
           
           void compute_shell_pair_derivative(
           const ECP &U, const GaussianShell &shellA, const GaussianShell &shellB,
           std::array<TwoIndex<double>, 9> &results) const;
           
           void compute_shell_pair_second_derivative(
           const ECP &U, const GaussianShell &shellA, const GaussianShell &shellB,
           std::array<TwoIndex<double>, 45> &results) const;
           
           void left_shell_derivative(
           const ECP &U, const GaussianShell &shellA, const GaussianShell &shellB,
           std::array<TwoIndex<double>, 3> &results) const;
           
           void left_shell_second_derivative(
           const ECP &U, const GaussianShell &shellA, const GaussianShell &shellB,
           std::array<TwoIndex<double>, 6> &results) const;
           
           void mixed_second_derivative(
           const ECP &U, const GaussianShell &shellA, const GaussianShell &shellB,
           std::array<TwoIndex<double>, 9> &results) const;
           
       };
   
   }
   #endif
