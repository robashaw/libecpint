
.. _program_listing_file__Users_robertshaw_devfiles_libecpint_include_libecpint_radial.hpp:

Program Listing for File radial.hpp
===================================

|exhale_lsh| :ref:`Return to documentation for file <file__Users_robertshaw_devfiles_libecpint_include_libecpint_radial.hpp>` (``/Users/robertshaw/devfiles/libecpint/include/libecpint/radial.hpp``)

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
   
   #ifndef RADIAL_HEAD
   #define RADIAL_HEAD
   
   #include <vector>
   #include "multiarr.hpp"
   #include "gaussquad.hpp"
   #include "ecp.hpp"
   #include "bessel.hpp"
   #include "gshell.hpp"
   
   namespace libecpint {
   
       constexpr double MIN_EXP = 0.002;
       class RadialIntegral
       {
       private:
           GCQuadrature bigGrid;
           GCQuadrature smallGrid;
           GCQuadrature primGrid; 
           BesselFunction bessie;
   
           double tolerance;
       
           static double integrand(double r, const double *p, int ix);
   
           void buildBessel(const std::vector<double> &r, int nr, int maxL, TwoIndex<double> &values, double weight = 1.0) const;
       
           double calcKij(double Na, double Nb, double zeta_a, double zeta_b, double R2) const;
       
           void buildU(const ECP &U, const int l, const int N, const GCQuadrature &grid, double *Utab) const;
       
           void buildF(
           const GaussianShell &shell, double A, int lstart, int lend,
           const std::vector<double> &r, int nr, int start, int end,
           TwoIndex<double> &F) const;
       
           int integrate(int maxL, int gridSize, const TwoIndex<double> &intValues, GCQuadrature &grid, std::vector<double> &values,
                   int start, int end, int offset = 0, int skip = 1) const;
           
           void compute_base_integrals(int N_min, int N_max, double p, double o_root_p, double P1,
                                       double P2, double P1_2, double P2_2, double X1, double X2,
                                       double oP1, double oP2, double* values) const;
                                       
           std::pair<double, bool> integrate_small(
               int N, int l1, int l2, double n, double a, double b, double A, double B) const;
           
       public:
           RadialIntegral();
       
           void init(int maxL, double tol = 1e-15, int small = 256, int large = 1024);
   
           struct Parameters {
         TwoIndex<double> p, P, P2, K;
           };
           Parameters buildParameters(
               const GaussianShell &shellA, const GaussianShell &shellB, const ShellPairData &data) const;
       
           void type1(int maxL, int N, int offset,
                  const ECP &U, const GaussianShell &shellA, const GaussianShell &shellB,
                  const ShellPairData &data, const Parameters & parameters, TwoIndex<double> &values) const;
       
           void type2(int lam, int l1start, int l1end, int l2start, int l2end, int N,
                  const ECP &U, const GaussianShell &shellA, const GaussianShell &shellB,
                  const ShellPairData &data, const Parameters & parameters, TwoIndex<double> &values) const;
   
           void type2(
           const std::vector<Triple> &triples, int nbase, int lam,
           const ECP &U, const GaussianShell &shellA, const GaussianShell &shellB,
           double A, double B, ThreeIndex<double> &radials) const;
       
           double estimate_type2(int N, int l1, int l2, double n, double a, double b, double A, double B) const;
       };
   
   }
   
   #endif
