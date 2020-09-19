
.. _program_listing_file__Users_robertshaw_devfiles_libecpint_new_include_libecpint_radial.hpp:

Program Listing for File radial.hpp
===================================

|exhale_lsh| :ref:`Return to documentation for file <file__Users_robertshaw_devfiles_libecpint_new_include_libecpint_radial.hpp>` (``/Users/robertshaw/devfiles/libecpint_new/include/libecpint/radial.hpp``)

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
   
       class RadialIntegral
       {
       private:
           GCQuadrature bigGrid;
           GCQuadrature smallGrid;
           GCQuadrature primGrid; 
           BesselFunction bessie;
       
           TwoIndex<double> p, P, P2, K;
       
           double tolerance;
       
           static double integrand(double r, double *p, int ix);
   
           void buildBessel(std::vector<double> &r, int nr, int maxL, TwoIndex<double> &values, double weight = 1.0);
       
           double calcKij(double Na, double Nb, double zeta_a, double zeta_b, double R2) const;
       
           void buildU(ECP &U, int l, int N, GCQuadrature &grid, double *Utab);
       
           void buildF(GaussianShell &shell, double A, int lstart, int lend, std::vector<double> &r, int nr, int start, int end, TwoIndex<double> &F);
       
           int integrate(int maxL, int gridSize, TwoIndex<double> &intValues, GCQuadrature &grid, std::vector<double> &values, int offset = 0, int skip = 1);
           
           void compute_base_integrals(int N_min, int N_max, double p, double o_root_p, double P1,
                                       double P2, double P1_2, double P2_2, double X1, double X2,
                                       double oP1, double oP2, double* values); 
                                       
           std::pair<double, bool> integrate_small(int N, int l1, int l2, double n, double a, double b, double A, double B);
           
       public:
           RadialIntegral();
       
           void init(int maxL, double tol = 1e-15, int small = 256, int large = 1024);
       
           void buildParameters(GaussianShell &shellA, GaussianShell &shellB, ShellPairData &data);
       
           void type1(int maxL, int N, int offset, ECP &U, GaussianShell &shellA, GaussianShell &shellB, ShellPairData &data, TwoIndex<double> &values);
       
           void type2(int lam, int l1start, int l1end, int l2start, int l2end, int N, ECP &U, GaussianShell &shellA, GaussianShell &shellB, ShellPairData &data, TwoIndex<double> &values);    
   
           void type2(std::vector<Triple> &triples, int nbase, int lam, ECP &U, GaussianShell &shellA, GaussianShell &shellB, double A, double B, ThreeIndex<double> &radials); 
       };
   
   }
   
   #endif
