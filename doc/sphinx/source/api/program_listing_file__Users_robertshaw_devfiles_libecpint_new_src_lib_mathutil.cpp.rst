
.. _program_listing_file__Users_robertshaw_devfiles_libecpint_new_src_lib_mathutil.cpp:

Program Listing for File mathutil.cpp
=====================================

|exhale_lsh| :ref:`Return to documentation for file <file__Users_robertshaw_devfiles_libecpint_new_src_lib_mathutil.cpp>` (``/Users/robertshaw/devfiles/libecpint_new/src/lib/mathutil.cpp``)

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
   
   #include "mathutil.hpp"
   #include <iostream>
   
   namespace libecpint {
       
       double FAC[MAX_FAC];
       double DFAC[MAX_DFAC];
       
       void initFactorials() {
       #ifndef FAC_INIT
       #define FAC_INIT
               FAC[0] = 1.0;
               DFAC[0] = 1.0;
               DFAC[1] = 1.0;
           
               for (int i = 1; i < MAX_FAC; i++)  FAC[i] = double(i) * FAC[i-1]; 
               for (int i = 2; i < MAX_DFAC; i++) DFAC[i] = double(i) * DFAC[i-2];
       #endif
       }
       
       // Compute all the real spherical harmonics Slm(theta, phi) for l,m up to lmax
       // x = std::cos (theta)
       TwoIndex<double> realSphericalHarmonics(int lmax, double x, double phi){
           TwoIndex<double> rshValues(lmax+1, 2*lmax+1, 0.0);
   
           if (lmax > 0) {
               // First calculate the associated Legendre polynomials, Plm(std::cos theta), ustd::sing the recursion relation
               // (l-m)Plm = x(2l - 1)P{l-1}m - (l+m-1)P{l-2}m
               // along with the zeroth order term
               // Pmm = (-1)^m (2m-1)!!(1-x^2)^{m/2}
               double x2 = x * x;
               double Plm[lmax+1][lmax+1]; 
               // First get all Pmm terms
               Plm[0][0] = 1.0;
               double sox2 = std::sqrt(std::max(0.0, 1.0 - x2));
               double ox2m = 1.0;
               for (int m = 1; m <= lmax; m++) {
                   ox2m *= -sox2;
                   Plm[m][m] = ox2m * DFAC[2*m-1]; 
               }
           
               // Then increment l for each m
               Plm[1][0] = x;
               Plm[0][1] = 0.0;
               for (int l = 2; l <= lmax; l++) {
                   ox2m = x * (2*l - 1);
                   for (int m = 0; m < l; m++) {
                       Plm[l][m] = ox2m * Plm[l-1][m] - (l + m - 1)*Plm[l-2][m];
                       Plm[l][m] /= ((double) (l -m));
                   }
                   Plm[l-1][l] = 0.0;
               }
           
               // Now we compute the spherical harmonics via
               // Slm(theta, phi) = Clm * Plm(std::cos(theta)) * std::cos(m * phi), m > 0
               // Sl{-m}(theta, phi) = Clm * Plm(std::cos(theta)) * std::sin(m * phi)
               // Sl0(theta, phi) = std::sqrt(2) * Cl0 * Pl0(std::cos(theta))
               // where Clm^2 = (2l + 1)*(l - m)! / (8*pi * (l+m)!)
               double osq4pi = 1.0 / std::sqrt(4.0 * M_PI); 
               int sign;
               for (int l = 0; l <= lmax; l++) {
                   rshValues(l, l) = osq4pi * std::sqrt(2.0 * l + 1.0) * Plm[l][0];
                   sign = -1;
                   for (int m = 1; m <= l; m++) {
                       ox2m = (2.0 * l + 1.0) * FAC[l-m] / FAC[l+m];
                       ox2m = sign * osq4pi * std::sqrt(2.0 * ox2m) * Plm[l][m];
                       rshValues(l, l+m) = ox2m * std::cos(m * phi);
                       rshValues(l, l-m) = ox2m * std::sin(m * phi);
                       sign *= -1;
                   }
               }
           
           } else {
               rshValues(0, 0) = 1.0 / std::sqrt(4.0 * M_PI);
           }
           
           return rshValues;
       }
       
       double frobenius_norm(TwoIndex<double>& mat) {
           return std::sqrt(std::inner_product(mat.data.begin(), mat.data.end(), mat.data.begin(), 0.0));
       }
       
   }
