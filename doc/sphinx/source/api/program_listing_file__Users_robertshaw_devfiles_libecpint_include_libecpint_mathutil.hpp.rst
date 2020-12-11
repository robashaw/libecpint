
.. _program_listing_file__Users_robertshaw_devfiles_libecpint_include_libecpint_mathutil.hpp:

Program Listing for File mathutil.hpp
=====================================

|exhale_lsh| :ref:`Return to documentation for file <file__Users_robertshaw_devfiles_libecpint_include_libecpint_mathutil.hpp>` (``/Users/robertshaw/devfiles/libecpint/include/libecpint/mathutil.hpp``)

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
   
   #ifndef MATHUTILHEADERDEF
   #define MATHUTILHEADERDEF
   
   #include <vector>
   #include <numeric>
   #include <cmath>
   #include "multiarr.hpp"
   
   #ifndef M_PI
   #define M_PI 3.14159265358979323846
   #endif
   
   #ifndef M_EULER
   #define M_EULER 2.71828182845904523536
   #endif
   
   #define MAX_FAC 100 
   #define MAX_DFAC 200 
   
   namespace libecpint {
       
       constexpr double ROOT_PI = 1.772453850905516; 
       constexpr double SINH_1 = 1.1752011936;
       constexpr int MAX_POW = 20;
   
       extern double FAC[MAX_FAC];     
       extern double DFAC[MAX_DFAC];   
       
       double pow_m2(double);
       double pow_m1(double);
       double pow_0(double);
       double pow_1(double);
       double pow_2(double);
       double pow_3(double);
       double pow_4(double);
       double pow_5(double);
       double pow_6(double);
       double pow_7(double);
       double pow_8(double);
       double pow_9(double);
       double pow_10(double);
       double pow_11(double);
       double pow_12(double);
       double pow_13(double);
       double pow_14(double);
       double pow_15(double);
       double pow_16(double);
       double pow_17(double);
       double pow_18(double);
       double pow_19(double);
       double pow_20(double);
       
       static double (*FAST_POW[23])(double) {pow_0, pow_1, pow_2, pow_3, pow_4, pow_5,
                                         pow_6, pow_7, pow_8, pow_9, pow_10, pow_11,
                                         pow_12, pow_13, pow_14, pow_15, pow_16, pow_17,
                                         pow_18, pow_19, pow_20, pow_m1, pow_m2};
       
       
       const double GAMMA[30] = {
           1.7724538509055,
           1.0,
           0.88622692545275,
           1.0,
           1.3293403881791,
           2.0,
           3.3233509704478,
           6.0,
           11.631728396567,
           24.0,
           52.342777784553,
           120.0,
           287.88527781504,
           720.0,
           1871.2543057978,
           5040.0,
           14034.407293483,
           40320.0,
           1.1929246199461e5,
           3.62880e5,
           1.1332783889488e6,
           3.628800e6,
           1.1899423083962e7,
           3.9916800e7,
           1.3684336546556e8, 
           4.79001600e8, 
           1.7105420683196e9, 
           6.227020800e9, 
           2.3092317922314e10,
           8.7178291200e10
       };
       
       TwoIndex<double> realSphericalHarmonics(int lmax, double x, double phi);  
       
       double frobenius_norm(const TwoIndex<double>& mat);
       
       void initFactorials(); 
   }
   
   #endif
