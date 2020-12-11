
.. _program_listing_file__Users_robertshaw_devfiles_libecpint_include_libecpint_bessel.hpp:

Program Listing for File bessel.hpp
===================================

|exhale_lsh| :ref:`Return to documentation for file <file__Users_robertshaw_devfiles_libecpint_include_libecpint_bessel.hpp>` (``/Users/robertshaw/devfiles/libecpint/include/libecpint/bessel.hpp``)

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
   
   #ifndef BESSEL_FUNCTION_HEAD
   #define BESSEL_FUNCTION_HEAD
   
   #include <vector>
   
   namespace libecpint {
   
       constexpr double SMALL = 1.0E-7; 
       constexpr int TAYLOR_CUT = 5; 
   
       class BesselFunction 
       {
       private:
           int lMax; 
           int N; 
           int order; 
           double scale; 
       
           double **K; 
           double ***dK; 
           double *C; 
       
           int tabulate(double accuracy);
       
       public:
           BesselFunction();
           
           BesselFunction(int lMax, int N, int order, double accuracy);
           
           ~BesselFunction();
       
           void init(int lMax, int N, int order, double accuracy);
       
           void calculate(double z, int maxL, std::vector<double> &values) const;
           
           double calculate(double z, int L) const;
           
           double upper_bound(double z, int L) const;
       };
   
   }
   #endif
