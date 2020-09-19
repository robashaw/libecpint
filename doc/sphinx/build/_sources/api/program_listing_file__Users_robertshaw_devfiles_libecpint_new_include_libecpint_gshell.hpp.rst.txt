
.. _program_listing_file__Users_robertshaw_devfiles_libecpint_new_include_libecpint_gshell.hpp:

Program Listing for File gshell.hpp
===================================

|exhale_lsh| :ref:`Return to documentation for file <file__Users_robertshaw_devfiles_libecpint_new_include_libecpint_gshell.hpp>` (``/Users/robertshaw/devfiles/libecpint_new/include/libecpint/gshell.hpp``)

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
   
   #ifndef GSHELL_HEAD
   #define GSHELL_HEAD
   
   #include <vector>
   
   namespace libecpint {
   
       struct GaussianShell {
           std::vector<double> exps; 
           std::vector<double> coeffs; 
           
           double* centerVec; 
           
           int l; 
       
           GaussianShell(double* A, int l);
           
           void addPrim(double exp, double c);
           
           int nprimitive() const { return exps.size(); }
           
           int ncartesian() const { return ((l+1)*(l+2))/2; }
           
           double* center() const { return centerVec; };
           
           double exp(int i) const { return exps[i]; }
           
           double coef(int i) const { return coeffs[i]; }
           
           int am() const { return l; }
           
           GaussianShell copy() const {
               GaussianShell result(centerVec, l);
               result.exps = exps;
               result.coeffs = coeffs;
               return result;
           }
       };
   
       struct ShellPairData {
           int LA;         
           int LB;         
           int maxLBasis;  
           int ncartA;     
           int ncartB;     
           double A[3];    
           double B[3];    
           double A2;      
           double Am;      
           double B2;      
           double Bm;      
           double RAB2;    
           double RABm;    
       };
   
   }
   
   #endif
