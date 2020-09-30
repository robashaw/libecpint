
.. _program_listing_file__Users_robertshaw_devfiles_libecpint_new_include_libecpint_angular.hpp:

Program Listing for File angular.hpp
====================================

|exhale_lsh| :ref:`Return to documentation for file <file__Users_robertshaw_devfiles_libecpint_new_include_libecpint_angular.hpp>` (``/Users/robertshaw/devfiles/libecpint_new/include/libecpint/angular.hpp``)

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
   
   #ifndef ANGULAR_HEAD
   #define ANGULAR_HEAD
   
   #include <vector>
   #include "multiarr.hpp"
   
   namespace libecpint {
       
       class AngularIntegral 
       {
       private: 
           int LB; 
           int LE; 
   
           int wDim; 
           int maxL; 
       
           FiveIndex<double> W; 
           SevenIndex<double> omega; 
       
           double calcG(int l, int m) const;
           double calcH1(int i, int j, int l, int m) const;
           double calcH2(int i, int j, int k, int m) const;
       
           ThreeIndex<double> Pijk(int maxI) const; 
       
           void makeW( FiveIndex<double> &U);
           void makeOmega(FiveIndex<double> &U);
       
       public:
       
           ThreeIndex<double> uklm(int lam, int mu) const;
       
           FiveIndex<double> makeU();
       
           AngularIntegral(); 
           
           AngularIntegral(int LB, int LE); 
           
           void init(int LB, int LE);
           
           void compute();
       
           void clear();
       
           double getIntegral(int k, int l, int m, int lam, int mu) const; 
           
           double getIntegral(int k, int l, int m, int lam, int mu, int rho, int sigma) const;
           
           int* getOmegaMults() { return omega.mults; }
           int* getOmegaDims() { return omega.dims; }
           std::vector<double>& getOmegaData() { return omega.data; }
           
           
           bool isZero(int k, int l, int m, int lam, int mu, double tolerance) const;
           
           bool isZero(int k, int l, int m, int lam, int mu, int rho, int sigma, double tolerance) const;  
       };
   
   }
   #endif
