
.. _program_listing_file__Users_robertshaw_devfiles_libecpint_include_libecpint_gaussquad.hpp:

Program Listing for File gaussquad.hpp
======================================

|exhale_lsh| :ref:`Return to documentation for file <file__Users_robertshaw_devfiles_libecpint_include_libecpint_gaussquad.hpp>` (``/Users/robertshaw/devfiles/libecpint/include/libecpint/gaussquad.hpp``)

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
   
   #ifndef GC_QUAD_HEAD
   #define GC_QUAD_HEAD
   
   #include <functional>
   #include <vector>
   
   namespace libecpint {
   
       enum GCTYPE {
           ONEPOINT, 
           TWOPOINT  
       };
   
       class GCQuadrature {
       private:
           int maxN; 
           int M;  
       
           std::vector<double> x; 
           std::vector<double> w; 
    
           double I; 
       
           GCTYPE t; 
       
           double sumTerms(const std::function<double(double, const double*, int)> &f,
                     const double *p, int limit, int start, int end, int shift, int skip) const;
   
       public:
           
           GCQuadrature();
           
           GCQuadrature(const GCQuadrature &other);
       
           void initGrid(int points, GCTYPE t);
       
           std::pair<double, bool> integrate(
               std::function<double(double, const double*, int)> &f,
               const double *params, double tolerance, int start, int end) const;
       
           void transformZeroInf();
           
           void transformRMinMax(double z, double p);  // Transfromation from [-1, 1] to [rmin, rmax] from Flores06
           void untransformRMinMax(double z, double p);
       
           int getN() const { return maxN; }
       
       std::vector<double>& getX() { return x; }
       const std::vector<double>& getX() const { return x; }
       };
   }
   
   #endif
