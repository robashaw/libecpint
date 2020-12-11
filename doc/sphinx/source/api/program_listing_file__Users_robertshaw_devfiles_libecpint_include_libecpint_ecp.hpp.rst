
.. _program_listing_file__Users_robertshaw_devfiles_libecpint_include_libecpint_ecp.hpp:

Program Listing for File ecp.hpp
================================

|exhale_lsh| :ref:`Return to documentation for file <file__Users_robertshaw_devfiles_libecpint_include_libecpint_ecp.hpp>` (``/Users/robertshaw/devfiles/libecpint/include/libecpint/ecp.hpp``)

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
   
   #ifndef ECP_HEAD
   #define ECP_HEAD
   
   #include <vector>
   #include <array>
   #include <map>
   #include <string>
   #include "config.hpp"
   
   namespace libecpint {
       
       const std::string atom_names[109] = {"h", "he", "li", "be", "b", "c", "n",
           "o", "f", "ne", "na", "mg", "al", "si", "p", "s", "cl", "ar",
           "k", "ca", "sc", "ti", "v", "cr", "mn", "fe", "co", "ni", "cu",
           "zn", "ga", "ge", "as", "se", "br", "kr", "rb", "sr", "y", "zr",
           "nb", "mo", "tc", "ru", "rh", "pd", "ag", "cd", "in", "sn", "sb",
           "te", "i", "xe", "cs", "ba", "la", "ce", "pr", "nd", "pm", "sm",
           "eu", "gd", "tb", "dy", "ho", "er", "tm", "yb", "lu", "hf", "ta",
           "w", "re", "os", "ir", "pt", "au", "hg", "tl", "pb", "bi", "po",
           "at", "rn", "fr", "ra", "ac", "th", "pa", "u", "np", "pu", "am",
           "cm", "bk", "cf", "es", "fm", "md", "no", "lr", "rf", "db", "sg",
           "bh", "hs", "mt" };
       
       struct GaussianECP {
           int n; 
           int l; 
           double a; 
           double d; 
       
           GaussianECP(); 
           
           GaussianECP(const int n, const int l, const double a, const double d);
           
           GaussianECP(const GaussianECP& other);
       };
   
       struct ECP {
           std::vector<GaussianECP> gaussians; 
           int N; 
           int L; 
           int atom_id; 
           double min_exp; 
           double min_exp_l[LIBECPINT_MAX_L+1]; 
           int    l_starts[LIBECPINT_MAX_L+2]; 
           
           std::array<double, 3> center_; 
       
           ECP();
           
           ECP(const double *_center);
           
           ECP(const ECP &other);
       
           void addPrimitive(const int n, const int l, const double a, const double d, const bool needSort = true);
           
           const double* center() const { return &center_[0]; }
           
           void setPos(const double x, const double y, const double z);
           
           void sort(); 
           
           GaussianECP& getGaussian(int i) { return gaussians[i]; }
       const GaussianECP& getGaussian(int i) const { return gaussians[i]; }
           int getN() const { return N; }
           
       bool noType1() const;
       
           double evaluate(const double r, const int l) const;
     
           int getL() const { return L; }
       
       };
   
       class ECPBasis {
       private:
           std::vector<ECP> basis;    
           std::vector<int> atomList; 
           int N; 
           int maxL; 
       
       public:
           ECPBasis(); 
           
           std::map<int, int> core_electrons;
           
           void addECP(const ECP &U, const int atom);
           
           ECP& getECP(const int i);
       const ECP& getECP(const int i) const;
           int getECPCore(const int q) const;
           
           int getAtom(int i) const { return atomList[i]; }
           
           int getMaxL() const { return maxL; }
           
           int getN() const { return N; }
           
           void addECP_from_file(
               const int q, const std::array<double, 3> & coords, const std::string & filename);
       };
   
   }
   
   #endif
