
.. _program_listing_file__Users_robertshaw_devfiles_libecpint_new_include_libecpint_qgen.hpp:

Program Listing for File qgen.hpp
=================================

|exhale_lsh| :ref:`Return to documentation for file <file__Users_robertshaw_devfiles_libecpint_new_include_libecpint_qgen.hpp>` (``/Users/robertshaw/devfiles/libecpint_new/include/libecpint/qgen.hpp``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   /* 
   *      Copyright (c) 2020 Robert Shaw
   *         This file is a part of Libecpint.
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
   
   #ifndef QGEN_HEAD
   #define QGEN_HEAD
   
   #include "multiarr.hpp"
   #include "radial.hpp"
   #include "angular.hpp"
   #include "gshell.hpp"
   #include "ecp.hpp"
   
   namespace libecpint {
       namespace qgen {
   
           void rolled_up(int lam, int LA, int LB, ThreeIndex<double>& radials, FiveIndex<double>& CA, FiveIndex<double>& CB, TwoIndex<double>& SA, TwoIndex<double>& SB, AngularIntegral& angints, ThreeIndex<double>& values);
   
           void rolled_up_special(int lam, int LA, int LB, ThreeIndex<double>& radials, FiveIndex<double>& CB, TwoIndex<double>& SB, AngularIntegral& angints, ThreeIndex<double>& values);
           
       void Q0_0_0(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q0_0_1(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q0_0_2(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q0_0_3(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q0_0_4(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q0_0_5(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q0_1_0(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q0_1_1(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q0_1_2(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q0_1_3(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q0_1_4(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q0_1_5(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q1_1_0(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q1_1_1(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q1_1_2(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q1_1_3(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q1_1_4(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q1_1_5(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q0_2_0(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q0_2_1(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q0_2_2(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q0_2_3(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q0_2_4(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q0_2_5(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q1_2_0(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q1_2_1(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q1_2_2(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q1_2_3(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q1_2_4(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q1_2_5(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q2_2_0(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q2_2_1(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q2_2_2(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q2_2_3(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q2_2_4(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q2_2_5(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q0_3_0(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q0_3_1(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q0_3_2(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q0_3_3(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q0_3_4(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q0_3_5(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q1_3_0(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q1_3_1(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q1_3_2(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q1_3_3(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q1_3_4(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q1_3_5(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q2_3_0(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q2_3_1(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q2_3_2(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q2_3_3(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q2_3_4(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q2_3_5(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q3_3_0(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q3_3_1(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q3_3_2(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q3_3_3(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q3_3_4(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q3_3_5(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q0_4_0(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q0_4_1(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q0_4_2(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q0_4_3(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q0_4_4(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q0_4_5(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q1_4_0(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q1_4_1(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q1_4_2(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q1_4_3(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q1_4_4(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q1_4_5(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q2_4_0(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q2_4_1(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q2_4_2(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q2_4_3(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q2_4_4(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q2_4_5(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q3_4_0(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q3_4_1(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q3_4_2(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q3_4_3(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q3_4_4(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q3_4_5(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q4_4_0(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q4_4_1(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q4_4_2(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q4_4_3(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q4_4_4(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q4_4_5(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q0_5_0(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q0_5_1(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q0_5_2(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q0_5_3(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q0_5_4(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q0_5_5(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q1_5_0(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q1_5_1(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q1_5_2(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q1_5_3(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q1_5_4(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q1_5_5(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q2_5_0(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q2_5_1(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q2_5_2(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q2_5_3(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q2_5_4(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q2_5_5(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q3_5_0(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q3_5_1(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q3_5_2(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q3_5_3(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q3_5_4(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q3_5_5(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q4_5_0(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q4_5_1(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q4_5_2(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q4_5_3(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q4_5_4(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q4_5_5(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q5_5_0(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q5_5_1(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q5_5_2(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q5_5_3(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q5_5_4(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
       void Q5_5_5(ECP&, GaussianShell&, GaussianShell&, FiveIndex<double>&, FiveIndex<double>&, TwoIndex<double>&, TwoIndex<double>&, double, double, RadialIntegral&, AngularIntegral&, ThreeIndex<double>&);
   
   }
   }
   #endif
