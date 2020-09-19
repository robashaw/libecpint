
.. _program_listing_file__Users_robertshaw_devfiles_libecpint_new_include_libecpint_multiarr.hpp:

Program Listing for File multiarr.hpp
=====================================

|exhale_lsh| :ref:`Return to documentation for file <file__Users_robertshaw_devfiles_libecpint_new_include_libecpint_multiarr.hpp>` (``/Users/robertshaw/devfiles/libecpint_new/include/libecpint/multiarr.hpp``)

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
   
   #ifndef MULTIARR_HEAD
   #define MULTIARR_HEAD
   
   #include <vector>
   #include <tuple> 
   #include <algorithm>
   #include <sstream>
   
   namespace libecpint {
   
       using Pair = std::tuple<int, int>;
       using Triple = std::tuple<int, int, int>;
       using Quintuple = std::tuple<int, int, int, int, int>; 
       using Heptuple = std::tuple<int, int, int, int, int, int, int>; 
       
       namespace aux{
           template<std::size_t...> struct seq{};
           template<std::size_t N, std::size_t... Is>
           struct gen_seq : gen_seq<N-1, N-1, Is...>{};
           template<std::size_t... Is>
           struct gen_seq<0, Is...> : seq<Is...>{};
           template<class Ch, class Tr, class Tuple, std::size_t... Is>
           void print_tuple(std::basic_ostream<Ch,Tr>& os, Tuple const& t, seq<Is...>){
               using swallow = int[];
               (void)swallow{0, (void(os << (Is == 0? "" : ", ") << std::get<Is>(t)), 0)...};
           }
       } 
       
       template<class Ch, class Tr, class... Args>
       auto operator<<(std::basic_ostream<Ch, Tr>& os, std::tuple<Args...> const& t)
           -> std::basic_ostream<Ch, Tr>&
       {
           os << "(";
           aux::print_tuple(os, t, aux::gen_seq<sizeof...(Args)>());
           return os << ")";
       }
       
       template<typename T>
       struct TwoIndex {
           int dims[2];
           std::vector<T> data;
           T& operator()(int i, int j) { return data[i * dims[1] + j]; }
           T operator()(int i, int j) const { return data[i * dims[1] + j]; }
           void assign(int dim1, int dim2, T value) {
               dims[0] = dim1; dims[1] = dim2;
               data.resize(dim1 * dim2);
               std::fill(data.begin(), data.end(), value);
           }
           TwoIndex<T> transpose() {
               TwoIndex<T> result(dims[1], dims[0]);
               for (int i = 0; i < dims[0]; i++) {
                   for (int j = 0; j < dims[1]; j++)
                       result.data[j * dims[0] + i] = data[i * dims[1] + j];
               }
               return result;
           }
           void add(const TwoIndex<T>& other) {
               std::transform (data.begin(), data.end(), other.data.begin(), data.begin(), std::plus<T>());
           }
           
           void multiply(T k) {
               std::transform(data.begin(), data.end(), data.begin(), [&k](T& c){return c*k;});
           }
           TwoIndex() { dims[0] = dims[1] = 0; }
           TwoIndex(int dim1, int dim2) {
               dims[0] = dim1; dims[1] = dim2;
               data.resize(dim1 * dim2);
           }
           TwoIndex(int dim1, int dim2, T value) { assign(dim1, dim2, value); }
           TwoIndex(const TwoIndex<T> &other) { 
               data = other.data;
               dims[0] = other.dims[0]; dims[1] = other.dims[1];
           }
       };
   
       template<typename T>
       struct ThreeIndex {
           int dims[3];
           std::vector<T> data;
           T& operator()(int i, int j, int k) { return data[i*dims[2]*dims[1] + j*dims[2] + k]; }
           T operator()(int i, int j, int k) const { return data[i*dims[2]*dims[1] + j*dims[2] + k]; }
           ThreeIndex(){ dims[0] = 0; dims[1] = 0; dims[2] = 0; }
           ThreeIndex(int dim1, int dim2, int dim3) {
               dims[0] = dim1; dims[1] = dim2; dims[2] = dim3;
               data.resize(dim1 * dim2 * dim3);
           }
           ThreeIndex(const ThreeIndex<T> &other)  { 
               data = other.data;
               for (int n = 0; n < 3; n++) dims[n] = other.dims[n]; 
           }
           void fill(T value) { std::fill(data.begin(), data.end(), value); }
       };
   
       template<typename T>
       struct FiveIndex {
           int dims[5];
           std::vector<T> data;
           T& operator()(int i, int j, int k, int l, int m) { 
               return data[m + dims[4] * (l + dims[3] * (k + dims[2] * (j + dims[1] * i)))]; 
           }
           T operator()(int i, int j, int k, int l, int m) const { 
               return data[m + dims[4] * (l + dims[3] * (k + dims[2] * (j + dims[1] * i)))];
           }
           FiveIndex() { dims[0] = dims[1] = dims[2] = dims[3] = dims[4] = 0; }
           FiveIndex(int dim1, int dim2, int dim3, int dim4, int dim5) {
               dims[0] = dim1; dims[1] = dim2; dims[2] = dim3; dims[3] = dim4; dims[4] = dim5;
               data.resize(dim1 * dim2 * dim3 * dim4 * dim5);
           }
           FiveIndex(const FiveIndex<T> &other) { 
               data = other.data;
               for (int n = 0; n < 5; n++) dims[n] = other.dims[n]; 
           }
       };
   
       template<typename T>
       struct SevenIndex {
           int dims[7];
           std::vector<T> data;
           T& operator()(int i, int j, int k, int l, int m, int n, int p) {
               return data[p + dims[6]*(n + dims[5] * (m + dims[4] * (l + dims[3] * (k + dims[2] * (j + dims[1] * i)))))];
           }
           T operator()(int i, int j, int k, int l, int m, int n, int p) const {
               return data[p + dims[6]*(n + dims[5] * (m + dims[4] * (l + dims[3] * (k + dims[2] * (j + dims[1] * i)))))];
           }
           SevenIndex() { dims[0] = dims[1] = dims[2] = dims[3] = dims[4] = dims[5] = dims[6] = 0; }
           SevenIndex(int dim1, int dim2, int dim3, int dim4, int dim5, int dim6, int dim7) {
               dims[0] = dim1; dims[1] = dim2; dims[2] = dim3; dims[3] = dim4; dims[4] = dim5; dims[5] = dim6; dims[6] = dim7;
               data.resize(dim1 * dim2 * dim3 * dim4 * dim5 * dim6 * dim7);
           }
           SevenIndex(const SevenIndex<T> &other) {
               data = other.data;
               for (int n = 0; n < 7; n++) dims[n] = other.dims[n];
           }
       };
   
   }
   #endif
