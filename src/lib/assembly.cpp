/*
 *      Copyright (c) 2020 Robert Shaw
 *		This file is a part of Libecpint.
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

#include "assembly.hpp"

namespace libecpint {

void pack_first_derivatives(const std::array<TwoIndex<double>, 3>& QA,
                            const std::array<TwoIndex<double>, 3>& QB, const bool onA,
                            const bool onB, const int ncartA, const int ncartB,
                            std::array<TwoIndex<double>, 9>& results) {
  if (!onA) {
    results[0] = QA[0];
    results[1] = QA[1];
    results[2] = QA[2];
    if (!onB) {
      results[3] = QB[0];
      results[4] = QB[1];
      results[5] = QB[2];
      for (int i = 6; i < 9; i++) results[i].assign(ncartA, ncartB, 0.0);
      for (int nA = 0; nA < ncartA; nA++) {
        for (int nB = 0; nB < ncartB; nB++) {
          results[6](nA, nB) = -1.0 * (results[0](nA, nB) + results[3](nA, nB));
          results[7](nA, nB) = -1.0 * (results[1](nA, nB) + results[4](nA, nB));
          results[8](nA, nB) = -1.0 * (results[2](nA, nB) + results[5](nA, nB));
        }
      }
    } else {
      results[3] = results[0];
      results[3].multiply(-1.0);
      results[4] = results[1];
      results[4].multiply(-1.0);
      results[5] = results[2];
      results[5].multiply(-1.0);
      for (int i = 6; i < 9; i++) results[i].assign(ncartA, ncartB, 0.0);
    }
  } else if (!onB) {
    results[3] = QB[0];
    results[4] = QB[1];
    results[5] = QB[2];
    results[0] = results[3];
    results[0].multiply(-1.0);
    results[1] = results[4];
    results[1].multiply(-1.0);
    results[2] = results[5];
    results[2].multiply(-1.0);
    for (int i = 6; i < 9; i++) results[i].assign(ncartA, ncartB, 0.0);
  } else {
    // else everything is zero
    for (auto& r : results) r.assign(ncartA, ncartB, 0.0);
  }
}

void pack_second_derivatives(const std::array<TwoIndex<double>, 6>& QAA,
                             const std::array<TwoIndex<double>, 9>& QAB,
                             const std::array<TwoIndex<double>, 6>& QBB, const bool onA,
                             const bool onB, const int ncartA, const int ncartB,
                             std::array<TwoIndex<double>, 45>& results) {
  for (auto& r : results) r.assign(ncartA, ncartB, 0.0);

  // Position of component j of a 3x3 block (xx, xy, xz, yx, yy, yz, zx, zy, zz) in a symmetric
  // block (xx, xy, xz, yy, yz, zz), and of its transpose in a 3x3 block
  int jaas[9] = {0, 1, 2, 1, 3, 4, 2, 4, 5};
  int jbbs[9] = {0, 3, 6, 1, 4, 7, 2, 5, 8};
  int jaa, jbb;
  if (!onA) {
    // AA (xx, xy, xz, yy, yz, zz)
    for (int i = 0; i < 6; i++) results[i] = QAA[i];
    if (!onB) {
      // AB (xx, xy, xz, yx, yy, yz, zx, zy, zz)
      for (int i = 6; i < 15; i++) results[i] = QAB[i - 6];
      // BB (xx, xy, xz, yy, yz, zz)
      for (int i = 24; i < 30; i++) results[i] = QBB[i - 24];
      for (int nA = 0; nA < ncartA; nA++) {
        for (int nB = 0; nB < ncartB; nB++) {
          for (int j = 0; j < 9; j++) {
            jaa = jaas[j];
            jbb = jbbs[j];
            // AC (xx, xy, xz, yx, yy, yz, zx, zy, zz)
            results[15 + j](nA, nB) = -1.0 * (QAA[jaa](nA, nB) + QAB[j](nA, nB));
            // BC (xx, xy, xz, yx, yy, yz, zx, zy, zz)
            results[30 + j](nA, nB) = -1.0 * (QBB[jaa](nA, nB) + QAB[jbb](nA, nB));
            // CC (xx, xy, xz, yy, yz, zz)
            results[39 + jaa](nA, nB) = -results[30 + j](nA, nB) - results[15 + j](nA, nB);
          }
        }
      }
    } else {
      // AB (xx, xy, xz, yx, yy, yz, zx, zy, zz)
      for (int i = 6; i < 15; i++) {
        results[i] = QAA[jaas[i - 6]];
        results[i].multiply(-1.0);
      }
      // BB (xx, xy, xz, yy, yz, zz)
      for (int i = 24; i < 30; i++) results[i] = QAA[i - 24];
    }
  } else if (!onB) {
    // BB (xx, xy, xz, yy, yz, zz)
    for (int i = 24; i < 30; i++) results[i] = QBB[i - 24];
    // AB (xx, xy, xz, yx, yy, yz, zx, zy, zz)
    for (int i = 6; i < 15; i++) {
      results[i] = QBB[jaas[i - 6]];
      results[i].multiply(-1.0);
    }
    // AA (xx, xy, xz, yy, yz, zz)
    for (int i = 0; i < 6; i++) results[i] = QBB[i];
  }
}

std::vector<double> unpack_hessian(const std::vector<double>& packed, const int natoms) {
  const int n3 = 3 * natoms;
  const int diag[3][3] = {{0, 1, 2}, {1, 3, 4}, {2, 4, 5}};  // xx, xy, xz, yy, yz, zz
  std::vector<double> hessian(n3 * n3);
  for (int a = 0; a < natoms; a++)
    for (int b = a; b < natoms; b++)
      for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) {
          const double h = a == b ? packed[H_START(a, a, natoms) + 3 + diag[i][j]]
                                  : packed[H_START(a, b, natoms) + 3 * i + j];
          hessian[(3 * a + i) * n3 + 3 * b + j] = h;
          hessian[(3 * b + j) * n3 + 3 * a + i] = h;
        }
  return hessian;
}

}  // namespace libecpint
