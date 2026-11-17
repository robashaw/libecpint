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

#include "qgen.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>

#include "mathutil.hpp"

namespace libecpint {
namespace qgen {
void rolled_up(const int lam, const int LA, const int LB, const ThreeIndex<double>& radials,
               const FiveIndex<double>& CA, const FiveIndex<double>& CB, const TwoIndex<double>& SA,
               const TwoIndex<double>& SB, const AngularIntegral& angint,
               ThreeIndex<double>& values) {
  double prefac = 16.0 * M_PI * M_PI;
  int L = LA + LB;

  int z1, z2;
  int w_l1, w_l2;
  double C, val;
  int w_size = 2 * lam + 1;

  // The angular integrals contracted with the spherical harmonics, w1_contr for shell A and
  // w2_contr for shell B, depend only on the (alpha_x, alpha_y, alpha_z) or (beta_x, beta_y,
  // beta_z) triple, so they are tabulated once for every triple rather than per shell function.
  // By inversion symmetry the angular integrals vanish (exactly, as stored) unless lam + l +
  // |alpha| is even, so only those l are computed; the others are left zero and not read below.
  auto tabulate = [&](const int Lmax, const TwoIndex<double>& S, std::vector<double>& tab) {
    const int stride = w_size * (lam + Lmax + 1);
    tab.assign((Lmax + 1) * (Lmax + 1) * (Lmax + 1) * stride, 0.0);
    for (int ax = 0; ax <= Lmax; ax++)
      for (int ay = 0; ax + ay <= Lmax; ay++)
        for (int az = 0; ax + ay + az <= Lmax; az++) {
          double* contr = &tab[((ax * (Lmax + 1) + ay) * (Lmax + 1) + az) * stride];
          for (int l = (lam + ax + ay + az) % 2; l <= lam + ax + ay + az; l += 2) {
            const int w_l = l * w_size + lam;
            for (int mu = -lam; mu <= lam; mu++) {
              // om[mu1] = Omega(alpha, lam, mu, l, mu1)
              const double* om = angint.getOmegaRow(ax, ay, az, lam, mu) + l * l + l;
              const double* S_l = &S.data[l * S.dims[1] + l];  // S_l[m] = S(l, l + m)
              // Reflections x -> -x and y -> -y: real harmonic S_lm has parity |m| + [m < 0] and
              // [m < 0] respectively, so the angular integral vanishes (exactly, as stored) unless
              // [mu1 < 0] = ay + [mu < 0] and mu1 = ax + |mu| + [mu < 0] + [mu1 < 0] (mod 2)
              const int neg = (ay + (mu < 0)) % 2, par = (ax + std::abs(mu) + (mu < 0) + neg) % 2;
              const int lo = neg ? -l : 0, hi = neg ? -1 : l;
              contr[w_l + mu] = 0.0;
              for (int mu1 = lo + ((lo + par) & 1); mu1 <= hi; mu1 += 2)
                contr[w_l + mu] += S_l[mu1] * om[mu1];
            }
          }
        }
    return stride;
  };
  std::vector<double> w1_tab, w2_tab;
  const int w1_stride = tabulate(LA, SA, w1_tab);
  const int w2_stride = tabulate(LB, SB, w2_tab);

  // G(alpha, beta, mu) = prefac sum_{lam1, lam2} radials(|alpha| + |beta|, lam1, lam2)
  //                      w1_contr(alpha)[lam1, mu] w2_contr(beta)[lam2, mu]
  // depends only on the triples alpha and beta, not on the shell functions they come from, so it
  // is computed once for each pair of triples (|alpha| <= LA, |beta| <= LB), and then
  //   values(na, nb, mu) = sum_{alpha <= a(na), beta <= b(nb)} CA(na, alpha) CB(nb, beta) G
  auto compact = [](const int Lmax, std::vector<int>& index) {  // box index -> compact index
    index.assign((Lmax + 1) * (Lmax + 1) * (Lmax + 1), -1);
    int n = 0;
    for (int x = 0; x <= Lmax; x++)
      for (int y = 0; x + y <= Lmax; y++)
        for (int z = 0; x + y + z <= Lmax; z++) index[(x * (Lmax + 1) + y) * (Lmax + 1) + z] = n++;
    return n;
  };
  std::vector<int> ia_of, ib_of;
  const int nta = compact(LA, ia_of);
  const int ntb = compact(LB, ib_of);
  std::vector<double> G(nta * ntb * w_size, 0.0);
  for (int ax = 0; ax <= LA; ax++)
    for (int ay = 0; ax + ay <= LA; ay++)
      for (int az = 0; ax + ay + az <= LA; az++) {
        const int alpha = ax + ay + az, boxa = (ax * (LA + 1) + ay) * (LA + 1) + az;
        const double* w1_contr = &w1_tab[boxa * w1_stride];
        for (int bx = 0; bx <= LB; bx++)
          for (int by = 0; bx + by <= LB; by++)
            for (int bz = 0; bx + by + bz <= LB; bz++) {
              const int beta = bx + by + bz, boxb = (bx * (LB + 1) + by) * (LB + 1) + bz;
              const double* w2_contr = &w2_tab[boxb * w2_stride];
              double* g = &G[(ia_of[boxa] * ntb + ib_of[boxb]) * w_size + lam];
              const int N = alpha + beta;
              // lam1 and lam2 of the parity of lam + alpha and lam + beta, as w1_contr and
              // w2_contr vanish otherwise
              for (int lam1 = (lam + alpha) % 2; lam1 <= lam + alpha; lam1 += 2) {
                w_l1 = lam1 * w_size + lam;
                for (int lam2 = (lam1 + N) % 2; lam2 <= lam + beta; lam2 += 2) {
                  w_l2 = lam2 * w_size + lam;
                  val = prefac * radials(N, lam1, lam2);
                  for (int mu = -lam; mu <= lam; mu++)
                    g[mu] += val * w1_contr[w_l1 + mu] * w2_contr[w_l2 + mu];
                }
              }
            }
      }

  // Loop over cartesian shell functions in alpha order, e.g. {xx xy, xz, yy, yz, zz} for l=2
  int na = 0;  // Rows are shellA
  for (int x1 = LA; x1 >= 0; x1--) {
    for (int r1 = LA - x1; r1 >= 0; r1--) {
      z1 = LA - x1 - r1;

      int nb = 0;  // Cols are shellB
      for (int x2 = LB; x2 >= 0; x2--) {
        for (int y2 = LB - x2; y2 >= 0; y2--) {
          z2 = LB - x2 - y2;

          for (int alpha_x = 0; alpha_x <= x1; alpha_x++) {
            for (int alpha_y = 0; alpha_y <= r1; alpha_y++) {
              for (int alpha_z = 0; alpha_z <= z1; alpha_z++) {
                const double CA_val = CA(0, na, alpha_x, alpha_y, alpha_z);
                if (CA_val == 0.0) continue;
                const double* ga =
                    &G[ia_of[(alpha_x * (LA + 1) + alpha_y) * (LA + 1) + alpha_z] * ntb * w_size];

                for (int beta_x = 0; beta_x <= x2; beta_x++) {
                  for (int beta_y = 0; beta_y <= y2; beta_y++) {
                    for (int beta_z = 0; beta_z <= z2; beta_z++) {
                      C = CA_val * CB(0, nb, beta_x, beta_y, beta_z);
                      if (C == 0.0) continue;
                      const double* g =
                          &ga[ib_of[(beta_x * (LB + 1) + beta_y) * (LB + 1) + beta_z] * w_size];
                      for (int mu = 0; mu < w_size; mu++) values(na, nb, mu) += C * g[mu];
                    }
                  }
                }
              }
            }
          }

          nb++;
        }
      }

      na++;
    }
  }
}

void rolled_up_special(const int lam, const int LA, const int LB, const ThreeIndex<double>& radials,
                       const FiveIndex<double>& CB, const TwoIndex<double>& SB,
                       const AngularIntegral& angint, ThreeIndex<double>& values) {
  double prefac = 8.0 * M_PI * std::sqrt(M_PI);
  int L = LA + LB;

  int z1, z2;
  double C, val1, val2;

  // Loop over cartesian shell functions in alpha order, e.g. {xx xy, xz, yy, yz, zz} for l=2
  int na = 0;  // Rows are shellA
  for (int x1 = LA; x1 >= 0; x1--) {
    for (int r1 = LA - x1; r1 >= 0; r1--) {
      z1 = LA - x1 - r1;

      int nb = 0;  // Cols are shellB
      for (int x2 = LB; x2 >= 0; x2--) {
        for (int y2 = LB - x2; y2 >= 0; y2--) {
          z2 = LB - x2 - y2;

          // Begin full ECP integral expansion
          int alpha = x1 + r1 + z1;

          for (int beta_x = 0; beta_x <= x2; beta_x++) {
            for (int beta_y = 0; beta_y <= y2; beta_y++) {
              for (int beta_z = 0; beta_z <= z2; beta_z++) {
                int beta = beta_x + beta_y + beta_z;
                int N = alpha + beta;
                C = CB(0, nb, beta_x, beta_y, beta_z);

                if (std::abs(C) > 1e-15) {
                  int lam2start = N % 2;
                  for (int lam2 = lam2start; lam2 <= lam + beta; lam2 += 2) {
                    val1 = prefac * C * radials(N, 0, lam2);

                    for (int mu2 = -lam2; mu2 <= lam2; mu2++) {
                      val2 = val1 * SB(lam2, lam2 + mu2);
                      for (int mu = -lam; mu <= lam; mu++) {
                        values(na, nb, lam + mu) +=
                            val2 * angint.getIntegral(x1, r1, z1, lam, mu, 0, 0) *
                            angint.getIntegral(beta_x, beta_y, beta_z, lam, mu, lam2, mu2);
                      }
                    }
                  }
                }
              }
            }
          }
          nb++;
        }
      }
      na++;
    }
  }
}
}  // namespace qgen
}  // namespace libecpint
