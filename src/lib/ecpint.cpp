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

#include "ecpint.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <set>
#include <stdexcept>
#include <string>

#include "assembly.hpp"
#include "mathutil.hpp"
#include "qgen.hpp"
#include "screening.hpp"

namespace libecpint {

// The largest power of r in the integrals is in the both-on-ECP type 2 term, FAST_POW[N + 1]
// with N = 2 + LA + LB + (n - 2), for ECP terms r^(n - 2) with n <= 4 (ECP::addPrimitive), which
// also indexes GAMMA[N]. Guard against running off the end of the hand-coded power table (entries
// beyond FAST_POW_MAX_EXPONENT are 1/z, 1/z^2, so this would be silently wrong) or GAMMA.
static_assert(5 + 2 * LIBECPINT_MAX_L <= FAST_POW_MAX_EXPONENT,
              "LIBECPINT_MAX_L is too large for the FAST_POW power table; extend mathutil's "
              "pow_* routines.");
static_assert(4 + 2 * LIBECPINT_MAX_L < sizeof(GAMMA) / sizeof(GAMMA[0]),
              "LIBECPINT_MAX_L is too large for the GAMMA table in mathutil.hpp.");

ECPIntegral::ECPIntegral(int maxLB, int maxLU, int deriv, double thresh, unsigned smallGrid,
                         unsigned bigGrid) {
  // Make sure library can perform requested integrals
  if (maxLB + deriv > LIBECPINT_MAX_L)
    throw std::invalid_argument("libecpint: basis angular momentum " + std::to_string(maxLB) +
                                " with derivatives of order " + std::to_string(deriv) +
                                " exceeds the maximum of " + std::to_string(LIBECPINT_MAX_L));
  if (maxLU > LIBECPINT_MAX_L)
    throw std::invalid_argument("libecpint: ECP angular momentum " + std::to_string(maxLU) +
                                " exceeds the maximum of " + std::to_string(LIBECPINT_MAX_L));

  // Initialise singletons
  initFactorials();
  zero = nonzero = skipped = 0;

  // Screening tolerance follows the requested integral threshold so that callers asking
  // for tight integrals are not silently denied small (e.g. high-exponent l>0) contributions.
  tolerance = thresh;

  // Initialise angular and radial integrators
  angInts.init(maxLB + deriv, maxLU);
  angInts.compute();
  type2LB = maxLB + deriv;
  type2LU = maxLU;
  buildType2Radials();
#ifdef DEBUG
  std::cout << "Initializing ECP radial integrator: thresh = " << thresh << "  grids: " << smallGrid
            << "/" << bigGrid << std::endl;
#endif
  radInts.init(2 * (maxLB + deriv) + maxLU, thresh, smallGrid, bigGrid);
};

void ECPIntegral::buildType2Radials() {
  // A radial integral (N, l1, l2) is needed for Q(LA, LB, lam) if, for some power triples alpha
  // of shell A (|alpha| <= LA) and beta of shell B (|beta| <= LB) in the binomial expansions
  // about the ECP, with N = |alpha| + |beta|, the angular factor
  //   sum_mu [sum_mu1 Omega(alpha, lam, mu, l1, mu1)] [sum_mu2 Omega(beta, lam, mu, l2, mu2)]
  // does not vanish (as the code generator, generate.cpp, selected them). The set for
  // (|alpha|, |beta|) = (a, b) is found once, and the set for (LA, LB) is the union over a <= LA,
  // b <= LB.
  const double prefac = 16.0 * M_PI * M_PI;
  const int L = type2LB;
  type2Radials.assign((L + 1) * (L + 1) * type2LU, Type2Radials{});
  std::vector<std::array<int, 3>> powers;  // every power triple with |alpha| <= L
  for (int n = 0; n <= L; n++)
    for (int x = n; x >= 0; x--)
      for (int y = n - x; y >= 0; y--) powers.push_back({x, y, n - x - y});

  for (int lam = 0; lam < type2LU; lam++) {
    const int w_size = 2 * lam + 1;
    // contr[t][l * w_size + lam + mu] = sum_mu1 Omega(powers[t], lam, mu, l, mu1)
    std::vector<std::vector<double>> contr(powers.size());
    for (size_t t = 0; t < powers.size(); t++) {
      const auto& p = powers[t];
      const int lmax = lam + p[0] + p[1] + p[2];
      contr[t].assign(w_size * (lmax + 1), 0.0);
      for (int l = 0; l <= lmax; l++)
        for (int mu = -lam; mu <= lam; mu++) {
          const double* om = angInts.getOmegaRow(p[0], p[1], p[2], lam, mu) + l * l + l;
          double& c = contr[t][l * w_size + lam + mu];
          for (int mu1 = -l; mu1 <= l; mu1++) c += om[mu1];
        }
    }

    std::vector<std::set<Triple>> by_ab((L + 1) * (L + 1));
    for (size_t ta = 0; ta < powers.size(); ta++) {
      const int a = powers[ta][0] + powers[ta][1] + powers[ta][2];
      for (size_t tb = 0; tb < powers.size(); tb++) {
        const int b = powers[tb][0] + powers[tb][1] + powers[tb][2], N = a + b;
        std::set<Triple>& found = by_ab[a * (L + 1) + b];
        for (int l1 = 0; l1 <= lam + a; l1++)
          for (int l2 = (l1 + N) % 2; l2 <= lam + b; l2 += 2) {
            double ang = 0.0;
            for (int mu = -lam; mu <= lam; mu++)
              ang += prefac * contr[ta][l1 * w_size + lam + mu] * contr[tb][l2 * w_size + lam + mu];
            if (std::fabs(ang) > 1e-15) found.insert(Triple{N, l1, l2});
          }
      }
    }

    for (int LB = 0; LB <= L; LB++)
      for (int LA = 0; LA <= LB; LA++) {
        std::set<Triple> needed;
        for (int a = 0; a <= LA; a++)
          for (int b = 0; b <= LB; b++)
            needed.insert(by_ab[a * (L + 1) + b].begin(), by_ab[a * (L + 1) + b].end());
        Type2Radials& r = type2Radials[(LA * (L + 1) + LB) * type2LU + lam];
        if (!needed.empty()) {
          const Triple& last = *needed.crbegin();
          r.nbase = std::max(0, std::get<0>(last) + std::get<1>(last) - 1);
        }
        for (const Triple& t : needed) {
          if (std::get<1>(t) <= std::get<2>(t))
            r.A.push_back(t);
          else
            r.B.push_back(Triple{std::get<0>(t), std::get<2>(t), std::get<1>(t)});
        }
      }
  }
}

double ECPIntegral::calcC(const int a, const int m, const double A) const {
  double value = 1.0 - 2 * ((a - m) % 2);
  value *= std::pow(A, a - m);
  value *= FAC[a] / (FAC[m] * FAC[a - m]);
  return value;
}

void ECPIntegral::makeC(FiveIndex<double>& C, const int L, const double* A) const {
  int z;
  double Ck, Cl;
  int na = 0;
  for (int x = L; x >= 0; x--) {
    for (int y = L - x; y >= 0; y--) {
      z = L - x - y;

      for (int k = 0; k <= x; k++) {
        Ck = calcC(x, k, A[0]);
        for (int l = 0; l <= y; l++) {
          Cl = calcC(y, l, A[1]);
          for (int m = 0; m <= z; m++) C(0, na, k, l, m) = Ck * Cl * calcC(z, m, A[2]);
        }
      }
      na++;
    }
  }
}

void ECPIntegral::type1(const ECP& U, const GaussianShell& shellA, const GaussianShell& shellB,
                        const ShellPairData& data, const FiveIndex<double>& CA,
                        const FiveIndex<double>& CB, const RadialIntegral::Parameters& parameters,
                        TwoIndex<double>& values) const {
  int LA = data.LA;
  int LB = data.LB;
  int maxLBasis = data.maxLBasis;

  // Build radial integrals
  int L = LA + LB;
  ThreeIndex<double> radials;
  radInts.type1(L, U, shellA, shellB, data, parameters, radials);

  // Unpack positions
  double Ax = data.A[0];
  double Ay = data.A[1];
  double Az = data.A[2];
  double Bx = data.B[0];
  double By = data.B[1];
  double Bz = data.B[2];

  // Calculate chi_ab for all ab in shells
  int z1, z2, lparity, mparity, msign, ix, k, l, m;
  double C;
  int na = 0, nb = 0;
  for (int x1 = LA; x1 >= 0; x1--) {
    for (int y1 = LA - x1; y1 >= 0; y1--) {
      z1 = LA - x1 - y1;
      nb = 0;

      for (int x2 = LB; x2 >= 0; x2--) {
        for (int y2 = LB - x2; y2 >= 0; y2--) {
          z2 = LB - x2 - y2;

          for (int k1 = 0; k1 <= x1; k1++) {
            for (int k2 = 0; k2 <= x2; k2++) {
              k = k1 + k2;

              for (int l1 = 0; l1 <= y1; l1++) {
                for (int l2 = 0; l2 <= y2; l2++) {
                  l = l1 + l2;

                  for (int m1 = 0; m1 <= z1; m1++) {
                    for (int m2 = 0; m2 <= z2; m2++) {
                      m = m1 + m2;
                      C = CA(0, na, k1, l1, m1) * CB(0, nb, k2, l2, m2);
                      if (fabs(C) > 1e-14) {
                        // Build radial integrals
                        ix = k + l + m;

                        // Certain terms can be neglected as the angular integrals will always be
                        // zero See Flores06 appendix for details.
                        lparity = ix % 2;
                        msign = 1 - 2 * (l % 2);
                        mparity = (lparity + m) % 2;

                        for (int lam = lparity; lam <= ix; lam += 2) {
                          for (int mu = mparity; mu <= lam; mu += 2)
                            values(na, nb) += C * angInts.getIntegral(k, l, m, lam, msign * mu) *
                                              radials(ix, lam, lam + msign * mu);
                        }
                      }
                    }
                  }
                }
              }
            }
          }

          values(na, nb) *= 4.0 * M_PI;
          nb++;
        }
      }

      na++;
    }
  }
}

void ECPIntegral::type2(const int lam, const ECP& U, const GaussianShell& shellA,
                        const GaussianShell& shellB, const ShellPairData& data,
                        const FiveIndex<double>& CA, const FiveIndex<double>& CB,
                        const RadialIntegral::Parameters& parameters,
                        ThreeIndex<double>& values) const {
  // Unpack some data for convenience
  int LA = data.LA;
  int LB = data.LB;
  int L = LA + LB;
  int maxLBasis = data.maxLBasis;

  double Am = data.Am;
  double Bm = data.Bm;

  if (data.A_on_ecp && data.B_on_ecp) {
    // Both on ECP, simplest case - see Shaw2017 supplementary material
    double prefactor = 4.0 * M_PI;
    int npA = shellA.nprimitive();
    int npB = shellB.nprimitive();
    int npC = U.getN();

    double zA, zB, zC, dA, dB, dC, p;
    int nC, z1, z2;

    // The radial factor here depends only on lam (via g.l) and the total angular
    // momentum LA + LB, not on the individual Cartesian components, so compute it once
    // rather than redundantly inside the (na, nb) loops.
    double value = 0.0;
    for (int c = 0; c < npC; c++) {
      const GaussianECP& g = U.getGaussian(c);
      if (g.l == lam) {
        zC = g.a;
        dC = g.d;
        nC = g.n;

        for (int a = 0; a < npA; a++) {
          zA = shellA.exp(a);
          dA = shellA.coef(a);

          for (int b = 0; b < npB; b++) {
            zB = shellB.exp(b);
            dB = shellB.coef(b);

            p = zA + zB + zC;

            double o_root_p = 1.0 / sqrt(p);
            int N = 2 + LA + LB + nC;
            value += 0.5 * dA * dB * dC * GAMMA[N] * FAST_POW[N + 1](o_root_p);
          }
        }
      }
    }

    int na = 0;
    for (int x1 = LA; x1 >= 0; x1--) {
      for (int r1 = LA - x1; r1 >= 0; r1--) {
        z1 = LA - x1 - r1;

        int nb = 0;
        for (int x2 = LB; x2 >= 0; x2--) {
          for (int y2 = LB - x2; y2 >= 0; y2--) {
            z2 = LB - x2 - y2;

            for (int mu = -lam; mu <= lam; mu++) {
              double angular = prefactor * angInts.getIntegral(x1, r1, z1, lam, mu, 0, 0) *
                               angInts.getIntegral(x2, y2, z2, lam, mu, 0, 0);
              values(na, nb, lam + mu) = angular * value;
            }
            nb++;
          }
        }

        na++;
      }
    }

  } else {
    // At least one of the shells is not on the ECP, so spherical harmonics will be required

    // cos and sin of the polar angles from the components (sin exactly 0 on the z axis)
    double xA = Am > 0 ? data.A[2] / Am : 0.0;
    double xB = Bm > 0 ? data.B[2] / Bm : 0.0;
    double sA = Am > 0 ? std::sqrt(data.A[0] * data.A[0] + data.A[1] * data.A[1]) / Am : 1.0;
    double sB = Bm > 0 ? std::sqrt(data.B[0] * data.B[0] + data.B[1] * data.B[1]) / Bm : 1.0;
    double phiA = atan2(data.A[1], data.A[0]);
    double phiB = atan2(data.B[1], data.B[0]);
    TwoIndex<double> SA, SB;
    realSphericalHarmonics(lam + LA, xA, sA, phiA, SA);
    realSphericalHarmonics(lam + LB, xB, sB, phiB, SB);

    // With one shell on the ECP, its radial factor is sum_a c_a exp(-a r^2) (of the Bessel
    // functions only K_0 = 1 survives at zero distance), so the radials needed, G(N, 0, l2) for
    // l2 = N mod 2, ..., lam + L_off, are integrated as for two off-centre shells, with the
    // on-centre shell first at distance 0
    auto on_centre_radials = [&](const GaussianShell& on, const GaussianShell& off, const int Loff,
                                 const double Doff, ThreeIndex<double>& radials) {
      std::vector<Triple> triples;
      for (int N = 0; N <= L; N++)
        for (int l2 = N % 2; l2 <= lam + Loff; l2 += 2) triples.push_back(Triple{N, 0, l2});
      radials.fill(0.0);
      radInts.type2(triples, 0, lam, U, on, off, 0.0, Doff, radials);
    };

    if (data.A_on_ecp) {
      ThreeIndex<double> radials(L + 1, 1, lam + LB + 1);
      std::fill(values.data.begin(), values.data.end(), 0.0);
      on_centre_radials(shellA, shellB, LB, Bm, radials);

      // a significant number of terms can be neglected a priori - see Shaw2017 supplementary
      // material.
      qgen::rolled_up_special(lam, LA, LB, radials, CB, SB, angInts, values);

    } else if (data.B_on_ecp) {
      // Same as above with A and B reversed
      ThreeIndex<double> radials(L + 1, 1, lam + LA + 1);
      ThreeIndex<double> tmpValues(values.dims[1], values.dims[0], values.dims[2]);
      std::fill(tmpValues.data.begin(), tmpValues.data.end(), 0.0);
      on_centre_radials(shellB, shellA, LA, Am, radials);

      // a significant number of terms can be neglected a priori - see Shaw2017 supplementary
      // material.
      qgen::rolled_up_special(lam, LB, LA, radials, CA, SA, angInts, tmpValues);
      // transcribe back into values
      for (int na = 0; na < values.dims[0]; na++)
        for (int nb = 0; nb < values.dims[1]; nb++)
          for (int nc = 0; nc < values.dims[2]; nc++) values(na, nb, nc) = tmpValues(nb, na, nc);
    } else {
      // Neither is on the ECP: an unrolled kernel for the smallest cases, otherwise the radial
      // integrals that the angular integrals select (buildType2Radials) contracted with them by
      // rolled_up. Both need LA <= LB, but symmetry means we can just swap the arguments if LB >
      // LA.
      auto integrate = [&](const GaussianShell& sA, const GaussianShell& sB, const int la,
                           const int lb, const FiveIndex<double>& cA, const FiveIndex<double>& cB,
                           const TwoIndex<double>& hA, const TwoIndex<double>& hB, const double am,
                           const double bm, ThreeIndex<double>& out) {
        if (const qgen::Kernel kernel = qgen::unrolled(la, lb, lam)) {
          kernel(U, sA, sB, cA, cB, hA, hB, am, bm, radInts, angInts, parameters, out);
          return;
        }
        assert(lb <= type2LB && lam < type2LU);
        const Type2Radials& t = type2Radials[(la * (type2LB + 1) + lb) * type2LU + lam];
        ThreeIndex<double> radials(lam + la + lb + 1, lam + la + 1, lam + lb + 1);
        radInts.type2(t.A, t.nbase, lam, U, sA, sB, am, bm, radials);
        // The radials with l1 > l2 are computed with the shells swapped
        ThreeIndex<double> radials_B(lam + la + lb + 1, lam + lb + 1, lam + la + 1);
        radInts.type2(t.B, t.nbase, lam, U, sB, sA, bm, am, radials_B);
        for (const Triple& r : t.B)
          radials(std::get<0>(r), std::get<2>(r), std::get<1>(r)) =
              radials_B(std::get<0>(r), std::get<1>(r), std::get<2>(r));
        qgen::rolled_up(lam, la, lb, radials, cA, cB, hA, hB, angInts, out);
      };
      if (LA <= LB)
        integrate(shellA, shellB, LA, LB, CA, CB, SA, SB, Am, Bm, values);
      else {
        ThreeIndex<double> temp_values(data.ncartB, data.ncartA, 2 * U.getL() + 1);
        integrate(shellB, shellA, LB, LA, CB, CA, SB, SA, Bm, Am, temp_values);
        for (int na = 0; na < data.ncartA; na++)
          for (int nb = 0; nb < data.ncartB; nb++)
            for (int nu = 0; nu < 2 * U.getL() + 1; nu++)
              values(na, nb, nu) = temp_values(nb, na, nu);
      }
    }
  }
}

void ECPIntegral::estimate_type2(const ECP& U, const GaussianShell& shellA,
                                 const GaussianShell& shellB, const ShellPairData& data,
                                 double* results) const {
  screening::channel_bounds(U, shellA, shellB, data, results);
}

void ECPIntegral::compute_shell_pair(const ECP& U, const GaussianShell& shellA,
                                     const GaussianShell& shellB, TwoIndex<double>& values,
                                     const int shiftA, const int shiftB) const {
  ShellPairData data;

  // Shift A and B to be relative to U
  const double* C = U.center();
  data.A[0] = shellA.center()[0] - C[0];
  data.A[1] = shellA.center()[1] - C[1];
  data.A[2] = shellA.center()[2] - C[2];
  data.B[0] = shellB.center()[0] - C[0];
  data.B[1] = shellB.center()[1] - C[1];
  data.B[2] = shellB.center()[2] - C[2];

  // Construct data that will be reused everywhere, and takes account of derivative shifts
  data.LA = shellA.am() + shiftA;
  data.LB = shellB.am() + shiftB;
  data.maxLBasis = data.LA > data.LB ? data.LA : data.LB;
  data.ncartA = (data.LA + 1) * (data.LA + 2) / 2;
  data.ncartB = (data.LB + 1) * (data.LB + 2) / 2;

  data.A2 = data.A[0] * data.A[0] + data.A[1] * data.A[1] + data.A[2] * data.A[2];
  data.Am = sqrt(data.A2);
  data.A_on_ecp = (data.Am < 1e-6);
  data.B2 = data.B[0] * data.B[0] + data.B[1] * data.B[1] + data.B[2] * data.B[2];
  data.Bm = sqrt(data.B2);
  data.B_on_ecp = (data.Bm < 1e-6);
  double RAB[3] = {data.A[0] - data.B[0], data.A[1] - data.B[1], data.A[2] - data.B[2]};
  data.RAB2 = RAB[0] * RAB[0] + RAB[1] * RAB[1] + RAB[2] * RAB[2];
  data.RABm = sqrt(data.RAB2);

  // Prepare the radial integrator
  const auto radIntParameters = radInts.buildParameters(shellA, shellB, data);

  // Construct coefficients
  FiveIndex<double> CA(1, data.ncartA, data.LA + 1, data.LA + 1, data.LA + 1);
  FiveIndex<double> CB(1, data.ncartB, data.LB + 1, data.LB + 1, data.LB + 1);
  makeC(CA, data.LA, data.A);
  makeC(CB, data.LB, data.B);

  // Channels whose contribution is bounded below CHANNEL_SCREEN * tolerance (estimate_type2) are
  // skipped. The factor allows for the number of blocks, each possibly missing skipped channels,
  // summed into one value of the molecular integrals or derivatives: one per ECP for integrals;
  // for first derivatives QA and QB (and C = -(QA + QB)) from two calls on raised shells each, so
  // ~10; for second derivatives QAA, QAB and QBB from 3-4 calls each, combined with weights up to
  // 4 (e.g. QAA + 2 QAB + QBB), so a few tens. Skipping at the tolerance itself gave second
  // derivatives off by 2e-14 (Ag3). 1e-2 is set by the second derivatives and kept for every order
  // (compute_shell_pair does not know which it serves): the cost is within timing noise.
  constexpr double CHANNEL_SCREEN = 1e-2;
  std::vector<double> screens(U.getL() + 1);
  estimate_type2(U, shellA, shellB, data, screens.data());
  const double channel_tol = CHANNEL_SCREEN * tolerance;

  // Calculate type1 integrals, if necessary
  values.assign(data.ncartA, data.ncartB, 0.0);
  if (!U.noType1() && screens[U.getL()] > channel_tol)
    type1(U, shellA, shellB, data, CA, CB, radIntParameters, values);

  std::vector<int> l_list;
  for (int l = 0; l < U.getL(); l++)
    if (screens[l] > channel_tol) l_list.push_back(l);

  // Now all the type2 integrals
  ThreeIndex<double> t2vals(data.ncartA, data.ncartB, 2 * U.getL() + 1);
  for (int l : l_list) {
    t2vals.fill(0.0);
    type2(l, U, shellA, shellB, data, CA, CB, radIntParameters, t2vals);

    for (int m = -l; m <= l; m++) {
      for (int na = 0; na < data.ncartA; na++) {
        for (int nb = 0; nb < data.ncartB; nb++) {
          values(na, nb) += t2vals(na, nb, l + m);
        }
      }
    }
  }
}

void ECPIntegral::left_shell_derivative(const ECP& U, const GaussianShell& shellA,
                                        const GaussianShell& shellB,
                                        std::array<TwoIndex<double>, 3>& results) const {
  int LA = shellA.am();
  int LB = shellB.am();

  int ncartB = (LB + 1) * (LB + 2) / 2;
  int ncartA = (LA + 1) * (LA + 2) / 2;
  int ncartA_minus = LA * (LA + 1) / 2;
  TwoIndex<double> Q_minus, Q_plus;

  for (auto& r : results) r.assign(ncartA, ncartB, 0.0);

  if (LA != 0) compute_shell_pair(U, shellA, shellB, Q_minus, -1, 0);

  // hack in the exponents to the coefficients
  GaussianShell tempA = shellA.copy();
  for (int i = 0; i < tempA.nprimitive(); i++) tempA.coeffs[i] *= tempA.exps[i];
  compute_shell_pair(U, tempA, shellB, Q_plus, 1, 0);

  // Now compile the derivatives
  if (LA != 0) {
    int nA = 0;
    int nA_minus, nA_plus;
    for (int k = LA; k >= 0; k--) {
      for (int l = LA - k; l >= 0; l--) {
        int m = LA - k - l;

        for (int nB = 0; nB < ncartB; nB++) {
          nA_plus = N_INDEX(l, m);
          nA_minus = std::min(nA_plus, Q_minus.dims[0] - 1);
          results[0](nA, nB) = -k * Q_minus(nA_minus, nB) + 2.0 * Q_plus(nA_plus, nB);

          nA_minus = l > 0 ? N_INDEX(l - 1, m) : 0;
          nA_plus = N_INDEX(l + 1, m);
          results[1](nA, nB) = -l * Q_minus(nA_minus, nB) + 2.0 * Q_plus(nA_plus, nB);

          nA_minus = m > 0 ? N_INDEX(l, m - 1) : 0;
          nA_plus = N_INDEX(l, m + 1);
          results[2](nA, nB) = -m * Q_minus(nA_minus, nB) + 2.0 * Q_plus(nA_plus, nB);
        }
        nA += 1;
      }
    }
  } else {
    for (int nB = 0; nB < ncartB; nB++) {
      results[0](0, nB) = 2.0 * Q_plus(0, nB);
      results[1](0, nB) = 2.0 * Q_plus(1, nB);
      results[2](0, nB) = 2.0 * Q_plus(2, nB);
    }
  }
}

void ECPIntegral::left_shell_second_derivative(const ECP& U, const GaussianShell& shellA,
                                               const GaussianShell& shellB,
                                               std::array<TwoIndex<double>, 6>& results) const {
  int LA = shellA.am();
  int LB = shellB.am();

  int ncartB = (LB + 1) * (LB + 2) / 2;
  int ncartA = (LA + 1) * (LA + 2) / 2;
  int ncartA_minus = std::max(1, (LA - 1) * (LA) / 2);
  TwoIndex<double> Q_minus, Q_plus, Q_0;

  for (auto& r : results) r.assign(ncartA, ncartB, 0.0);

  if (LA > 1)
    compute_shell_pair(U, shellA, shellB, Q_minus, -2, 0);
  else
    Q_minus.assign(ncartA_minus, ncartB, 0.0);

  // hack in the exponents to the coefficients
  GaussianShell tempA = shellA.copy();
  for (int i = 0; i < tempA.nprimitive(); i++) tempA.coeffs[i] *= tempA.exps[i];
  compute_shell_pair(U, tempA, shellB, Q_0, 0, 0);

  // and for the l+2
  for (int i = 0; i < tempA.nprimitive(); i++) tempA.coeffs[i] *= tempA.exps[i];
  compute_shell_pair(U, tempA, shellB, Q_plus, 2, 0);

  // Now compile the derivatives
  int nA = 0;
  int nA_mm, nA_pp, nA_mp, nA_pm;
  for (int k = LA; k >= 0; k--) {
    for (int l = LA - k; l >= 0; l--) {
      int m = LA - k - l;

      for (int nB = 0; nB < ncartB; nB++) {
        nA_mp = nA_pp = N_INDEX(l, m);  // dxx
        nA_mm = std::min(nA_mp, Q_minus.dims[0] - 1);
        results[0](nA, nB) = k * (k - 1) * Q_minus(nA_mm, nB) - 2.0 * (2 * k + 1) * Q_0(nA_mp, nB) +
                             4.0 * Q_plus(nA_pp, nB);

        nA_pm = l > 0 ? N_INDEX(l - 1, m) : 0;
        nA_mm = k > 0 ? nA_pm : 0;  // dxy
        nA_pp = N_INDEX(l + 1, m);
        nA_mp = k > 0 ? nA_pp : 0;
        results[1](nA, nB) = k * l * Q_minus(nA_mm, nB) - 2.0 * k * Q_0(nA_mp, nB) -
                             2.0 * l * Q_0(nA_pm, nB) + 4.0 * Q_plus(nA_pp, nB);

        nA_pm = m > 0 ? N_INDEX(l, m - 1) : 0;
        nA_mm = k > 0 ? nA_pm : 0;  // dxz
        nA_pp = N_INDEX(l, m + 1);
        nA_mp = k > 0 ? nA_pp : 0;
        results[2](nA, nB) = k * m * Q_minus(nA_mm, nB) - 2.0 * k * Q_0(nA_mp, nB) -
                             2.0 * m * Q_0(nA_pm, nB) + 4.0 * Q_plus(nA_pp, nB);

        nA_mm = l > 1 ? N_INDEX(l - 2, m) : 0;  // dyy
        nA_mp = N_INDEX(l, m);
        nA_pp = N_INDEX(l + 2, m);
        results[3](nA, nB) = l * (l - 1) * Q_minus(nA_mm, nB) - 2.0 * (2 * l + 1) * Q_0(nA_mp, nB) +
                             4.0 * Q_plus(nA_pp, nB);

        nA_mm = l * m > 0 ? N_INDEX(l - 1, m - 1) : 0;  // dyz
        nA_mp = l > 0 ? N_INDEX(l - 1, m + 1) : 0;
        nA_pm = m > 0 ? N_INDEX(l + 1, m - 1) : 0;
        nA_pp = N_INDEX(l + 1, m + 1);
        results[4](nA, nB) = l * m * Q_minus(nA_mm, nB) - 2.0 * l * Q_0(nA_mp, nB) -
                             2.0 * m * Q_0(nA_pm, nB) + 4.0 * Q_plus(nA_pp, nB);

        nA_mm = m > 1 ? N_INDEX(l, m - 2) : 0;  // dzz
        nA_mp = N_INDEX(l, m);
        nA_pp = N_INDEX(l, m + 2);
        results[5](nA, nB) = m * (m - 1) * Q_minus(nA_mm, nB) - 2.0 * (2 * m + 1) * Q_0(nA_mp, nB) +
                             4.0 * Q_plus(nA_pp, nB);
      }
      nA += 1;
    }
  }
}

void ECPIntegral::mixed_second_derivative(const ECP& U, const GaussianShell& shellA,
                                          const GaussianShell& shellB,
                                          std::array<TwoIndex<double>, 9>& results) const {
  int LA = shellA.am();
  int LB = shellB.am();

  int ncartB = (LB + 1) * (LB + 2) / 2;
  int ncartA = (LA + 1) * (LA + 2) / 2;
  int ncartB_minus = std::max(1, (LB) * (LB + 1) / 2);
  int ncartA_minus = std::max(1, (LA) * (LA + 1) / 2);
  int ncartB_plus = (LB + 2) * (LB + 3) / 2;
  int ncartA_plus = (LA + 2) * (LA + 3) / 2;
  TwoIndex<double> Q_mm, Q_mp, Q_pm, Q_pp;

  for (auto& r : results) r.assign(ncartA, ncartB, 0.0);

  GaussianShell tempA = shellA.copy();
  for (int i = 0; i < tempA.nprimitive(); i++) tempA.coeffs[i] *= tempA.exps[i];
  GaussianShell tempB = shellB.copy();
  for (int i = 0; i < tempB.nprimitive(); i++) tempB.coeffs[i] *= tempB.exps[i];

  if (LA > 0) {
    if (LB > 0) {
      compute_shell_pair(U, shellA, shellB, Q_mm, -1, -1);
      compute_shell_pair(U, tempA, shellB, Q_pm, 1, -1);
    } else {
      Q_mm.assign(ncartA_minus, ncartB_minus, 0.0);
      Q_pm.assign(ncartA_plus, ncartB_minus, 0.0);
    }
    compute_shell_pair(U, shellA, tempB, Q_mp, -1, 1);
  } else if (LB > 0) {
    compute_shell_pair(U, tempA, shellB, Q_pm, 1, -1);
    Q_mm.assign(ncartA_minus, ncartB_minus, 0.0);
    Q_mp.assign(ncartA_minus, ncartB_plus, 0.0);
  } else {
    Q_mm.assign(ncartA_minus, ncartB_minus, 0.0);
    Q_mp.assign(ncartA_minus, ncartB_plus, 0.0);
    Q_pm.assign(ncartA_plus, ncartB_minus, 0.0);
  }
  compute_shell_pair(U, tempA, tempB, Q_pp, 1, 1);

  // Now compile the derivatives
  int nA = 0;
  int nB = 0;
  int nA_m[3], nA_p[3], nB_m[3], nB_p[3], AL[3], BL[3];
  for (int ka = LA; ka >= 0; ka--) {
    for (int la = LA - ka; la >= 0; la--) {
      int ma = LA - ka - la;
      AL[0] = ka;
      AL[1] = la;
      AL[2] = ma;
      nA_p[0] = N_INDEX(la, ma);
      nA_m[0] = std::min(nA_p[0], Q_mm.dims[0] - 1);
      nA_m[1] = la > 0 ? N_INDEX(la - 1, ma) : 0;
      nA_m[2] = ma > 0 ? N_INDEX(la, ma - 1) : 0;
      nA_p[1] = N_INDEX(la + 1, ma);
      nA_p[2] = N_INDEX(la, ma + 1);

      nB = 0;
      for (int kb = LB; kb >= 0; kb--) {
        for (int lb = LB - kb; lb >= 0; lb--) {
          int mb = LB - kb - lb;
          nB_p[0] = N_INDEX(lb, mb);
          nB_m[0] = std::min(nB_p[0], Q_mm.dims[1] - 1);
          nB_m[1] = lb > 0 ? N_INDEX(lb - 1, mb) : 0;
          nB_m[2] = mb > 0 ? N_INDEX(lb, mb - 1) : 0;
          nB_p[1] = N_INDEX(lb + 1, mb);
          nB_p[2] = N_INDEX(lb, mb + 1);
          BL[0] = kb;
          BL[1] = lb;
          BL[2] = mb;

          for (int p = 0; p < 3; p++) {
            for (int q = 0; q < 3; q++) {
              results[3 * p + q](nA, nB) =
                  AL[p] * BL[q] * Q_mm(nA_m[p], nB_m[q]) - 2.0 * BL[q] * Q_pm(nA_p[p], nB_m[q]) -
                  2.0 * AL[p] * Q_mp(nA_m[p], nB_p[q]) + 4.0 * Q_pp(nA_p[p], nB_p[q]);
            }
          }

          nB += 1;
        }
      }
      nA += 1;
    }
  }
}

void ECPIntegral::compute_shell_pair_derivative(const ECP& U, const GaussianShell& shellA,
                                                const GaussianShell& shellB,
                                                std::array<TwoIndex<double>, 9>& results) const {
  const bool onA = on_ecp(shellA, U), onB = on_ecp(shellB, U);

  // Calculate shell derivatives (QB as d/dB <A|U|B>, from d/dB <B|U|A>)
  std::array<TwoIndex<double>, 3> QA, QB;
  if (!onA) left_shell_derivative(U, shellA, shellB, QA);
  if (!onB) {
    left_shell_derivative(U, shellB, shellA, QB);
    for (auto& q : QB) q = q.transpose();
  }

  pack_first_derivatives(QA, QB, onA, onB, shellA.ncartesian(), shellB.ncartesian(), results);
}

void ECPIntegral::compute_shell_pair_derivative_Bfield(
    const ECP& U, const GaussianShell& shellA, const GaussianShell& shellB,
    std::array<TwoIndex<double>, 3>& results) const {
  const int LA = shellA.am(), LB = shellB.am();
  for (auto& r : results) r.assign(shellA.ncartesian(), shellB.ncartesian(), 0.0);

  // (A - C) x r = A x C + (A - C) x (r - A), and r - A raises shell A (likewise for B), so
  //   results[k] = [(A - B) x C]_k <A|U|B> + [(A - C) x <A+1|]_k U|B> - <A|U [(B - C) x |B+1>]_k.
  // The term of a shell on the ECP vanishes. When A and B share a centre off the ECP the block is
  // not zero: the bra and ket factors act on different angular coordinates of the semi-local
  // projectors, so they do not cancel.
  const bool onA = on_ecp(shellA, U), onB = on_ecp(shellB, U);
  if (onA && onB) return;

  const double *A = shellA.center(), *B = shellB.center(), *C = U.center();
  const double AmC[3] = {A[0] - C[0], A[1] - C[1], A[2] - C[2]};
  const double BmC[3] = {B[0] - C[0], B[1] - C[1], B[2] - C[2]};
  const double AmB[3] = {A[0] - B[0], A[1] - B[1], A[2] - B[2]};
  const double ABxC[3] = {AmB[1] * C[2] - AmB[2] * C[1], AmB[2] * C[0] - AmB[0] * C[2],
                          AmB[0] * C[1] - AmB[1] * C[0]};

  TwoIndex<double> Q, QAp, QBp;
  compute_shell_pair(U, shellA, shellB, Q);
  if (!onA) compute_shell_pair(U, shellA, shellB, QAp, 1, 0);
  if (!onB) compute_shell_pair(U, shellA, shellB, QBp, 0, 1);

  int nA = 0;
  for (int kA = LA; kA >= 0; kA--) {
    for (int lA = LA - kA; lA >= 0; lA--, nA++) {
      const int mA = LA - kA - lA;
      // index of the component raised in x, y, z in the shell of angular momentum LA + 1
      const int rA[3] = {N_INDEX(lA, mA), N_INDEX(lA + 1, mA), N_INDEX(lA, mA + 1)};
      int nB = 0;
      for (int kB = LB; kB >= 0; kB--) {
        for (int lB = LB - kB; lB >= 0; lB--, nB++) {
          const int mB = LB - kB - lB;
          const int rB[3] = {N_INDEX(lB, mB), N_INDEX(lB + 1, mB), N_INDEX(lB, mB + 1)};
          for (int k = 0; k < 3; k++) {
            const int p = (k + 1) % 3, q = (k + 2) % 3;  // (D x V)_k = D_p V_q - D_q V_p
            double v = ABxC[k] * Q(nA, nB);
            if (!onA) v += AmC[p] * QAp(rA[q], nB) - AmC[q] * QAp(rA[p], nB);
            if (!onB) v -= BmC[p] * QBp(nA, rB[q]) - BmC[q] * QBp(nA, rB[p]);
            results[k](nA, nB) = v;
          }
        }
      }
    }
  }
}

void ECPIntegral::compute_shell_pair_second_derivative(
    const ECP& U, const GaussianShell& shellA, const GaussianShell& shellB,
    std::array<TwoIndex<double>, 45>& results) const {
  const bool onA = on_ecp(shellA, U), onB = on_ecp(shellB, U);

  // Calculate shell derivatives (QBB as d2/dB2 <A|U|B>, from d2/dB2 <B|U|A>)
  std::array<TwoIndex<double>, 6> QAA, QBB;
  std::array<TwoIndex<double>, 9> QAB;
  if (!onA) left_shell_second_derivative(U, shellA, shellB, QAA);
  if (!onB) {
    left_shell_second_derivative(U, shellB, shellA, QBB);
    for (auto& q : QBB) q = q.transpose();
  }
  if (!onA && !onB) mixed_second_derivative(U, shellA, shellB, QAB);

  pack_second_derivatives(QAA, QAB, QBB, onA, onB, shellA.ncartesian(), shellB.ncartesian(),
                          results);
}

}  // namespace libecpint
