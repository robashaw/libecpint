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

#include "radial.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <vector>

#include "gaussquad.hpp"
#include "mathutil.hpp"
#include "screening.hpp"
#include "shared_grid.hpp"

namespace libecpint {

RadialIntegral::RadialIntegral() {}

void RadialIntegral::init(int maxL, double tol, int /*small*/, int /*large*/) {
  // The Bessel function table is always converged to BESSEL_ACCURACY, independently of tol: tol is
  // a screening threshold (0 to skip nothing), and the series overflows if converged to 0
  bessie.init(maxL, 1600, 200, BESSEL_ACCURACY);

  tolerance = tol;
}

double RadialIntegral::calcKij(const double Na, const double Nb, const double zeta_a,
                               const double zeta_b, const double R2) const {
  double muij = zeta_a * zeta_b / (zeta_a + zeta_b);
  return Na * Nb * std::exp(-muij * R2);
}

RadialIntegral::Parameters RadialIntegral::buildParameters(const GaussianShell& shellA,
                                                           const GaussianShell& shellB,
                                                           const ShellPairData& data) const {
  int npA = shellA.nprimitive();
  int npB = shellB.nprimitive();

  // Initialise result object
  Parameters result;
  auto& p = result.p;
  auto& P = result.P;
  auto& P2 = result.P2;
  auto& K = result.K;

  p.assign(npA, npB, 0.0);
  P.assign(npA, npB, 0.0);
  P2.assign(npA, npB, 0.0);
  K.assign(npA, npB, 0.0);

  double Pvec[3];
  double zetaA, zetaB;
  for (int a = 0; a < npA; a++) {
    zetaA = shellA.exp(a);

    for (int b = 0; b < npB; b++) {
      zetaB = shellB.exp(b);

      p(a, b) = zetaA + zetaB;
      for (int n = 0; n < 3; n++) Pvec[n] = (zetaA * data.A[n] + zetaB * data.B[n]) / p(a, b);

      P2(a, b) = Pvec[0] * Pvec[0] + Pvec[1] * Pvec[1] + Pvec[2] * Pvec[2];
      P(a, b) = std::sqrt(P2(a, b));
      K(a, b) = calcKij(1.0, 1.0, zetaA, zetaB, data.RAB2);
    }
  }
  return result;
}

using namespace screening;
using namespace detail;

namespace {
// Nodes r and weights w of the composite Gauss-Legendre rule of the shared grid on a window of an
// integrand with Gaussian exponent zt: equal panels no wider than SG_PANEL_WIDTH / sqrt(zt), each
// with SG_PANEL_POINTS points, which integrate r^k exp(-zt (r - p)^2) to ~1e-14 relative for
// k <= 40 (see shared_grid.hpp). A single panel over the window (>= 16 / sqrt(zt) wide) gave
// radials off by up to 6e-13 relative, which the angular sums amplify for shells far from the ECP.
// Returns the number of panels.
int window_rule(const Window& window, const double zt, std::vector<double>& r,
                std::vector<double>& w) {
  const auto& gl = gauss_legendre<SG_PANEL_POINTS>();
  const double length = window.rmax - window.rmin;
  const int np = std::max(1, static_cast<int>(std::ceil(length * std::sqrt(zt) / SG_PANEL_WIDTH)));
  const double half = 0.5 * length / np;
  r.resize(np * SG_PANEL_POINTS);
  w.resize(np * SG_PANEL_POINTS);
  for (int p = 0; p < np; p++)
    for (int g = 0; g < SG_PANEL_POINTS; g++) {
      r[p * SG_PANEL_POINTS + g] = window.rmin + half * (2 * p + 1 + gl.x[g]);
      w[p * SG_PANEL_POINTS + g] = half * gl.w[g];
    }
  return np;
}
}  // namespace

void RadialIntegral::type2(const std::vector<Triple>& triples, const int /*nbase*/, const int lam,
                           const ECP& U, const GaussianShell& shellA, const GaussianShell& shellB,
                           const double A, const double B, ThreeIndex<double>& radials) const {
  // All radials are integrated numerically: primitive triples with both basis exponents <=
  // SG_TIGHT on the shared grid, the others each with the same composite rule (window_rule) on a
  // window bracketing the integrand. The closed forms in terms of base integrals (Shaw2017) that
  // were used where available lost accuracy to cancellation, for diffuse primitives and also for a
  // tight primitive paired with a moderate one (up to 2% relative).
  if (triples.empty()) return;
  const int npA = shellA.nprimitive();
  const int npB = shellB.nprimitive();
  int lmaxA = 0, lmaxB = 0, Nmin = std::get<0>(triples[0]), Nmax = Nmin;
  for (const Triple& t : triples) {
    lmaxA = std::max(lmaxA, std::get<1>(t));
    lmaxB = std::max(lmaxB, std::get<2>(t));
    Nmin = std::min(Nmin, std::get<0>(t));
    Nmax = std::max(Nmax, std::get<0>(t));
  }
  // Per primitive triple, wrkq holds weight * exponential factor * r^k at each node rq, for
  // k - u.n - 2 in [Nmin, Nmax], and Kqa, Kqb the scaled Bessel functions of each order.
  const int nk = Nmax - Nmin + 1;
  // Scratch, per thread (compute_shell_pair may run concurrently) and reused across calls; every
  // element read is written first in this call
  thread_local std::vector<double> rq, wq, wrkq, Kqa, Kqb, bes;
  thread_local std::vector<size_t> nlast;
  bes.resize(std::max(lmaxA, lmaxB) + 1);

  // Threshold for each primitive triple's bound (Nmax = LA + LB)
  const double thresh = triple_threshold(tolerance, A, B, Nmax, npA * npB * count_channel(U, lam));

  // |coef| type2_bound < thresh, tested first with the cheap bound and the largest |coef|
  double cA = 0.0, cB = 0.0, cU = 0.0;
  for (int na = 0; na < npA; na++) cA = std::max(cA, std::abs(shellA.coef(na)));
  for (int nb = 0; nb < npB; nb++) cB = std::max(cB, std::abs(shellB.coef(nb)));
  for (const auto& u : U.gaussians)
    if (u.l == lam) cU = std::max(cU, std::abs(u.d));
  const double log_cmax = std::log(cA * cB * cU / thresh);
  auto negligible = [&](const double coef, const int kmin, const int kmax, const double n,
                        const double a, const double b) {
    return log_cmax + type2_log_bound_cheap(kmax, n, a, b, A, B) < -1e-6 ||
           std::abs(coef) * type2_bound(kmin, kmax, n, a, b, A, B) < thresh;
  };

  // Primitive triples with both basis exponents <= SG_TIGHT, on the shared grid, if any of them is
  // not provably negligible (most diffuse first, as those usually are not)
  bool shared = false;
  double sumA = 0.0, sumB = 0.0;  // sum of |c| over those primitives, bounding |F|
  for (int na = 0; na < npA; na++)
    if (shellA.exp(na) <= SG_TIGHT) sumA += std::abs(shellA.coef(na));
  for (int nb = 0; nb < npB; nb++)
    if (shellB.exp(nb) <= SG_TIGHT) sumB += std::abs(shellB.coef(nb));
  for (const auto& u : U.gaussians) {
    if (u.l != lam || shared) continue;
    for (int na = npA - 1; na >= 0 && !shared; na--) {
      if (shellA.exp(na) > SG_TIGHT) continue;
      for (int nb = npB - 1; nb >= 0 && !shared; nb--) {
        if (shellB.exp(nb) > SG_TIGHT) continue;
        shared = !negligible(shellA.coef(na) * shellB.coef(nb) * u.d, Nmin + u.n + 2,
                             Nmax + u.n + 2, u.a, shellA.exp(na), shellB.exp(nb));
      }
    }
  }
  if (shared) {
    const auto grid = shared_grid(U);
    const auto FA = shared_factor(*grid, shellA, A, lmaxA, bessie);
    const auto FB = shared_factor(*grid, shellB, B, lmaxB, bessie);
    const auto WN = shared_weights(*grid, lam, Nmax);
    const size_t n = grid->r.size();
    // The number of nodes before the integrand (bounded by |W_N| sumA sumB) becomes negligible
    nlast.resize(nk);
    for (int k = 0; k < nk; k++) {
      const double* S = &WN->S[(Nmin + k) * (n + 1)];
      nlast[k] = std::partition_point(
                     S, S + n, [&](const double s) { return s * sumA * sumB >= 1e-6 * thresh; }) -
                 S;
    }
    for (const Triple& t : triples) {
      const int N = std::get<0>(t), l1 = std::get<1>(t), l2 = std::get<2>(t);
      const double* w = &WN->W[N * n];
      const double* Fa = &FA->F[l1 * n];
      const double* Fb = &FB->F[l2 * n];
      // Four partial sums, so that the loop vectorises (the order is fixed, so results are
      // still deterministic)
      const size_t m = nlast[N - Nmin];
      double s0 = 0.0, s1 = 0.0, s2 = 0.0, s3 = 0.0;
      size_t i = 0;
      for (; i + 4 <= m; i += 4) {
        s0 += w[i] * Fa[i] * Fb[i];
        s1 += w[i + 1] * Fa[i + 1] * Fb[i + 1];
        s2 += w[i + 2] * Fa[i + 2] * Fb[i + 2];
        s3 += w[i + 3] * Fa[i + 3] * Fb[i + 3];
      }
      for (; i < m; i++) s0 += w[i] * Fa[i] * Fb[i];
      radials(N, l1, l2) += (s0 + s1) + (s2 + s3);
    }
  }

  // Remaining primitive triples, each on its own window: those with a basis exponent > SG_TIGHT.
  // Loop over primitives in ECP, only considering correct ang. momentum
  thread_local std::vector<int> tightB;
  tightB.clear();
  for (int nb = 0; nb < npB; nb++)
    if (shellB.exp(nb) > SG_TIGHT) tightB.push_back(nb);
  for (const auto& u : U.gaussians) {
    if (u.l != lam) continue;
    const int kmin = Nmin + u.n + 2, kmax = Nmax + u.n + 2;
    // Loop over primitives in orbital basis shells: every nb if a is tight, else the tight nb
    for (int na = 0; na < npA; na++) {
      const double a = shellA.exp(na);
      const double da = shellA.coef(na);
      const bool all = a > SG_TIGHT;
      const int countB = all ? npB : static_cast<int>(tightB.size());
      for (int ib = 0; ib < countB; ib++) {
        const int nb = all ? ib : tightB[ib];
        const double b = shellB.exp(nb);
        const double db = shellB.coef(nb);
        const double coef = da * db * u.d;
        // Skip primitive triples whose every radial is provably negligible (e.g. a tight
        // primitive far from the ECP)
        if (negligible(coef, kmin, kmax, u.a, a, b)) continue;

        const double xA = 2.0 * a * A, xB = 2.0 * b * B;  // Bessel arguments are xA r, xB r
        const Window window = integrand_window(kmin, kmax, lmaxA, lmaxB, u.a, a, b, A, B);
        const int npanel = window_rule(window, u.a + a + b, rq, wq);
        const int nq = rq.size();
        wrkq.resize(nk * nq);
        Kqa.resize((lmaxA + 1) * nq);
        Kqb.resize((lmaxB + 1) * nq);
        for (int g = 0; g < nq; g++) {
          const double r = rq[g], rA = r - A, rB = r - B;
          double wrk =
              wq[g] * std::exp(-u.a * r * r - a * rA * rA - b * rB * rB) * FAST_POW[kmin](r);
          for (int n = 0; n < nk; n++, wrk *= r) wrkq[n * nq + g] = wrk;
          // A shell on the ECP (x = 0) has K_0 = 1 and K_l = 0 for l > 0; calculate only writes
          // K_0 then, so the other orders are set here
          auto bessel = [&](const double x, const int lmax, double* K) {
            if (x > 0.0) {
              bessie.calculate(x * r, lmax, bes);
              for (int l = 0; l <= lmax; l++) K[l * nq + g] = bes[l];
            } else
              for (int l = 0; l <= lmax; l++) K[l * nq + g] = l == 0 ? 1.0 : 0.0;
          };
          bessel(xA, lmaxA, Kqa.data());
          bessel(xB, lmaxB, Kqb.data());
        }

        for (const Triple& t : triples) {
          const int N = std::get<0>(t), l1 = std::get<1>(t), l2 = std::get<2>(t);
          const double* wrk = &wrkq[(N - Nmin) * nq];
          const double* Ka = &Kqa[l1 * nq];
          const double* Kb = &Kqb[l2 * nq];
          double sum = 0.0;
          for (int g = 0; g < nq; g++) sum += wrk[g] * Ka[g] * Kb[g];
#ifndef NDEBUG
          // Debug builds: the (positive) integrand must be negligible at the window ends, and the
          // rule must agree with a higher-order one on the same panels. Both are measured
          // against the radial with l1 = l2 = 0 and the same k, which bounds all l (K_l <= K_0):
          // relative to itself, a tiny high-l radial carries the absolute error of the tabulated
          // Bessel functions at small argument.
          {
            const int k = N + u.n + 2;
            double scale = 0.0;
            for (int g = 0; g < nq; g++) scale += wrk[g] * Kqa[g] * Kqb[g];
            auto f = [&](const double r) {
              const double rA = r - A, rB = r - B;
              return FAST_POW[k](r) * std::exp(-u.a * r * r - a * rA * rA - b * rB * rB) *
                     (xA > 0.0 ? bessie.calculate(xA * r, l1) : double(l1 == 0)) *
                     (xB > 0.0 ? bessie.calculate(xB * r, l2) : double(l2 == 0));
            };
            const auto& ref = gauss_legendre<96>();
            const double half = 0.5 * (window.rmax - window.rmin) / npanel;
            double I = 0.0;
            for (int p = 0; p < npanel; p++)
              for (int g = 0; g < 96; g++)
                I += half * ref.w[g] * f(window.rmin + half * (2 * p + 1 + ref.x[g]));
            const double ends = (window.rmin > 0.0 ? f(window.rmin) : 0.0) + f(window.rmax);
            assert(ends * (window.rmax - window.rmin) <= 1e-12 * scale + tolerance &&
                   "type 2 radial integrand not negligible at window end");
            assert(std::abs(sum - I) <= 1e-11 * scale + tolerance &&
                   "fixed type 2 radial quadrature disagrees with 96-point rule");
          }
#endif
          radials(N, l1, l2) += coef * sum;
        }
      }
    }
  }
}
void RadialIntegral::type1(const int L, const ECP& U, const GaussianShell& shellA,
                           const GaussianShell& shellB, const ShellPairData& data,
                           const Parameters& parameters, ThreeIndex<double>& radials) const {
  // radials(N, lam, lam + mu) = sum_ab d_a d_b K_ab Y_lam,mu(P_ab)
  //                             * int r^(N+2) U_L(r) exp(-p (r - P)^2) K_lam(2 p P r) dr
  // for N <= L, lam <= N with lam = N mod 2, where p = a + b and P = |a A + b B| / p. For each
  // primitive pair and local ECP primitive this is the type 2 integrand with a single Gaussian
  // (p, P), so it is integrated with the same rule on the same window, for all N at once.
  radials = ThreeIndex<double>(L + 1, L + 1, 2 * L + 1);
  radials.fill(0.0);
  const int Lloc = U.getL();
  // Threshold for each primitive triple's bound
  const double thresh =
      triple_threshold(tolerance, data.Am, data.Bm, L,
                       shellA.nprimitive() * shellB.nprimitive() * count_channel(U, Lloc));
  std::vector<double> rq, wq, wrk, K, bes(L + 1), R((L + 1) * (L + 1));
  TwoIndex<double> harmonics;
  for (int na = 0; na < shellA.nprimitive(); na++) {
    for (int nb = 0; nb < shellB.nprimitive(); nb++) {
      const double p = parameters.p(na, nb), P = parameters.P(na, nb);
      const double cab = shellA.coef(na) * shellB.coef(nb) * parameters.K(na, nb);
      if (cab == 0.0) continue;
      std::fill(R.begin(), R.end(), 0.0);
      bool any = false;
      for (const auto& u : U.gaussians) {
        if (u.l != Lloc || u.d == 0.0) continue;
        const int kmin = u.n + 2, kmax = L + u.n + 2;
        if (std::abs(cab * u.d) * type2_bound(kmin, kmax, u.a, p, 0.0, P, 0.0) < thresh) continue;
        any = true;
        const auto window = integrand_window(kmin, kmax, L, 0, u.a, p, 0.0, P, 0.0);
        window_rule(window, u.a + p, rq, wq);
        const int nq = rq.size();
        wrk.resize((L + 1) * nq);
        K.resize((L + 1) * nq);
        const double x = 2.0 * p * P;
        for (int g = 0; g < nq; g++) {
          const double r = rq[g], rP = r - P;
          double e = wq[g] * std::exp(-u.a * r * r - p * rP * rP) * FAST_POW[kmin](r);
          for (int N = 0; N <= L; N++, e *= r) wrk[N * nq + g] = e;
          if (x > 0.0) {
            bessie.calculate(x * r, L, bes);
            for (int l = 0; l <= L; l++) K[l * nq + g] = bes[l];
          } else
            for (int l = 0; l <= L; l++) K[l * nq + g] = l == 0 ? 1.0 : 0.0;
        }
        for (int N = 0; N <= L; N++)
          for (int lam = N % 2; lam <= N; lam += 2) {
            double sum = 0.0;
            for (int g = 0; g < nq; g++) sum += wrk[N * nq + g] * K[lam * nq + g];
            R[N * (L + 1) + lam] += u.d * sum;
          }
      }
      if (!any) continue;

      // Direction of the product centre, with sin(theta) from the transverse components so that
      // it is exactly 0 on the z axis
      const double za = shellA.exp(na), zb = shellB.exp(nb);
      double Pv[3];
      for (int i = 0; i < 3; i++) Pv[i] = (za * data.A[i] + zb * data.B[i]) / p;
      const double Pt = std::sqrt(Pv[0] * Pv[0] + Pv[1] * Pv[1]);
      const double Pm = std::sqrt(Pt * Pt + Pv[2] * Pv[2]);
      const bool origin = Pm < 1e-12;
      realSphericalHarmonics(L, origin ? 0.0 : Pv[2] / Pm, origin ? 1.0 : Pt / Pm,
                             std::atan2(Pv[1], Pv[0]), harmonics);
      for (int N = 0; N <= L; N++)
        for (int lam = N % 2; lam <= N; lam += 2)
          for (int mu = -lam; mu <= lam; mu++)
            radials(N, lam, lam + mu) += cab * harmonics(lam, lam + mu) * R[N * (L + 1) + lam];
    }
  }
}
}  // namespace libecpint
