/*
 *      Copyright (c) 2020 Robert Shaw
 *		This file was generated as a part of Libecpint.
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

#include <algorithm>
#include <cassert>
#include <cmath>
#include <utility>

#include "mathutil.hpp"
#include "radial.hpp"

namespace libecpint {

namespace {
// Region [rmin, rmax] outside which r^k exp(-n r^2 - a (r - A)^2 - b (r - B)^2) K_l1(2aAr)
// K_l2(2bBr), with k in [kmin, kmax], l1 <= l1max and l2 <= l2max, is negligible. The Gaussian
// factor is exactly exp(-zt (r - pt)^2) up to a constant, with zt = n + a + b and pt = (aA + bB) /
// zt. The scaled Bessel functions K_l(z) = exp(-z) i_l(z) have logarithmic derivatives (in r)
// between -1/r and l/r, so the remaining factor can only move the maximum to between the
// stationary points of exp(-zt (r - pt)^2) r^m for m = kmin - 2 and m = kmax + l1max + l2max.
// The window extends 7 and 9 Gaussian widths beyond those (as transformRMinMax does about pt).
//
// A finite window is needed because the integrand can be far narrower than the spacing of a grid
// on [0, inf) (width 1/sqrt(zt) ~ 0.02 for a tight core s primitive), and the adaptive quadrature
// then reports convergence to zero.
std::pair<double, double> integrand_window(const int kmin, const int kmax, const int l1max,
                                           const int l2max, const double n, const double a,
                                           const double b, const double A, const double B) {
  const double zt = n + a + b;
  const double pt = (a * A + b * B) / zt;
  const double width = 1.0 / std::sqrt(zt);
  auto stationary = [&](const double m) {
    const double disc = 0.25 * pt * pt + 0.5 * m / zt;
    return disc > 0.0 ? 0.5 * pt + std::sqrt(disc) : 0.0;
  };
  const double rmin = stationary(kmin - 2) - 7.0 * width;
  const double rmax = stationary(kmax + l1max + l2max) + 9.0 * width;
  return {std::max(rmin, 0.0), rmax};
}

// N-point Gauss-Legendre rule on [-1, 1]. Numerical type 2 radials use a fixed rule on the window
// above rather than an adaptive one: there the integrand is a Gaussian of at most ~20 widths times
// factors varying on the scale of r, so the rule converges exponentially (against adaptive
// quadrature on real shell pairs: 24 points 1e-6, 32 points 3e-11, 48 points 2e-14 relative).
// Fixed nodes also make the result a smooth function of the geometry, which an adaptive stopping
// rule is not.
template <int N>
struct GaussLegendre {
  double x[N], w[N];
  GaussLegendre() {
    for (int i = 0; i < N; i++) {
      double z = std::cos(M_PI * (i + 0.75) / (N + 0.5)), z1, dp;
      do {
        double p1 = 1.0, p2 = 0.0;
        for (int j = 1; j <= N; j++) {
          const double p3 = p2;
          p2 = p1;
          p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3) / j;
        }
        dp = N * (z * p1 - p2) / (z * z - 1.0);
        z1 = z;
        z = z1 - p1 / dp;
      } while (std::abs(z - z1) > 1e-15);
      x[i] = z;
      w[i] = 2.0 / ((1.0 - z * z) * dp * dp);
    }
  }
};
template <int N>
const GaussLegendre<N>& gauss_legendre() {
  static const GaussLegendre<N> gl;
  return gl;
}
constexpr int TYPE2_QUAD_POINTS = 48;

// Rigorous upper bound on |G(k, l1, l2)| for all k in [kmin, kmax] and all l1, l2. Writing
// n r^2 + a (r - A)^2 + b (r - B)^2 = zt (r - pt)^2 + D, and as 0 <= K_l <= 1,
//   |G| <= e^-D int_0^inf r^k e^(-zt (r - pt)^2) dr
//       <= e^-D max_r [r^k e^(-zt (r - pt)^2 / 2)] sqrt(2 pi / zt),
// with r^k <= max(r^kmin, r^kmax).
double type2_bound(const int kmin, const int kmax, const double n, const double a, const double b,
                   const double A, const double B) {
  const double zt = n + a + b;
  const double pt = (a * A + b * B) / zt;
  const double D = a * A * A + b * B * B - zt * pt * pt;
  auto peak = [&](const int k) {
    const double rs = 0.5 * pt + std::sqrt(0.25 * pt * pt + k / zt);
    const double d = rs - pt;
    return std::pow(rs, k) * std::exp(-0.5 * zt * d * d);
  };
  return std::exp(-D) * std::max(peak(kmin), peak(kmax)) * std::sqrt(2.0 * M_PI / zt);
}
}  // namespace

void RadialIntegral::type2(const std::vector<Triple>& triples, const int /*nbase*/, const int lam,
                           const ECP& U, const GaussianShell& shellA, const GaussianShell& shellB,
                           const double A, const double B, ThreeIndex<double>& radials) const {
  // All radials are integrated numerically with a fixed Gauss-Legendre rule on a window
  // bracketing the integrand. The closed forms in terms of base integrals (Shaw2017) that were
  // used where available lost accuracy to cancellation, for diffuse primitives and also for a
  // tight primitive paired with a moderate one (up to 2% relative), and with the factors shared
  // by all triples tabulated once per primitive triple, quadrature costs little more.
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
  // Per primitive triple, wexpq holds weight * exponential factor at each node, rkq the powers
  // r^k (k - u.n - 2 in [Nmin, Nmax]), and Kqa, Kqb the scaled Bessel functions of each order.
  const auto& gl = gauss_legendre<TYPE2_QUAD_POINTS>();
  constexpr int nq = TYPE2_QUAD_POINTS;
  const int nk = Nmax - Nmin + 1;
  std::vector<double> wexpq(nq), rkq(nk * nq), Kqa((lmaxA + 1) * nq), Kqb((lmaxB + 1) * nq);
  std::vector<double> bes(std::max(lmaxA, lmaxB) + 1);

  // Loop over primitives in ECP, only considering correct ang. momentum
  for (const auto& u : U.gaussians) {
    if (u.l != lam) continue;
    const int kmin = Nmin + u.n + 2, kmax = Nmax + u.n + 2;
    // Loop over primitives in orbital basis shells
    for (int na = 0; na < npA; na++) {
      const double a = shellA.exp(na);
      const double da = shellA.coef(na);
      for (int nb = 0; nb < npB; nb++) {
        const double b = shellB.exp(nb);
        const double db = shellB.coef(nb);
        const double coef = da * db * u.d;
        if (std::abs(coef) < tolerance) continue;
        // Skip primitive triples whose every radial is provably negligible (e.g. a tight
        // primitive far from the ECP)
        if (std::abs(coef) * type2_bound(kmin, kmax, u.a, a, b, A, B) < tolerance) continue;

        const double xA = 2.0 * a * A, xB = 2.0 * b * B;  // Bessel arguments are xA r, xB r
        const std::pair<double, double> window =
            integrand_window(kmin, kmax, lmaxA, lmaxB, u.a, a, b, A, B);
        const double half = 0.5 * (window.second - window.first);
        const double mid = 0.5 * (window.second + window.first);
        for (int g = 0; g < nq; g++) {
          const double r = half * gl.x[g] + mid, rA = r - A, rB = r - B;
          wexpq[g] = half * gl.w[g] * std::exp(-u.a * r * r - a * rA * rA - b * rB * rB);
          double rk = FAST_POW[kmin](r);
          for (int n = 0; n < nk; n++, rk *= r) rkq[n * nq + g] = rk;
          // r > 0 and xA, xB > 0, so every order up to lmax is written
          bessie.calculate(xA * r, lmaxA, bes);
          for (int l = 0; l <= lmaxA; l++) Kqa[l * nq + g] = bes[l];
          bessie.calculate(xB * r, lmaxB, bes);
          for (int l = 0; l <= lmaxB; l++) Kqb[l * nq + g] = bes[l];
        }

        for (const Triple& t : triples) {
          const int N = std::get<0>(t), l1 = std::get<1>(t), l2 = std::get<2>(t);
          const double* rk = &rkq[(N - Nmin) * nq];
          const double* Ka = &Kqa[l1 * nq];
          const double* Kb = &Kqb[l2 * nq];
          double sum = 0.0;
          for (int g = 0; g < nq; g++) sum += rk[g] * wexpq[g] * Ka[g] * Kb[g];
#ifndef NDEBUG
          // Debug builds: the (positive) integrand must be negligible at the window ends, and the
          // fixed rule must agree with a higher-order one on the same window. Both are measured
          // against the radial with l1 = l2 = 0 and the same k, which bounds all l (K_l <= K_0):
          // relative to itself, a tiny high-l radial carries the absolute error of the tabulated
          // Bessel functions at small argument.
          {
            const int k = N + u.n + 2;
            double scale = 0.0;
            for (int g = 0; g < nq; g++) scale += rk[g] * wexpq[g] * Kqa[g] * Kqb[g];
            auto f = [&](const double r) {
              const double rA = r - A, rB = r - B;
              return FAST_POW[k](r) * std::exp(-u.a * r * r - a * rA * rA - b * rB * rB) *
                     bessie.calculate(xA * r, l1) * bessie.calculate(xB * r, l2);
            };
            const auto& ref = gauss_legendre<96>();
            double I = 0.0;
            for (int g = 0; g < 96; g++) I += half * ref.w[g] * f(half * ref.x[g] + mid);
            const double ends = (window.first > 0.0 ? f(window.first) : 0.0) + f(window.second);
            assert(ends * 2.0 * half <= 1e-12 * scale + tolerance &&
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
}  // namespace libecpint
