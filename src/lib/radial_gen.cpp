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

#include "mathutil.hpp"
#include "radial.hpp"
#ifdef USING_CERF
#include "cerf.h"
#else
#include "Faddeeva.hpp"
#endif
#include <algorithm>
#include <cmath>
#include <iostream>
#include <utility>

namespace libecpint {

namespace {
// Smallest min(x, y) = min(aA, bB), indexed by l1 + l2, for which the closed forms below are used.
// Below it they lose accuracy to cancellation (measured against independent quadrature: e.g. with
// l1 + l2 = 4, relative errors 5e-6 at min(x, y) = 0.16 and 0.5 at 0.04), and the radial integral
// is computed numerically instead. Chosen so the closed forms stay within ~1e-10 relative.
constexpr int CLOSED_FORM_NL = 9;
constexpr double CLOSED_FORM_MIN_XY[CLOSED_FORM_NL] = {0.0, 0.04, 0.25, 0.5, 0.6,
                                                       0.6, 0.6,  1.0,  1.0};

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
}  // namespace

void RadialIntegral::compute_base_integrals(const int N_min, const int N_max, const double p,
                                            const double o_root_p, const double P1, const double P2,
                                            const double P1_2, const double P2_2, const double X1,
                                            const double X2, const double oP1, const double oP2,
                                            double* values) const {
  // Recursively construct the base integrals in order F2, G3, F4, G5, etc... as described in
  // Shaw2017

  int imax = N_max / 2;
  int imin = (N_min + 1) / 2;
  int gmax = (N_max - 1) / 2;
  int gmin = N_min / 2;

  double P1_2k = 1.0;
  double P2_2k = 1.0;

  for (int k = 2; k < imin; k++) {
    P1_2k *= P1_2;
    P2_2k *= P2_2;
  }

  double ck, dk, ek, val;
  double C0 = o_root_p * ROOT_PI;
  for (int n = imin; n <= imax; n++) {
    ck = C0;
    dk = P1_2k * X1;
    ek = P2_2k * X2;
    val = ck * (dk - ek);

    for (int k = n - 1; k > 1; k--) {
      ck *= 2 * k * (2 * k - 1) * (n - k - 0.5) / ((2 * n - 2 * k) * (2 * n - 2 * k - 1) * p);
      dk *= oP1;
      ek *= oP2;
      val += ck * (dk - ek);
    }

    if (n > 1) {
      ck *= 2 * (n - 1.5) / ((2 * n - 2) * (2 * n - 3) * p);
      val += ck * (X1 - X2);
    }

    values[2 * n - N_min] = val;

    P1_2k *= P1_2;
    P2_2k *= P2_2;
  }

  P1_2k = P1;
  P2_2k = P2;
  for (int k = 1; k < gmin; k++) {
    P1_2k *= P1_2;
    P2_2k *= P2_2;
  }

  for (int n = gmin; n <= gmax; n++) {
    ck = C0;
    dk = P1_2k * X1;
    ek = P2_2k * X2;
    val = ck * (dk - ek);

    for (int k = n - 1; k > 0; k--) {
      ck *= 2 * k * (2 * k + 1) * (n - k - 0.5) / ((2 * n - 2 * k) * (2 * n - 1 - 2 * k) * p);
      dk *= oP1;
      ek *= oP2;
      val += ck * (dk - ek);
    }

    values[2 * n + 1 - N_min] = val;

    P1_2k *= P1_2;
    P2_2k *= P2_2;
  }
}

std::pair<double, bool> RadialIntegral::integrate_small(const int N, const int l1, const int l2,
                                                        const double n, const double a,
                                                        const double b, const double A,
                                                        const double B, const double rmin,
                                                        const double rmax) const {
  // Fallback for type2() when the small grid on [rmin, rmax] does not converge: the large grid
  // on the same window.
  const double aA = 2.0 * a * A;
  const double bB = 2.0 * b * B;

  thread_local std::vector<double> Ftab;
  auto tabulate = [&](const std::vector<double>& gridPoints, const int gridSize) {
    Ftab.resize(gridSize);
    for (int i = 0; i < gridSize; i++) {
      const double z = gridPoints[i];
      const double zA = z - A;
      const double zB = z - B;
      Ftab[i] = FAST_POW[N](z) * exp(-n * z * z - a * zA * zA - b * zB * zB) *
                bessie.calculate(aA * z, l1) * bessie.calculate(bB * z, l2);
    }
  };

  // Reused per-thread scratch to avoid reallocating the grid on every call
  thread_local GCQuadrature transformedGrid;
  transformedGrid = bigGrid;
  transformedGrid.transformInterval(rmin, rmax);
  tabulate(transformedGrid.getX(), transformedGrid.getN());
  return transformedGrid.integrate(Ftab.data(), tolerance, 0, transformedGrid.getN() - 1);
}

void RadialIntegral::type2(const std::vector<Triple>& triples, const int nbase, const int lam,
                           const ECP& U, const GaussianShell& shellA, const GaussianShell& shellB,
                           const double A, const double B, ThreeIndex<double>& radials) const {
  int npA = shellA.nprimitive();
  int npB = shellB.nprimitive();

  // Base integrals are fully recomputed for each primitive triple, so a single reused
  // buffer suffices - this avoids a heap allocation in the innermost loop.
  std::vector<double> base_integrals(nbase + 2);

  // Radial integrals without a usable closed form are integrated numerically. Per primitive pair,
  // the exponential factor and every Bessel order on the quadrature grid are shared by all such
  // triples, so they are tabulated once (lazily) rather than per triple.
  if (triples.empty()) return;
  int lmaxA = 0, lmaxB = 0, Nmin = std::get<0>(triples[0]), Nmax = Nmin;
  for (const Triple& t : triples) {
    lmaxA = std::max(lmaxA, std::get<1>(t));
    lmaxB = std::max(lmaxB, std::get<2>(t));
    Nmin = std::min(Nmin, std::get<0>(t));
    Nmax = std::max(Nmax, std::get<0>(t));
  }
  // The small grid, mapped per primitive triple onto the window where the integrand lives
  GCQuadrature qgrid = primGrid;
  const int nq = qgrid.getN();
  const std::vector<double>& rq = qgrid.getX();
  std::vector<double> expq(nq), Kqa((lmaxA + 1) * nq), Kqb((lmaxB + 1) * nq), Fq(nq);
  std::vector<double> bes(std::max(lmaxA, lmaxB) + 1);

  // Loop over primitives in ECP, only considering correct ang. momentum
  for (const auto& u : U.gaussians) {
    if (u.l == lam) {
      // Loop over primitives in orbital basis shellß
      for (int na = 0; na < npA; na++) {
        double a = shellA.exp(na);
        double da = shellA.coef(na);

        for (int nb = 0; nb < npB; nb++) {
          double b = shellB.exp(nb);
          double db = shellB.coef(nb);

          // Construct values that will be reused across all radial integrals
          double p = u.a + a + b;
          double x = a * A;
          double y = b * B;

          bool tables_built = false;
          std::pair<double, double> window;
          auto quadrature = [&](const int kk, const int ii, const int jj) {
            if (!tables_built) {
              window =
                  integrand_window(Nmin + u.n + 2, Nmax + u.n + 2, lmaxA, lmaxB, u.a, a, b, A, B);
              qgrid = primGrid;
              qgrid.transformInterval(window.first, window.second);
              for (int g = 0; g < nq; g++) {
                const double r = rq[g], rA = r - A, rB = r - B;
                expq[g] = std::exp(-u.a * r * r - a * rA * rA - b * rB * rB);
                std::fill(bes.begin(), bes.end(), 0.0);
                bessie.calculate(2.0 * x * r, lmaxA, bes);
                for (int l = 0; l <= lmaxA; l++) Kqa[l * nq + g] = bes[l];
                std::fill(bes.begin(), bes.end(), 0.0);
                bessie.calculate(2.0 * y * r, lmaxB, bes);
                for (int l = 0; l <= lmaxB; l++) Kqb[l * nq + g] = bes[l];
              }
              tables_built = true;
            }
            for (int g = 0; g < nq; g++)
              Fq[g] = FAST_POW[kk](rq[g]) * expq[g] * Kqa[ii * nq + g] * Kqb[jj * nq + g];
            std::pair<double, bool> q = qgrid.integrate(Fq.data(), tolerance, 0, nq - 1);
            if (!q.second)  // large grid on the same window
              q = integrate_small(kk, ii, jj, u.a, a, b, A, B, window.first, window.second);
            if (!q.second) std::cout << "Quadrature failed" << std::endl;
            return q.first;
          };

          double P1 = (x + y) / p;
          double P2 = (y - x) / p;
          double P1_2 = P1 * P1;
          double P2_2 = P2 * P2;
          double oP1 = 1.0 / P1_2;
          double oP2 = std::abs(P2) < 1e-7 ? 0.0 : 1.0 / P2_2;
          double root_p = sqrt(p);
          double o_root_p = 1.0 / root_p;
          double aAbB = a * A * A + b * B * B;
          double Kab = 1.0 / (16.0 * x * y);
          double X1 = exp(p * P1_2 - aAbB) * Kab;
          double X2 = exp(p * P2_2 - aAbB) * Kab;

          double x2 = x * x;
          double y2 = y * y;
          double p2 = p * p;

          double result = 0.0;

          // G1A, G1B may not be required, but it seems to be quicker to calculate than to check if
          // needed
#ifdef USING_CERF
          double daw1 = X1 * dawson(root_p * P1);
          double daw2 = X2 * dawson(root_p * P2);
#else
          double daw1 = X1 * Faddeeva::Dawson(root_p * P1);
          double daw2 = X2 * Faddeeva::Dawson(root_p * P2);
#endif
          double G1B = 2.0 * ROOT_PI * (daw1 - daw2);
          double G1A = 2.0 * ROOT_PI * (daw1 + daw2);
          double H2 = ROOT_PI * (X1 + X2) * o_root_p;

          // Compute base integrals
          double* values = base_integrals.data();
          compute_base_integrals(2, 3 + nbase, p, o_root_p, P1, P2, P1_2, P2_2, X1, X2, oP1, oP2,
                                 values);

          // Loop over all radial integrals required, divert to generated code
          for (const Triple& triple : triples) {
            int i = std::get<1>(triple);
            int j = std::get<2>(triple);
            int k = std::get<0>(triple) + u.n + 2;

            int ijk = i * 10000 + j * 100 + k;
            double result = 0.0;
            if (i + j < CLOSED_FORM_NL && std::min(x, y) < CLOSED_FORM_MIN_XY[i + j]) {
              // The closed forms are polynomials in 1/x, 1/y and cancel catastrophically for a
              // diffuse primitive (small x = aA or y = bB), so integrate numerically instead.
              result = quadrature(k, i, j);
            } else if (a * b > MIN_EXP) {  // && b > MIN_EXP) {
              switch (ijk) {
                case 2: {
                  result = (1) * values[0];
                  break;
                }

                case 4: {
                  result += (1) * values[2];
                  break;
                }

                case 6: {
                  result += (1) * values[4];
                  break;
                }

                case 8: {
                  result += (1) * values[6];
                  break;
                }

                case 10: {
                  result += (1) * values[8];
                  break;
                }

                case 12: {
                  result += (1) * values[10];
                  break;
                }

                case 101: {
                  result = (p / y) * values[0];
                  result += (-x / y) * G1A;
                  break;
                }

                case 103: {
                  result = (-1 / (2 * y)) * values[0];
                  result += (1) * values[1];
                  break;
                }

                case 105: {
                  result += (-1 / (2 * y)) * values[2];
                  result += (1) * values[3];
                  break;
                }

                case 107: {
                  result += (-1 / (2 * y)) * values[4];
                  result += (1) * values[5];
                  break;
                }

                case 109: {
                  result += (-1 / (2 * y)) * values[6];
                  result += (1) * values[7];
                  break;
                }

                case 111: {
                  result += (-1 / (2 * y)) * values[8];
                  result += (1) * values[9];
                  break;
                }

                case 10102: {
                  result = (-(p / 2 + y2) / (x * y)) * values[0];
                  result += (p / x) * values[1];
                  break;
                }

                case 10104: {
                  result = (1 / (2 * x * y)) * values[0];
                  result += (-(p / 2 + y2) / (x * y)) * values[2];
                  result += (-1 / x) * values[1];
                  result += (p / x) * values[3];
                  break;
                }

                case 10106: {
                  result += (1 / (x * y)) * values[2];
                  result += (-(p / 2 + y2) / (x * y)) * values[4];
                  result += (-2 / x) * values[3];
                  result += (p / x) * values[5];
                  break;
                }

                case 10108: {
                  result += (3 / (2 * x * y)) * values[4];
                  result += (-(p / 2 + y2) / (x * y)) * values[6];
                  result += (-3 / x) * values[5];
                  result += (p / x) * values[7];
                  break;
                }

                case 10110: {
                  result += (4 / (2 * x * y)) * values[6];
                  result += (-(p / 2 + y2) / (x * y)) * values[8];
                  result += (-4 / x) * values[7];
                  result += (p / x) * values[9];
                  break;
                }

                case 202: {
                  result = (-3 * p / (2 * y2) + 1) * values[0];
                  result += (3 * x / (2 * y2)) * G1A;
                  break;
                }

                case 204: {
                  result = (3 / (4 * y2)) * values[0];
                  result += (1) * values[2];
                  result += (-3 / (2 * y)) * values[1];
                  break;
                }

                case 206: {
                  result += (3 / (4 * y2)) * values[2];
                  result += (1) * values[4];
                  result += (-3 / (2 * y)) * values[3];
                  break;
                }

                case 208: {
                  result += (3 / (4 * y2)) * values[4];
                  result += (1) * values[6];
                  result += (-3 / (2 * y)) * values[5];
                  break;
                }

                case 210: {
                  result += (3 / (4 * y2)) * values[6];
                  result += (1) * values[8];
                  result += (-3 / (2 * y)) * values[7];
                  break;
                }

                case 10201: {
                  result = (-p * (p + 2 * x2) / (2 * x * y2)) * values[0];
                  result += (x2 / y2) * G1A;
                  result += (p / y) * H2;
                  break;
                }

                case 10203: {
                  result = ((3 * p + 2 * y2) / (4 * x * y2)) * values[0];
                  result += (p / x) * values[2];
                  result += (-(3 * p / 2 + y2) / (x * y)) * values[1];
                  break;
                }

                case 10205: {
                  result = (-3 / (4 * x * y2)) * values[0];
                  result += ((3 * p - 2 * y2) / (4 * x * y2)) * values[2];
                  result += (p / x) * values[4];
                  result += (3 / (2 * x * y)) * values[1];
                  result += (-(3 * p / 2 + y2) / (x * y)) * values[3];
                  break;
                }

                case 10207: {
                  result += (-3 / (2 * x * y2)) * values[2];
                  result += (3 * (p - 2 * y2) / (4 * x * y2)) * values[4];
                  result += (p / x) * values[6];
                  result += (3 / (x * y)) * values[3];
                  result += (-(3 * p / 2 + y2) / (x * y)) * values[5];
                  break;
                }

                case 10209: {
                  result += (-9 / (4 * x * y2)) * values[4];
                  result += ((3 * p - 10 * y2) / (4 * x * y2)) * values[6];
                  result += (p / x) * values[8];
                  result += (9 / (2 * x * y)) * values[5];
                  result += (-(3 * p / 2 + y2) / (x * y)) * values[7];
                  break;
                }

                case 20202: {
                  result = ((3 * p2 + 4 * y2 * (p + y2)) / (4 * x2 * y2)) * values[0];
                  result += (p2 / x2) * values[2];
                  result += (-p * (3 * p + 4 * y2) / (2 * x2 * y)) * values[1];
                  break;
                }

                case 20204: {
                  result = (-(3 * p / 2 + y2) / (x2 * y2)) * values[0];
                  result += ((3 * p2 + 4 * y2 * (-p + y2)) / (4 * x2 * y2)) * values[2];
                  result += (p2 / x2) * values[4];
                  result += ((3 * p + 2 * y2) / (x2 * y)) * values[1];
                  result += (-p * (3 * p + 4 * y2) / (2 * x2 * y)) * values[3];
                  break;
                }

                case 20206: {
                  result = (3 / (2 * x2 * y2)) * values[0];
                  result += (-3 * p / (x2 * y2)) * values[2];
                  result += ((3 * p2 + 4 * y2 * (-3 * p + y2)) / (4 * x2 * y2)) * values[4];
                  result += (p2 / x2) * values[6];
                  result += (-3 / (x2 * y)) * values[1];
                  result += (2 * (3 * p + 2 * y2) / (x2 * y)) * values[3];
                  result += (-p * (3 * p + 4 * y2) / (2 * x2 * y)) * values[5];
                  break;
                }

                case 20208: {
                  result += (9 / (2 * x2 * y2)) * values[2];
                  result += (3 * (-3 * p + 2 * y2) / (2 * x2 * y2)) * values[4];
                  result += ((3 * p2 + 4 * y2 * (-5 * p + y2)) / (4 * x2 * y2)) * values[6];
                  result += (p2 / x2) * values[8];
                  result += (-9 / (x2 * y)) * values[3];
                  result += (3 * (3 * p + 2 * y2) / (x2 * y)) * values[5];
                  result += (-p * (3 * p + 4 * y2) / (2 * x2 * y)) * values[7];
                  break;
                }

                case 301: {
                  result = (p * (-5 * p + 5 * x2 + 2 * y2) / (2 * (y2 * y))) * values[0];
                  result += (x * (15 * p - 10 * x2 + 6 * y2) / (4 * (y2 * y))) * G1A;
                  result += (-5 * p * x / (2 * y2)) * H2;
                  break;
                }

                case 303: {
                  result = (15 * p / (4 * (y2 * y)) - 3 / y) * values[0];
                  result += (1) * values[1];
                  result += (-15 * x / (4 * (y2 * y))) * G1A;
                  break;
                }

                case 305: {
                  result = (-15 / (8 * (y2 * y))) * values[0];
                  result += (-3 / y) * values[2];
                  result += (15 / (4 * y2)) * values[1];
                  result += (1) * values[3];
                  break;
                }

                case 307: {
                  result += (-15 / (8 * (y2 * y))) * values[2];
                  result += (-3 / y) * values[4];
                  result += (15 / (4 * y2)) * values[3];
                  result += (1) * values[5];
                  break;
                }

                case 309: {
                  result += (-15 / (8 * (y2 * y))) * values[4];
                  result += (-3 / y) * values[6];
                  result += (15 / (4 * y2)) * values[5];
                  result += (1) * values[7];
                  break;
                }

                case 10302: {
                  result =
                      ((5 * p2 + 10 * p * x2 - 2 * p * y2 - 4 * (y2 * y2)) / (4 * x * (y2 * y))) *
                      values[0];
                  result += (p / x) * values[1];
                  result += (-5 * x2 / (2 * (y2 * y))) * G1A;
                  result += (-5 * p / (2 * y2)) * H2;
                  break;
                }

                case 10304: {
                  result = (-(15 * p + 6 * y2) / (8 * x * (y2 * y))) * values[0];
                  result += (-(3 * p + y2) / (x * y)) * values[2];
                  result += (3 * (5 * p + 2 * y2) / (4 * x * y2)) * values[1];
                  result += (p / x) * values[3];
                  break;
                }

                case 10306: {
                  result = (15 / (8 * x * (y2 * y))) * values[0];
                  result += (3 * (-5 * p + 6 * y2) / (8 * x * (y2 * y))) * values[2];
                  result += (-(3 * p + y2) / (x * y)) * values[4];
                  result += (-15 / (4 * x * y2)) * values[1];
                  result += ((15 * p + 2 * y2) / (4 * x * y2)) * values[3];
                  result += (p / x) * values[5];
                  break;
                }

                case 10308: {
                  result += (15 / (4 * x * (y2 * y))) * values[2];
                  result += (3 * (-5 * p + 14 * y2) / (8 * x * (y2 * y))) * values[4];
                  result += (-(3 * p + y2) / (x * y)) * values[6];
                  result += (-15 / (2 * x * y2)) * values[3];
                  result += ((15 * p - 2 * y2) / (4 * x * y2)) * values[5];
                  result += (p / x) * values[7];
                  break;
                }

                case 20301: {
                  result =
                      (p * (3 * p2 + 2 * p * x2 + 4 * (x2 * x2) + 4 * x2 * y2 - 4 * (y2 * y2)) /
                       (4 * x2 * (y2 * y))) *
                      values[0];
                  result += (p2 / x2) * values[1];
                  result += (-(x2 * x) / (y2 * y)) * G1A;
                  result += (-p * (3 * p + 2 * x2 + 2 * y2) / (2 * x * y2)) * H2;
                  break;
                }

                case 20303: {
                  result =
                      (-(15 * p2 + 12 * p * y2 + 4 * (y2 * y2)) / (8 * x2 * (y2 * y))) * values[0];
                  result += (-p * (3 * p + 2 * y2) / (x2 * y)) * values[2];
                  result += ((15 * p2 + 4 * y2 * (3 * p + y2)) / (4 * x2 * y2)) * values[1];
                  result += (p2 / x2) * values[3];
                  break;
                }

                case 20305: {
                  result = (3 * (5 * p + 2 * y2) / (4 * x2 * (y2 * y))) * values[0];
                  result += (3 * (-5 * p2 + 12 * p * y2 + 4 * (y2 * y2)) / (8 * x2 * (y2 * y))) *
                            values[2];
                  result += (-p * (3 * p + 2 * y2) / (x2 * y)) * values[4];
                  result += (-(15 * p + 6 * y2) / (2 * x2 * y2)) * values[1];
                  result += ((15 * p2 + 4 * y2 * (p + y2)) / (4 * x2 * y2)) * values[3];
                  result += (p2 / x2) * values[5];
                  break;
                }

                case 20307: {
                  result = (-15 / (4 * x2 * (y2 * y))) * values[0];
                  result += (3 * (5 * p - 2 * y2) / (2 * x2 * (y2 * y))) * values[2];
                  result +=
                      ((-15 * p2 + 84 * p * y2 + 28 * (y2 * y2)) / (8 * x2 * (y2 * y))) * values[4];
                  result += (-p * (3 * p + 2 * y2) / (x2 * y)) * values[6];
                  result += (15 / (2 * x2 * y2)) * values[1];
                  result += (-(15 * p + 4 * y2) / (x2 * y2)) * values[3];
                  result += ((15 * p2 + 4 * y2 * (-p + y2)) / (4 * x2 * y2)) * values[5];
                  result += (p2 / x2) * values[7];
                  break;
                }

                case 30302: {
                  result =
                      (-(15 * (p2 * p) + 18 * p2 * y2 + 12 * p * (y2 * y2) + 8 * (y2 * y2 * y2)) /
                       (8 * (x2 * x) * (y2 * y))) *
                      values[0];
                  result += (-3 * p2 * (p + y2) / ((x2 * x) * y)) * values[2];
                  result += (3 * p * (5 * p2 + 6 * p * y2 + 4 * (y2 * y2)) / (4 * (x2 * x) * y2)) *
                            values[1];
                  result += ((p2 * p) / (x2 * x)) * values[3];
                  break;
                }

                case 30304: {
                  result =
                      (3 * (15 * p2 + 12 * p * y2 + 4 * (y2 * y2)) / (8 * (x2 * x) * (y2 * y))) *
                      values[0];
                  result += ((-15 * (p2 * p) + 54 * p2 * y2 + 4 * (y2 * y2) * (9 * p - 2 * y2)) /
                             (8 * (x2 * x) * (y2 * y))) *
                            values[2];
                  result += (-3 * p2 * (p + y2) / ((x2 * x) * y)) * values[4];
                  result += (-(45 * p2 + 12 * y2 * (3 * p + y2)) / (4 * (x2 * x) * y2)) * values[1];
                  result +=
                      (3 * p * (5 * p2 + 2 * y2 * (p + 2 * y2)) / (4 * (x2 * x) * y2)) * values[3];
                  result += ((p2 * p) / (x2 * x)) * values[5];
                  break;
                }

                case 30306: {
                  result = (-(45 * p + 18 * y2) / (4 * (x2 * x) * (y2 * y))) * values[0];
                  result +=
                      (3 * (15 * p2 - 12 * p * y2 - 4 * (y2 * y2)) / (4 * (x2 * x) * (y2 * y))) *
                      values[2];
                  result += ((-15 * (p2 * p) + 126 * p2 * y2 + 4 * (y2 * y2) * (21 * p - 2 * y2)) /
                             (8 * (x2 * x) * (y2 * y))) *
                            values[4];
                  result += (-3 * p2 * (p + y2) / ((x2 * x) * y)) * values[6];
                  result += (9 * (5 * p + 2 * y2) / (2 * (x2 * x) * y2)) * values[1];
                  result += (-(45 * p2 + 12 * y2 * (2 * p + y2)) / (2 * (x2 * x) * y2)) * values[3];
                  result +=
                      (3 * p * (5 * p2 + 2 * y2 * (-p + 2 * y2)) / (4 * (x2 * x) * y2)) * values[5];
                  result += ((p2 * p) / (x2 * x)) * values[7];
                  break;
                }

                case 402: {
                  result = ((-5 * p * (-7 * p + 7 * x2 + 4 * y2) / 4 + (y2 * y2)) / (y2 * y2)) *
                           values[0];
                  result += (5 * x * (-21 * p + 14 * x2 - 6 * y2) / (8 * (y2 * y2))) * G1A;
                  result += (35 * p * x / (4 * (y2 * y))) * H2;
                  break;
                }

                case 404: {
                  result = (15 * (-7 * p + 6 * y2) / (8 * (y2 * y2))) * values[0];
                  result += (1) * values[2];
                  result += (-5 / y) * values[1];
                  result += (105 * x / (8 * (y2 * y2))) * G1A;
                  break;
                }

                case 406: {
                  result = (105 / (16 * (y2 * y2))) * values[0];
                  result += (45 / (4 * y2)) * values[2];
                  result += (1) * values[4];
                  result += (-105 / (8 * (y2 * y))) * values[1];
                  result += (-5 / y) * values[3];
                  break;
                }

                case 408: {
                  result += (105 / (16 * (y2 * y2))) * values[2];
                  result += (45 / (4 * y2)) * values[4];
                  result += (1) * values[6];
                  result += (-105 / (8 * (y2 * y))) * values[3];
                  result += (-5 / y) * values[5];
                  break;
                }

                case 10401: {
                  result =
                      (-p * (-7 * p2 - 28 * p * x2 + 2 * p * y2 + 14 * (x2 * x2) + 4 * x2 * y2) /
                       (4 * x * (y2 * y2))) *
                      values[0];
                  result += (x2 * (-35 * p + 14 * x2 - 10 * y2) / (4 * (y2 * y2))) * G1A;
                  result += (p * (-7 * p + 7 * x2 + 2 * y2) / (2 * (y2 * y))) * H2;
                  break;
                }

                case 10403: {
                  result = (-(35 * p2 + 70 * p * x2 - 20 * p * y2 - 32 * (y2 * y2)) /
                            (8 * x * (y2 * y2))) *
                           values[0];
                  result += (p / x) * values[2];
                  result += (-(5 * p + y2) / (x * y)) * values[1];
                  result += (35 * x2 / (4 * (y2 * y2))) * G1A;
                  result += (35 * p / (4 * (y2 * y))) * H2;
                  break;
                }

                case 10405: {
                  result = (15 * (7 * p + 2 * y2) / (16 * x * (y2 * y2))) * values[0];
                  result += (45 * p / (4 * x * y2) + 3 / x) * values[2];
                  result += (p / x) * values[4];
                  result += (-(105 * p + 30 * y2) / (8 * x * (y2 * y))) * values[1];
                  result += (-(5 * p + y2) / (x * y)) * values[3];
                  break;
                }

                case 10407: {
                  result = (-105 / (16 * x * (y2 * y2))) * values[0];
                  result += (15 * (7 * p - 10 * y2) / (16 * x * (y2 * y2))) * values[2];
                  result += (45 * p / (4 * x * y2) + 2 / x) * values[4];
                  result += (p / x) * values[6];
                  result += (105 / (8 * x * (y2 * y))) * values[1];
                  result += (5 * (-21 * p + 2 * y2) / (8 * x * (y2 * y))) * values[3];
                  result += (-(5 * p + y2) / (x * y)) * values[5];
                  break;
                }

                case 20402: {
                  result = (-(21 * (p2 * p) + 14 * p2 * x2 - 6 * p2 * y2 + 28 * p * (x2 * x2) +
                              28 * p * x2 * y2 - 36 * p * (y2 * y2) - 8 * (y2 * y2 * y2)) /
                            (8 * x2 * (y2 * y2))) *
                           values[0];
                  result += (p2 / x2) * values[2];
                  result += (-p * (5 * p + 2 * y2) / (x2 * y)) * values[1];
                  result += (7 * (x2 * x) / (2 * (y2 * y2))) * G1A;
                  result += (7 * p * (3 * p + 2 * x2 + 2 * y2) / (4 * x * (y2 * y))) * H2;
                  break;
                }

                case 20404: {
                  result = (3 * (35 * p2 + 20 * p * y2 + 4 * (y2 * y2)) / (16 * x2 * (y2 * y2))) *
                           values[0];
                  result += ((45 * p2 + 4 * y2 * (6 * p + y2)) / (4 * x2 * y2)) * values[2];
                  result += (p2 / x2) * values[4];
                  result += (-(105 * p2 + 60 * p * y2 + 12 * (y2 * y2)) / (8 * x2 * (y2 * y))) *
                            values[1];
                  result += (-p * (5 * p + 2 * y2) / (x2 * y)) * values[3];
                  break;
                }

                case 20406: {
                  result = (-(105 * p + 30 * y2) / (8 * x2 * (y2 * y2))) * values[0];
                  result +=
                      (3 * (35 * p2 - 100 * p * y2 - 28 * (y2 * y2)) / (16 * x2 * (y2 * y2))) *
                      values[2];
                  result += ((45 * p2 + 4 * y2 * (4 * p + y2)) / (4 * x2 * y2)) * values[4];
                  result += (p2 / x2) * values[6];
                  result += (15 * (7 * p + 2 * y2) / (4 * x2 * (y2 * y))) * values[1];
                  result +=
                      ((-105 * p2 + 20 * p * y2 + 4 * (y2 * y2)) / (8 * x2 * (y2 * y))) * values[3];
                  result += (-p * (5 * p + 2 * y2) / (x2 * y)) * values[5];
                  break;
                }

                case 30401: {
                  result = (-p *
                            (15 * (p2 * p) + 6 * p2 * x2 + 4 * p * (x2 * x2) + 24 * p * x2 * y2 -
                             36 * p * (y2 * y2) + 8 * (x2 * x2 * x2) + 8 * (x2 * x2) * y2 +
                             8 * x2 * (y2 * y2) - 16 * (y2 * y2 * y2)) /
                            (8 * (x2 * x) * (y2 * y2))) *
                           values[0];
                  result += ((p2 * p) / (x2 * x)) * values[2];
                  result += (-p2 * (5 * p + 3 * y2) / ((x2 * x) * y)) * values[1];
                  result += ((x2 * x2) / (y2 * y2)) * G1A;
                  result += (p *
                             (15 * p2 + 6 * p * x2 + 20 * p * y2 + 4 * (x2 * x2) + 4 * x2 * y2 +
                              4 * (y2 * y2)) /
                             (4 * x2 * (y2 * y))) *
                            H2;
                  break;
                }

                case 30403: {
                  result =
                      ((105 * (p2 * p) + 90 * p2 * y2 + 36 * p * (y2 * y2) + 8 * (y2 * y2 * y2)) /
                       (16 * (x2 * x) * (y2 * y2))) *
                      values[0];
                  result +=
                      (3 * p * (15 * p2 + 4 * y2 * (3 * p + y2)) / (4 * (x2 * x) * y2)) * values[2];
                  result += ((p2 * p) / (x2 * x)) * values[4];
                  result +=
                      (-(105 * (p2 * p) + 90 * p2 * y2 + 36 * p * (y2 * y2) + 8 * (y2 * y2 * y2)) /
                       (8 * (x2 * x) * (y2 * y))) *
                      values[1];
                  result += (-p2 * (5 * p + 3 * y2) / ((x2 * x) * y)) * values[3];
                  break;
                }

                case 30405: {
                  result =
                      (-(315 * p2 + 180 * p * y2 + 36 * (y2 * y2)) / (16 * (x2 * x) * (y2 * y2))) *
                      values[0];
                  result += (-(-105 * (p2 * p) + 450 * p2 * y2 + 252 * p * (y2 * y2) +
                               40 * (y2 * y2 * y2)) /
                             (16 * (x2 * x) * (y2 * y2))) *
                            values[2];
                  result +=
                      (3 * p * (15 * p2 + 4 * y2 * (2 * p + y2)) / (4 * (x2 * x) * y2)) * values[4];
                  result += ((p2 * p) / (x2 * x)) * values[6];
                  result +=
                      (9 * (35 * p2 + 20 * p * y2 + 4 * (y2 * y2)) / (8 * (x2 * x) * (y2 * y))) *
                      values[1];
                  result += ((-105 * (p2 * p) + 30 * p2 * y2 + 4 * (y2 * y2) * (3 * p - 2 * y2)) /
                             (8 * (x2 * x) * (y2 * y))) *
                            values[3];
                  result += (-p2 * (5 * p + 3 * y2) / ((x2 * x) * y)) * values[5];
                  break;
                }

                case 40402: {
                  result = ((105 * (p2 * p2) + 120 * (p2 * p) * y2 + 72 * p2 * (y2 * y2) +
                             32 * p * (y2 * y2 * y2) + 16 * (y2 * y2 * y2 * y2)) /
                            (16 * (x2 * x2) * (y2 * y2))) *
                           values[0];
                  result += (3 * p2 * (15 * p2 + 8 * y2 * (2 * p + y2)) / (4 * (x2 * x2) * y2)) *
                            values[2];
                  result += ((p2 * p2) / (x2 * x2)) * values[4];
                  result +=
                      (-p *
                       (105 * (p2 * p) + 120 * p2 * y2 + 72 * p * (y2 * y2) + 32 * (y2 * y2 * y2)) /
                       (8 * (x2 * x2) * (y2 * y))) *
                      values[1];
                  result += (-(p2 * p) * (5 * p + 4 * y2) / ((x2 * x2) * y)) * values[3];
                  break;
                }

                case 40404: {
                  result =
                      (-(105 * (p2 * p) + 90 * p2 * y2 + 36 * p * (y2 * y2) + 8 * (y2 * y2 * y2)) /
                       (4 * (x2 * x2) * (y2 * y2))) *
                      values[0];
                  result += ((105 * (p2 * p2) - 600 * (p2 * p) * y2 - 504 * p2 * (y2 * y2) -
                              160 * p * (y2 * y2 * y2) + 16 * (y2 * y2 * y2 * y2)) /
                             (16 * (x2 * x2) * (y2 * y2))) *
                            values[2];
                  result += (p2 * (45 * p2 + 32 * p * y2 + 24 * (y2 * y2)) / (4 * (x2 * x2) * y2)) *
                            values[4];
                  result += ((p2 * p2) / (x2 * x2)) * values[6];
                  result +=
                      ((105 * (p2 * p) + 90 * p2 * y2 + 36 * p * (y2 * y2) + 8 * (y2 * y2 * y2)) /
                       (2 * (x2 * x2) * (y2 * y))) *
                      values[1];
                  result += (p * (-105 * (p2 * p) + 40 * p2 * y2 + (y2 * y2) * (24 * p - 32 * y2)) /
                             (8 * (x2 * x2) * (y2 * y))) *
                            values[3];
                  result += (-(p2 * p) * (5 * p + 4 * y2) / ((x2 * x2) * y)) * values[5];
                  break;
                }

                default: {
                  if (estimate_type2(k, i, j, u.a, a, b, A, B) > tolerance) {
                    result = quadrature(k, i, j);
                  }
                }
              }
            } else {
              if (estimate_type2(k, i, j, u.a, a, b, A, B) > tolerance) {
                result = quadrature(k, i, j);
              }
            }

            radials(k - 2 - u.n, i, j) += da * db * u.d * result;
          }
        }
      }
    }
  }
}
}  // namespace libecpint
