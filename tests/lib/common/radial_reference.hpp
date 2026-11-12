/*
 * Independent reference values for type-2 (semi-local) ECP radial integrals, for regression tests.
 *
 * Deliberately shares no code with the library: modified spherical Bessel functions are evaluated
 * from their power series / exact finite sums, and the radial integral by composite Gauss-Legendre
 * quadrature over a range bounded using the r^k exp(-u r^2) envelope.
 *
 * Definitions (as in RadialIntegral::type2 for shells A, B both off the ECP centre):
 *   K_l(z) = exp(-z) i_l(z)
 *   G(k, l1, l2) = int_0^inf r^k exp(-u r^2 - a (r - A)^2 - b (r - B)^2) K_l1(2 a A r) K_l2(2 b B r) dr
 * where u is the ECP exponent, a and b the shell exponents and A, B the distances of the shell
 * centres from the ECP centre.
 */
#ifndef LIBECPINT_TESTS_RADIAL_REFERENCE_HPP
#define LIBECPINT_TESTS_RADIAL_REFERENCE_HPP

#include <cmath>
#include <vector>

namespace radref {

/// exp(-z) i_l(z), l >= 0, z >= 0
inline double K(const int l, const double z) {
  if (z < 1e-12) return l == 0 ? 1.0 : 0.0;
  if (z < 8.0) {
    // i_l(z) = z^l sum_k (z^2/2)^k / (k! (2l+2k+1)!!); all terms positive
    double dfac = 1.0;
    for (int m = 1; m <= 2 * l + 1; m += 2) dfac *= m;
    double term = std::pow(z, l) / dfac, sum = term;
    for (int k = 1; k < 400; k++) {
      term *= (0.5 * z * z) / (k * (2.0 * l + 2.0 * k + 1.0));
      sum += term;
      if (term < 1e-18 * sum) break;
    }
    return std::exp(-z) * sum;
  }
  // exact: i_l(z) = [e^z R_l(-z) - (-1)^l e^-z R_l(z)] / (2z), R_l(x) = sum_k (l+k)!/(k!(l-k)!) (2x)^-k
  double rm = 0.0, rp = 0.0, c = 1.0;
  for (int k = 0; k <= l; k++) {
    if (k > 0) c *= double((l - k + 1) * (l + k)) / k;
    rm += c * std::pow(-2.0 * z, -k);
    rp += c * std::pow(2.0 * z, -k);
  }
  const double sgn = (l % 2) ? -1.0 : 1.0;
  return (rm - sgn * std::exp(-2.0 * z) * rp) / (2.0 * z);
}

/// 16-point Gauss-Legendre nodes/weights on [-1, 1]
inline void gauss_legendre16(std::vector<double>& x, std::vector<double>& w) {
  const int n = 16;
  x.assign(n, 0.0);
  w.assign(n, 0.0);
  for (int i = 0; i < n; i++) {
    double t = std::cos(M_PI * (i + 0.75) / (n + 0.5)), dp = 0.0;
    for (int it = 0; it < 100; it++) {
      double p0 = 1.0, p1 = t;
      for (int k = 2; k <= n; k++) {
        double p2 = ((2 * k - 1) * t * p1 - (k - 1) * p0) / k;
        p0 = p1;
        p1 = p2;
      }
      dp = n * (t * p1 - p0) / (t * t - 1.0);
      const double dt = p1 / dp;
      t -= dt;
      if (std::abs(dt) < 1e-16) break;
    }
    x[i] = t;
    w[i] = 2.0 / ((1.0 - t * t) * dp * dp);
  }
}

/// G(k, l1, l2) as defined above, by composite Gauss-Legendre quadrature
inline double G(const int k, const int l1, const int l2, const double u, const double a, const double b,
                const double A, const double B) {
  // integrand <= r^k exp(-u r^2); choose rmax where that bound is negligible (< ~1e-32)
  double rmax = std::sqrt(80.0 / u) + 1.0;
  for (int it = 0; it < 20; it++) rmax = std::sqrt((80.0 + k * std::log(std::max(rmax, 1.0))) / u) + 1.0;
  std::vector<double> x, w;
  gauss_legendre16(x, w);
  const int panels = 2000;
  const double h = rmax / panels;
  double sum = 0.0;
  for (int p = 0; p < panels; p++) {
    const double mid = (p + 0.5) * h;
    for (int i = 0; i < 16; i++) {
      const double r = mid + 0.5 * h * x[i];
      const double f = std::pow(r, k) * std::exp(-u * r * r - a * (r - A) * (r - A) - b * (r - B) * (r - B)) *
                       K(l1, 2.0 * a * A * r) * K(l2, 2.0 * b * B * r);
      sum += 0.5 * h * w[i] * f;
    }
  }
  return sum;
}

}  // namespace radref

#endif
