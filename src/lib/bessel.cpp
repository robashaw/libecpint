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

#include "bessel.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>

#include "mathutil.hpp"

namespace libecpint {

// Constructor
BesselFunction::BesselFunction() {}
BesselFunction::BesselFunction(const int _lMax, const int _N, const int _order,
                               const double accuracy) {
  init(_lMax, _N, _order, accuracy);
}

void BesselFunction::init(const int _lMax, const int _N, const int _order, const double accuracy) {
  // Check parameters
  lMax = _lMax > -1 ? _lMax : 0;
  N = _N > 0 ? _N : 1;
  order = _order > 0 ? _order : 1;
  scale = N / 16.0;

  int lDim = lMax + TAYLOR_CUT + 1;
  K_stride = lDim;
  dK_stride1 = lDim;
  dK_stride0 = (TAYLOR_CUT + 1) * lDim;

  K.assign((N + 1) * lDim, 0.0);
  C.assign(lMax + TAYLOR_CUT, 0.0);
  dK.assign((N + 1) * (TAYLOR_CUT + 1) * lDim, 0.0);

  tabulate(accuracy);
}

// Tabulate the bessel function values
int BesselFunction::tabulate(const double accuracy) {
  int retval = 0;  // 0 for success, -1 for not converged
  // Series expansion for bessel function, K, is given by:
  // K_l(z) ~ z^l sum_{j=0 to infty} F_j(z) / (2j + 2l + 1)!!
  // where F_j(z) = e^(-z) * (z^2/2)^j / j!
  int lmax = lMax + TAYLOR_CUT;

  std::vector<double> F(order + 1);  // F_j above

  K[0] = 1.0;
  double z, z2;  // z and z^2 / 2
  double ratio;  // F_j(z) / (2j+1)!!
  for (int i = 0; i <= N; i++) {
    // Calculate K(z) at equally spaced points z = 16/N to 16
    z = i / (N / 16.0);
    z2 = z * z / 2.0;

    F[0] = exp(-z);
    ratio = F[0] / DFAC[0];
    K[i * K_stride + 0] = ratio;

    // Series expansion for K_0(z)
    int l = order;
    int j;
    for (j = 1; j <= l; j++) {
      if (ratio < accuracy) {
        // Reached convergence
        break;
      }

      F[j] = F[j - 1] * z2 / ((double)j);
      ratio = F[j] / DFAC[2 * j + 1];
      K[i * K_stride + 0] += ratio;
    }
    // if ( ratio > accuracy ) { retval = -1; break; } // Not converged

    // Calculate K_l from K_0
    z2 = z;
    for (l = 1; l <= lmax; l++) {
      ratio = 0;
      for (int m = 0; m < j; m++) ratio += F[m] / DFAC[2 * l + 2 * m + 1];
      K[i * K_stride + l] = z2 * ratio;
      z2 *= z;
    }
  }

  // Determine coefficients for derivative recurrence
  for (int i = 1; i < lmax; i++) C[i] = i / (2.0 * i + 1.0);

  // Determine the necessary derivatives from
  // K_l^(n+1) = C_l K_(l-1)^(n) + (C_l + 1/(2l+1))K_(l+1)^(n) - K_l^(n)
  for (int ix = 0; ix < N + 1; ix++) {
    int dK_ix = ix * dK_stride0;
    // Copy K values into dK
    for (int l = 0; l <= lMax + TAYLOR_CUT; l++) dK[dK_ix + l] = K[ix * K_stride + l];

    // Then the rest
    for (int n = 1; n < TAYLOR_CUT + 1; n++) {
      int dK_ix_n = dK_ix + n * dK_stride1;
      int dK_ix_nm1 = dK_ix + (n - 1) * dK_stride1;
      dK[dK_ix_n + 0] = dK[dK_ix_nm1 + 1] - dK[dK_ix_nm1 + 0];
      for (int l = 1; l <= lMax + TAYLOR_CUT - n; l++)
        dK[dK_ix_n + l] = C[l] * dK[dK_ix_nm1 + l - 1] +
                          (C[l] + 1.0 / (2.0 * l + 1.0)) * dK[dK_ix_nm1 + l + 1] -
                          dK[dK_ix_nm1 + l];
    }
  }

  // Store K^(n) / n!, so that the Taylor series is a polynomial in (z - z0) evaluated without
  // divisions
  for (int ix = 0; ix < N + 1; ix++) {
    double ofac = 1.0;
    for (int n = 1; n < TAYLOR_CUT + 1; n++) {
      ofac /= n;
      for (int l = 0; l < dK_stride1; l++) dK[ix * dK_stride0 + n * dK_stride1 + l] *= ofac;
    }
  }

  return retval;
}

// Get an upper bound for M_l(z)
double BesselFunction::upper_bound(const double z, const int L) const {
  // find nearest point (on left) in tabulated values
  int ix = std::floor(N * z / 16.0);
  int minix = L > 0 ? 1 : 0;
  ix = std::min(N, std::max(minix, ix));
  int lx = std::min(L, lMax);
  return K[ix * K_stride + lx];
}

// Calculate modified spherical Bessel function K_l(z), weighted with an exponential factor e^(-z)
// for l = 0 to lMax. This restricts K(z) to the interval [0,1].
void BesselFunction::calculate(const double z, int maxL, std::vector<double>& values) const {
  if (lMax < maxL) {
    std::cout << "Asked for " << maxL << " but only initialised to maximum L = " << lMax << "\n";
    maxL = lMax;
  }

  // Set K_0(z) = 1.0, and K_l(z) = 0.0 (for l != 0) if z <= 0
  if (z <= 0) values[0] = 1.0;
  // Zeroth order case
  // K_l(z) ~ (1-z)*z^l / (2l + 1)!!
  else if (z < SMALL) {
    values[0] = 1.0 - z;
    for (int l = 1; l <= maxL; l++) values[l] = values[l - 1] * z / (2.0 * l + 1.0);
  }
  // Large z case
  // K_l(z) = [R_l(-z) - (-1)^l e^{-2z} R_l(z)]/(2z), where R_l(z) = sum_{k=0 to l}
  // (l+k)!/[k!(l-k)!] (2z)^{-k}, by the recurrence of the scaled modified spherical Bessel
  // functions, K_{l+1} = K_{l-1} - (2l+1)/z K_l, upwards from exact K_0 and K_1 while 2l+1 <= z.
  // Beyond that, upwards amplifies errors by (2l+1)/z per order (1e-14 at l = 16, 1e-8 at l = 30
  // just above z = 16; dropping the e^{-2z} terms of K_0 and K_1 gives 2e-13 at l = 14), so the
  // higher orders come from the recurrence downwards, which is stable there (Miller's algorithm):
  // started from K_{top+1} = 0 far enough above lMax that the error this makes is damped by
  // prod ((2l+1)/z)^2 > 1e20, and normalised to the last order from upwards. top depends on lMax,
  // not maxL, so that each K_l is the same whatever maxL is asked for (shared_factor reuses
  // values computed for a larger maxL; with top from maxL, they differed in the last bit).
  else if (z > 16.0) {
    const double oz = 1.0 / z, e2 = std::exp(-2.0 * z);
    values[0] = 0.5 * oz * (1.0 - e2);
    if (maxL > 0) values[1] = 0.5 * oz * ((1.0 - oz) + e2 * (1.0 + oz));
    const int lup = std::min(maxL, static_cast<int>(0.5 * (z - 1.0)));
    for (int l = 1; l < lup; l++) values[l + 1] = values[l - 1] - (2 * l + 1) * oz * values[l];
    if (lup < maxL) {
      int top = lMax;
      for (double damp = 1.0; damp < 1e20; top++) damp *= std::pow((2 * top + 3) * oz, 2);
      double kp = 0.0, k = 1.0;  // K_{l+1} and K_l, up to a common factor
      for (int l = top; l > lup; l--) {
        if (l <= maxL) values[l] = k;
        const double km = kp + (2 * l + 1) * oz * k;
        kp = k;
        k = km;
      }
      const double scale = values[lup] / k;
      for (int l = lup + 1; l <= maxL; l++) values[l] *= scale;
    }
  }
  // SMALL < z < 16
  // Use Taylor series around pretabulated values in class
  // 5 terms is usually sufficient for machine accuracy
  else {
    // Index of abscissa z in table, and K(z) ~ sum_{n=0 to 5} [K^(n)(z0) / n!] (z-z0)^n by Horner
    const int ix = static_cast<int>(z * scale + 0.5);
    const double dz = z - ix / scale;  // z - z0
    const double* d = &dK[ix * dK_stride0];
    const int s1 = dK_stride1;
    for (int l = 0; l <= maxL; l++) {
      double t = d[TAYLOR_CUT * s1 + l];
      for (int n = TAYLOR_CUT - 1; n >= 0; n--) t = t * dz + d[n * s1 + l];
      values[l] = t;
    }
  }
}

// Calculate a modified spherical bessel function value at a point for only a single L
// method the same as in calculate for multiple L, but with efficiencies
double BesselFunction::calculate(const double z, const int L) const {
  double value = 0.0;

  if (z <= 0)
    value = 1.0;
  else if (z < SMALL) {
    // K_L(z) ~ (1 - z) z^L / (2L+1)!! = (1 - z) prod_{k=1}^{L} z/(2k+1)
    value = 1.0 - z;
    for (int k = 1; k < L + 1; k++) value *= z / (2.0 * k + 1.0);
  } else if (z > 16.0) {
    // as for all orders: the closed form R_L(-z)/(2z) cancels badly for L near z
    std::vector<double> values(L + 1);
    calculate(z, L, values);
    value = values[L];
  } else {
    const int ix = static_cast<int>(z * scale + 0.5);
    const double dz = z - ix / scale;  // z - z0
    const double* d = &dK[ix * dK_stride0 + L];
    value = d[TAYLOR_CUT * dK_stride1];
    for (int n = TAYLOR_CUT - 1; n >= 0; n--) value = value * dz + d[n * dK_stride1];
  }

  return value;
}
}  // namespace libecpint
