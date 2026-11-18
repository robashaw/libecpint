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

#include "screening.hpp"

#include <cmath>
#include <limits>

#include "mathutil.hpp"

namespace libecpint {
namespace screening {

void channel_bounds(const ECP& U, const GaussianShell& shellA, const GaussianShell& shellB,
                    const ShellPairData& data, double* results) {
  // See the formula in screening.hpp; variable names follow it
  const double alpha = shellA.min_exp, beta = shellB.min_exp;
  // a_bound = (LA / (2 e sigma_a))^(LA/2) sum_i |c_i| alpha_i^(-LA/2), and likewise b_bound; the
  // sums do not depend on l
  auto coef_sum = [](const GaussianShell& shell, const int L) {
    double sum = 0.0;
    for (int i = 0; i < shell.nprimitive(); i++)
      sum += std::abs(shell.coef(i)) * (L == 0 ? 1.0 : std::pow(shell.exp(i), -0.5 * L));
    return sum;
  };
  const double sum_a = coef_sum(shellA, data.LA), sum_b = coef_sum(shellB, data.LB);
  for (int l = 0; l <= U.getL(); l++) {
    const double eta = U.min_exp_l[l];
    const double an = alpha + eta, bn = beta + eta;
    // sigma in (0, 1]: for sigma > 1 the leftover Gaussian exp(-(1 - sigma) alpha_i x^2) grows,
    // and bounding it by the smallest exponent would no longer be an upper bound
    double sigma_a = data.A2 < 1e-6
                         ? 0.5 * an / alpha
                         : 0.5 * data.LA * an * an / (alpha * (eta * eta * data.A2 + data.LA * an));
    double sigma_b = data.B2 < 1e-6
                         ? 0.5 * bn / beta
                         : 0.5 * data.LB * bn * bn / (beta * (eta * eta * data.B2 + data.LB * bn));
    sigma_a = std::min(sigma_a, 1.0);
    sigma_b = std::min(sigma_b, 1.0);
    const double at = (1.0 - sigma_a) * alpha, bt = (1.0 - sigma_b) * beta;

    const double a_bound =
        sum_a * (data.LA == 0 ? 1.0 : std::pow(data.LA / (2.0 * M_EULER * sigma_a), 0.5 * data.LA));
    const double b_bound =
        sum_b * (data.LB == 0 ? 1.0 : std::pow(data.LB / (2.0 * M_EULER * sigma_b), 0.5 * data.LB));

    const double xp = at * at * data.A2 + bt * bt * data.B2;
    const double T0 = 2.0 * at * bt * data.Am * data.Bm;
    double ab_bound = 0.0;
    for (int k = U.l_starts[l]; k < U.l_starts[l + 1]; k++) {
      const GaussianECP& g = U.getGaussian(k);
      const double zeta = at + bt + g.a, T = T0 / zeta;
      // log S(T), S(T) >= sinh(T)/T; and T coth(T)
      const double logS = T > 1.0 ? T - std::log(2.0 * T) : std::log(SINH_1);
      const double Tc = T > 1e-4 ? T / std::tanh(T) : 1.0 + T * T / 3.0;
      // F_n >= I_n / I_2 for the ECP power r^(n-2) (g.n holds n - 2)
      double F;
      switch (g.n) {
        case 0:
          F = 1.0;
          break;  // n = 2, exact
        case -2:
          F = 2.0 * zeta * Tc;
          break;  // n = 0
        case -1:
          F = std::sqrt(2.0 * zeta * Tc);
          break;  // n = 1
        case 2:
          F = (0.5 + xp / zeta + Tc) / zeta;
          break;  // n = 4, exact
        case 1:
          F = std::sqrt((0.5 + xp / zeta + Tc) / zeta);
          break;  // n = 3
        default:
          F = std::numeric_limits<double>::infinity();  // no bound: never screen
      }
      // one exponential, so that no factor underflows or overflows on its own
      ab_bound += std::abs(g.d) * FAST_POW[3](std::sqrt(M_PI / g.a)) * F *
                  std::exp(xp / zeta - at * data.A2 - bt * data.B2 + logS);
    }
    results[l] = (2 * l + 1) * (2 * l + 1) * a_bound * b_bound * ab_bound;
  }
}

double shell_bound(const int la, const double alpha, const double A2, const double eta) {
  double sigma;
  if (A2 < 1e-6) {
    sigma = 0.5 * (1.0 + eta / alpha);
  } else {
    sigma = 1.0 / (2.0 * alpha * (eta * eta * A2 + la * (alpha + eta)));
    sigma = sigma * la * (alpha + eta) * (alpha + eta);
  }
  sigma = std::min(sigma, 1.0);  // the bound needs sigma <= 1

  double atilde = (1.0 - sigma) * alpha;
  double Na = la / (2 * M_EULER * alpha * sigma);
  Na = FAST_POW[la](std::sqrt(Na));
  double result = atilde * eta * A2 / (atilde + eta);
  result = std::exp(-result) * Na;
  return result;
}

double shell_ecp_bound(const GaussianShell& shell, const ECP& U, const int deriv) {
  double csum = 0.0;
  for (int i = 0; i < shell.nprimitive(); i++)
    csum += std::abs(shell.coef(i)) * std::pow(2.0 * shell.exp(i) + shell.l + 1.0, deriv);
  // sum_k |d_k| exp(-eta_k r^2) phi_A phi_B <= sqrt(sum_k |d_k| (phi_A e_k)^2)
  // sqrt(sum_k |d_k| (phi_B e_k)^2) with e_k = exp(-eta_k r^2 / 2) (Cauchy-Schwarz), so each shell
  // takes half of each ECP primitive's Gaussian
  const double A2 = shell_ecp_distance2(shell, U);
  double sum = 0.0;
  for (const auto& g : U.gaussians)
    if (g.d != 0.0)
      sum +=
          std::abs(g.d) * std::pow(shell_bound(shell.l + deriv, shell.min_exp, A2, 0.5 * g.a), 2);
  return csum * std::sqrt(sum);
}

std::vector<double> shell_ecp_bounds(const std::vector<GaussianShell>& shells, const ECPBasis& ecps,
                                     const int deriv) {
  const int necp = ecps.getN();
  std::vector<double> bounds(shells.size() * necp);
  for (size_t s = 0; s < shells.size(); s++)
    for (int i = 0; i < necp; i++)
      bounds[s * necp + i] = shell_ecp_bound(shells[s], ecps.getECP(i), deriv);
  return bounds;
}

std::vector<double> shell_ecp_bfield_bounds(const std::vector<GaussianShell>& shells,
                                            const ECPBasis& ecps) {
  const int necp = ecps.getN();
  std::vector<double> bounds(shells.size() * necp);
  for (size_t s = 0; s < shells.size(); s++)
    for (int i = 0; i < necp; i++) {
      const ECP& U = ecps.getECP(i);
      const double *R = shells[s].center(), *C = U.center();
      const double RxC[3] = {R[1] * C[2] - R[2] * C[1], R[2] * C[0] - R[0] * C[2],
                             R[0] * C[1] - R[1] * C[0]};
      const double cross = std::sqrt(RxC[0] * RxC[0] + RxC[1] * RxC[1] + RxC[2] * RxC[2]);
      bounds[s * necp + i] =
          (1.0 + cross) * shell_ecp_bound(shells[s], U, 0) +
          std::sqrt(shell_ecp_distance2(shells[s], U)) * shell_ecp_bound(shells[s], U, 1);
    }
  return bounds;
}

std::vector<std::pair<int, double>> screen_shell_ecps(const std::vector<double>& bounds,
                                                      const int s, const int necp,
                                                      const double thresh) {
  std::vector<std::pair<int, double>> result;
  for (int i = 0; i < necp; i++)
    if (bounds[s * necp + i] > thresh) result.push_back({i, bounds[s * necp + i]});
  return result;
}

double shell_tolerance(const double tolerance, const int maxLB, const double min_alpha) {
  double denominator =
      FAST_POW[maxLB + 3]((maxLB + 3.0) / min_alpha) * FAST_POW[3](M_PI / (2 * maxLB + 3.0));
  denominator /= FAST_POW[maxLB](2.0 * M_EULER);
  return tolerance / std::sqrt(denominator);
}

}  // namespace screening
}  // namespace libecpint
