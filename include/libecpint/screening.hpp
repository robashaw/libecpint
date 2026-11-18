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

#ifndef SCREENING_HEAD
#define SCREENING_HEAD

#include <algorithm>
#include <cmath>
#include <utility>
#include <vector>

#include "ecp.hpp"
#include "gshell.hpp"
#include "mathutil.hpp"

namespace libecpint {

/**
 * Bounds and estimates used to skip negligible work, at every level: whole shells against an ECP
 * (shell_ecp_bound), shell pairs against an ECP channel (channel_bounds), and primitive triples in
 * the radial integrals (type2_bound). Kept apart from the integration code so that their
 * correctness (a bound must bound) can be tested separately from that of the integrals.
 *
 * Notation throughout: a radial integrand
 *   G(k, l1, l2) = int r^k exp(-n r^2 - a (r - A)^2 - b (r - B)^2) K_l1(2aAr) K_l2(2bBr) dr,
 * with n an ECP exponent, a and b basis exponents, A and B the distances of the shells from the
 * ECP, and K_l(z) = exp(-z) i_l(z) the scaled modified spherical Bessel functions (0 <= K_l <= 1).
 */
namespace screening {

/// An interval [rmin, rmax] of the radial coordinate
struct Window {
  double rmin, rmax;
};

/**
 * Region [rmin, rmax] outside which the integrand of G(k, l1, l2), for k in [kmin, kmax],
 * l1 <= l1max and l2 <= l2max, is negligible:
 *
 *   rmin = max(0, r_s(kmin - 2) - 7 w),  rmax = r_s(kmax + l1max + l2max) + 9 w,
 *   r_s(m) = pt/2 + sqrt(pt^2/4 + m / (2 zt))  (0 if the root is complex),  w = 1/sqrt(zt),
 *   zt = n + a + b,  pt = (aA + bB) / zt.
 *
 * The Gaussian factor is exactly exp(-zt (r - pt)^2) up to a constant. The scaled Bessel functions
 * have logarithmic derivatives (in r) between -1/r and l/r, so the remaining factor can only move
 * the maximum to between the stationary points r_s(m) of exp(-zt (r - pt)^2) r^m for m = kmin - 2
 * and m = kmax + l1max + l2max. The window extends 7 and 9 Gaussian widths beyond those.
 *
 * A finite window is needed because the integrand can be far narrower than the spacing of a grid
 * on [0, inf) (width 1/sqrt(zt) ~ 0.02 for a tight core s primitive).
 */
inline Window integrand_window(const int kmin, const int kmax, const int l1max, const int l2max,
                               const double n, const double a, const double b, const double A,
                               const double B) {
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

/**
 * Rigorous upper bound on |G(k, l1, l2)| for all k in [kmin, kmax] and all l1, l2:
 *
 *   bound = exp(-D) max_{k in {kmin, kmax}} [r_k^k exp(-zt (r_k - pt)^2 / 2)] sqrt(2 pi / zt)
 *   zt = n + a + b,  pt = (aA + bB) / zt,  D = a A^2 + b B^2 - zt pt^2  (>= 0)
 *   r_k = pt/2 + sqrt(pt^2/4 + k/zt)   (the maximiser of r^k exp(-zt (r - pt)^2 / 2))
 *
 * Steps: n r^2 + a (r - A)^2 + b (r - B)^2 = zt (r - pt)^2 + D exactly; 0 <= K_l <= 1;
 *   int_0^inf r^k e^(-zt (r - pt)^2) dr <= max_r [r^k e^(-zt (r - pt)^2 / 2)] int e^(-zt (r -
 *   pt)^2 / 2) dr;
 * and r^k <= max(r^kmin, r^kmax) for k in [kmin, kmax] (pointwise, as the maximum of a convex
 * function of k on an interval is at an end).
 */
inline double type2_bound(const int kmin, const int kmax, const double n, const double a,
                          const double b, const double A, const double B) {
  const double zt = n + a + b;
  const double pt = (a * A + b * B) / zt;
  const double D = a * A * A + b * B * B - zt * pt * pt;
  auto log_peak = [&](const int k) {
    const double rs = 0.5 * pt + std::sqrt(0.25 * pt * pt + k / zt);
    const double d = rs - pt;
    return (k > 0 ? k * std::log(rs) : 0.0) - 0.5 * zt * d * d;
  };
  return std::exp(std::max(log_peak(kmin), log_peak(kmax)) - D) * std::sqrt(2.0 * M_PI / zt);
}

/**
 * Cheap upper bound on log(type2_bound(kmin, kmax, n, a, b, A, B)), with one sqrt and no log or
 * exp, to discard most negligible primitive triples before the exact bound:
 *
 *   log bound <= kmax max(0, pt + sqrt(kmax / zt) - 1) + (2 pi / zt - 1) / 2 - D
 *
 * from log x <= x - 1, dropping -zt (r_k - pt)^2 / 2 <= 0, and r_k <= pt + sqrt(kmax / zt) for
 * all k <= kmax (zt, pt, D and r_k as in type2_bound).
 */
inline double type2_log_bound_cheap(const int kmax, const double n, const double a, const double b,
                                    const double A, const double B) {
  const double zt = n + a + b;
  const double pt = (a * A + b * B) / zt;
  const double D = a * A * A + b * B * B - zt * pt * pt;
  return kmax * std::max(0.0, pt + std::sqrt(kmax / zt) - 1.0) + 0.5 * (2.0 * M_PI / zt - 1.0) - D;
}

/**
 * Threshold for each primitive triple's bound (type2_bound times |coefficients|) in the type 1 and
 * type 2 radial integrals:
 *
 *   threshold = tolerance / ((2 + 2 max(A, B))^L ntriples),  L = LA + LB
 *
 * The radials are multiplied downstream by products of binomial factors in the components of A
 * and B, whose sums are prod_c (1 + |A_c|)^(a_c) <= (1 + |A|)^LA (and likewise for B); the extra
 * factor 2^L is headroom for the angular sums (a generous heuristic, not a proven bound).
 * Dividing by the number of triples keeps the neglected ones together below the tolerance.
 */
inline double triple_threshold(const double tolerance, const double A, const double B, const int L,
                               const int ntriples) {
  return tolerance / (std::pow(2.0 + 2.0 * std::max(A, B), L) * std::max(1, ntriples));
}

/// @return the number of primitives with non-zero coefficient in channel l of U
inline int count_channel(const ECP& U, const int l) {
  int n = 0;
  for (const auto& u : U.gaussians) n += u.l == l && u.d != 0.0;
  return n;
}

/**
 * Upper bound on the contribution of each ECP channel l (0..L, L = local) to the integrals over
 * a shell pair, used to skip channels below the tolerance. It bounds the full integral (basis
 * function polynomials included), so it holds after the downstream angular assembly.
 *
 * Notation: A, B are the distances of the shells from the ECP, LA, LB their angular momenta;
 * alpha, beta their smallest exponents, alpha_i, beta_j all exponents, c_i, c_j coefficients;
 * eta the smallest exponent in channel l; ECP primitives d_k r^(n_k - 2) exp(-eta_k r^2).
 *
 *   result_l = (2l + 1)^2 a_bound b_bound sum_k |d_k| (pi/eta_k)^(3/2) F_{n_k}
 *                 exp(xp/zeta_k - at A^2 - bt B^2) S(T_k)
 *
 *   a_bound = sum_i |c_i| (LA / (2 e alpha_i sigma_a))^(LA/2)   (1 per primitive if LA = 0)
 *   sigma_a = min(1, LA (alpha + eta)^2 / (2 alpha (eta^2 A^2 + LA (alpha + eta))))  [A > 0]
 *           = min(1, (alpha + eta) / (2 alpha))                                       [A = 0]
 *   at = (1 - sigma_a) alpha,  and likewise b_bound, sigma_b, bt with B, LB, beta
 *   zeta_k = at + bt + eta_k,  xp = at^2 A^2 + bt^2 B^2,  T_k = 2 at bt A B / zeta_k
 *   S(T) = exp(T) / (2T) for T > 1, sinh(1) otherwise   (S(T) >= sinh(T)/T)
 *
 * Steps: |x|^LA exp(-alpha_i x^2) <= (LA/(2 e alpha_i sigma))^(LA/2) exp(-(1 - sigma) alpha_i
 * x^2) for any sigma in (0, 1], and then exp(-(1 - sigma) alpha_i x^2) <= exp(-(1 - sigma)
 * alpha x^2), which needs sigma <= 1 (hence the clamp; sigma minimises the bound otherwise).
 * The projector is bounded by Schwarz, and the angular integrals of the remaining Gaussians
 * exactly, leaving per ECP primitive the radial integral
 *   I_n = 1/(kA kB) int_0^inf r^(n-2) exp(-zeta r^2) sinh(kA r) sinh(kB r) dr,
 *   kA = 2 at A, kB = 2 bt B,
 * which for n = 2 is exactly (sqrt(pi)/4) zeta^(-3/2) exp(xp/zeta) sinh(T)/T. For other n,
 * F_n >= I_n / I_2 (with Tc = T coth T >= 1):
 *   n = 2: F = 1;  n = 0: F = 2 zeta Tc  (sinh x <= x cosh x on both sinh factors);
 *   n = 1: F = sqrt(2 zeta Tc)  (Cauchy-Schwarz between n = 0 and n = 2);
 *   n = 4: F = (1/2 + xp/zeta + Tc) / zeta  (exact: I_4 = -dI_2/dzeta);
 *   n = 3: F = sqrt(F_4)  (Cauchy-Schwarz);  any other n: no bound (infinite, never skipped).
 * All exponentials are combined into one so that none underflows on its own.
 *
 * @param data - the shell pair, relative to the ECP centre (as built by compute_shell_pair)
 * @param results - array of U.getL() + 1 bounds, one per channel
 */
void channel_bounds(const ECP& U, const GaussianShell& shellA, const GaussianShell& shellB,
                    const ShellPairData& data, double* results);

/**
 * Shape factor bounding a primitive of a shell against an ECP, used for shell screening:
 *
 *   shell_bound = (la / (2 e alpha sigma))^(la/2) exp(-at eta A^2 / (at + eta)),
 *   at = (1 - sigma) alpha,
 *   sigma = min(1, la (alpha + eta)^2 / (2 alpha (eta^2 A^2 + la (alpha + eta))))   [A > 0]
 *         = min(1, (alpha + eta) / (2 alpha))                                        [A = 0]
 *
 * It bounds max_r |r - A|^la exp(-alpha |r - A|^2) exp(-eta r^2) for a shell at distance A from
 * the ECP: |x|^la exp(-alpha x^2) <= (la / (2 e alpha sigma))^(la/2) exp(-(1 - sigma) alpha x^2)
 * for sigma in (0, 1], and max_r exp(-at |r - A|^2 - eta r^2) = exp(-at eta A^2 / (at + eta)).
 * sigma must not exceed 1: the leftover Gaussian would then grow and its maximum be unbounded.
 * It is a shape factor only: the coefficients are included by shell_ecp_bound.
 * @param la - angular momentum of the shell
 * @param alpha - smallest exponent of the shell
 * @param A2 - squared distance of the shell from the ECP
 * @param eta - smallest exponent of the ECP
 */
double shell_bound(int la, double alpha, double A2, double eta);

/// @return the squared distance between the centres of shell and U
inline double shell_ecp_distance2(const GaussianShell& shell, const ECP& U) {
  double dx = shell.center()[0] - U.center_[0];
  double dy = shell.center()[1] - U.center_[1];
  double dz = shell.center()[2] - U.center_[2];
  return dx * dx + dy * dy + dz * dz;
}

/**
 * Estimate of the size of a shell's contribution to integrals over an ECP (or to their derivatives
 * of order deriv), such that the product for two shells estimates their integral:
 *
 *   shell_ecp_bound = shell_bound(l + deriv, alpha, A^2, eta / 2) sum_i |c_i| (2 a_i + l + 1)^deriv
 *                     sqrt(sum_k |d_k|)
 *
 * The integrand of a shell pair is bounded by |phi_A exp(-eta r^2 / 2)| |phi_B exp(-eta r^2 / 2)|
 * sum_k |d_k| r^n_k, so each shell takes half of the ECP's most diffuse Gaussian (taking all of
 * it for each shell underestimated pairs on an atom far from the ECP by ~1e3: an Au f shell with
 * an s shell on the same atom, 4.9 bohr from the other Au's ECP, estimated at 4e-14 against
 * integrals of 1.6e-12). shell_bound only describes the shape, so it is scaled by the size of the
 * contraction coefficients and (as a square root, for each of the two shells) of the ECP
 * coefficients, which can be ~1e3. Derivative integrals are integrals over shells of angular
 * momentum up to l + deriv, so for those the shape is taken at l + deriv; differentiating a
 * primitive gives terms with coefficients 2 a c and l c, and twice gives coefficients summing to at
 * most (2 a + l)^2 + 2 a, so (2 a + l + 1)^deriv bounds the coefficient sum for deriv <= 2. This is
 * an estimate, not a proven bound.
 */
double shell_ecp_bound(const GaussianShell& shell, const ECP& U, int deriv);

/// shell_ecp_bound for every shell and ECP, bounds[s * necp + i]
std::vector<double> shell_ecp_bounds(const std::vector<GaussianShell>& shells, const ECPBasis& ecps,
                                     int deriv);

/**
 * As shell_ecp_bounds, for the magnetic field derivatives (compute_shell_pair_derivative_Bfield),
 *
 *   bound = (1 + |R x C|) shell_ecp_bound(deriv = 0) + |R - C| shell_ecp_bound(deriv = 1),
 *
 * for a shell at R and an ECP at C. A block is [(A - B) x C] <A|U|B> plus (A - C) x <A+1|U|B> and
 * (B - C) x <A|U|B+1>, and the product of the bounds of A and B covers each of these terms: the
 * raised shells are bounded by the first-derivative bounds, which take the shape at l + 1 and
 * coefficients at least as large. As for shell_ecp_bound, this is an estimate, not a proven bound.
 */
std::vector<double> shell_ecp_bfield_bounds(const std::vector<GaussianShell>& shells,
                                            const ECPBasis& ecps);

/// @return the ECPs (index, bound) for which the bound of shell s exceeds thresh
std::vector<std::pair<int, double>> screen_shell_ecps(const std::vector<double>& bounds, int s,
                                                      int necp, double thresh);

/**
 * Threshold on shell_ecp_bound below which a whole shell is skipped, given the shell-pair
 * tolerance:
 *
 *   shell_tolerance = tolerance / sqrt(((L + 3) / alpha)^(L + 3) (pi / (2L + 3))^3 / (2 e)^L)
 *
 * with L the largest angular momentum and alpha the smallest exponent in the basis.
 */
double shell_tolerance(double tolerance, int maxLB, double min_alpha);

}  // namespace screening
}  // namespace libecpint

#endif
