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
 * The shell and channel bounds are built from the shell-bounding Gaussians of McKenzie et al.
 * (S. C. McKenzie, E. Epifanovsky, G. M. J. Barca, A. T. B. Gilbert and P. M. W. Gill,
 * J. Phys. Chem. A 122 (2018), doi:10.1021/acs.jpca.7b12679); each function says where it follows
 * or departs from them.
 *
 * Notation throughout: a radial integrand
 * <!-- The formula is kept on one line: exhale, which renders this namespace description, cannot
 * handle formulae over several lines. -->
 * \f$G(k,l_1,l_2) = \int_0^\infty r^k e^{-nr^2-a(r-A)^2-b(r-B)^2} K_{l_1}(2aAr) K_{l_2}(2bBr) dr\f$
 * with \f$n\f$ an ECP exponent, \f$a\f$ and \f$b\f$ basis exponents, \f$A\f$ and \f$B\f$ the
 * distances of the shells from the ECP, and \f$K_l(z) = e^{-z} i_l(z)\f$ the scaled modified
 * spherical Bessel functions (\f$0 \le K_l \le 1\f$).
 */
namespace screening {

/// An interval [rmin, rmax] of the radial coordinate
struct Window {
  double rmin, rmax;
};

/**
 * Region \f$[r_\mathrm{min}, r_\mathrm{max}]\f$ outside which the integrand of
 * \f$G(k, l_1, l_2)\f$, for \f$k \in [k_\mathrm{min}, k_\mathrm{max}]\f$,
 * \f$l_1 \le l_{1,\mathrm{max}}\f$ and \f$l_2 \le l_{2,\mathrm{max}}\f$, is negligible:
 * \f[
 *   \begin{aligned}
 *   r_\mathrm{min} &= \max(0,\ r_s(k_\mathrm{min} - 2) - 7w), \qquad
 *   r_\mathrm{max} = r_s(k_\mathrm{max} + l_{1,\mathrm{max}} + l_{2,\mathrm{max}}) + 9w, \\
 *   r_s(m) &= \frac{p_t}{2} + \sqrt{\frac{p_t^2}{4} + \frac{m}{2 \zeta_t}}
 *     \quad (0 \text{ if the root is complex}), \qquad w = \frac{1}{\sqrt{\zeta_t}}, \\
 *   \zeta_t &= n + a + b, \qquad p_t = \frac{aA + bB}{\zeta_t}.
 *   \end{aligned}
 * \f]
 *
 * The Gaussian factor is exactly \f$e^{-\zeta_t (r - p_t)^2}\f$ up to a constant. The scaled Bessel
 * functions have logarithmic derivatives (in \f$r\f$) between \f$-1/r\f$ and \f$l/r\f$, so the
 * remaining factor can only move the maximum to between the stationary points \f$r_s(m)\f$ of
 * \f$e^{-\zeta_t (r - p_t)^2} r^m\f$ for \f$m = k_\mathrm{min} - 2\f$ and
 * \f$m = k_\mathrm{max} + l_{1,\mathrm{max}} + l_{2,\mathrm{max}}\f$. The window extends 7 and 9
 * Gaussian widths beyond those.
 *
 * A finite window is needed because the integrand can be far narrower than the spacing of a grid
 * on \f$[0, \infty)\f$ (width \f$1/\sqrt{\zeta_t} \approx 0.02\f$ for a tight core s primitive).
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
 * Rigorous upper bound on \f$|G(k, l_1, l_2)|\f$ for all \f$k \in [k_\mathrm{min},
 * k_\mathrm{max}]\f$ and all \f$l_1, l_2\f$:
 * \f[
 *   \begin{aligned}
 *   \mathrm{bound} &= e^{-D} \max_{k \in \{k_\mathrm{min}, k_\mathrm{max}\}}
 *     \left[ r_k^k\, e^{-\zeta_t (r_k - p_t)^2 / 2} \right] \sqrt{\frac{2\pi}{\zeta_t}}, \\
 *   \zeta_t &= n + a + b, \qquad p_t = \frac{aA + bB}{\zeta_t}, \qquad
 *   D = a A^2 + b B^2 - \zeta_t p_t^2 \ \ (\ge 0), \\
 *   r_k &= \frac{p_t}{2} + \sqrt{\frac{p_t^2}{4} + \frac{k}{\zeta_t}}
 *     \quad \text{(the maximiser of } r^k e^{-\zeta_t (r - p_t)^2 / 2} \text{)}.
 *   \end{aligned}
 * \f]
 *
 * Steps: \f$n r^2 + a (r - A)^2 + b (r - B)^2 = \zeta_t (r - p_t)^2 + D\f$ exactly;
 * \f$0 \le K_l \le 1\f$;
 * \f[
 *   \int_0^\infty r^k e^{-\zeta_t (r - p_t)^2} dr \le
 *   \max_r \left[ r^k e^{-\zeta_t (r - p_t)^2 / 2} \right] \int e^{-\zeta_t (r - p_t)^2 / 2} dr;
 * \f]
 * and \f$r^k \le \max(r^{k_\mathrm{min}}, r^{k_\mathrm{max}})\f$ for
 * \f$k \in [k_\mathrm{min}, k_\mathrm{max}]\f$ (pointwise, as the maximum of a convex function of
 * \f$k\f$ on an interval is at an end).
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
 * \f[
 *   \log \mathrm{bound} \le k_\mathrm{max} \max\left(0,\ p_t + \sqrt{k_\mathrm{max} / \zeta_t} - 1
 *   \right) + \frac{1}{2} \left( \frac{2\pi}{\zeta_t} - 1 \right) - D
 * \f]
 * from \f$\log x \le x - 1\f$, dropping \f$-\zeta_t (r_k - p_t)^2 / 2 \le 0\f$, and
 * \f$r_k \le p_t + \sqrt{k_\mathrm{max} / \zeta_t}\f$ for all \f$k \le k_\mathrm{max}\f$
 * (\f$\zeta_t\f$, \f$p_t\f$, \f$D\f$ and \f$r_k\f$ as in type2_bound).
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
 * \f[
 *   \mathrm{threshold} = \frac{\mathrm{tolerance}}{(2 + 2 \max(A, B))^L\, n_\mathrm{triples}},
 *   \qquad L = L_A + L_B
 * \f]
 *
 * The radials are multiplied downstream by products of binomial factors in the components of
 * \f$\mathbf{A}\f$ and \f$\mathbf{B}\f$, whose sums are
 * \f$\prod_c (1 + |A_c|)^{a_c} \le (1 + |\mathbf{A}|)^{L_A}\f$ (and likewise for \f$\mathbf{B}\f$);
 * the extra factor \f$2^L\f$ is headroom for the angular sums (a generous heuristic, not a proven
 * bound).
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
 * Upper bound on the contribution of each ECP channel \f$l\f$ (0..L, L = local) to the integrals
 * over a shell pair, used to skip channels below the tolerance. It bounds the full integral (basis
 * function polynomials included), so it holds after the downstream angular assembly.
 *
 * Notation: \f$A, B\f$ are the distances of the shells from the ECP, \f$L_A, L_B\f$ their angular
 * momenta; \f$\alpha, \beta\f$ their smallest exponents, \f$\alpha_i, \beta_j\f$ all exponents,
 * \f$c_i, c_j\f$ coefficients; \f$\eta\f$ the smallest exponent in channel \f$l\f$; ECP primitives
 * \f$d_k r^{n_k - 2} e^{-\eta_k r^2}\f$.
 * \f[
 *   \mathrm{result}_l = (2l + 1)^2\, a_\mathrm{bound}\, b_\mathrm{bound} \sum_k |d_k|
 *     \left(\frac{\pi}{\eta_k}\right)^{3/2} F_{n_k}\,
 *     \exp\left(\frac{x_p}{\zeta_k} - \tilde{a} A^2 - \tilde{b} B^2\right) S(T_k)
 * \f]
 * \f[
 *   \begin{aligned}
 *   a_\mathrm{bound} &= \sum_i |c_i| \left(\frac{L_A}{2 e \alpha_i \sigma_a}\right)^{L_A/2}
 *     \quad (1 \text{ per primitive if } L_A = 0), \\
 *   \sigma_a &= \begin{cases}
 *     \min\left(1,\ \dfrac{L_A (\alpha + \eta)^2}{2 \alpha (\eta^2 A^2 + L_A (\alpha + \eta))}
 *       \right) & A > 0 \\
 *     \min\left(1,\ \dfrac{\alpha + \eta}{2 \alpha}\right) & A = 0
 *   \end{cases} \\
 *   \tilde{a} &= (1 - \sigma_a) \alpha, \quad \text{and likewise } b_\mathrm{bound},\ \sigma_b,\
 *     \tilde{b} \text{ with } B,\ L_B,\ \beta, \\
 *   \zeta_k &= \tilde{a} + \tilde{b} + \eta_k, \qquad x_p = \tilde{a}^2 A^2 + \tilde{b}^2 B^2,
 *     \qquad T_k = \frac{2 \tilde{a} \tilde{b} A B}{\zeta_k}, \\
 *   S(T) &= \begin{cases} e^T / (2T) & T > 1 \\ \sinh 1 & \text{otherwise} \end{cases}
 *     \qquad (S(T) \ge \sinh(T) / T).
 *   \end{aligned}
 * \f]
 *
 * Steps: \f$|x|^{L_A} e^{-\alpha_i x^2} \le (L_A / (2 e \alpha_i \sigma))^{L_A/2}
 * e^{-(1 - \sigma) \alpha_i x^2}\f$ for any \f$\sigma \in (0, 1]\f$, and then
 * \f$e^{-(1 - \sigma) \alpha_i x^2} \le e^{-(1 - \sigma) \alpha x^2}\f$, which needs
 * \f$\sigma \le 1\f$ (hence the clamp; \f$\sigma\f$ minimises the bound otherwise).
 * The projector is bounded by Schwarz, and the angular integrals of the remaining Gaussians
 * exactly, leaving per ECP primitive the radial integral
 * \f[
 *   I_n = \frac{1}{k_A k_B} \int_0^\infty r^{n-2} e^{-\zeta r^2} \sinh(k_A r) \sinh(k_B r)\, dr,
 *   \qquad k_A = 2 \tilde{a} A, \quad k_B = 2 \tilde{b} B,
 * \f]
 * which for \f$n = 2\f$ is exactly
 * \f$\frac{\sqrt{\pi}}{4} \zeta^{-3/2} e^{x_p / \zeta} \sinh(T) / T\f$. For other \f$n\f$,
 * \f$F_n \ge I_n / I_2\f$ (with \f$T_c = T \coth T \ge 1\f$):
 * \f[
 *   F_n = \begin{cases}
 *     1 & n = 2 \\
 *     2 \zeta T_c & n = 0 \quad (\sinh x \le x \cosh x \text{ on both sinh factors}) \\
 *     \sqrt{2 \zeta T_c} & n = 1 \quad \text{(Cauchy-Schwarz between } n = 0 \text{ and } n = 2) \\
 *     (1/2 + x_p / \zeta + T_c) / \zeta & n = 4 \quad (\text{exact: } I_4 = -dI_2 / d\zeta) \\
 *     \sqrt{F_4} & n = 3 \quad \text{(Cauchy-Schwarz)}
 *   \end{cases}
 * \f]
 * and for any other \f$n\f$ there is no bound (infinite, never skipped).
 * All exponentials are combined into one so that none underflows on its own.
 *
 * This is the contracted three-centre projected bound \f$\langle B_3 \rangle_\ell\f$ of McKenzie
 * et al. (their eqs. 44-46), with three changes: the factors \f$F_n\f$, which extend it from
 * \f$r^0\f$ terms (our n = 2; they reconstruct other powers as sums of \f$r^0\f$ terms) to the
 * powers of real ECPs; \f$\sigma\f$ taken with the smallest exponent of each channel rather than of
 * the whole ECP; and the clamp \f$\sigma \le 1\f$ (see shell_bound).
 *
 * @param data - the shell pair, relative to the ECP centre (as built by compute_shell_pair)
 * @param results - array of U.getL() + 1 bounds, one per channel
 */
void channel_bounds(const ECP& U, const GaussianShell& shellA, const GaussianShell& shellB,
                    const ShellPairData& data, double* results);

/**
 * Shape factor bounding a primitive of a shell against an ECP, used for shell screening:
 * \f[
 *   \begin{aligned}
 *   \mathrm{shell\_bound} &= \left(\frac{l_a}{2 e \alpha \sigma}\right)^{l_a/2}
 *     \exp\left(-\frac{\tilde{a} \eta A^2}{\tilde{a} + \eta}\right),
 *     \qquad \tilde{a} = (1 - \sigma) \alpha, \\
 *   \sigma &= \begin{cases}
 *     \min\left(1,\ \dfrac{l_a (\alpha + \eta)^2}{2 \alpha (\eta^2 A^2 + l_a (\alpha + \eta))}
 *       \right) & A > 0 \\
 *     \min\left(1,\ \dfrac{\alpha + \eta}{2 \alpha}\right) & A = 0
 *   \end{cases}
 *   \end{aligned}
 * \f]
 *
 * It bounds \f$\max_r |\mathbf{r} - \mathbf{A}|^{l_a} e^{-\alpha |\mathbf{r} - \mathbf{A}|^2}
 * e^{-\eta r^2}\f$ for a shell at distance \f$A\f$ from the ECP:
 * \f$|x|^{l_a} e^{-\alpha x^2} \le (l_a / (2 e \alpha \sigma))^{l_a/2}
 * e^{-(1 - \sigma) \alpha x^2}\f$
 * for \f$\sigma \in (0, 1]\f$, and
 * \f$\max_r e^{-\tilde{a} |\mathbf{r} - \mathbf{A}|^2 - \eta r^2} =
 * e^{-\tilde{a} \eta A^2 / (\tilde{a} + \eta)}\f$.
 * \f$\sigma\f$ must not exceed 1: the leftover Gaussian would then grow and its maximum be
 * unbounded.
 * It is a shape factor only: the coefficients are included by shell_ecp_bound.
 *
 * This is the shell-bounding Gaussian of McKenzie et al. against the ECP Gaussian,
 * \f$[\mathring{a}]_L\f$ (their eqs. 24, 25 and 28), with their approximation to the optimal
 * \f$\sigma\f$ (eq. 31), clamped here: their eq. 31 exceeds 1 when \f$\eta > \alpha\f$ near the
 * ECP.
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
 * \f[
 *   \mathrm{shell\_ecp\_bound} = \mathrm{shell\_bound}(l + o,\ \alpha,\ A^2,\ \eta / 2)
 *     \sum_i |c_i| (2 a_i + l + 1)^{o} \sqrt{\textstyle\sum_k |d_k|}
 * \f]
 * with \f$o\f$ the derivative order (deriv).
 *
 * The integrand of a shell pair is bounded by
 * \f$|\phi_A e^{-\eta r^2 / 2}|\, |\phi_B e^{-\eta r^2 / 2}| \sum_k |d_k| r^{n_k}\f$,
 * so each shell takes half of the ECP's most diffuse Gaussian (taking all of
 * it for each shell underestimated pairs on an atom far from the ECP by ~1e3: an Au f shell with
 * an s shell on the same atom, 4.9 bohr from the other Au's ECP, estimated at 4e-14 against
 * integrals of 1.6e-12). shell_bound only describes the shape, so it is scaled by the size of the
 * contraction coefficients and (as a square root, for each of the two shells) of the ECP
 * coefficients, which can be ~1e3. Derivative integrals are integrals over shells of angular
 * momentum up to \f$l + o\f$, so for those the shape is taken at \f$l + o\f$; differentiating a
 * primitive gives terms with coefficients \f$2ac\f$ and \f$lc\f$, and twice gives coefficients
 * summing to at most \f$(2a + l)^2 + 2a\f$, so \f$(2a + l + 1)^o\f$ bounds the coefficient sum for
 * \f$o \le 2\f$. This is an estimate, not a proven bound.
 *
 * It is not from McKenzie et al.: their pair bound evaluates the exact three-centre Gaussian
 * product, which does not factorise into one value per shell, whereas this one is tabulated once
 * per (shell, ECP) and the pair test is a product of two table entries.
 */
double shell_ecp_bound(const GaussianShell& shell, const ECP& U, int deriv);

/// shell_ecp_bound for every shell and ECP, bounds[s * necp + i]
std::vector<double> shell_ecp_bounds(const std::vector<GaussianShell>& shells, const ECPBasis& ecps,
                                     int deriv);

/**
 * As shell_ecp_bounds, for the magnetic field derivatives (compute_shell_pair_derivative_Bfield),
 *
 * \f[
 *   \mathrm{bound} = (1 + |\mathbf{R} \times \mathbf{C}|)\, \mathrm{shell\_ecp\_bound}(o = 0)
 *     + |\mathbf{R} - \mathbf{C}|\, \mathrm{shell\_ecp\_bound}(o = 1)
 * \f]
 * for a shell at \f$\mathbf{R}\f$ and an ECP at \f$\mathbf{C}\f$. A block is
 * \f$[(\mathbf{A} - \mathbf{B}) \times \mathbf{C}] \langle A|U|B \rangle\f$ plus
 * \f$(\mathbf{A} - \mathbf{C}) \times \langle A+1|U|B \rangle\f$ and
 * \f$(\mathbf{B} - \mathbf{C}) \times \langle A|U|B+1 \rangle\f$, and the product of the bounds of
 * A and B covers each of these terms: the
 * raised shells are bounded by the first-derivative bounds, which take the shape at l + 1 and
 * coefficients at least as large. As for shell_ecp_bound, this is an estimate, not a proven bound.
 */
std::vector<double> shell_ecp_bfield_bounds(const std::vector<GaussianShell>& shells,
                                            const ECPBasis& ecps);

/**
 * The largest bound of any shell against each ECP,
 * \f$\mathrm{partner}_i = \max_s \mathrm{bounds}[s \cdot n_\mathrm{ecp} + i]\f$, for the shell test
 * of screen_shell_ecps.
 */
std::vector<double> partner_bounds(const std::vector<double>& bounds, int necp);

/**
 * The ECPs (index, bound) against which shell s is not screened out:
 * \f[
 *   \mathrm{bounds}[s \cdot n_\mathrm{ecp} + i] \times \mathrm{partner}_i > \mathrm{tolerance}.
 * \f]
 * A shell pair is screened out when the product of its two bounds is at most the tolerance, so a
 * shell skipped here would have every one of its pairs screened out too: the shell test only saves
 * work, and changes no result. (This takes the place of the shell test of McKenzie et al. (2018),
 * which compares a shell's shape factor with \f$\tau / Q\f$, Q bounding any partner shell over
 * all space. That bound does not apply to shell_ecp_bound, which includes the coefficients and
 * derivative factors, and the largest actual partner bound is both exact and cheap.)
 */
std::vector<std::pair<int, double>> screen_shell_ecps(const std::vector<double>& bounds, int s,
                                                      int necp, const std::vector<double>& partner,
                                                      double tolerance);

}  // namespace screening
}  // namespace libecpint

#endif
