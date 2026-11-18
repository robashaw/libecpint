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

#ifndef ECPINT_HEAD
#define ECPINT_HEAD

#include <array>
#include <vector>

#include "angular.hpp"
#include "bessel.hpp"
#include "config.hpp"
#include "ecp.hpp"
#include "gaussquad.hpp"
#include "gshell.hpp"
#include "multiarr.hpp"
#include "radial.hpp"

namespace libecpint {

/// returns the index of the basis function with y,z angular momentum components l,m
#define N_INDEX(l, m) (((l + m) * (l + m + 1)) / 2 + m)

/**
 * \class ECPIntegral
 * \brief Calculates ECP integrals.
 *
 * Given an ECP basis, and orbital bases, this will calculate the ECP integrals over all ECP
 * centers.
 *
 * REFERENCES:
 * (Shaw2017) R. A. Shaw, J. G. Hill, J. Chem. Phys. 147 (2017), 074108
 * (Flores06) R. Flores-Moreno et al., J. Comput. Chem. 27 (2006), 1009
 * (MM81) L. E. McMurchie and E. R. Davidson, J. Comp. Phys. 44 (1981), 289 - 301
 */
class ECPIntegral {
 private:
  RadialIntegral radInts;   ///< The interface to the radial integral calculation
  AngularIntegral angInts;  ///< The angular integrals, which can be reused over all ECP centers

  /// Screening tolerance: channel contributions whose upper bound falls below 1e-2 times this,
  /// and primitive triples whose bound falls below their share of it, are skipped. Set from the
  /// constructor `thresh` so the caller can tighten it; a fixed value here would silently drop
  /// small but genuine l>0 integrals.
  double tolerance;

  /// Worker functions for calculating binomial expansion coefficients
  double calcC(int a, int m, double A) const;

  /// The radial integrals (N, l1, l2) needed for the type 2 integrals over shells of angular
  /// momenta LA <= LB in ECP channel lam, those whose angular integrals do not vanish: A with
  /// l1 <= l2, B the others with l1 and l2 swapped; nbase as in RadialIntegral::type2
  struct Type2Radials {
    std::vector<Triple> A, B;
    int nbase = 0;
  };
  int type2LB;  ///< largest basis angular momentum in type2Radials (maxLB + deriv)
  int type2LU;  ///< number of ECP channels in type2Radials (maxLU)
  /// Type2Radials for every (LA <= LB, lam), at (LA * (type2LB + 1) + LB) * type2LU + lam, found
  /// once at construction from the angular integrals
  std::vector<Type2Radials> type2Radials;
  /// Finds type2Radials, from the angular integrals
  void buildType2Radials();

 public:
  int skipped, zero, nonzero;

  /**
   * Constructs the coefficients in the binomial expansion (see REF. Shaw2017)
   * @param C - reference to a FiveIndex array to store the results in
   * @param L - maximum angular momentum to go up to in expansion
   * @param A - xyz coordinates for the center to calculate over
   */
  void makeC(FiveIndex<double>& C, int L, const double* A) const;

  /**
   * Creates an ECP integrator, initialising the radial and angular parts for subsequent
   * calculations.
   * @param maxLB - the maximum angular momentum in the orbital basis
   * @param maxLU - the maximum angular momentum in the ECP basis
   * @param deriv - the maximum order of derivative to be calculated (TODO: derivs currently being
   * implemented)
   * @param thresh - the tolerance for convergence of radial integrals (defaults to 1e-15)
   * @param small - the maximum number of quadrature points for the small radial integration grid
   * (default 256, minimum recommended)
   * @param large - the maximum number of quadrature points for the large radial integration grid
   * (default 1024, minimum recommended)
   * @throws std::invalid_argument if maxLB + deriv or maxLU exceeds LIBECPINT_MAX_L
   */
  ECPIntegral(int maxLB, int maxLU, int deriv = 0, double thresh = 1e-15, unsigned smallGrid = 256,
              unsigned bigGrid = 1024);

  /**
   * Calculates the type 1 integrals for the given ECP center over the given shell pair, using
   * quadrature
   * @param U - reference to ECP
   * @param shellA - the first basis shell (rows in values)
   * @param shellB - the second basis shell (cols in values)
   * @param data - wrapper for data about shell pair
   * @param CA - binomial expansion coefficients for shellA, made with makeC
   * @param CB - binomial expansion coefficients for shellB, made with makeC
   * @param parameters - pre-calculated parameters for the radial integral
   * @param values - array in which results are returned
   */
  void type1(const ECP& U, const GaussianShell& shellA, const GaussianShell& shellB,
             const ShellPairData& data, const FiveIndex<double>& CA, const FiveIndex<double>& CB,
             const RadialIntegral::Parameters& parameters, TwoIndex<double>& values) const;

  /**
   * Calculates the type 2 integrals for the given ECP center over the given shell pair
   * @param l - angular momentum shell of ECP to calculate over
   * @param U - reference to ECP
   * @param shellA - the first basis shell (rows in values)
   * @param shellB - the second basis shell (cols in values)
   * @param data - wrapper for data about shell pair
   * @param CA - binomial expansion coefficients for shellA, made with makeC
   * @param CB - binomial expansion coefficients for shellB, made with makeC
   * @param parameters - pre-calculated parameters for the radial integral
   * @param values - array in which results are returned
   */
  void type2(int l, const ECP& U, const GaussianShell& shellA, const GaussianShell& shellB,
             const ShellPairData& data, const FiveIndex<double>& CA, const FiveIndex<double>& CB,
             const RadialIntegral::Parameters& parameters, ThreeIndex<double>& values) const;

  /**
   * Upper bound on the contribution of each ECP channel l (0..L, L = local) to the integrals over
   * a shell pair, used to skip channels below 1e-2 times the tolerance (compute_shell_pair); see
   * screening::channel_bounds in
   * screening.hpp for the formula.
   *
   * @param results - array of U.getL() + 1 bounds, one per channel
   */
  void estimate_type2(const ECP& U, const GaussianShell& shellA, const GaussianShell& shellB,
                      const ShellPairData& data, double* results) const;

  /**
   * Computes the overall ECP integrals over the given ECP center and shell pair. This is the lower
   * level API, where you want finer control over the calculation. Results are returned with rows
   * corresponding to shellA and cols to shellB, with the Cartesian functions in alpha order e.g.
   * {xxx, xxy, xxz, xyy, xyz, xzz, yyy, yyz, yzz, zzz} l = 3
   *
   * @param U - reference to the ECP to calculate the integral over
   * @param shellA - the first basis shell (rows in values)
   * @param shellB - the second basis shell (cols in values)
   * @param values - reference to TwoIndex array where the results will be stored
   */
  void compute_shell_pair(const ECP& U, const GaussianShell& shellA, const GaussianShell& shellB,
                          TwoIndex<double>& values, int shiftA = 0, int shiftB = 0) const;

  /**
   * Computes the overall ECP integral first derivatives over the given ECP center, C, and shell
   * pair (A | B) The results are placed in order [Ax, Ay, Az, Bx, By, Bz, Cx, Cy, Cz] and are
   * calculated so that each component can always be added to the relevant total derivative. E.g. if
   * A = B, then the contribution to the total derivative for that coordinate on the x axis will be
   * Ax + Bx. The order for each derivative matrices matches that specified in compute_shell_pair
   *
   * @param U - reference to the ECP
   * @param shellA - the first basis shell (rows in values)
   * @param shellB - the second basis shell (cols in values)
   * @param results - reference to array of 9 TwoIndex arrays where the results will be stored
   */
  void compute_shell_pair_derivative(const ECP& U, const GaussianShell& shellA,
                                     const GaussianShell& shellB,
                                     std::array<TwoIndex<double>, 9>& results) const;

  /**
   * Computes the overall ECP integral second derivatives over the given ECP center, C, and shell
   * pair (A | B) The results are placed in order [AA, AB, AC, BB, BC, CC] with components [xx, xy,
   * xz, yy, yz, zz] for AA, BB, and CC, and [xx, xy, xz, yx, yy, yz, zx, zy, zz] for AB, AC, and
   * BC. As for the first derivatives, the components are calculated such that they can usually be
   * added to the relevant total derivative. However, this is more complicated than for first
   * derivatives, especially in the instance where A=B. It's recommended to look at the
   * compute_second_derivatives interface in api.cpp for how to handle this. The order for each
   * derivative matrices matches that specified in compute_shell_pair
   *
   * @param U - reference to the ECP
   * @param shellA - the first basis shell (rows in values)
   * @param shellB - the second basis shell (cols in values)
   * @param results - reference to array of 45 TwoIndex arrays where the results will be stored
   */
  void compute_shell_pair_second_derivative(const ECP& U, const GaussianShell& shellA,
                                            const GaussianShell& shellB,
                                            std::array<TwoIndex<double>, 45>& results) const;

  /**
   * Worker function to calculate the derivative of the integral <A | C | B> with respect to A.
   * This is given as <d_q A(l_q) | C | B> = l_q*<A(l_q-1) | C | B> - 2*mu*<A(l_q + 1) | C | B>
   * where l_q is the angular momentum component of A in the q coordinate, and mu is the exponent of
   * A.
   *
   * @param U - reference to the ECP
   * @param shellA - the first basis shell (rows in values)
   * @param shellB - the second basis shell (cols in values)
   * @param results - reference to array of 3 TwoIndex arrays for the [x, y, z] derivatives
   */
  void left_shell_derivative(const ECP& U, const GaussianShell& shellA, const GaussianShell& shellB,
                             std::array<TwoIndex<double>, 3>& results) const;

  /**
   * Worker function to calculate the second derivatives of the integral <A | C | B> with respect to
   * AA. This is given as <d_p d_q A(l_p, l_q) | C | B> = l_p*l_q*<A(l_p-1, l_q-1) | C | B> -
   * 2*mu*l_p*<A(l_p-1, l_q+1) | C | B>
   * - 2*mu*l_q*<A(l_p+1, l_q-1) | C | B> + 4*mu^2*<A(l_p+1, l_q+1) | C | B >
   *
   * @param U - reference to the ECP
   * @param shellA - the first basis shell (rows in values)
   * @param shellB - the second basis shell (cols in values)
   * @param results - reference to array of 6 TwoIndex arrays for the [xx, xy, xz, yy, yz, zz]
   * derivatives
   */
  void left_shell_second_derivative(const ECP& U, const GaussianShell& shellA,
                                    const GaussianShell& shellB,
                                    std::array<TwoIndex<double>, 6>& results) const;

  /**
   * Worker function to calculate the second derivatives of the integral <A | C | B> with respect to
   * AB. This is given as <d_p A(l_p) | C | d_q B(l_q)> = l_p*l_q*<A(l_p-1) | C | B(l_q-1)> -
   * 2*mu_B*l_p*<A(l_p-1) | C | B(l_q+1)>
   * - 2*mu_A*l_q*<A(l_p+1) | C | B(l_q-1)> + 4*mu_A*mu_B*<A(l_p+1) | C | B(l_q+1) >
   *
   * @param U - reference to the ECP
   * @param shellA - the first basis shell (rows in values)
   * @param shellB - the second basis shell (cols in values)
   * @param results - reference to array of 9 TwoIndex arrays for the [xx, xy, xz, yx, yy, yz, zx,
   * zy, zz] derivatives
   */
  void mixed_second_derivative(const ECP& U, const GaussianShell& shellA,
                               const GaussianShell& shellB,
                               std::array<TwoIndex<double>, 9>& results) const;
  /**
   * Computes the first derivatives with respect to a uniform external magnetic field of the ECP
   * integral over the ECP centre, C, and the shell pair (A | B) of gauge-including atomic orbitals
   * (London orbitals, chi_A exp(-i/2 (B x A).r)), with the ECP's projectors gauge-transformed to
   * its own centre. For k = x, y, z,
   *
   *   results[k] = <A| [(A - C) x r]_k U |B> - <A| U [(B - C) x r']_k |B>,
   *
   * with r (r') the position, relative to the global origin, on the bra (ket) side of the
   * semi-local projectors. The derivative is i/2 results[k], so each matrix is real and
   * antisymmetric (results[k] for (B | A) is minus the transpose), and depends on the origin as
   * GIAO derivatives do. A block with both shells on the ECP is zero. The order of each matrix
   * matches that of compute_shell_pair.
   *
   * The shells are raised by one in angular momentum, so this needs an ECPIntegral constructed
   * with deriv >= 1.
   *
   * @param U - reference to the ECP
   * @param shellA - the first basis shell (rows in values)
   * @param shellB - the second basis shell (cols in values)
   * @param results - reference to array of 3 TwoIndex arrays for the [x, y, z] derivatives
   */
  void compute_shell_pair_derivative_Bfield(const ECP& U, const GaussianShell& shellA,
                                            const GaussianShell& shellB,
                                            std::array<TwoIndex<double>, 3>& results) const;
};
}  // namespace libecpint
#endif
