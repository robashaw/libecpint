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

#ifndef RADIAL_HEAD
#define RADIAL_HEAD

#include <vector>

#include "bessel.hpp"
#include "ecp.hpp"
#include "gshell.hpp"
#include "multiarr.hpp"

namespace libecpint {

constexpr double MIN_EXP = 0.002;
/// Accuracy to which the tabulated Bessel functions are converged, whatever the tolerance
constexpr double BESSEL_ACCURACY = 1e-15;
/**
 * \class RadialIntegral
 * \brief Abstracts the calculation of radial integrals for ECP integration.
 *
 * This should not be used directly, and is owned by ECPIntegral.
 * It calculates the type 1 and type 2 radial integrals by fixed Gauss-Legendre quadratures.
 */
class RadialIntegral {
 private:
  /// Modified spherical Bessel function of the first kind
  BesselFunction bessie;

  /// Tolerance for the integrals; contributions provably below it are neglected
  double tolerance;

  double calcKij(double Na, double Nb, double zeta_a, double zeta_b, double R2) const;

 public:
  /// Default constructor creates an empty object
  RadialIntegral();

  /**
   * Initialises the object, in turn intialising the BesselFunction
   * @param maxL - the maximum angular momentum of integral needed
   * @param tol - the tolerance for the integrals (defaults to 1e-15): contributions provably below
   * it are skipped, and 0 skips nothing. It does not affect the Bessel functions, which are always
   * tabulated to BESSEL_ACCURACY.
   * @param small - unused (formerly the size of the small adaptive quadrature grid)
   * @param large - unused (formerly the size of the large adaptive quadrature grid)
   */
  void init(int maxL, double tol = 1e-15, int small = 256, int large = 1024);

  /// struct to store all parameters needed in both type 1 and 2 integrations
  struct Parameters {
    /// Matrices of parameters needed in both type 1 and 2 integrations
    TwoIndex<double> p, P, P2, K;
  };
  /**
   * Given two GaussianShells, builds the parameters needed by both kind of integral.
   * @param shellA - the first GaussianShell
   * @param shellB - the second GaussianShell
   * @param data - the data container for the shell pair
   * @returns the parameters needed in both type 1 and 2 integrations
   */
  Parameters buildParameters(const GaussianShell& shellA, const GaussianShell& shellB,
                             const ShellPairData& data) const;

  /**
   * Calculates all type 1 radial integrals, with the angular dependence on the Gaussian product
   * centres included, for powers N = 0, ..., L, by a fixed Gauss-Legendre quadrature per primitive
   * pair and local ECP primitive.
   * @param L - the maximum power N (LA + LB)
   * @param U - the ECP to be integrated over
   * @param shellA - the first GaussianShell
   * @param shellB - the second GaussianShell
   * @param data - the data container for the shell pair
   * @param parameters - pre-calculated parameters for the radial integral
   * @param radials - returned, indexed as (N, lam, lam + mu) for lam <= N with lam = N mod 2
   */
  void type1(int L, const ECP& U, const GaussianShell& shellA, const GaussianShell& shellB,
             const ShellPairData& data, const Parameters& parameters,
             ThreeIndex<double>& radials) const;

  /**
   * Calculates the requested type 2 radial integrals
   * @code{.unparsed}
   *   G(N, l1, l2) = sum c_a c_b d_u int r^(N+2+n_u) exp(-u r^2 - a (r - A)^2 - b (r - B)^2)
   *                  K_l1(2 a A r) K_l2(2 b B r) dr,
   * @endcode
   * over the primitives a of shellA, b of shellB and u of the ECP channel lam. Primitive triples
   * with both basis exponents <= 30 are integrated on a composite Gauss-Legendre grid shared by
   * all shell pairs (per thread, cached), the others each by the same composite rule on a window
   * bracketing the integrand. Triples provably below the tolerance are skipped. A or B may
   * be 0 (a shell on the ECP), in which case only l = 0 is non-zero for that shell.
   *
   * @param triples - vector of triples of form {N, l1, l2} of all required radial integrals
   * @param nbase - the number of base integrals for the closed forms (currently unused)
   * @param lam - the ECP shell angular momentum to be calculated over
   * @param U - the ECP to be integrated over
   * @param shellA - the first GaussianShell
   * @param shellB - the second GaussianShell
   * @param A - the magnitude of the distance of shellA from the ECP
   * @param B - the magnitude of the distance of shellB from the ECP
   * @param radials - the array to return the integrals in, indexed as (N, l1, l2)
   */
  void type2(const std::vector<Triple>& triples, int nbase, int lam, const ECP& U,
             const GaussianShell& shellA, const GaussianShell& shellB, double A, double B,
             ThreeIndex<double>& radials) const;
};

}  // namespace libecpint

#endif
