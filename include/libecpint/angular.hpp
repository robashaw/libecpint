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

#ifndef ANGULAR_HEAD
#define ANGULAR_HEAD

#include <vector>

#include "multiarr.hpp"

namespace libecpint {

/**
 * \class AngularIntegral
 * \brief Calculates and stores the angular integrals needed for ECP integration.
 *
 * This should not usually be created directly, it is instead owned by an ECPIntegral object,
 * so that integrals can be performed over multiple ECP centers without duplicating work.
 */
class AngularIntegral {
 private:
  int LB;  ///< Maximum angular momentum of orbital basis
  int LE;  ///< Maximum angular momentum of ECP basis

  int maxL;  ///< Largest power sum k + l + m and angular momentum in W: max(2 LB, LB + LE)
  int nrs;   ///< Number of (rho, sigma) pairs in each row of omega: (LB + LE + 1)^2

  /// Type 1 angular integrals W(k, l, m, lam, mu), for k + l + m <= maxL and lam <= maxL, at
  /// simplex(k, l, m) * (maxL + 1)^2 + lam^2 + lam + mu. Only the entries read are stored (a dense
  /// array over every index took gigabytes at high angular momentum).
  std::vector<double> W;
  /// Type 2 angular integrals Omega(k, l, m, lam, mu, rho, sigma), for k + l + m <= LB, lam <= LE
  /// (the ECP channel) and rho <= LB + LE (the basis side), at
  /// (simplex(k, l, m) * (LE + 1)^2 + lam^2 + lam + mu) * nrs + rho^2 + rho + sigma
  std::vector<double> omega;

  /// Index of (k, l, m) among all power triples, ordered by k + l + m, so that the triples with
  /// k + l + m <= n come first for every n
  static int simplex(const int k, const int l, const int m) {
    const int n = k + l + m, j = l + m;
    return n * (n + 1) * (n + 2) / 6 + j * (j + 1) / 2 + m;
  }

  /// Worker functions for calculating terms in the USP to spherical transformation coefficients
  double calcG(int l, int m) const;
  double calcH1(int i, int j, int l, int m) const;
  double calcH2(int i, int j, int k, int m) const;

  /**
   * Calculates the polynomial integrals, int (x^i y^j z^k dOm) where dOm is the solid angle
   * differential
   * @param maxI - the maximum power, i, to determine
   * @return ThreeIndex of polynomials P(i, j, k), where strictly i >= j >= k
   */
  ThreeIndex<double> Pijk(int maxI) const;

  /**
   * Builds the type 1 angular integrals
   * @param U - the USP to spherical transformation coefficients
   */
  void makeW(const FiveIndex<double>& U);
  /**
   * Builds the type 2 angular integrals
   * @param U - the USP to spherical transformation coefficients
   */
  void makeOmega(const FiveIndex<double>& U);

 public:
  /**
   * Calculates all possible USP to spherical transformation coefficients for a given angular
   * momentum
   * @param lam - the angular momentum
   * @param mu - the subshell
   * @return ThreeIndex of the values U_lam,mu(k, l, m)
   */
  ThreeIndex<double> uklm(int lam, int mu) const;

  /**
   * Builds the USP to spherical transformation coefficients for use in calculating the type 1 and 2
   * integrals
   * @return FiveIndex of the coefficients U(lam, lam+mu, k, l, m)
   */
  FiveIndex<double> makeU() const;

  /// Default constructor creates empty object
  AngularIntegral();

  /**
   *  Specified constructor calls init with given arguments
   * @param LB - the maximum angular momentum of the orbital basis
   * @param LE - the maximum angular momentum of the ECP basis
   */
  AngularIntegral(int LB, int LE);

  /**
   * Initialises the object, must be called before anything else if default constructor was used.
   * @param LB - the maximum angular momentum of the orbital basis
   * @param LE - the maximum angular momentum of the ECP basis
   */
  void init(int LB, int LE);

  /**
   * Computes the type 1 and 2 angular integrals
   */
  void compute();

  /// TODO: Clears the W and omega arrays
  void clear();

  /**
   * Returns the type 1 angular integral W(k, l, m, lam, mu)
   * @param k - x index
   * @param l - y index
   * @param m - z index
   * @param lam - angular momentum
   * @param mu - subshell
   * @return value of type 1 angular integral
   */
  double getIntegral(const int k, const int l, const int m, const int lam, const int mu) const {
    return W[simplex(k, l, m) * (maxL + 1) * (maxL + 1) + lam * lam + lam + mu];
  }

  /**
   * Returns the type 2 angular integral Omega(k, l, m, lam, mu, rho, sigma)
   * @param k - x index
   * @param l - y index
   * @param m - z index
   * @param lam - angular momentum of current ECP shell (<= LE)
   * @param mu - subshell of lam
   * @param rho - angular momentum of current basis shell
   * @param sigma - subshell of rho
   * @return value of type 2 angular integral
   */
  double getIntegral(const int k, const int l, const int m, const int lam, const int mu,
                     const int rho, const int sigma) const {
    return getOmegaRow(k, l, m, lam, mu)[rho * rho + rho + sigma];
  }

  /// The type 2 angular integrals Omega(k, l, m, lam, mu, rho, sigma) for every (rho, sigma), at
  /// rho^2 + rho + sigma
  const double* getOmegaRow(const int k, const int l, const int m, const int lam,
                            const int mu) const {
    return &omega[(simplex(k, l, m) * (LE + 1) * (LE + 1) + lam * lam + lam + mu) * nrs];
  }

  /// is W(k, l, m, lam, mu) zero to within a given tolerance?
  bool isZero(int k, int l, int m, int lam, int mu, double tolerance) const;

  /// is Omega(k, l, m, lam, mu, rho, sigma) zero to within a given tolerance?
  bool isZero(int k, int l, int m, int lam, int mu, int rho, int sigma, double tolerance) const;
};

}  // namespace libecpint
#endif
