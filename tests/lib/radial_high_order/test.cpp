/*
 * Regression test for type-2 radial integrals beyond the closed forms, and the ECP integral
 * derivatives that depend on them.
 *
 * For shells A, B both off the ECP centre, radial integrals G(k, l1, l2) without a closed form
 * (e.g. l1 or l2 >= 5, which occur once lambda + l >= 5, and for derivatives of f/g shells) fall
 * back to a numerical integration. That integrated on a fixed window around the Gaussian-product
 * centre (aA + bB)/(u + a + b), which misses most of the integrand once high-order Bessel factors
 * shift its weight, yet reported convergence. Integrals were wrong by percent to several hundred
 * percent, e.g. first derivatives of f/g shells with a d projector (def2 ECPs with def2-TZVP and
 * larger bases).
 *
 * Checks
 *  (a) G(k, l1, l2) beyond the closed forms against independent quadrature (radial_reference.hpp);
 *  (b) d/dR sum_ij <a_i|U|b_j> for an (l, R | l, R) shell pair against 5-point finite differences
 *      of the integrals, for l = 2..4 with s, p, d projectors.
 */
#include "ecp.hpp"
#include "gshell.hpp"
#include "ecpint.hpp"
#include "mathutil.hpp"
#include "multiarr.hpp"
#include "radial.hpp"
#include "../common/radial_reference.hpp"
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>

using namespace libecpint;

int main() {
  initFactorials();
  std::cout << std::scientific << std::setprecision(3);
  int nfail = 0;

  // (a) radial integrals beyond the closed forms
  {
    RadialIntegral rad;
    rad.init(20, 1e-15, 256, 1024);
    const double reltol = 1e-8;
    // u, a, b, A, B
    const std::vector<std::array<double, 5>> sets = {
        {1.2, 0.8, 0.6, 2.1, 1.7}, {3.0, 0.3, 1.1, 1.0, 2.5}, {0.5, 1.5, 0.4, 3.0, 0.8}, {1.2, 0.8, 0.6, 2.185, 1.758}};
    // (l1, l2, k): no closed form for max(l1, l2) > 4, for l1 > l2, or for k > 12 - l1 - l2
    const std::vector<std::array<int, 3>> triples = {{5, 0, 5}, {0, 5, 7}, {5, 5, 2}, {5, 4, 9}, {4, 5, 5},
                                                     {6, 3, 5}, {2, 6, 8}, {6, 6, 4}, {7, 5, 6}, {1, 0, 3},
                                                     {1, 1, 12}, {0, 0, 14}, {2, 2, 10}};
    double worst = 0.0;
    for (const auto& t : triples) {
      const int l1 = t[0], l2 = t[1], k = t[2];
      for (const auto& s : sets) {
        const double u = s[0], a = s[1], b = s[2], A = s[3], B = s[4];
        const int nin = (k % 2) ? 1 : 2;
        const int N = k - nin;
        double O[3] = {0.0, 0.0, 0.0};
        ECP U(O);
        U.addPrimitive(nin, 0, u, 1.0, true);
        std::array<double, 3> c0 = {0.0, 0.0, 0.0};
        GaussianShell sa(c0, 0), sb(c0, 0);
        sa.addPrim(a, 1.0);
        sb.addPrim(b, 1.0);
        const int nbase = 20;
        ThreeIndex<double> radials(nbase, 10, 10);
        std::vector<Triple> tl = {Triple{N, l1, l2}};
        rad.type2(tl, nbase, 0, U, sa, sb, A, B, radials);
        const double ref = radref::G(k, l1, l2, u, a, b, A, B);
        const double rel = std::abs(radials(N, l1, l2) - ref) / std::abs(ref);
        worst = std::max(worst, rel);
        if (rel > reltol) {
          nfail++;
          std::cout << "FAIL radial (l1, l2, k) = (" << l1 << ", " << l2 << ", " << k << ")  value "
                    << radials(N, l1, l2) << " ref " << ref << " rel " << rel << std::endl;
        }
      }
    }
    std::cout << "(a) high-order radial integrals: worst relative error " << worst << " (tolerance " << reltol
              << ")" << std::endl;
  }

  // (b) first derivatives of off-centre f/g shell pairs vs finite differences
  {
    const double reltol = 1e-6;  // limited by the finite-difference reference
    double worst = 0.0;
    double O[3] = {0.0, 0.0, 0.0};
    const std::array<double, 3> R = {0.31, 0.52, 4.10};
    auto sumV = [&](const ECPIntegral& eng, const ECP& U, int l, std::array<double, 3> C) {
      GaussianShell s(C, l);
      s.addPrim(0.8, 1.0);
      TwoIndex<double> v;
      eng.compute_shell_pair(U, s, s, v);
      double sum = 0.0;
      for (double x : v.data) sum += x;
      return sum;
    };
    for (int l = 2; l <= 4; l++) {
      for (int lam = 0; lam <= 2; lam++) {
        ECP U(O);
        U.addPrimitive(2, lam + 1, 1.0, 0.0, false);  // zero local channel above the projector
        U.addPrimitive(2, lam, 1.2, 10.0, true);
        ECPIntegral eng(l, lam + 1, 1, 1e-15);
        GaussianShell s(R, l);
        s.addPrim(0.8, 1.0);
        std::array<TwoIndex<double>, 9> d;
        eng.compute_shell_pair_derivative(U, s, s, d);
        const double h = 1e-3;
        for (int k = 0; k < 3; k++) {
          double anal = 0.0;
          for (double x : d[k].data) anal += x;
          for (double x : d[3 + k].data) anal += x;
          double e[4];
          const int st[4] = {-2, -1, 1, 2};
          for (int i = 0; i < 4; i++) {
            auto C = R;
            C[k] += st[i] * h;
            e[i] = sumV(eng, U, l, C);
          }
          const double fd = (e[0] - 8 * e[1] + 8 * e[2] - e[3]) / (12 * h);
          const double rel = std::abs(anal - fd) / std::abs(fd);
          worst = std::max(worst, rel);
          if (rel > reltol) {
            nfail++;
            std::cout << "FAIL derivative l = " << l << " lambda = " << lam << " xyz = " << k << "  analytic " << anal
                      << " fd " << fd << " rel " << rel << std::endl;
          }
        }
      }
    }
    std::cout << "(b) derivatives vs finite differences: worst relative error " << worst << " (tolerance " << reltol
              << ")" << std::endl;
  }

  if (nfail) {
    std::cout << "FAIL: " << nfail << " checks outside tolerance" << std::endl;
    return 1;
  }
  std::cout << "PASS" << std::endl;
  return 0;
}
