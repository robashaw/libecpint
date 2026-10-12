/*
 * Regression test for the closed-form type-2 radial integrals in radial_gen.cpp.
 *
 * RadialIntegral::type2(triples, ...) evaluates G(k, l1, l2) (see radial_reference.hpp) for
 * l1 <= l2 <= 4 by generated closed forms in terms of base integrals. One of them,
 * (l1, l2, k) = (1, 1, 10), had wrong coefficients, giving errors up to ~40x; it is reached by
 * high-L shell pairs (e.g. l = 5) off the ECP centre.
 *
 * Checks every closed-form triple, l1 <= l2 <= 4, k = l1 + l2 (mod 2), 1 <= k <= 12 - l1 - l2,
 * for several parameter sets against independent quadrature.
 */
#include "ecp.hpp"
#include "gshell.hpp"
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
  RadialIntegral rad;
  rad.init(20, 1e-15, 256, 1024);

  // Wrong formulas show up as O(1) errors. (Small a*A or b*B, where the closed forms would cancel,
  // is integrated numerically instead; see RadialSmallArgument.)
  const double reltol = 1e-9;
  // u, a, b, A, B
  const std::vector<std::array<double, 5>> sets = {
      {1.2, 0.8, 0.6, 2.1, 1.7}, {3.0, 0.3, 1.1, 1.0, 2.5}, {0.5, 1.5, 0.4, 3.0, 0.8},
      {2.0, 2.0, 2.0, 1.4, 1.4}, {1.2, 0.8, 0.6, 2.185, 1.758}};

  double worst = 0.0;
  int nfail = 0, ntot = 0;
  std::cout << std::scientific << std::setprecision(3);
  for (int l1 = 0; l1 <= 4; l1++) {
    for (int l2 = l1; l2 <= 4; l2++) {
      for (int k = 1; k <= 12 - l1 - l2; k++) {
        if ((k + l1 + l2) % 2) continue;
        for (const auto& s : sets) {
          const double u = s[0], a = s[1], b = s[2], A = s[3], B = s[4];
          const int nin = (k % 2) ? 1 : 2;  // ECP power r^(nin-2); radial power k = N + (nin - 2) + 2
          const int N = k - nin;
          double O[3] = {0.0, 0.0, 0.0};
          ECP U(O);
          U.addPrimitive(nin, 0, u, 1.0, true);
          std::array<double, 3> c0 = {0.0, 0.0, 0.0};
          GaussianShell sa(c0, 0), sb(c0, 0);
          sa.addPrim(a, 1.0);
          sb.addPrim(b, 1.0);

          const int nbase = 16;
          ThreeIndex<double> radials(nbase, 8, 8);
          std::vector<Triple> t = {Triple{N, l1, l2}};
          rad.type2(t, nbase, 0, U, sa, sb, A, B, radials);

          const double ref = radref::G(k, l1, l2, u, a, b, A, B);
          const double rel = std::abs(radials(N, l1, l2) - ref) / std::abs(ref);
          worst = std::max(worst, rel);
          ntot++;
          if (rel > reltol) {
            nfail++;
            std::cout << "FAIL (l1, l2, k) = (" << l1 << ", " << l2 << ", " << k << ")  u a b A B = " << u << " "
                      << a << " " << b << " " << A << " " << B << "  value " << radials(N, l1, l2) << " ref "
                      << ref << " rel " << rel << std::endl;
          }
        }
      }
    }
  }

  std::cout << ntot << " closed-form radial integrals checked; worst relative error " << worst << " (tolerance "
            << reltol << ")" << std::endl;
  if (nfail) {
    std::cout << "FAIL: " << nfail << " values outside tolerance" << std::endl;
    return 1;
  }
  std::cout << "PASS" << std::endl;
  return 0;
}
