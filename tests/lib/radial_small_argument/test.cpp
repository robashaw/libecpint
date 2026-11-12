/*
 * Regression test for type-2 radial integrals with a diffuse primitive.
 *
 * The closed forms for G(k, l1, l2), l1 <= l2 <= 4, are polynomials in 1/x and 1/y with x = aA and
 * y = bB applied to base integrals. For a diffuse primitive (small x or y, e.g. exponent 0.005 at
 * 1.6 bohr gives x = 0.008) the terms cancel catastrophically: relative errors reached 2e4 at
 * x = 0.02 and 1e9 at x = 0.005 for l1 + l2 >= 4. Such integrals are now computed numerically.
 *
 * Checks every closed-form-range triple (l1 <= l2 <= 4, k = l1 + l2 (mod 2), 1 <= k <= 12 - l1 - l2)
 * with x, and separately y, between 0.0125 and 0.4 against independent quadrature
 * (radial_reference.hpp).
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

  const double reltol = 1e-8;
  const double abstol = 1e-15;  // the library's absolute integration threshold; some radials are ~1e-12
  const double u = 2.0, A = 1.6, bhealthy = 1.0, B = 2.0;  // healthy partner: y = bB = 2
  const double xs[] = {0.0125, 0.025, 0.05, 0.1, 0.2, 0.4};

  double worst = 0.0;
  int nfail = 0, ntot = 0;
  std::cout << std::scientific << std::setprecision(3);
  for (int l1 = 0; l1 <= 4; l1++) {
    for (int l2 = l1; l2 <= 4; l2++) {
      for (int k = 1; k <= 12 - l1 - l2; k++) {
        if ((k + l1 + l2) % 2) continue;
        const int nin = (k % 2) ? 1 : 2;  // ECP power r^(nin-2); radial power k = N + (nin - 2) + 2
        const int N = k - nin;
        for (double xv : xs) {
          for (int mirror = 0; mirror < 2; mirror++) {
            // mirror = 0: small x (shell a diffuse); mirror = 1: small y (shell b diffuse)
            const double a = mirror ? bhealthy * B / A : xv / A;
            const double b = mirror ? xv / B : bhealthy;
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
            if (std::abs(radials(N, l1, l2) - ref) > abstol) worst = std::max(worst, rel);
            ntot++;
            if (rel > reltol && std::abs(radials(N, l1, l2) - ref) > abstol) {
              nfail++;
              std::cout << "FAIL (l1, l2, k) = (" << l1 << ", " << l2 << ", " << k << ")  "
                        << (mirror ? "y = bB = " : "x = aA = ") << xv << "  value " << radials(N, l1, l2)
                        << " ref " << ref << " rel " << rel << std::endl;
            }
          }
        }
      }
    }
  }

  std::cout << ntot << " radial integrals with a diffuse primitive checked; worst relative error " << worst
            << " (tolerance " << reltol << ", or absolute " << abstol << ")" << std::endl;
  if (nfail) {
    std::cout << "FAIL: " << nfail << " values outside tolerance" << std::endl;
    return 1;
  }
  std::cout << "PASS" << std::endl;
  return 0;
}
