/*
 * Regression test for type-2 radial integrals whose integrand is narrow and far from the ECP.
 *
 * With a tight primitive off the ECP centre (e.g. an O 1s exponent 2266 at 2.1 bohr from La), the
 * integrand of G(k, l1, l2) is a peak of width ~1/sqrt(a) = 0.02 bohr around r = A. On the
 * adaptive grid over [0, inf) used for numerical radials, the peak fell between abscissae, every
 * sample was ~0 and the quadrature reported convergence to zero. Paired with a diffuse primitive
 * (which routes the radial to quadrature), whole ECP matrix elements came out wrong by ~1e-4 and
 * gradients by up to 1e-1 for La2O3/def2-SVP. Numerical radials are now integrated over a
 * finite window bracketing the integrand.
 *
 * Errors are measured relative to the largest radial over (l1, l2, k) for the same primitives,
 * the scale on which they enter the integrals, and checked against independent quadrature
 * (radial_reference.hpp). The adaptive quadrature stops on an absolute criterion (threshold
 * 1e-15), which leaves errors up to ~1e-14 in radials that are themselves ~1e-13; relative to the
 * scale these stay below 1e-5. Before the fix, errors were up to 1.0 (stock: 0.3) of the scale.
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
#include <vector>

using namespace libecpint;

int main() {
  initFactorials();
  RadialIntegral rad;
  rad.init(20, 1e-15, 256, 1024);

  const double tol = 1e-5;
  const double u = 3.3099, A = 2.11;  // La def2 ECP s exponent; La-O distance in La2O3
  const double as[] = {2266.18, 340.87, 21.4796};  // O def2-SVP 1s primitives
  const double bs[] = {0.0251, 0.276415};          // diffuse partners
  const double Bs[] = {2.11, 1.89};                // same centre, other centre

  double worst = 0.0;
  int nfail = 0, ntot = 0;
  std::cout << std::scientific << std::setprecision(3);
  for (double a : as) {
    for (double b : bs) {
      for (double B : Bs) {
        std::vector<std::array<int, 3>> idx;
        std::vector<double> val, ref;
        for (int l1 = 0; l1 <= 3; l1++) {
          for (int l2 = 0; l2 <= 4; l2++) {
            for (int k = 2; k <= 9; k++) {
              double O[3] = {0.0, 0.0, 0.0};
              ECP U(O);
              U.addPrimitive(2, 0, u, 1.0, true);  // r^0 exp(-u r^2): radial power k = N + 2
              std::array<double, 3> c0 = {0.0, 0.0, 0.0};
              GaussianShell sa(c0, 0), sb(c0, 0);
              sa.addPrim(a, 1.0);
              sb.addPrim(b, 1.0);
              const int N = k - 2, nbase = 16;
              ThreeIndex<double> radials(nbase, 8, 8);
              std::vector<Triple> t = {Triple{N, l1, l2}};
              rad.type2(t, nbase, 0, U, sa, sb, A, B, radials);
              idx.push_back({l1, l2, k});
              val.push_back(radials(N, l1, l2));
              ref.push_back(radref::G(k, l1, l2, u, a, b, A, B));
            }
          }
        }
        double scale = 0.0;
        for (double r : ref) scale = std::max(scale, std::abs(r));
        for (size_t n = 0; n < val.size(); n++) {
          const double err = std::abs(val[n] - ref[n]) / scale;
          worst = std::max(worst, err);
          ntot++;
          if (err > tol) {
            nfail++;
            std::cout << "FAIL a = " << a << " b = " << b << " B = " << B << " (l1, l2, k) = (" << idx[n][0]
                      << ", " << idx[n][1] << ", " << idx[n][2] << ")  value " << val[n] << " ref " << ref[n]
                      << " err/scale " << err << std::endl;
          }
        }
      }
    }
  }

  std::cout << ntot << " radial integrals with a tight off-centre primitive checked; worst error "
            << worst << " of scale (tolerance " << tol << ")" << std::endl;
  if (nfail) {
    std::cout << "FAIL: " << nfail << " values outside tolerance" << std::endl;
    return 1;
  }
  std::cout << "PASS" << std::endl;
  return 0;
}
