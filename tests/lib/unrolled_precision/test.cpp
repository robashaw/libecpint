/*
 * Regression test for the precision of the unrolled (generated) type-2 integral kernels.
 *
 * For LA, LB <= LIBECPINT_MAX_UNROL the generator writes the angular prefactors into the kernel
 * source as decimal literals. These were written with the default stream precision of 6
 * significant digits (e.g. 16 pi^2 -> 157.914), so every off-centre s/p type-2 integral was wrong
 * at the 1e-6 - 1e-5 relative level.
 *
 * Checks, for shells A, B both off the ECP centre and a single s projector
 * U = d exp(-u r^2) |s><s| (n = 2), against independent references (radial_reference.hpp):
 *   <s|U|s>   = 4 pi d G(2, 0, 0)
 *   <p_x|U|s> = (1 / 2a) d/dA_x <s|U|s>      (unnormalised primitives, Gaussian derivative identity)
 */
#include "ecp.hpp"
#include "gshell.hpp"
#include "ecpint.hpp"
#include "multiarr.hpp"
#include "../common/radial_reference.hpp"
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>

using namespace libecpint;

int main() {
  double O[3] = {0.0, 0.0, 0.0};
  const double reltol = 1e-10;
  // u, d, a, b, A(xyz), B(xyz)
  const std::vector<std::array<double, 10>> cases = {
      {1.2, 10.0, 0.8, 0.6, 0.31, 0.52, 2.10, -0.70, 0.20, 1.60},
      {3.0, -4.0, 0.3, 1.1, 0.70, 0.90, 0.40, 1.20, -1.10, 1.90},
      {0.5, 2.5, 1.5, 0.4, -1.30, 0.20, 2.60, 0.40, 0.60, -0.50},
  };

  double worst = 0.0;
  std::cout << std::scientific << std::setprecision(3);
  for (const auto& c : cases) {
    const double u = c[0], d = c[1], a = c[2], b = c[3];
    const std::array<double, 3> RA = {c[4], c[5], c[6]}, RB = {c[7], c[8], c[9]};

    ECP U(O);
    U.addPrimitive(2, 1, 1.0, 0.0, false);  // zero local (p-ul) channel so that l = 0 is a projector
    U.addPrimitive(2, 0, u, d, true);
    ECPIntegral eng(1, 1, 0, 1e-15);

    auto vss_ref = [&](std::array<double, 3> A) {
      const double An = std::sqrt(A[0] * A[0] + A[1] * A[1] + A[2] * A[2]);
      const double Bn = std::sqrt(RB[0] * RB[0] + RB[1] * RB[1] + RB[2] * RB[2]);
      return 4.0 * M_PI * d * radref::G(2, 0, 0, u, a, b, An, Bn);
    };

    GaussianShell sA(RA, 0), pA(RA, 1), sB(RB, 0);
    sA.addPrim(a, 1.0);
    pA.addPrim(a, 1.0);
    sB.addPrim(b, 1.0);

    TwoIndex<double> vss, vps;
    eng.compute_shell_pair(U, sA, sB, vss);
    eng.compute_shell_pair(U, pA, sB, vps);

    const double ss_ref = vss_ref(RA);
    const double h = 1e-3;
    double e[4];
    const int st[4] = {-2, -1, 1, 2};
    for (int i = 0; i < 4; i++) {
      auto A = RA;
      A[0] += st[i] * h;
      e[i] = vss_ref(A);
    }
    const double ps_ref = (e[0] - 8 * e[1] + 8 * e[2] - e[3]) / (12 * h) / (2.0 * a);

    const double rel_ss = std::abs(vss(0, 0) - ss_ref) / std::abs(ss_ref);
    const double rel_ps = std::abs(vps(0, 0) - ps_ref) / std::abs(ps_ref);
    std::cout << "<s|U|s> " << vss(0, 0) << " ref " << ss_ref << " rel " << rel_ss << "   <px|U|s> "
              << vps(0, 0) << " ref " << ps_ref << " rel " << rel_ps << std::endl;
    worst = std::max({worst, rel_ss, rel_ps});
  }

  std::cout << "worst relative error " << worst << " (tolerance " << reltol << ")" << std::endl;
  if (worst > reltol) {
    std::cout << "FAIL" << std::endl;
    return 1;
  }
  std::cout << "PASS" << std::endl;
  return 0;
}
