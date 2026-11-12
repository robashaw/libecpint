/*
 * Regression test for spurious early convergence of the adaptive Gauss-Chebyshev quadrature.
 *
 * The one-point (Perez92) scheme stopped on the first level where its convergence test held. On a
 * coarse level that can happen by coincidence, and the result then jumps between very different
 * values under tiny changes of the integrand. In La2O3/def2-SVP (large-core La ECPs, compressed
 * geometry) single ECP matrix elements jumped by up to 7% (type 1 integral of an O d shell with the
 * O 1s contraction) when the geometry moved by 1e-14 bohr, putting ~4e-7 Eh steps into the energy
 * and making finite-difference gradients wrong by ~5e-4. The two-point (Perez93) scheme, then used
 * for windowed type 2 radials, did the same at the ~1e-8 level.
 *
 * Checks two such shell pairs (with the La ECP) at the geometries where the jumps occurred against
 * the polynomial interpolating the same integrals at nearby geometries, i.e. that the integrals
 * are smooth there; no reference values are needed.
 */
#include "ecp.hpp"
#include "ecpint.hpp"
#include "gshell.hpp"
#include "multiarr.hpp"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace libecpint;

namespace {
// La2O3 (bohr) and the displacement direction along which the jumps were found
const double X0[5][3] = {{0.0, 0.944863, -0.3286479},
                         {0.0, -0.944863, -0.3286479},
                         {0.0, 1.88973, 1.5610781},
                         {0.0, 0.0, 1.5610781},
                         {0.0, -1.88973, 1.5610781}};
const double DZ[5] = {0.0, 0.8660254037844386, -0.28867513459481287, -0.28867513459481287,
                      -0.28867513459481287};

// La def2 ECP (46 core electrons): n, l, exponent, coefficient
const double LA_ECP[10][4] = {{2, 3, 4.0286, -36.010016},  {2, 0, 3.3099, 91.932177},
                              {2, 0, 1.655, -3.788764},    {2, 0, 4.0286, 36.010016},
                              {2, 1, 2.8368, 63.759486},   {2, 1, 1.4184, -0.647958},
                              {2, 1, 4.0286, 36.010016},   {2, 2, 2.0213, 36.116173},
                              {2, 2, 1.0107, 0.219114},    {2, 2, 4.0286, 36.010016}};

struct Shell {
  int centre, l;
  std::vector<double> exps, coefs;
};
struct Case {
  const char* what;
  Shell a, b;
  double t;  // displacement at which the jump occurred
};

// def2-SVP shells, coefficients as passed by Psi4 (normalized)
const Case CASES[] = {
    {"type 1: O d x O 1s",
     {3, 2, {1.2}, {2.2645224825069943}},
     {2,
      0,
      {2266.1767785, 340.87010191, 77.363135167, 21.47964494, 6.6589433124},
      {-1.2615859724901695, -2.2748500427418206, -3.3479623286497575, -3.3299980530970594,
       -1.3203828288991972}},
     -0.004997862070479},
    {"type 2: La f x O p",
     {1, 3, {0.25683}, {0.069128541039011268}},
     {3, 1, {17.721504317, 3.863550544, 1.0480920883}, {3.2460107641875786, 2.5734959172661265, 1.1208372892488978}},
     0.011181024626829},
};

std::vector<double> integrals(const ECPIntegral& eng, const Case& c, const double t) {
  double X[5][3];
  for (int i = 0; i < 5; i++) {
    X[i][0] = X0[i][0];
    X[i][1] = X0[i][1];
    X[i][2] = X0[i][2] + t * DZ[i];
  }
  ECP U(X[0]);  // ECP on the first La
  for (int k = 0; k < 10; k++)
    U.addPrimitive((int)LA_ECP[k][0], (int)LA_ECP[k][1], LA_ECP[k][2], LA_ECP[k][3], k == 9);
  GaussianShell sa(X[c.a.centre], c.a.l), sb(X[c.b.centre], c.b.l);
  for (size_t i = 0; i < c.a.exps.size(); i++) sa.addPrim(c.a.exps[i], c.a.coefs[i]);
  for (size_t i = 0; i < c.b.exps.size(); i++) sb.addPrim(c.b.exps[i], c.b.coefs[i]);
  TwoIndex<double> v;
  eng.compute_shell_pair(U, sa, sb, v);
  return v.data;
}
}  // namespace

int main() {
  ECPIntegral eng(3, 3, 0);
  const double h = 1e-3, tol = 1e-9;
  const int ks[6] = {-3, -2, -1, 1, 2, 3};
  int nfail = 0;
  std::cout << std::scientific << std::setprecision(3);
  for (const Case& c : CASES) {
    const std::vector<double> v0 = integrals(eng, c, c.t);
    std::vector<std::vector<double>> vk;
    for (int k : ks) vk.push_back(integrals(eng, c, c.t + k * h));
    double scale = 0.0, worst = 0.0;
    for (double x : v0) scale = std::max(scale, std::abs(x));
    for (size_t e = 0; e < v0.size(); e++) {
      // Lagrange interpolation through the six neighbouring geometries, evaluated at c.t
      double interp = 0.0;
      for (int i = 0; i < 6; i++) {
        double w = 1.0;
        for (int j = 0; j < 6; j++)
          if (j != i) w *= (0.0 - ks[j]) / double(ks[i] - ks[j]);
        interp += w * vk[i][e];
      }
      worst = std::max(worst, std::abs(v0[e] - interp) / scale);
    }
    std::cout << c.what << ": max |integral - interpolation| = " << worst
              << " of the largest element (tolerance " << tol << ")" << std::endl;
    if (worst > tol) nfail++;
  }
  if (nfail) {
    std::cout << "FAIL" << std::endl;
    return 1;
  }
  std::cout << "PASS" << std::endl;
  return 0;
}
