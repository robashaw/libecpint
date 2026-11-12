/*
 * Sweep test for type-2 radial integrals with shells A and B both off the ECP centre.
 *
 * These radials are integrated with a fixed Gauss-Legendre rule on a window bracketing the
 * integrand, and primitive triples whose radials are provably negligible are skipped. Both are
 * only correct if the window contains the integrand and the rule resolves it, for every primitive
 * triple and every triple set the generated code can pass. This checks type2() over a broad
 * sweep - ECP exponents and powers, shell exponents from very diffuse to core-tight, near and far,
 * equal and unequal distances - passing the whole set of (N, l1, l2) for a shell pair at once, as
 * the generated code does, since that gives the widest windows.
 *
 * Errors are measured relative to the largest radial over the set for the same primitives, the
 * scale on which they enter the integrals, against independent quadrature (radial_reference.hpp).
 * The closed forms previously used for l1 <= l2 <= 4 failed this sweep (up to 1e-3 of the scale,
 * 2% relative) for a core-tight primitive paired with a moderate one, as did a one-point
 * screening estimate (skipped radials up to 0.26 of the scale); the fixed rule is within ~4e-12.
 */
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

#include "../common/radial_reference.hpp"
#include "ecp.hpp"
#include "gshell.hpp"
#include "mathutil.hpp"
#include "multiarr.hpp"
#include "radial.hpp"

using namespace libecpint;

int main() {
  initFactorials();
  RadialIntegral rad;
  rad.init(20, 1e-15, 256, 1024);

  const double reltol = 1e-10;
  const double abstol = 2e-15;  // the library skips contributions below 1e-15
  // (LA, LB, lambda): d|d with an f projector, g|g with a d projector, and h|h (first derivatives
  // of g shells) with an f projector
  const int shells[][3] = {{2, 2, 3}, {4, 4, 2}, {5, 5, 3}};
  // ECP exponent and power n (r^(n-2) exp(-u r^2))
  const double us[] = {0.25, 3.3, 40.0};
  const int ns[] = {2, 1, 0};
  const double exps[] = {0.005, 0.07, 1.0, 20.0, 2000.0};
  const double geoms[][2] = {{0.8, 0.8}, {2.1, 1.9}, {1.5, 6.0}};  // (A, B) in bohr

  // Reference quadrature nodes: 16-point Gauss-Legendre on each of `panels` panels of [0, rmax]
  std::vector<double> x16, w16;
  radref::gauss_legendre16(x16, w16);
  const int panels = 1000;

  double worst = 0.0;
  int nfail = 0, ntot = 0;
  std::cout << std::scientific << std::setprecision(3);
  for (const auto& sh : shells) {
    const int LA = sh[0], LB = sh[1], lam = sh[2];
    std::vector<Triple> triples;
    for (int N = 0; N <= LA + LB; N++)
      for (int l1 = 0; l1 <= lam + LA; l1++)
        for (int l2 = (l1 + N) % 2; l2 <= lam + LB; l2 += 2) triples.push_back(Triple{N, l1, l2});
    // nbase as in the generator (from the lexicographically largest triple)
    int nbase = 0;
    for (const Triple& t : triples) nbase = std::max(nbase, std::get<0>(t) + std::get<1>(t) - 1);
    const int lmaxA = lam + LA, lmaxB = lam + LB, kmaxN = LA + LB;

    for (int iu = 0; iu < 3; iu++) {
      const double u = us[iu];
      const int nin = ns[iu];
      for (double a : exps) {
        for (double b : exps) {
          for (const auto& g : geoms) {
            const double A = g[0], B = g[1];
            double O[3] = {0.0, 0.0, 0.0};
            ECP U(O);
            U.addPrimitive(nin, lam, u, 1.0, true);
            std::array<double, 3> c0 = {0.0, 0.0, 0.0};
            GaussianShell sa(c0, LA), sb(c0, LB);
            sa.addPrim(a, 1.0);
            sb.addPrim(b, 1.0);
            ThreeIndex<double> radials(nbase + 2, lmaxA + 1, lmaxB + 1);
            radials.fill(0.0);
            rad.type2(triples, nbase, lam, U, sa, sb, A, B, radials);

            // Reference: tabulate the exponential and Bessel factors once on the reference nodes,
            // skipping nodes where the integrand underflows
            const int kmax = kmaxN + nin;
            double rmax = std::sqrt(80.0 / u) + 1.0;
            for (int it = 0; it < 20; it++)
              rmax = std::sqrt((80.0 + kmax * std::log(std::max(rmax, 1.0))) / u) + 1.0;
            const double h = rmax / panels;
            std::vector<double> wrk, Ka, Kb;  // wrk[q * (kmax + 1) + k] = weight * e * r^k
            for (int p = 0; p < panels; p++) {
              for (int i = 0; i < 16; i++) {
                const double r = (p + 0.5) * h + 0.5 * h * x16[i];
                const double e =
                    std::exp(-u * r * r - a * (r - A) * (r - A) - b * (r - B) * (r - B));
                if (e < 1e-250) continue;
                double wk = 0.5 * h * w16[i] * e;
                for (int k = 0; k <= kmax; k++, wk *= r) wrk.push_back(wk);
                for (int l = 0; l <= lmaxA; l++) Ka.push_back(radref::K(l, 2.0 * a * A * r));
                for (int l = 0; l <= lmaxB; l++) Kb.push_back(radref::K(l, 2.0 * b * B * r));
              }
            }
            std::vector<double> val, ref;
            for (const Triple& t : triples) {
              const int N = std::get<0>(t), l1 = std::get<1>(t), l2 = std::get<2>(t);
              const int k = N + nin;  // r^k with k = N + (nin - 2) + 2
              double sum = 0.0;
              const size_t nref = wrk.size() / (kmax + 1);
              for (size_t q = 0; q < nref; q++)
                sum +=
                    wrk[q * (kmax + 1) + k] * Ka[q * (lmaxA + 1) + l1] * Kb[q * (lmaxB + 1) + l2];
              val.push_back(radials(N, l1, l2));
              ref.push_back(sum);
            }

            double scale = 0.0;
            for (double r : ref) scale = std::max(scale, std::abs(r));
            for (size_t n = 0; n < val.size(); n++) {
              const double err = std::abs(val[n] - ref[n]);
              ntot++;
              if (scale > abstol) worst = std::max(worst, err / scale);
              if (err > reltol * scale + abstol) {
                nfail++;
                if (nfail <= 20)
                  std::cout << "FAIL (LA, LB, lam) = (" << LA << ", " << LB << ", " << lam
                            << ") u = " << u << " n = " << nin << " a = " << a << " b = " << b
                            << " A = " << A << " B = " << B << " (N, l1, l2) = ("
                            << std::get<0>(triples[n]) << ", " << std::get<1>(triples[n]) << ", "
                            << std::get<2>(triples[n]) << ")  value " << val[n] << " ref " << ref[n]
                            << " err/scale " << err / scale << std::endl;
              }
            }
          }
        }
      }
    }
  }

  std::cout << ntot << " off-centre type 2 radial integrals checked; worst error " << worst
            << " of scale (tolerance " << reltol << ", or absolute " << abstol << ")" << std::endl;
  if (nfail) {
    std::cout << "FAIL: " << nfail << " values outside tolerance" << std::endl;
    return 1;
  }
  std::cout << "PASS" << std::endl;
  return 0;
}
