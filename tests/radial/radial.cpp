/*
 * Accuracy of the type 2 radial integrals (RadialIntegral::type2) against independent quadrature
 * (tests/support/radial_reference.hpp), over parameter sets that each reproduce a past failure:
 *
 *   G(k, l1, l2) = int_0^inf r^k exp(-u r^2 - a (r - A)^2 - b (r - B)^2) K_l1(2aAr) K_l2(2bBr) dr
 *
 * for one ECP primitive (exponent u) and one primitive on each of two shells at distances A, B
 * from the ECP. Each test prints its worst error, so that a drift in accuracy shows before it
 * fails.
 */
#include "radial.hpp"

#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <cstdio>
#include <vector>

#include "../support/radial_reference.hpp"
#include "ecp.hpp"
#include "ecpint.hpp"
#include "gshell.hpp"
#include "mathutil.hpp"
#include "multiarr.hpp"

using namespace libecpint;

namespace {

const RadialIntegral& radial_integrator() {
  static const RadialIntegral* rad = [] {
    auto* r = new RadialIntegral;
    r->init(20, 1e-15, 256, 1024);
    return r;
  }();
  return *rad;
}

/// The library's G(k, l1, l2) for one primitive triple, with the ECP power r^(n - 2), n = 1 or 2,
/// chosen to match the parity of k (k = N + n for the library's N)
double library_G(const int k, const int l1, const int l2, const double u, const double a,
                 const double b, const double A, const double B) {
  const int n = (k % 2) ? 1 : 2, N = k - n;
  double O[3] = {0.0, 0.0, 0.0};
  ECP U(O);
  U.addPrimitive(n, 0, u, 1.0, true);
  const std::array<double, 3> c0 = {0.0, 0.0, 0.0};
  GaussianShell sa(c0, 0), sb(c0, 0);
  sa.addPrim(a, 1.0);
  sb.addPrim(b, 1.0);
  ThreeIndex<double> radials(N + 1, l1 + 1, l2 + 1);
  radials.fill(0.0);
  radial_integrator().type2({Triple{N, l1, l2}}, 0, 0, U, sa, sb, A, B, radials);
  return radials(N, l1, l2);
}

void info(const char* what, const double worst) {
  std::printf("[ info ] %s: worst error %.2e\n", what, worst);
}

}  // namespace

// Every (l1, l2, k) with l1 <= l2 <= 4, k = l1 + l2 (mod 2), 1 <= k <= 12 - l1 - l2: the range
// once covered by generated closed forms, one of which, (1, 1, 10), had wrong coefficients
// (errors up to ~40x, reached by high-L shell pairs off the ECP centre). Measured: 3e-14.
TEST(Type2, LowOrders) {
  const double reltol = 1e-12;
  // u, a, b, A, B
  const std::vector<std::array<double, 5>> sets = {{1.2, 0.8, 0.6, 2.1, 1.7},
                                                   {3.0, 0.3, 1.1, 1.0, 2.5},
                                                   {0.5, 1.5, 0.4, 3.0, 0.8},
                                                   {2.0, 2.0, 2.0, 1.4, 1.4},
                                                   {1.2, 0.8, 0.6, 2.185, 1.758}};
  double worst = 0.0;
  for (int l1 = 0; l1 <= 4; l1++)
    for (int l2 = l1; l2 <= 4; l2++)
      for (int k = 1; k <= 12 - l1 - l2; k++) {
        if ((k + l1 + l2) % 2) continue;
        for (const auto& s : sets) {
          const double ref = radref::G(k, l1, l2, s[0], s[1], s[2], s[3], s[4]);
          const double rel =
              std::abs(library_G(k, l1, l2, s[0], s[1], s[2], s[3], s[4]) - ref) / std::abs(ref);
          worst = std::max(worst, rel);
          EXPECT_LE(rel, reltol) << "(l1, l2, k) = (" << l1 << ", " << l2 << ", " << k
                                 << "), u a b A B = " << s[0] << " " << s[1] << " " << s[2] << " "
                                 << s[3] << " " << s[4];
        }
      }
  info("low orders (relative)", worst);
}

// A diffuse primitive gives small Bessel arguments x = aA or y = bB (exponent 0.005 at 1.6 bohr
// gives x = 0.008). The closed forms once used here, polynomials in 1/x and 1/y, cancelled
// catastrophically: relative errors of 2e4 at x = 0.02 and 1e9 at x = 0.005 for l1 + l2 >= 4.
// Measured: 1e-20 absolute on radials of ~1e-8 with high-order Bessel functions of small argument
// (the absolute precision of the Bessel table), otherwise below 1e-14 relative.
TEST(Type2, DiffusePrimitive) {
  const double reltol = 1e-12, abstol = 1e-18;
  const double u = 2.0, A = 1.6, b_partner = 1.0, B = 2.0;  // the partner has y = bB = 2
  double worst = 0.0;
  for (int l1 = 0; l1 <= 4; l1++)
    for (int l2 = l1; l2 <= 4; l2++)
      for (int k = 1; k <= 12 - l1 - l2; k++) {
        if ((k + l1 + l2) % 2) continue;
        for (const double x : {0.0125, 0.025, 0.05, 0.1, 0.2, 0.4})
          for (const bool small_y : {false, true}) {
            const double a = small_y ? b_partner * B / A : x / A;
            const double b = small_y ? x / B : b_partner;
            const double ref = radref::G(k, l1, l2, u, a, b, A, B);
            const double err = std::abs(library_G(k, l1, l2, u, a, b, A, B) - ref);
            if (err > abstol) worst = std::max(worst, err / std::abs(ref));
            EXPECT_LE(err, reltol * std::abs(ref) + abstol)
                << "(l1, l2, k) = (" << l1 << ", " << l2 << ", " << k << "), "
                << (small_y ? "y = bB = " : "x = aA = ") << x;
          }
      }
  info("diffuse primitive (relative, above 1e-18 absolute)", worst);
}

// Orders beyond the old closed forms (max(l1, l2) > 4, l1 > l2, or k > 12 - l1 - l2), as for
// lambda + l >= 5 and derivatives of f/g shells. Their numerical integration once used a fixed
// window around the Gaussian-product centre, which misses most of the integrand once high-order
// Bessel factors shift its weight, while reporting convergence: errors of percent to several
// hundred percent. Measured: 2e-14.
TEST(Type2, HighOrders) {
  const double reltol = 1e-12;
  const std::vector<std::array<double, 5>> sets = {{1.2, 0.8, 0.6, 2.1, 1.7},
                                                   {3.0, 0.3, 1.1, 1.0, 2.5},
                                                   {0.5, 1.5, 0.4, 3.0, 0.8},
                                                   {1.2, 0.8, 0.6, 2.185, 1.758}};
  // (l1, l2, k)
  const std::vector<std::array<int, 3>> triples = {
      {5, 0, 5}, {0, 5, 7}, {5, 5, 2}, {5, 4, 9},  {4, 5, 5},  {6, 3, 5}, {2, 6, 8},
      {6, 6, 4}, {7, 5, 6}, {1, 0, 3}, {1, 1, 12}, {0, 0, 14}, {2, 2, 10}};
  double worst = 0.0;
  for (const auto& t : triples)
    for (const auto& s : sets) {
      const double ref = radref::G(t[2], t[0], t[1], s[0], s[1], s[2], s[3], s[4]);
      const double rel =
          std::abs(library_G(t[2], t[0], t[1], s[0], s[1], s[2], s[3], s[4]) - ref) / std::abs(ref);
      worst = std::max(worst, rel);
      EXPECT_LE(rel, reltol) << "(l1, l2, k) = (" << t[0] << ", " << t[1] << ", " << t[2] << ")";
    }
  info("high orders (relative)", worst);
}

// The integrals that the high orders above feed: d/dR sum_ij <a_i|U|b_j> for an (l, R | l, R)
// shell pair, l = 2..4, with s, p and d projectors, against 5-point finite differences. With the
// old fixed window, first derivatives of f/g shells with a d projector (def2 ECPs with def2-TZVP
// and larger) were wrong by percent or more. Measured: 7e-8, the finite differences' accuracy.
TEST(Type2, HighOrderDerivatives) {
  const double reltol = 1e-6;  // limited by the finite differences
  double worst = 0.0;
  double O[3] = {0.0, 0.0, 0.0};
  const std::array<double, 3> R = {0.31, 0.52, 4.10};
  auto sum_values = [&](const ECPIntegral& eng, const ECP& U, const int l,
                        const std::array<double, 3>& C) {
    GaussianShell s(C, l);
    s.addPrim(0.8, 1.0);
    TwoIndex<double> v;
    eng.compute_shell_pair(U, s, s, v);
    double sum = 0.0;
    for (const double x : v.data) sum += x;
    return sum;
  };
  for (int l = 2; l <= 4; l++)
    for (int lam = 0; lam <= 2; lam++) {
      ECP U(O);
      U.addPrimitive(2, lam + 1, 1.0, 0.0, false);  // zero local channel above the projector
      U.addPrimitive(2, lam, 1.2, 10.0, true);
      const ECPIntegral eng(l, lam + 1, 1, 1e-15);
      GaussianShell s(R, l);
      s.addPrim(0.8, 1.0);
      std::array<TwoIndex<double>, 9> d;
      eng.compute_shell_pair_derivative(U, s, s, d);
      const double h = 1e-3;
      for (int k = 0; k < 3; k++) {
        double analytic = 0.0;
        for (const double x : d[k].data) analytic += x;
        for (const double x : d[3 + k].data) analytic += x;
        double e[4];
        const int steps[4] = {-2, -1, 1, 2};
        for (int i = 0; i < 4; i++) {
          auto C = R;
          C[k] += steps[i] * h;
          e[i] = sum_values(eng, U, l, C);
        }
        const double fd = (e[0] - 8 * e[1] + 8 * e[2] - e[3]) / (12 * h);
        const double rel = std::abs(analytic - fd) / std::abs(fd);
        worst = std::max(worst, rel);
        EXPECT_LE(rel, reltol) << "l = " << l << ", lambda = " << lam << ", xyz = " << k;
      }
    }
  info("high-order derivatives vs finite differences (relative)", worst);
}

// A tight primitive off the ECP centre (an O 1s exponent 2266 at 2.1 bohr from La) makes the
// integrand a peak of width ~0.02 bohr around r = A. On the adaptive grid once used, the peak fell
// between abscissae and the quadrature converged to zero: whole matrix elements of La2O3/def2-SVP
// wrong by ~1e-4, gradients by up to 1e-1. Errors relative to the largest radial over (l1, l2, k)
// for the same primitives, the scale on which they enter the integrals. Measured: 1e-14.
TEST(Type2, NarrowIntegrand) {
  const double tol = 1e-12;
  const double u = 3.3099, A = 2.11;  // La def2 ECP s exponent; La-O distance in La2O3
  double worst = 0.0;
  for (const double a : {2266.18, 340.87, 21.4796})  // O def2-SVP 1s primitives
    for (const double b : {0.0251, 0.276415})        // diffuse partners
      for (const double B : {2.11, 1.89}) {          // same centre, other centre
        std::vector<std::array<int, 3>> idx;
        std::vector<double> val, ref;
        for (int l1 = 0; l1 <= 3; l1++)
          for (int l2 = 0; l2 <= 4; l2++)
            for (int k = 2; k <= 9; k++) {
              idx.push_back({l1, l2, k});
              val.push_back(library_G(k, l1, l2, u, a, b, A, B));
              ref.push_back(radref::G(k, l1, l2, u, a, b, A, B));
            }
        double scale = 0.0;
        for (const double r : ref) scale = std::max(scale, std::abs(r));
        for (size_t n = 0; n < val.size(); n++) {
          const double err = std::abs(val[n] - ref[n]) / scale;
          worst = std::max(worst, err);
          EXPECT_LE(err, tol) << "a = " << a << ", b = " << b << ", B = " << B
                              << ", (l1, l2, k) = (" << idx[n][0] << ", " << idx[n][1] << ", "
                              << idx[n][2] << ")";
        }
      }
  info("narrow integrand (of scale)", worst);
}

// A broad sweep of type2() for shells both off the ECP: ECP exponents and powers, shell exponents
// from very diffuse to core-tight, near and far, equal and unequal distances, passing the whole
// set of (N, l1, l2) for a shell pair at once as the generated code does (the widest windows).
// The windows must contain the integrand, the rule resolve it, and the skipped primitive triples
// be negligible. The closed forms once used failed this sweep (up to 1e-3 of the scale) for a
// core-tight primitive with a moderate one, as did a one-point screening estimate (skipped radials
// up to 0.26 of the scale). Errors relative to the largest radial of the set. Measured: 2e-14.
TEST(Type2, ShellPairSweep) {
  const double reltol = 1e-12;
  const double abstol = 2e-15;  // the library skips contributions below 1e-15
  // (LA, LB, lambda): d|d with an f projector, g|g with a d projector, and h|h (first derivatives
  // of g shells) with an f projector
  const int shells[][3] = {{2, 2, 3}, {4, 4, 2}, {5, 5, 3}};
  const double us[] = {0.25, 3.3, 40.0};
  const int ns[] = {2, 1, 0};  // ECP power r^(n - 2)
  const double exps[] = {0.005, 0.07, 1.0, 20.0, 2000.0};
  const double geoms[][2] = {{0.8, 0.8}, {2.1, 1.9}, {1.5, 6.0}};  // (A, B) in bohr

  // Reference: 16-point Gauss-Legendre on each of 1000 panels of [0, rmax], with the exponential
  // and Bessel factors tabulated once per primitive triple
  std::vector<double> x16, w16;
  radref::gauss_legendre16(x16, w16);
  const int panels = 1000;

  double worst = 0.0;
  for (const auto& sh : shells) {
    const int LA = sh[0], LB = sh[1], lam = sh[2];
    std::vector<Triple> triples;
    for (int N = 0; N <= LA + LB; N++)
      for (int l1 = 0; l1 <= lam + LA; l1++)
        for (int l2 = (l1 + N) % 2; l2 <= lam + LB; l2 += 2) triples.push_back(Triple{N, l1, l2});
    const int lmaxA = lam + LA, lmaxB = lam + LB, kmaxN = LA + LB;
    for (int iu = 0; iu < 3; iu++) {
      const double u = us[iu];
      const int n = ns[iu];
      for (const double a : exps)
        for (const double b : exps)
          for (const auto& g : geoms) {
            const double A = g[0], B = g[1];
            double O[3] = {0.0, 0.0, 0.0};
            ECP U(O);
            U.addPrimitive(n, lam, u, 1.0, true);
            const std::array<double, 3> c0 = {0.0, 0.0, 0.0};
            GaussianShell sa(c0, LA), sb(c0, LB);
            sa.addPrim(a, 1.0);
            sb.addPrim(b, 1.0);
            ThreeIndex<double> radials(kmaxN + 1, lmaxA + 1, lmaxB + 1);
            radials.fill(0.0);
            radial_integrator().type2(triples, 0, lam, U, sa, sb, A, B, radials);

            const int kmax = kmaxN + n;
            double rmax = std::sqrt(80.0 / u) + 1.0;
            for (int it = 0; it < 20; it++)
              rmax = std::sqrt((80.0 + kmax * std::log(std::max(rmax, 1.0))) / u) + 1.0;
            const double h = rmax / panels;
            std::vector<double> wrk, Ka, Kb;  // wrk[q * (kmax + 1) + k] = weight * e * r^k
            for (int p = 0; p < panels; p++)
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
            const size_t nref = wrk.size() / (kmax + 1);
            std::vector<double> val, ref;
            for (const Triple& t : triples) {
              const int N = std::get<0>(t), l1 = std::get<1>(t), l2 = std::get<2>(t);
              double sum = 0.0;
              for (size_t q = 0; q < nref; q++)
                sum += wrk[q * (kmax + 1) + N + n] * Ka[q * (lmaxA + 1) + l1] *
                       Kb[q * (lmaxB + 1) + l2];
              val.push_back(radials(N, l1, l2));
              ref.push_back(sum);
            }
            double scale = 0.0;
            for (const double r : ref) scale = std::max(scale, std::abs(r));
            for (size_t m = 0; m < val.size(); m++) {
              const double err = std::abs(val[m] - ref[m]);
              if (scale > abstol) worst = std::max(worst, err / scale);
              EXPECT_LE(err, reltol * scale + abstol)
                  << "(LA, LB, lam) = (" << LA << ", " << LB << ", " << lam << "), u = " << u
                  << ", n = " << n << ", a = " << a << ", b = " << b << ", A = " << A
                  << ", B = " << B << ", (N, l1, l2) = (" << std::get<0>(triples[m]) << ", "
                  << std::get<1>(triples[m]) << ", " << std::get<2>(triples[m]) << ")";
            }
          }
    }
  }
  info("shell pair sweep (of scale)", worst);
}

// Primitive triples with a basis exponent above SG_TIGHT, integrated each on its own window, with
// every radial checked relative to itself: the type 2 integrals combine them with coefficients up
// to ~(1 + d sqrt(2a))^l larger than the result (see ARCHITECTURE.md), so a small relative error
// in one radial matters even if it is small on the scale of the set. A single 48-point panel over
// each window gave up to 6e-13 (now panels as on the shared grid). Measured: 1.4e-14.
TEST(Type2, TightPrimitives) {
  const double reltol = 1e-13;
  // a, b, A, B, u: tight pairs on the same and different centres, near and far, and one shell on
  // the ECP
  const std::vector<std::array<double, 5>> sets = {
      {50, 50, 2, 2, 1},     {50, 1.4, 2, 2, 1},   {200, 200, 1.5, 1.5, 3},
      {50, 50, 0.3, 0.3, 1}, {40, 8, 2, 0.7, 0.5}, {1000, 1000, 1, 1, 2},
      {50, 50, 2, 2, 10},    {60, 60, 0, 2, 1},    {2266, 1.2, 2.113, 3.5, 1}};
  // Unscreened, so that no radial is skipped for being small in absolute terms
  RadialIntegral rad;
  rad.init(12, 0.0, 256, 1024);
  const int Nmax = 6;
  double worst = 0.0;
  for (const auto& s : sets) {
    const double a = s[0], b = s[1], A = s[2], B = s[3], u = s[4];
    const int l1max = A > 0 ? 6 : 0, l2max = B > 0 ? 6 : 0;
    double O[3] = {0.0, 0.0, 0.0};
    ECP U(O);
    U.addPrimitive(2, 0, u, 1.0, true);  // r^0 exp(-u r^2): radial power k = N + 2
    const std::array<double, 3> c0 = {0.0, 0.0, 0.0};
    GaussianShell sa(c0, 0), sb(c0, 0);
    sa.addPrim(a, 1.0);
    sb.addPrim(b, 1.0);
    std::vector<Triple> triples;
    for (int N = 0; N <= Nmax; N++)
      for (int l1 = 0; l1 <= l1max; l1++)
        for (int l2 = 0; l2 <= l2max; l2++) triples.push_back(Triple{N, l1, l2});
    ThreeIndex<double> radials(Nmax + 1, l1max + 1, l2max + 1);
    radials.fill(0.0);
    rad.type2(triples, 0, 0, U, sa, sb, A, B, radials);
    for (const Triple& t : triples) {
      const int N = std::get<0>(t), l1 = std::get<1>(t), l2 = std::get<2>(t);
      const double ref = radref::G(N + 2, l1, l2, u, a, b, A, B);
      const double rel = std::abs(radials(N, l1, l2) - ref) / std::abs(ref);
      worst = std::max(worst, rel);
      EXPECT_LE(rel, reltol) << "a = " << a << ", b = " << b << ", A = " << A << ", B = " << B
                             << ", u = " << u << ", (N, l1, l2) = (" << N << ", " << l1 << ", "
                             << l2 << ")";
    }
  }
  info("tight primitives (relative)", worst);
}

// Smoothness of whole shell-pair integrals under tiny geometry changes, at geometries where an
// adaptive Gauss-Chebyshev quadrature once converged early by coincidence: in La2O3/def2-SVP
// (large-core La ECP, compressed geometry) single matrix elements jumped by up to 7% (type 1, an
// O d shell with the O 1s contraction) when the geometry moved by 1e-14 bohr, making
// finite-difference gradients wrong by ~5e-4; windowed type 2 radials did the same at ~1e-8.
// Checks the integrals at those geometries against the polynomial interpolating them at six
// nearby geometries (no reference values needed), relative to the largest element. Measured: 5e-15.
namespace smoothness {

// La2O3 (bohr) and the displacement direction along which the jumps were found
const double X0[5][3] = {{0.0, 0.944863, -0.3286479},
                         {0.0, -0.944863, -0.3286479},
                         {0.0, 1.88973, 1.5610781},
                         {0.0, 0.0, 1.5610781},
                         {0.0, -1.88973, 1.5610781}};
const double DZ[5] = {0.0, 0.8660254037844386, -0.28867513459481287, -0.28867513459481287,
                      -0.28867513459481287};

// La def2 ECP (46 core electrons): n, l, exponent, coefficient
const double LA_ECP[10][4] = {{2, 3, 4.0286, -36.010016}, {2, 0, 3.3099, 91.932177},
                              {2, 0, 1.655, -3.788764},   {2, 0, 4.0286, 36.010016},
                              {2, 1, 2.8368, 63.759486},  {2, 1, 1.4184, -0.647958},
                              {2, 1, 4.0286, 36.010016},  {2, 2, 2.0213, 36.116173},
                              {2, 2, 1.0107, 0.219114},   {2, 2, 4.0286, 36.010016}};

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
     {3,
      1,
      {17.721504317, 3.863550544, 1.0480920883},
      {3.2460107641875786, 2.5734959172661265, 1.1208372892488978}},
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
    U.addPrimitive(int(LA_ECP[k][0]), int(LA_ECP[k][1]), LA_ECP[k][2], LA_ECP[k][3], k == 9);
  GaussianShell sa(X[c.a.centre], c.a.l), sb(X[c.b.centre], c.b.l);
  for (size_t i = 0; i < c.a.exps.size(); i++) sa.addPrim(c.a.exps[i], c.a.coefs[i]);
  for (size_t i = 0; i < c.b.exps.size(); i++) sb.addPrim(c.b.exps[i], c.b.coefs[i]);
  TwoIndex<double> v;
  eng.compute_shell_pair(U, sa, sb, v);
  return v.data;
}

}  // namespace smoothness

TEST(Smoothness, La2O3) {
  using namespace smoothness;
  const ECPIntegral eng(3, 3, 0);
  const double h = 1e-3, tol = 1e-13;
  const int ks[6] = {-3, -2, -1, 1, 2, 3};
  for (const Case& c : CASES) {
    const std::vector<double> v0 = integrals(eng, c, c.t);
    std::vector<std::vector<double>> vk;
    for (const int k : ks) vk.push_back(integrals(eng, c, c.t + k * h));
    double scale = 0.0, worst = 0.0;
    for (const double x : v0) scale = std::max(scale, std::abs(x));
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
    EXPECT_LE(worst, tol) << c.what;
    info(c.what, worst);
  }
}

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  initFactorials();
  return RUN_ALL_TESTS();
}
