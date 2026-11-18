/*
 * ECP integrals with all centres coincident, which have closed forms, checked against both the
 * library and the test reference engine (tests/support/integral_reference.hpp). The second
 * validates the engine that the suite's reference values come from, on cases where the exact
 * answer is known.
 *
 *   RefAtom     s, p and d shells against three ECPs, against exact values from sympy
 *               (generate_ref_atom.py, writing ref_atom.txt: angular monomial integrals over the
 *               sphere and half-integer Gamma functions, nothing shared with either engine)
 *   ScreenScan  s and p shells against a one-term local ECP, over exponents where the type 1
 *               integrals once collapsed, against one-line closed forms
 *
 * Each test prints its worst errors, so that a drift in accuracy shows before it fails.
 */
#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

#include "../support/integral_reference.hpp"
#include "ecp.hpp"
#include "ecpint.hpp"
#include "gshell.hpp"
#include "mathutil.hpp"
#include "multiarr.hpp"

using namespace libecpint;

namespace {

struct ShellSpec {
  int l;
  std::vector<std::array<double, 2>> prims;  // (exponent, coefficient)
};

/// Integrals over every pair (ia, ib), ia >= ib, of shells at the origin against one ECP, from the
/// library and from the reference engine, in the library's order (row-major blocks)
void integrals(const std::vector<intref::EcpPrimitive>& ecp, const std::vector<ShellSpec>& shells,
               const std::vector<std::array<int, 2>>& pairs, std::vector<double>& lib,
               std::vector<double>& ref) {
  double O[3] = {0.0, 0.0, 0.0};
  const std::array<double, 3> origin = {0.0, 0.0, 0.0};
  ECP U(O);
  int maxLU = 0, maxLB = 0;
  for (const auto& p : ecp) {
    U.addPrimitive(p.n, p.l, p.a, p.d);
    maxLU = std::max(maxLU, p.l);
  }
  U.sort();
  std::vector<GaussianShell> lshells;
  std::vector<intref::Shell> rshells;
  for (const auto& sh : shells) {
    GaussianShell g(origin, sh.l);
    intref::Shell r{origin, sh.l, {}, {}, 0};
    for (const auto& p : sh.prims) {
      g.addPrim(p[0], p[1]);
      r.exps.push_back(p[0]);
      r.coefs.push_back(p[1]);
    }
    lshells.push_back(g);
    rshells.push_back(r);
    maxLB = std::max(maxLB, sh.l);
  }
  const ECPIntegral engine(maxLB, maxLU, 0, 0.0);  // unscreened
  const intref::Engine reference(rshells, {intref::Ecp{origin, ecp, 0}}, 0);
  for (const auto& p : pairs) {
    TwoIndex<double> v;
    engine.compute_shell_pair(U, lshells[p[0]], lshells[p[1]], v);
    for (int ca = 0; ca < v.dims[0]; ca++)
      for (int cb = 0; cb < v.dims[1]; cb++) {
        lib.push_back(v(ca, cb));
        ref.push_back(reference.get(0, p[0], p[1], 0, 0, ca, cb));
      }
  }
}

/// Worst error of values against exact ones: relative where exact != 0, absolute where it is 0
struct Errors {
  double rel = 0.0, abs = 0.0;
};
Errors errors(const std::vector<double>& values, const std::vector<double>& exact) {
  Errors e;
  for (size_t i = 0; i < exact.size(); i++) {
    const double d = std::abs(values[i] - exact[i]);
    if (exact[i] != 0.0)
      e.rel = std::max(e.rel, d / std::abs(exact[i]));
    else
      e.abs = std::max(e.abs, d);
  }
  return e;
}

}  // namespace

// Scenarios as in generate_ref_atom.py (keep the two in step):
//   A: moderate exponents, mixed channels (general sanity);
//   B: a steep local channel (exponents up to 1e5) across s/p/d shells and diffuse AOs, where the
//      type 1 grid, scaled by the AO exponents alone, once collapsed, as did l > 0 screening;
//   C: a diffuse AO pair against a standard local primitive (exponent ratio ~250), the ordinary
//      augmented-basis case that screening once zeroed.
// Measured: 7e-15 for both (the library's was 5e-13 while the GAMMA table in mathutil.hpp held
// 14-digit half-integer values, amplified 18x where the local and semi-local parts cancel).
TEST(Analytic, RefAtom) {
  const double reltol = 1e-13, abstol = 1e-15;
  std::vector<double> exact;
  {
    std::ifstream in(std::string(ANALYTIC_DIR) + "/ref_atom.txt");
    ASSERT_TRUE(in.is_open()) << "cannot open ref_atom.txt";
    for (double x; in >> x;) exact.push_back(x);
  }
  const std::vector<std::array<int, 2>> pairs = {{0, 0}, {1, 0}, {1, 1}, {2, 0}, {2, 1}, {2, 2}};
  std::vector<double> lib, ref;
  integrals(
      {{2, 0, 1.5, 2.0}, {3, 0, 0.7, 0.4}, {2, 1, 1.2, -1.0}, {2, 2, 0.9, 0.5}, {4, 2, 0.6, 0.25}},
      {{0, {{1.1, 1.0}}}, {1, {{0.8, 1.0}, {2.5, 0.7}}}, {2, {{0.7, 1.0}}}}, pairs, lib, ref);
  integrals({{2, 0, 1000.0, 1.0}, {2, 0, 100000.0, 0.5}},
            {{0, {{0.05, 1.0}}}, {1, {{0.05, 1.0}}}, {2, {{0.5, 1.0}}}}, pairs, lib, ref);
  integrals({{2, 0, 10.0, 1.0}}, {{0, {{0.02, 1.0}}}, {1, {{0.02, 1.0}}}}, {{0, 0}, {1, 0}, {1, 1}},
            lib, ref);
  ASSERT_EQ(exact.size(), lib.size()) << "ref_atom.txt does not match the scenarios";

  for (const auto* v : {&lib, &ref}) {
    const char* who = v == &lib ? "library" : "reference engine";
    for (size_t i = 0; i < exact.size(); i++) {
      if (exact[i] != 0.0)
        EXPECT_LE(std::abs((*v)[i] - exact[i]), reltol * std::abs(exact[i]))
            << who << ", value " << i << ": " << (*v)[i] << " against " << exact[i];
      else
        EXPECT_LE(std::abs((*v)[i]), abstol) << who << ", value " << i << " should be 0";
    }
    const Errors e = errors(*v, exact);
    std::printf("[ info ] ref_atom, %s: worst relative %.2e, absolute (zeros) %.2e\n", who, e.rel,
                e.abs);
  }
}

// On-centre s and p_x shells (one primitive, exponent alpha, coefficient 1) against a local ECP
// r^(n - 2) exp(-g r^2), with c = 2 alpha + g:
//   <s|U|s>     n = 1: 2 pi / c          n = 2: pi^1.5 / c^1.5
//   <p_x|U|p_x> n = 1: 2 pi / (3 c^2)    n = 2: pi^1.5 / (2 c^2.5)
// These once collapsed to ~0 for an ECP much steeper than the AO pair (g >> 2 alpha), the type 1
// grid being scaled by the AO exponents alone, and for l >= 1 through a hard-coded screening
// tolerance. All values are checked relative to themselves, the integrals being unscreened.
// Measured: library 3e-15, reference engine 3e-16.
TEST(Analytic, ScreenScan) {
  const double reltol = 1e-13;
  double worst_lib = 0.0, worst_ref = 0.0;
  for (const int n : {1, 2})
    for (const double alpha : {1e-2, 1.0, 1e2, 1e4})
      for (const double g : {1e0, 1e2, 1e4, 1e6, 1e7}) {
        const double c = 2.0 * alpha + g;
        const double s_exact = n == 1 ? 2.0 * M_PI / c : std::pow(M_PI / c, 1.5);
        const double p_exact =
            n == 1 ? 2.0 * M_PI / (3.0 * c * c) : std::pow(M_PI, 1.5) / (2.0 * std::pow(c, 2.5));
        std::vector<double> lib, ref;
        integrals({{n, 0, g, 1.0}}, {{0, {{alpha, 1.0}}}, {1, {{alpha, 1.0}}}}, {{0, 0}, {1, 1}},
                  lib, ref);
        // lib/ref: <s|U|s>, then the p block, whose first element is <p_x|U|p_x>
        for (const auto& [exact, k] : {std::pair<double, int>{s_exact, 0}, {p_exact, 1}}) {
          const double el = std::abs(lib[k] - exact) / exact, er = std::abs(ref[k] - exact) / exact;
          worst_lib = std::max(worst_lib, el);
          worst_ref = std::max(worst_ref, er);
          EXPECT_LE(el, reltol) << "library, " << (k ? "p" : "s") << " n = " << n
                                << ", alpha = " << alpha << ", g = " << g;
          EXPECT_LE(er, reltol) << "reference engine, " << (k ? "p" : "s") << " n = " << n
                                << ", alpha = " << alpha << ", g = " << g;
        }
      }
  std::printf("[ info ] screen_scan: worst relative, library %.2e, reference engine %.2e\n",
              worst_lib, worst_ref);
}

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  initFactorials();
  return RUN_ALL_TESTS();
}
