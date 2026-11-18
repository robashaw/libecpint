/*
 * Magnetic field derivatives of single ECP blocks
 * (ECPIntegral::compute_shell_pair_derivative_Bfield) in the configurations with special cases,
 * against the reference engine's field mode (tests/support/integral_reference.hpp), which
 * multiplies the basis functions by their field factors pointwise rather than raising the shells as
 * the library does, and against exact properties of the derivatives:
 *
 *   SameCentreOffEcp  A and B on one centre off the ECP: not zero (the first implementation
 *                     returned zero for these)
 *   OneOnEcp          one shell on the ECP, whose term vanishes (and is skipped by the library)
 *   Distinct          A, B and C all distinct
 *   AllOnEcp          A = B = C: exactly zero
 *   Antisymmetry      the block for (B | A) is minus the transpose of that for (A | B)
 *   Translation       moving every centre by T adds [(A - B) x T]_k <A|U|B>: the origin dependence
 *                     of GIAO derivatives, fixing the convention of absolute positions
 *   AngularMomenta    every shell angular momentum the library allows (LIBECPINT_MAX_L - 1, the
 *                     shells being raised by one)
 *   OriginalInputs    the inputs of the regression tests that came with the implementation (a very
 *                     diffuse s shell against an f shell, and a shell on the ECP)
 *
 * Comparisons with the reference engine use the suite's tolerance (compare.hpp), with the
 * conditioning scales of reference_integrals. The whole-molecule derivatives, through
 * ECPIntegrator::compute_bfield_derivs, are checked by the suite (Suite.Stored, Suite.Live,
 * Suite.Screening), and their assembly by units/assembly. Each test prints its worst error, so that
 * a drift in accuracy shows before it fails.
 */
#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <cstdio>
#include <functional>
#include <string>
#include <vector>

#include "compare.hpp"
#include "ecp.hpp"
#include "ecpint.hpp"
#include "gshell.hpp"
#include "integral_reference.hpp"
#include "mathutil.hpp"
#include "multiarr.hpp"
#include "reference.hpp"

using namespace libecpint;
using V3 = std::array<double, 3>;

namespace {

struct ShellSpec {
  V3 centre;
  int l;
  std::vector<std::array<double, 2>> prims;  // (exponent, coefficient)
};

struct EcpSpec {
  V3 centre;
  std::vector<intref::EcpPrimitive> prims;  // n, l, a, d: d r^(n-2) exp(-a r^2) in channel l
};

// Off the origin, so that the [(A - B) x C] term is not zero
const V3 C0 = {0.1, -0.2, 1.4}, P1 = {-0.5, -0.5, -0.5}, P2 = {0.9, 0.4, -0.3};

// Semi-local channels 0-2 and a local channel 3 (type 1 integrals), with n = 1, 2, 3 and 4
std::vector<intref::EcpPrimitive> mixed_ecp() {
  return {{2, 0, 2.015, 0.438}, {1, 0, 0.7, 0.3},  {2, 1, 1.226, 0.418},  {4, 1, 0.9, -0.2},
          {2, 2, 3.416, 10.48}, {3, 2, 1.1, -1.5}, {2, 3, 3.574, -4.478}, {2, 3, 0.8, 0.5}};
}

/// Worst error over a set of blocks: as a fraction of the suite's tolerance, relative to the
/// largest reference value, and that value
struct Error {
  double ratio = 0.0, rel = 0.0, big = 0.0;
};

/// The library's and (if reference) the reference engine's blocks, with their conditioning
/// scales, for every pair (sa, sb), sa >= sb, of shells against one ECP
struct Blocks {
  std::vector<GaussianShell> shells;
  ECP U;
  std::vector<std::vector<std::array<TwoIndex<double>, 3>>> lib, ref;  // [sa][sb][k]
  std::vector<std::vector<double>> scale;                              // [sa][sb]

  Blocks(const EcpSpec& ecp, const std::vector<ShellSpec>& specs, const bool reference = true)
      : U(ecp.centre.data()) {
    int maxLU = 0, maxLB = 0;
    for (const auto& p : ecp.prims) {
      U.addPrimitive(p.n, p.l, p.a, p.d);
      maxLU = std::max(maxLU, p.l);
    }
    U.sort();
    std::vector<intref::Shell> rshells;
    for (const auto& sh : specs) {
      GaussianShell g(sh.centre, sh.l);
      intref::Shell r{sh.centre, sh.l, {}, {}, 0};
      for (const auto& p : sh.prims) {
        g.addPrim(p[0], p[1]);
        r.exps.push_back(p[0]);
        r.coefs.push_back(p[1]);
      }
      shells.push_back(g);
      rshells.push_back(r);
      maxLB = std::max(maxLB, sh.l);
    }
    const ECPIntegral engine(maxLB, maxLU, 1, 0.0);  // unscreened
    const int n = shells.size();
    lib.resize(n);
    for (int sa = 0; sa < n; sa++)
      for (int sb = 0; sb <= sa; sb++) {
        lib[sa].emplace_back();
        engine.compute_shell_pair_derivative_Bfield(U, shells[sa], shells[sb], lib[sa][sb]);
      }
    if (!reference) return;
    const intref::Engine field(rshells, {intref::Ecp{ecp.centre, ecp.prims, 0}}, 0, true);
    ref.resize(n);
    scale.resize(n);
    for (int sa = 0; sa < n; sa++)
      for (int sb = 0; sb <= sa; sb++) {
        ref[sa].emplace_back();
        scale[sa].push_back(
            testsupport::bfield_conditioning_scale(engine, U, shells[sa], shells[sb]));
        for (int k = 0; k < 3; k++) {
          TwoIndex<double>& r = ref[sa][sb][k];
          r.assign(shells[sa].ncartesian(), shells[sb].ncartesian(), 0.0);
          for (int ca = 0; ca < r.dims[0]; ca++)
            for (int cb = 0; cb < r.dims[1]; cb++)
              r(ca, cb) = field.get_bfield(0, sa, sb, k, ca, cb);
        }
      }
  }

  /// The worst error of the library over the pairs selected by use(sa, sb), against the tolerance
  /// DEFAULT_REL_TOL * max |reference| + CONDITIONING_TOL * scale over them
  Error error(const std::function<bool(int, int)>& use) const {
    Error e;
    const int n = lib.size();
    for (int sa = 0; sa < n; sa++)
      for (int sb = 0; sb <= sa; sb++)
        if (use(sa, sb))
          for (int k = 0; k < 3; k++)
            for (const double r : ref[sa][sb][k].data) e.big = std::max(e.big, std::abs(r));
    double err = 0.0;
    for (int sa = 0; sa < n; sa++)
      for (int sb = 0; sb <= sa; sb++) {
        if (!use(sa, sb)) continue;
        const double tol =
            testsupport::DEFAULT_REL_TOL * e.big + testsupport::CONDITIONING_TOL * scale[sa][sb];
        for (int k = 0; k < 3; k++)
          for (size_t i = 0; i < ref[sa][sb][k].data.size(); i++) {
            const double d = std::abs(lib[sa][sb][k].data[i] - ref[sa][sb][k].data[i]);
            err = std::max(err, d);
            e.ratio = std::max(e.ratio, d / tol);
          }
      }
    e.rel = err / e.big;
    return e;
  }
};

void report(const std::string& what, const Error& e) {
  std::printf("[ info ] %s: worst error %.1e of the tolerance, %.2e relative (largest %.2e)\n",
              what.c_str(), e.ratio, e.rel, e.big);
}

bool same(const V3& a, const V3& b) {
  return std::abs(a[0] - b[0]) + std::abs(a[1] - b[1]) + std::abs(a[2] - b[2]) < 1e-12;
}

// s, p and d shells on P1, p and d on the ECP, f on P2: every configuration of centres
const std::vector<ShellSpec> MIXED_SHELLS = {
    {P1, 0, {{3.1, 0.4}, {0.45, 0.7}}},
    {P1, 1, {{1.4, 0.6}, {0.3, 0.5}}},
    {P1, 2, {{0.9, 1.0}}},
    {C0, 1, {{2.2, 0.8}, {0.5, 0.4}}},
    {C0, 2, {{1.2, 1.0}}},
    {P2, 3, {{0.7, 1.0}, {0.25, 0.3}}},
};

const Blocks& mixed() {
  static const Blocks b({C0, mixed_ecp()}, MIXED_SHELLS);
  return b;
}

const V3& centre(const int s) { return MIXED_SHELLS[s].centre; }

}  // namespace

TEST(BField, SameCentreOffEcp) {
  const Error e = mixed().error([](const int sa, const int sb) {
    return same(centre(sa), centre(sb)) && !same(centre(sa), C0);
  });
  EXPECT_GT(e.big, 1e-2);  // far from zero
  EXPECT_LE(e.ratio, 1.0);
  report("same centre off the ECP", e);
}

TEST(BField, OneOnEcp) {
  const Error e = mixed().error(
      [](const int sa, const int sb) { return same(centre(sa), C0) != same(centre(sb), C0); });
  EXPECT_GT(e.big, 1e-2);
  EXPECT_LE(e.ratio, 1.0);
  report("one shell on the ECP", e);
}

TEST(BField, Distinct) {
  const Error e = mixed().error([](const int sa, const int sb) {
    return !same(centre(sa), centre(sb)) && !same(centre(sa), C0) && !same(centre(sb), C0);
  });
  EXPECT_GT(e.big, 1e-2);
  EXPECT_LE(e.ratio, 1.0);
  report("all centres distinct", e);
}

TEST(BField, AllOnEcp) {
  const Blocks& b = mixed();
  for (int sa : {3, 4})
    for (int sb = 3; sb <= sa; sb++)
      for (int k = 0; k < 3; k++) {
        for (const double v : b.lib[sa][sb][k].data) EXPECT_EQ(v, 0.0);
        for (const double v : b.ref[sa][sb][k].data) EXPECT_EQ(v, 0.0);
      }
}

TEST(BField, Antisymmetry) {
  const Blocks& b = mixed();
  const ECPIntegral engine(3, 3, 1, 0.0);
  double err = 0.0, big = 0.0;
  for (size_t sa = 0; sa < b.shells.size(); sa++)
    for (size_t sb = 0; sb <= sa; sb++) {
      std::array<TwoIndex<double>, 3> swapped;
      engine.compute_shell_pair_derivative_Bfield(b.U, b.shells[sb], b.shells[sa], swapped);
      for (int k = 0; k < 3; k++)
        for (int ca = 0; ca < b.shells[sa].ncartesian(); ca++)
          for (int cb = 0; cb < b.shells[sb].ncartesian(); cb++) {
            err = std::max(err, std::abs(swapped[k](cb, ca) + b.lib[sa][sb][k](ca, cb)));
            big = std::max(big, std::abs(b.lib[sa][sb][k](ca, cb)));
          }
    }
  EXPECT_LE(err, 1e-13 * big);
  std::printf("[ info ] antisymmetry: worst %.2e relative (largest %.2e)\n", err / big, big);
}

TEST(BField, Translation) {
  const V3 T = {0.7, -1.3, 2.1};
  auto shift = [&](const V3& x) { return V3{x[0] + T[0], x[1] + T[1], x[2] + T[2]}; };
  std::vector<ShellSpec> moved = MIXED_SHELLS;
  for (auto& s : moved) s.centre = shift(s.centre);
  const Blocks& b = mixed();
  const Blocks m({shift(C0), mixed_ecp()}, moved, false);
  const ECPIntegral engine(3, 3, 1, 0.0);
  double err = 0.0, big = 0.0;
  for (size_t sa = 0; sa < b.shells.size(); sa++)
    for (size_t sb = 0; sb <= sa; sb++) {
      TwoIndex<double> Q;
      engine.compute_shell_pair(b.U, b.shells[sa], b.shells[sb], Q);
      const double* A = b.shells[sa].center();
      const double* B = b.shells[sb].center();
      const double D[3] = {A[0] - B[0], A[1] - B[1], A[2] - B[2]};
      const double DxT[3] = {D[1] * T[2] - D[2] * T[1], D[2] * T[0] - D[0] * T[2],
                             D[0] * T[1] - D[1] * T[0]};
      for (int k = 0; k < 3; k++)
        for (size_t i = 0; i < Q.data.size(); i++) {
          const double expected = b.lib[sa][sb][k].data[i] + DxT[k] * Q.data[i];
          err = std::max(err, std::abs(m.lib[sa][sb][k].data[i] - expected));
          big = std::max(big, std::abs(expected));
        }
    }
  EXPECT_LE(err, 1e-13 * big);
  std::printf("[ info ] translation: worst %.2e relative (largest %.2e)\n", err / big, big);
}

// Each angular momentum on P2 with a p shell on P1, and with itself and an s shell on P2 (same
// centre). The ECP has no local part, whose quadrature in the reference engine is slow for high l
// pairs (the other tests cover it).
TEST(BField, AngularMomenta) {
  const std::vector<intref::EcpPrimitive> semilocal = {
      {2, 0, 2.015, 0.438}, {2, 1, 1.226, 0.418}, {2, 2, 3.416, 10.48}, {2, 3, 1.0, 0.0}};
  for (int l = 0; l < LIBECPINT_MAX_L; l++) {
    const Blocks b({C0, semilocal},
                   {{P1, 1, {{0.35, 0.5}}}, {P2, 0, {{0.6, 1.0}}}, {P2, l, {{1.3, 0.6}}}});
    const Error e = b.error([](const int sa, const int) { return sa == 2; });
    EXPECT_LE(e.ratio, 1.0) << "l = " << l;
    report("l = " + std::to_string(l), e);
  }
}

// The inputs of the two regression tests of the original implementation (bfield_test1/2), which
// stored the library's own output (including zeros for the same-centre blocks)
TEST(BField, OriginalInputs) {
  // a local channel with zero coefficient, so no type 1 integrals
  const std::vector<intref::EcpPrimitive> prims = {{2, 4, 1.0, 0.0},
                                                   {2, 0, 2.015185434, 0.437903605},
                                                   {2, 1, 1.225772391, 0.418499875},
                                                   {2, 2, 3.415601371, 10.482661239},
                                                   {2, 3, 3.574403408, -4.478284427}};
  const V3 A = {0.5, 1.5, 0.9}, B = {-0.5, -0.5, -0.5};
  const std::vector<std::array<double, 2>> primsB = {{0.02, 1.1}, {1.41, 0.9}, {8.71, 0.4}};
  // 1: a p shell on the ECP and a d shell off it
  const Blocks one({A, prims}, {{A, 1, {{0.15, 0.3}, {2.49, 0.9}, {10.11, 1.2}}}, {B, 2, primsB}});
  // 2: a very diffuse s shell and an f shell, both off the ECP
  const Blocks two({{0.0, 0.0, 1.4}, prims},
                   {{A, 0, {{0.005, 0.3}, {2.49, 0.9}, {10.11, 1.2}}}, {B, 3, primsB}});
  for (const Blocks* b : {&one, &two}) {
    const Error e = b->error([](int, int) { return true; });
    EXPECT_GT(e.big, 1e-2);
    EXPECT_LE(e.ratio, 1.0);
    report(b == &one ? "original test 1" : "original test 2", e);
  }
}

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  initFactorials();
  return RUN_ALL_TESTS();
}
