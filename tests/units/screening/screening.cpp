/*
 * Unit tests of the screening bounds and estimates in screening.hpp, separately from the accuracy
 * of the integrals: a bound must bound, and whatever is skipped must be below the tolerance.
 *
 *   Radial      type2_bound >= |G|, type2_log_bound_cheap >= log type2_bound, and the integrand
 *               is negligible outside integrand_window (G by the reference quadrature)
 *   ShellBound  shell_bound >= the function it bounds, sampled
 *   Systems     on every system of the manifest: channel_bounds >= the largest element of each
 *               channel's contribution, and every (shell pair, ECP) skipped by the default
 *               shell-pair screening has integrals (and derivatives) below the tolerance
 */
#include "screening.hpp"

#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <sstream>
#include <string>
#include <vector>

#include "../../support/radial_reference.hpp"
#include "ecpint.hpp"
#include "system.hpp"

using namespace libecpint;
using namespace libecpint::screening;

namespace {

// Parameters of a radial integral G(k, l1, l2) (notation of screening.hpp)
struct Radial {
  double n, a, b, A, B;
};

std::vector<Radial> radial_sweep() {
  std::vector<Radial> v;
  for (const double n : {0.5, 5.0, 50.0})
    for (const double a : {0.05, 1.0, 30.0})
      for (const double b : {0.1, 8.0})
        for (const double A : {0.0, 0.7, 3.0})
          for (const double B : {0.0, 1.5, 4.0}) v.push_back({n, a, b, A, B});
  return v;
}

const std::array<std::pair<int, int>, 3> K_RANGES = {{{0, 2}, {2, 6}, {1, 10}}};

std::string describe(const Radial& p, const int k, const int l1, const int l2) {
  std::ostringstream s;
  s << "k=" << k << " l1=" << l1 << " l2=" << l2 << " n=" << p.n << " a=" << p.a << " b=" << p.b
    << " A=" << p.A << " B=" << p.B;
  return s.str();
}

// The integrand of G, and its integral over [r0, r1] by composite Gauss-Legendre
double integrand(const int k, const int l1, const int l2, const Radial& p, const double r) {
  return std::pow(r, k) *
         std::exp(-p.n * r * r - p.a * (r - p.A) * (r - p.A) - p.b * (r - p.B) * (r - p.B)) *
         radref::K(l1, 2.0 * p.a * p.A * r) * radref::K(l2, 2.0 * p.b * p.B * r);
}

double integral(const int k, const int l1, const int l2, const Radial& p, const double r0,
                const double r1) {
  if (r1 <= r0) return 0.0;
  std::vector<double> x, w;
  radref::gauss_legendre16(x, w);
  const int panels = 400;
  const double h = (r1 - r0) / panels;
  double sum = 0.0;
  for (int q = 0; q < panels; q++)
    for (int i = 0; i < 16; i++)
      sum += 0.5 * h * w[i] * integrand(k, l1, l2, p, r0 + (q + 0.5 + 0.5 * x[i]) * h);
  return sum;
}

}  // namespace

TEST(Radial, Type2BoundBoundsIntegral) {
  // The bound does not depend on l1, l2 and K_l decreases with l, so l = 0 is the tightest case
  const std::array<std::pair<int, int>, 2> ls = {{{0, 0}, {2, 1}}};
  for (const Radial& p : radial_sweep())
    for (const auto& [kmin, kmax] : K_RANGES) {
      const double bound = type2_bound(kmin, kmax, p.n, p.a, p.b, p.A, p.B);
      for (const int k : {kmin, (kmin + kmax) / 2, kmax})
        for (const auto& [l1, l2] : ls) {
          const double g = radref::G(k, l1, l2, p.n, p.a, p.b, p.A, p.B);
          EXPECT_LE(std::abs(g), bound)
              << describe(p, k, l1, l2) << " kmin=" << kmin << " kmax=" << kmax;
        }
    }
}

TEST(Radial, CheapLogBoundBoundsLogBound) {
  for (const double n : {0.02, 0.5, 5.0, 50.0, 500.0})
    for (const double a : {0.01, 0.2, 1.0, 30.0, 1e4})
      for (const double b : {0.01, 0.3, 8.0, 1e3})
        for (const double A : {0.0, 0.3, 1.5, 4.0, 10.0})
          for (const double B : {0.0, 0.9, 3.0, 12.0})
            for (int kmax = 0; kmax <= 12; kmax++)
              for (int kmin = 0; kmin <= kmax; kmin++) {
                const double exact = std::log(type2_bound(kmin, kmax, n, a, b, A, B));
                const double cheap = type2_log_bound_cheap(kmax, n, a, b, A, B);
                EXPECT_GE(cheap, exact)
                    << describe({n, a, b, A, B}, kmax, 0, 0) << " kmin=" << kmin;
              }
}

TEST(Radial, IntegrandNegligibleOutsideWindow) {
  // The integrand is positive, so the parts outside the window are compared to the whole
  const double eps = 1e-15;
  for (const Radial& p : radial_sweep())
    for (const auto& [kmin, kmax] : K_RANGES)
      for (const int lmax : {0, 3}) {
        const Window win = integrand_window(kmin, kmax, lmax, lmax, p.n, p.a, p.b, p.A, p.B);
        double rend = std::sqrt(80.0 / p.n) + 1.0;  // as in radref::G
        for (int it = 0; it < 20; it++)
          rend = std::sqrt((80.0 + kmax * std::log(std::max(rend, 1.0))) / p.n) + 1.0;
        for (const int k : {kmin, kmax})
          for (const int l : {0, lmax}) {
            const double total = radref::G(k, l, l, p.n, p.a, p.b, p.A, p.B);
            const double outside =
                integral(k, l, l, p, 0.0, win.rmin) + integral(k, l, l, p, win.rmax, rend);
            EXPECT_LE(outside, eps * total)
                << describe(p, k, l, l) << " window [" << win.rmin << ", " << win.rmax << "]";
          }
      }
}

TEST(ShellBound, BoundsSampledFunction) {
  // shell_bound bounds |r - A|^la exp(-alpha |r - A|^2) exp(-eta r^2), with the ECP at the origin
  // and the shell at (0, 0, A); sampled in (z, rho), rho the distance from the axis
  for (int la = 0; la <= 5; la++)
    for (const double alpha : {0.03, 0.5, 10.0})
      for (const double A : {0.0, 0.5, 2.0, 6.0})
        for (const double eta : {0.1, 1.0, 20.0}) {
          const double bound = shell_bound(la, alpha, A * A, eta);
          double worst = 0.0;
          const double range = A + 3.0 * std::sqrt((la + 1.0) / alpha);
          for (int iz = -400; iz <= 400; iz++)
            for (int ir = 0; ir <= 100; ir++) {
              const double z = range * iz / 400.0, rho = range * ir / 100.0;
              const double d2 = rho * rho + (z - A) * (z - A);
              const double f = std::pow(d2, 0.5 * la) * std::exp(-alpha * d2) *
                               std::exp(-eta * (rho * rho + z * z));
              worst = std::max(worst, f);
            }
          EXPECT_GE(bound, worst) << "la=" << la << " alpha=" << alpha << " A=" << A
                                  << " eta=" << eta;
        }
}

// ---- On the systems of the manifest ----

namespace {

std::vector<std::string> system_names() {
  std::vector<std::string> names;
  for (const auto& spec : testsupport::read_manifest(testsupport::data_dir() + "/systems.txt"))
    names.push_back(spec.name);
  return names;
}

class Systems : public ::testing::TestWithParam<std::string> {
 protected:
  void SetUp() override { sys_ = testsupport::load_system(GetParam()); }
  testsupport::System sys_;
};

// The shell pair relative to the ECP, as compute_shell_pair builds it
ShellPairData pair_data(const ECP& U, const GaussianShell& A, const GaussianShell& B) {
  ShellPairData d;
  for (int c = 0; c < 3; c++) {
    d.A[c] = A.center()[c] - U.center()[c];
    d.B[c] = B.center()[c] - U.center()[c];
  }
  d.LA = A.am();
  d.LB = B.am();
  d.maxLBasis = std::max(d.LA, d.LB);
  d.ncartA = (d.LA + 1) * (d.LA + 2) / 2;
  d.ncartB = (d.LB + 1) * (d.LB + 2) / 2;
  d.A2 = d.A[0] * d.A[0] + d.A[1] * d.A[1] + d.A[2] * d.A[2];
  d.Am = std::sqrt(d.A2);
  d.A_on_ecp = d.Am < 1e-6;
  d.B2 = d.B[0] * d.B[0] + d.B[1] * d.B[1] + d.B[2] * d.B[2];
  d.Bm = std::sqrt(d.B2);
  d.B_on_ecp = d.Bm < 1e-6;
  const double R[3] = {d.A[0] - d.B[0], d.A[1] - d.B[1], d.A[2] - d.B[2]};
  d.RAB2 = R[0] * R[0] + R[1] * R[1] + R[2] * R[2];
  d.RABm = std::sqrt(d.RAB2);
  return d;
}

// U with the coefficients of every channel but l set to zero (GaussianECP stores n - 2)
ECP only_channel(const ECP& U, const int l) {
  ECP V(U.center());
  for (const auto& g : U.gaussians) V.addPrimitive(g.n + 2, g.l, g.a, g.l == l ? g.d : 0.0);
  return V;
}

double max_abs(const TwoIndex<double>& m) {
  double v = 0.0;
  for (const double x : m.data) v = std::max(v, std::abs(x));
  return v;
}

}  // namespace

TEST_P(Systems, ChannelBoundsBoundChannels) {
  libecpint::ECPIntegrator f;
  sys_.setup(f, 0);
  const ECPIntegral engine(f.maxLB, f.ecps.getMaxL(), 0, 0.0);
  int nviolations = 0;
  for (int e = 0; e < f.ecps.getN(); e++) {
    const ECP& U = f.ecps.getECP(e);
    std::vector<ECP> channels;
    for (int l = 0; l <= U.getL(); l++) channels.push_back(only_channel(U, l));
    std::vector<double> bounds(U.getL() + 1);
    for (size_t a = 0; a < f.shells.size(); a++)
      for (size_t b = 0; b <= a; b++) {
        const GaussianShell &A = f.shells[a], &B = f.shells[b];
        channel_bounds(U, A, B, pair_data(U, A, B), bounds.data());
        for (int l = 0; l <= U.getL(); l++) {
          if (count_channel(U, l) == 0) continue;
          TwoIndex<double> block;
          engine.compute_shell_pair(channels[l], A, B, block);
          const double actual = max_abs(block);
          if (actual > bounds[l] && nviolations++ < 10)
            ADD_FAILURE() << "ECP " << e << " channel " << l << ", shells " << a << " (l=" << A.am()
                          << ") and " << b << " (l=" << B.am() << "): bound " << bounds[l]
                          << " < largest element " << actual;
        }
      }
  }
  EXPECT_EQ(nviolations, 0);
}

TEST_P(Systems, SkippedShellPairsBelowTolerance) {
  // Mirrors the loops of assembly.hpp: for shells s2 <= s1, the pair is skipped against ECP i when
  // shell s1's bound is at most shell_tolerance, or the product of the bounds at most tolerance.
  // Every block skipped must then be at most tolerance (each skipped pair by itself).
  for (int deriv = 0; deriv <= sys_.spec.deriv; deriv++) {
    libecpint::ECPIntegrator f;
    sys_.setup(f, deriv);
    const int necp = f.ecps.getN();
    const std::vector<double> bounds = shell_ecp_bounds(f.shells, f.ecps, deriv);
    const ECPIntegral engine(f.maxLB, f.ecps.getMaxL(), deriv, 0.0);
    double worst = 0.0;
    int nskipped = 0, nviolations = 0;
    for (size_t s1 = 0; s1 < f.shells.size(); s1++)
      for (size_t s2 = 0; s2 <= s1; s2++)
        for (int i = 0; i < necp; i++) {
          const double b1 = bounds[s1 * necp + i], b2 = bounds[s2 * necp + i];
          if (b1 > f.shell_tolerance && b1 * b2 > f.tolerance) continue;
          nskipped++;
          const ECP& U = f.ecps.getECP(i);
          const GaussianShell &A = f.shells[s1], &B = f.shells[s2];
          double actual = 0.0;
          if (deriv == 0) {
            TwoIndex<double> block;
            engine.compute_shell_pair(U, A, B, block);
            actual = max_abs(block);
          } else if (deriv == 1) {
            std::array<TwoIndex<double>, 9> blocks;
            engine.compute_shell_pair_derivative(U, A, B, blocks);
            for (const auto& m : blocks) actual = std::max(actual, max_abs(m));
          } else {
            std::array<TwoIndex<double>, 45> blocks;
            engine.compute_shell_pair_second_derivative(U, A, B, blocks);
            for (const auto& m : blocks) actual = std::max(actual, max_abs(m));
          }
          worst = std::max(worst, actual);
          if (actual > f.tolerance && nviolations++ < 10)
            ADD_FAILURE() << "deriv " << deriv << ": skipped shells " << s1 << " (l=" << A.am()
                          << ") and " << s2 << " (l=" << B.am() << ") against ECP " << i
                          << " (bounds " << b1 << ", " << b2 << ") have elements up to " << actual;
        }
    EXPECT_EQ(nviolations, 0) << "deriv " << deriv << ": " << nskipped
                              << " (shell pair, ECP) skipped, largest element " << worst;
    std::printf("[ info ] %s deriv %d: %d skipped, largest element %.3e (tolerance %.1e)\n",
                GetParam().c_str(), deriv, nskipped, worst, f.tolerance);
  }
}

INSTANTIATE_TEST_SUITE_P(Manifest, Systems, ::testing::ValuesIn(system_names()),
                         [](const auto& info) { return info.param; });

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
