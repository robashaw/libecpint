/*
 * Unit tests of the shared radial grid (shared_grid.hpp): the LRU cache, the accuracy of the grid
 * and of the cached factors against direct evaluation, and the independence of the results from
 * what the caches already hold. Each thread has its own caches, so a new thread starts empty.
 */
#include "shared_grid.hpp"

#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <thread>
#include <vector>

#include "../../support/radial_reference.hpp"
#include "bessel.hpp"
#include "mathutil.hpp"
#include "radial.hpp"

using namespace libecpint;
using namespace libecpint::detail;

namespace {

// Runs f in a new thread, i.e. with empty caches
template <typename F>
void in_new_thread(F&& f) {
  std::thread t(std::forward<F>(f));
  t.join();
}

// An ECP with r^-2, r^-1 and r^0 terms and a range of exponents
ECP test_ecp() {
  ECP U(std::array<double, 3>{0.0, 0.0, 0.0}.data());
  U.addPrimitive(0, 0, 0.9, 3.0);
  U.addPrimitive(2, 0, 12.0, 40.0);
  U.addPrimitive(1, 1, 2.5, -6.0);
  U.addPrimitive(2, 1, 0.3, 0.8);
  U.addPrimitive(2, 2, 1.1, -2.0);
  U.addPrimitive(0, 2, 0.05, 0.4);
  return U;
}

// A shell with diffuse and tight primitives, at distance A from the ECP
GaussianShell test_shell(const double A) {
  GaussianShell s(std::array<double, 3>{0.0, 0.0, A}, 2);
  s.addPrim(150.0, 0.2);  // above SG_TIGHT, so not on the shared grid
  s.addPrim(20.0, 0.5);
  s.addPrim(1.3, 0.7);
  s.addPrim(0.04, 0.3);
  return s;
}

// Composite Gauss-Legendre integral of f over [0, rend]
template <typename F>
double integrate(F&& f, const double rend) {
  std::vector<double> x, w;
  radref::gauss_legendre16(x, w);
  const int panels = 4000;
  const double h = rend / panels;
  double sum = 0.0;
  for (int p = 0; p < panels; p++)
    for (int i = 0; i < 16; i++) sum += 0.5 * h * w[i] * f((p + 0.5 + 0.5 * x[i]) * h);
  return sum;
}

}  // namespace

TEST(LruCache, FindAndInsert) {
  LruCache<int> cache(1000);
  EXPECT_EQ(cache.find({1.0, 2.0}), nullptr);
  cache.insert({1.0, 2.0}, 7, 10);
  cache.insert({3.0}, 8, 10);
  ASSERT_NE(cache.find({1.0, 2.0}), nullptr);
  EXPECT_EQ(*cache.find({1.0, 2.0}), 7);
  EXPECT_EQ(*cache.find({3.0}), 8);
  EXPECT_EQ(cache.find({1.0}), nullptr);
  cache.insert({3.0}, 9, 10);  // replaces
  EXPECT_EQ(*cache.find({3.0}), 9);
  EXPECT_EQ(cache.size(), 2u);
}

TEST(LruCache, EvictsLeastRecentlyUsed) {
  LruCache<int> cache(100);
  cache.insert({1.0}, 1, 40);
  cache.insert({2.0}, 2, 40);
  cache.find({1.0});           // 2 is now the least recently used
  cache.insert({3.0}, 3, 40);  // over budget: evicts 2
  EXPECT_EQ(cache.size(), 2u);
  EXPECT_NE(cache.find({1.0}), nullptr);
  EXPECT_EQ(cache.find({2.0}), nullptr);
  EXPECT_NE(cache.find({3.0}), nullptr);
}

TEST(LruCache, KeepsOneItemOverBudget) {
  LruCache<int> cache(10);
  cache.insert({1.0}, 1, 40);
  EXPECT_NE(cache.find({1.0}), nullptr);
  cache.insert({2.0}, 2, 40);
  EXPECT_EQ(cache.size(), 1u);
  EXPECT_NE(cache.find({2.0}), nullptr);
}

TEST(SharedGrid, IntegratesGaussiansUpToTightLimit) {
  // sum_i Uw[lam][i] f(r_i) against int U_lam(r) f(r) dr, for f(r) = r^(N+2) exp(-z (r - p)^2)
  // with z up to 2 SG_TIGHT (a product of two primitives with exponents up to SG_TIGHT), and powers
  // up to r^40: N = LA + LB, plus up to r^(l1 + l2) from the Bessel functions near r = 0. The
  // powers narrow the integrand most at p = 0, the reason for the small p
  const ECP U = test_ecp();
  const auto grid = shared_grid(U);
  const double rend = std::sqrt(80.0 / U.min_exp) + 5.0;
  for (int lam = 0; lam <= U.getL(); lam++)
    for (const double z : {0.02, 1.0, 10.0, 2.0 * SG_TIGHT})
      for (const double p : {0.0, 0.1, 0.3, 1.0, 3.0})
        for (int N = 0; N <= 38; N += 2) {
          auto f = [&](const double r) {
            return std::pow(r, N + 2) * std::exp(-z * (r - p) * (r - p));
          };
          double sum = 0.0;
          for (size_t i = 0; i < grid->r.size(); i++) sum += grid->Uw[lam][i] * f(grid->r[i]);
          const double ref =
              integrate([&](const double r) { return U.evaluate(r, lam) * f(r); }, rend);
          const double scale =
              integrate([&](const double r) { return std::abs(U.evaluate(r, lam) * f(r)); }, rend);
          EXPECT_NEAR(sum, ref, 1e-13 * scale)
              << "lam=" << lam << " z=" << z << " p=" << p << " N=" << N;
        }
}

TEST(SharedFactor, MatchesDirectEvaluation) {
  // F_l(r) = sum over primitives with a <= SG_TIGHT of c exp(-a (r - A)^2) K_l(2 a A r), with the
  // library's K_l (their accuracy is tested in units/bessel)
  const ECP U = test_ecp();
  const BesselFunction bessie(6, 1600, 200, BESSEL_ACCURACY);
  const auto grid = shared_grid(U);
  const size_t n = grid->r.size();
  for (const double A : {0.0, 0.8, 2.5}) {
    const GaussianShell s = test_shell(A);
    const int lmax = 6;
    const auto F = shared_factor(*grid, s, A, lmax, bessie);
    std::vector<double> bes(lmax + 1);
    for (int l = 0; l <= lmax; l++) {
      double scale = 0.0;
      std::vector<double> ref(n, 0.0);
      for (size_t i = 0; i < n; i++) {
        const double r = grid->r[i];
        for (int k = 0; k < s.nprimitive(); k++)
          if (s.exp(k) <= SG_TIGHT) {
            bessie.calculate(2.0 * s.exp(k) * A * r, lmax, bes);
            ref[i] += s.coef(k) * std::exp(-s.exp(k) * (r - A) * (r - A)) * bes[l];
          }
        scale = std::max(scale, std::abs(ref[i]));
      }
      for (size_t i = 0; i < n; i++)
        ASSERT_NEAR(F->F[l * n + i], ref[i], 1e-15 * scale)
            << "A=" << A << " l=" << l << " r=" << grid->r[i];
    }
  }
}

TEST(SharedFactor, IndependentOfCachedLmax) {
  // A factor cached for a larger lmax is returned for a smaller one: its values must be identical
  const ECP U = test_ecp();
  const BesselFunction bessie(6, 1600, 200, BESSEL_ACCURACY);
  const GaussianShell s = test_shell(1.7);
  std::vector<double> fresh;
  in_new_thread([&]() { fresh = shared_factor(*shared_grid(U), s, 1.7, 2, bessie)->F; });
  in_new_thread([&]() {
    const auto grid = shared_grid(U);
    shared_factor(*grid, s, 1.7, 6, bessie);
    const auto cached = shared_factor(*grid, s, 1.7, 2, bessie);
    EXPECT_EQ(cached->lmax, 6);
    const size_t n = grid->r.size();
    ASSERT_EQ(fresh.size(), 3 * n);
    EXPECT_TRUE(std::equal(fresh.begin(), fresh.end(), cached->F.begin()));
  });
}

TEST(SharedWeights, DefinitionAndCachedNmax) {
  const ECP U = test_ecp();
  const int lam = 1;
  std::vector<double> fresh;
  in_new_thread([&]() { fresh = shared_weights(*shared_grid(U), lam, 3)->W; });
  in_new_thread([&]() {
    const auto grid = shared_grid(U);
    const size_t n = grid->r.size();
    const auto big = shared_weights(*grid, lam, 9);
    // W_N(r_i) = Uw[lam][i] r_i^(N + 2), and S its suffix maxima of |W_N|
    for (int N = 0; N <= 9; N++) {
      double smax = 0.0;
      for (size_t i = n; i-- > 0;) {
        const double w = grid->Uw[lam][i] * std::pow(grid->r[i], N + 2);
        EXPECT_NEAR(big->W[N * n + i], w, 1e-14 * std::abs(w) + 1e-300) << "N=" << N << " i=" << i;
        smax = std::max(smax, std::abs(big->W[N * n + i]));
        EXPECT_EQ(big->S[N * (n + 1) + i], smax) << "N=" << N << " i=" << i;
      }
    }
    const auto small = shared_weights(*grid, lam, 3);
    EXPECT_EQ(small, big);  // the cached one, with values for N <= 3 identical to a fresh build
    ASSERT_EQ(fresh.size(), 4 * n);
    EXPECT_TRUE(std::equal(fresh.begin(), fresh.end(), big->W.begin()));
  });
}

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  initFactorials();  // for BesselFunction, as ECPIntegral does
  return RUN_ALL_TESTS();
}
