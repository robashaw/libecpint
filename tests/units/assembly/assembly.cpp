/*
 * Unit tests of the bookkeeping in assembly.hpp, which turns per-(ECP, shell pair) results into
 * molecular integrals and derivatives. None of it affects accuracy, so it is tested separately from
 * the integrals: the packing functions exactly, on synthetic input; the assembly loops against
 * finite differences of the library's own assembled integrals (consistency, whatever their
 * accuracy), on systems covering every coincidence of the centres A, B and C.
 */
#include "assembly.hpp"

#include <array>
#include <cmath>
#include <vector>

#include "api.hpp"
#include "gtest/gtest.h"

using namespace libecpint;

namespace {

// A matrix whose elements identify it: value(m, i, j) = 1000 m + 10 i + j + 0.5
TwoIndex<double> tagged(const int m, const int na, const int nb) {
  TwoIndex<double> M(na, nb);
  for (int i = 0; i < na; i++)
    for (int j = 0; j < nb; j++) M(i, j) = 1000.0 * m + 10.0 * i + j + 0.5;
  return M;
}

const int NA = 3, NB = 6;  // p and d shells

// Index of (p, q) in a symmetric block (xx, xy, xz, yy, yz, zz)
int sym(const int p, const int q) {
  static const int ix[3][3] = {{0, 1, 2}, {1, 3, 4}, {2, 4, 5}};
  return ix[p][q];
}

// The full 9 x 9 Hessian of one element (i, j) over centres X, Y in {A, B, C} and coordinates p,
// q, from the 45 packed matrices [AA, AB, AC, BB, BC, CC]
double full(const std::array<TwoIndex<double>, 45>& r, const int X, const int p, const int Y,
            const int q, const int i, const int j) {
  if (X > Y) return full(r, Y, q, X, p, i, j);
  if (X == Y) {
    const int start[3] = {0, 24, 39};
    return r[start[X] + sym(p, q)](i, j);
  }
  const int start = X == 0 ? (Y == 1 ? 6 : 15) : 30;
  return r[start + 3 * p + q](i, j);
}

}  // namespace

TEST(PackFirst, AllOffCentre) {
  std::array<TwoIndex<double>, 3> QA, QB;
  for (int p = 0; p < 3; p++) {
    QA[p] = tagged(p, NA, NB);
    QB[p] = tagged(10 + p, NA, NB);
  }
  std::array<TwoIndex<double>, 9> r;
  pack_first_derivatives(QA, QB, false, false, NA, NB, r);
  for (int p = 0; p < 3; p++)
    for (int i = 0; i < NA; i++)
      for (int j = 0; j < NB; j++) {
        EXPECT_EQ(r[p](i, j), QA[p](i, j));
        EXPECT_EQ(r[3 + p](i, j), QB[p](i, j));
        EXPECT_EQ(r[6 + p](i, j), -(QA[p](i, j) + QB[p](i, j)));
      }
}

TEST(PackFirst, OnCentreCases) {
  std::array<TwoIndex<double>, 3> QA, QB;
  for (int p = 0; p < 3; p++) {
    QA[p] = tagged(p, NA, NB);
    QB[p] = tagged(10 + p, NA, NB);
  }
  for (const bool onA : {false, true})
    for (const bool onB : {false, true}) {
      std::array<TwoIndex<double>, 9> r;
      pack_first_derivatives(QA, QB, onA, onB, NA, NB, r);
      for (int p = 0; p < 3; p++)
        for (int i = 0; i < NA; i++)
          for (int j = 0; j < NB; j++) {
            // A shell on the ECP is moved by its other partner's derivative; C carries nothing
            const double a = onA ? (onB ? 0.0 : -QB[p](i, j)) : QA[p](i, j);
            const double b = onB ? (onA ? 0.0 : -QA[p](i, j)) : QB[p](i, j);
            EXPECT_EQ(r[p](i, j), a) << onA << onB;
            EXPECT_EQ(r[3 + p](i, j), b) << onA << onB;
            // translational invariance holds exactly in every case
            EXPECT_EQ(r[p](i, j) + r[3 + p](i, j) + r[6 + p](i, j), 0.0) << onA << onB;
          }
    }
}

TEST(PackSecond, AllOffCentre) {
  std::array<TwoIndex<double>, 6> QAA, QBB;
  std::array<TwoIndex<double>, 9> QAB;
  for (int k = 0; k < 6; k++) {
    QAA[k] = tagged(k, NA, NB);
    QBB[k] = tagged(20 + k, NA, NB);
  }
  for (int k = 0; k < 9; k++) QAB[k] = tagged(10 + k, NA, NB);
  std::array<TwoIndex<double>, 45> r;
  pack_second_derivatives(QAA, QAB, QBB, false, false, NA, NB, r);
  for (int i = 0; i < NA; i++)
    for (int j = 0; j < NB; j++)
      for (int p = 0; p < 3; p++)
        for (int q = 0; q < 3; q++) {
          EXPECT_EQ(full(r, 0, p, 0, q, i, j), QAA[sym(p, q)](i, j));
          EXPECT_EQ(full(r, 1, p, 1, q, i, j), QBB[sym(p, q)](i, j));
          EXPECT_EQ(full(r, 0, p, 1, q, i, j), QAB[3 * p + q](i, j));
          // d/dC = -(d/dA + d/dB), applied to each index
          EXPECT_EQ(full(r, 0, p, 2, q, i, j), -(QAA[sym(p, q)](i, j) + QAB[3 * p + q](i, j)));
          EXPECT_EQ(full(r, 1, p, 2, q, i, j), -(QBB[sym(p, q)](i, j) + QAB[3 * q + p](i, j)));
          for (int Y = 0; Y < 3; Y++) {
            double sum = 0.0;
            for (int X = 0; X < 3; X++) sum += full(r, X, p, Y, q, i, j);
            EXPECT_NEAR(sum, 0.0, 1e-9) << "invariance, Y = " << Y << " p = " << p << " q = " << q;
          }
        }
}

TEST(PackSecond, OnCentreCases) {
  std::array<TwoIndex<double>, 6> QAA, QBB;
  std::array<TwoIndex<double>, 9> QAB;
  for (int k = 0; k < 6; k++) {
    QAA[k] = tagged(k, NA, NB);
    QBB[k] = tagged(20 + k, NA, NB);
  }
  for (int k = 0; k < 9; k++) QAB[k] = tagged(10 + k, NA, NB);
  for (const bool onA : {false, true})
    for (const bool onB : {false, true}) {
      if (!onA && !onB) continue;
      std::array<TwoIndex<double>, 45> r;
      pack_second_derivatives(QAA, QAB, QBB, onA, onB, NA, NB, r);
      for (int i = 0; i < NA; i++)
        for (int j = 0; j < NB; j++)
          for (int p = 0; p < 3; p++)
            for (int q = 0; q < 3; q++) {
              // The shell off the ECP (if any) carries everything: Q = its second derivative
              const double Q =
                  onA && onB ? 0.0 : (onA ? QBB[sym(p, q)](i, j) : QAA[sym(p, q)](i, j));
              EXPECT_EQ(full(r, 0, p, 0, q, i, j), Q);
              EXPECT_EQ(full(r, 1, p, 1, q, i, j), Q);
              EXPECT_EQ(full(r, 0, p, 1, q, i, j), -Q);
              for (int Y = 0; Y < 3; Y++) {
                EXPECT_EQ(full(r, 2, p, Y, q, i, j), 0.0);
              }
            }
    }
}

TEST(UnpackHessian, RoundTrip) {
  const int natoms = 3, n3 = 9;
  std::vector<double> H(n3 * n3);
  for (int r = 0; r < n3; r++)
    for (int c = r; c < n3; c++) H[r * n3 + c] = H[c * n3 + r] = 1.0 + r + 0.1 * c;
  std::vector<double> packed(n3 * (n3 + 1) / 2, -1.0);
  for (int a = 0; a < natoms; a++)
    for (int b = a; b < natoms; b++)
      for (int p = 0; p < 3; p++)
        for (int q = 0; q < 3; q++) {
          if (a == b && q < p) continue;
          const int k =
              a == b ? H_START(a, a, natoms) + 3 + sym(p, q) : H_START(a, b, natoms) + 3 * p + q;
          packed[k] = H[(3 * a + p) * n3 + 3 * b + q];
        }
  for (const double v : packed) EXPECT_NE(v, -1.0) << "packed layout has a gap";
  EXPECT_EQ(unpack_hessian(packed, natoms), H);
}

TEST(OnEcp, Threshold) {
  ECP U(std::array<double, 3>{1.0, 2.0, 3.0}.data());
  GaussianShell on(std::array<double, 3>{1.0, 2.0, 3.0 + 5e-7}, 0);
  GaussianShell off(std::array<double, 3>{1.0, 2.0, 3.0 + 2e-6}, 0);
  EXPECT_TRUE(on_ecp(on, U));
  EXPECT_FALSE(on_ecp(off, U));
}

TEST(ShellPairScreen, NoneSkipsNothing) {
  const auto screen = ShellPairScreen::none(4, 3);
  for (int s = 0; s < 4; s++) EXPECT_EQ(screen.ecps_for(s, 3).size(), 3u);
  EXPECT_GT(screen.bounds[0] * screen.bounds[1], screen.tolerance);
}

// ---- Assembly loops, against finite differences of the assembled library integrals ----

namespace {

// A small system: atom 0 has shells and an ECP, atom 1 shells only, atom 2 an ECP only; so shell
// pairs and ECPs cover A = B = C, A = B != C, A = C != B, B = C != A and all distinct.
struct Tiny {
  static constexpr int lmax = 2;  // the largest shell angular momentum
  std::vector<std::array<double, 3>> atoms = {{0.1, 0.2, -0.3}, {1.3, -0.4, 0.9}, {-0.8, 1.1, 0.6}};
  void setup(ECPIntegrator& f, const int deriv) const {
    // shells: atom 0 s and p, atom 1 d
    std::vector<double> coords, exps = {3.0, 0.4, 1.1, 0.3, 0.9, 0.25},
                                coefs = {0.6, 0.5, 0.7, 0.4, 0.8, 0.3};
    const int ams[3] = {0, 1, 2}, lens[3] = {2, 2, 2}, shell_atom[3] = {0, 0, 1};
    for (const int a : shell_atom) coords.insert(coords.end(), atoms[a].begin(), atoms[a].end());
    f.set_gaussian_basis(3, coords.data(), exps.data(), coefs.data(), ams, lens);
    // ECPs on atoms 0 and 2: local and two semi-local channels
    std::vector<double> ecoords;
    for (const int a : {0, 2}) ecoords.insert(ecoords.end(), atoms[a].begin(), atoms[a].end());
    const double eexps[8] = {1.5, 2.0, 0.8, 1.2, 1.5, 2.0, 0.8, 1.2};
    const double ecoefs[8] = {-1.0, 3.0, 2.0, -0.5, -1.0, 3.0, 2.0, -0.5};
    const int eams[8] = {2, 0, 1, 1, 2, 0, 1, 1}, ens[8] = {2, 2, 2, 1, 2, 2, 2, 1};
    const int elens[2] = {4, 4};
    f.set_ecp_basis(2, ecoords.data(), eexps, ecoefs, eams, ens, elens);
    f.init(deriv);
  }
};

// Moves atom a of f by (dx, dy, dz): its shells and its ECP
void move(ECPIntegrator& f, const Tiny& t, const int a, const std::array<double, 3>& d) {
  for (auto& s : f.shells)
    if (s.atom_id == a)
      for (int k = 0; k < 3; k++) s.localCenter[k] = t.atoms[a][k] + d[k];
  for (int i = 0; i < f.ecps.getN(); i++) {
    ECP& U = f.ecps.getECP(i);
    if (U.atom_id == a) U.setPos(t.atoms[a][0] + d[0], t.atoms[a][1] + d[1], t.atoms[a][2] + d[2]);
  }
}

// Four-point central difference of quantity(f) with respect to coordinate p of atom a
template <typename Q>
std::vector<double> finite_difference(ECPIntegrator& f, const Tiny& t, const int a, const int p,
                                      Q quantity) {
  const double h = 1e-3, steps[4] = {-2.0, -1.0, 1.0, 2.0}, w[4] = {1.0, -8.0, 8.0, -1.0};
  std::vector<double> result;
  for (int k = 0; k < 4; k++) {
    std::array<double, 3> d = {0.0, 0.0, 0.0};
    d[p] = steps[k] * h;
    move(f, t, a, d);
    const std::vector<double> v = quantity(f);
    result.resize(v.size(), 0.0);
    for (size_t i = 0; i < v.size(); i++) result[i] += w[k] * v[i] / (12.0 * h);
  }
  move(f, t, a, {0.0, 0.0, 0.0});
  return result;
}

double max_abs(const std::vector<double>& v) {
  double m = 0.0;
  for (const double x : v) m = std::max(m, std::abs(x));
  return m;
}

}  // namespace

TEST(AssemblyLoops, FirstDerivativesMatchFiniteDifferences) {
  const Tiny t;
  ECPIntegrator f;
  t.setup(f, 1);
  ASSERT_EQ(f.natoms, 3);
  f.compute_first_derivs();
  for (int a = 0; a < 3; a++)
    for (int p = 0; p < 3; p++) {
      const auto fd = finite_difference(f, t, a, p, [](ECPIntegrator& g) {
        g.compute_integrals();
        return g.integrals.data;
      });
      const auto& d = f.first_derivs[3 * a + p].data;
      const double scale = max_abs(fd);
      for (size_t i = 0; i < d.size(); i++)
        EXPECT_NEAR(d[i], fd[i], 1e-8 * std::max(1.0, scale)) << "atom " << a << " coord " << p;
    }
}

TEST(AssemblyLoops, SecondDerivativesMatchFiniteDifferences) {
  const Tiny t;
  ECPIntegrator f;
  t.setup(f, 2);
  f.compute_second_derivs();
  const int natoms = 3;
  for (int b = 0; b < natoms; b++)
    for (int q = 0; q < 3; q++) {
      // derivatives of all first-derivative matrices with respect to coordinate q of atom b
      const auto fd = finite_difference(f, t, b, q, [](ECPIntegrator& g) {
        g.first_derivs.clear();
        g.compute_first_derivs();
        std::vector<double> v;
        for (const auto& M : g.first_derivs) v.insert(v.end(), M.data.begin(), M.data.end());
        return v;
      });
      const size_t n2 = f.ncart * f.ncart;
      for (int a = 0; a < natoms; a++)
        for (int p = 0; p < 3; p++) {
          const int row = 3 * a + p, col = 3 * b + q;
          const int lo = std::min(row, col), hi = std::max(row, col);
          const int A = lo / 3, B = hi / 3, P = lo % 3, Qc = hi % 3;
          const int k =
              A == B ? H_START(A, A, natoms) + 3 + sym(P, Qc) : H_START(A, B, natoms) + 3 * P + Qc;
          const auto& d = f.second_derivs[k].data;
          const double scale = max_abs(fd);
          for (size_t i = 0; i < n2; i++)
            EXPECT_NEAR(d[i], fd[row * n2 + i], 1e-7 * std::max(1.0, scale))
                << "d2/d(atom " << a << ", " << p << ") d(atom " << b << ", " << q << ")";
        }
    }
}

TEST(AssemblyLoops, GradientSumsToZero) {
  const Tiny t;
  ECPIntegrator f;
  t.setup(f, 1);
  f.compute_first_derivs();
  for (int p = 0; p < 3; p++)
    for (size_t i = 0; i < f.first_derivs[0].data.size(); i++) {
      double sum = 0.0, scale = 0.0;
      for (int a = 0; a < 3; a++) {
        sum += f.first_derivs[3 * a + p].data[i];
        scale = std::max(scale, std::abs(f.first_derivs[3 * a + p].data[i]));
      }
      EXPECT_NEAR(sum, 0.0, 1e-14 * std::max(1.0, scale));
    }
}

TEST(AssemblyLoops, ContractionsMatchDenseDerivatives) {
  const Tiny t;
  ECPIntegrator f;
  t.setup(f, 2);
  f.compute_integrals();
  f.compute_first_derivs();
  f.compute_second_derivs();
  const int n = f.ncart, n3 = 3 * f.natoms;
  std::vector<double> D(n * n);
  for (int k = 0; k < n * n; k++) D[k] = std::sin(1.0 + 0.37 * k);  // not symmetric
  auto contract = [&](const TwoIndex<double>& M) {
    double s = 0.0;
    for (int k = 0; k < n * n; k++) s += D[k] * M.data[k];
    return s;
  };
  const auto g = f.compute_gradient(D.data());
  for (int i = 0; i < n3; i++) EXPECT_NEAR(g[i], contract(f.first_derivs[i]), 1e-12);
  std::vector<double> packed;
  for (const auto& M : f.second_derivs) packed.push_back(contract(M));
  const auto expected = unpack_hessian(packed, f.natoms);
  const auto h = f.compute_hessian(D.data());
  for (int i = 0; i < n3 * n3; i++) EXPECT_NEAR(h[i], expected[i], 1e-12);
}

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
