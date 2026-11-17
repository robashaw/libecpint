/*
 * Unit tests of the derivative recurrences of ECPIntegral (left_shell_derivative,
 * left_shell_second_derivative, mixed_second_derivative), which build the derivatives of a shell
 * pair from integrals over shells of angular momentum l +- 1. Checked against finite differences of
 * compute_shell_pair with respect to the shell centres, i.e. for consistency with the integrals,
 * whatever their accuracy.
 */
#include <array>
#include <cmath>
#include <functional>

#include "ecpint.hpp"
#include "gtest/gtest.h"

using namespace libecpint;

namespace {

struct TestPair {
  ECP U;
  GaussianShell A, B;
  TestPair(const int la, const int lb)
      : U(std::array<double, 3>{0.2, -0.1, 0.3}.data()),
        A(std::array<double, 3>{0.9, 0.4, -0.5}, la),
        B(std::array<double, 3>{-0.6, 1.0, 0.7}, lb) {
    U.addPrimitive(2, 2, 1.4, -0.8);
    U.addPrimitive(2, 0, 2.2, 4.0);
    U.addPrimitive(2, 1, 0.9, 2.5);
    U.addPrimitive(1, 1, 1.6, -0.7, true);
    A.addPrim(2.5, 0.5);
    A.addPrim(0.35, 0.6);
    B.addPrim(1.2, 0.7);
    B.addPrim(0.2, 0.4);
  }
};

// Four-point central difference of f() with respect to x (restored afterwards)
TwoIndex<double> fd(double& x, const std::function<TwoIndex<double>()>& f) {
  const double h = 1e-3, x0 = x;
  const double steps[4] = {-2.0, -1.0, 1.0, 2.0}, w[4] = {1.0, -8.0, 8.0, -1.0};
  TwoIndex<double> result;
  for (int k = 0; k < 4; k++) {
    x = x0 + steps[k] * h;
    TwoIndex<double> v = f();
    if (k == 0) result.assign(v.dims[0], v.dims[1], 0.0);
    for (size_t i = 0; i < v.data.size(); i++) result.data[i] += w[k] * v.data[i] / (12.0 * h);
  }
  x = x0;
  return result;
}

void expect_close(const TwoIndex<double>& value, const TwoIndex<double>& ref, const double rel,
                  const std::string& what) {
  ASSERT_EQ(value.data.size(), ref.data.size()) << what;
  double scale = 1.0;
  for (const double r : ref.data) scale = std::max(scale, std::abs(r));
  for (size_t i = 0; i < ref.data.size(); i++)
    EXPECT_NEAR(value.data[i], ref.data[i], rel * scale) << what << " element " << i;
}

class Recurrences : public ::testing::TestWithParam<std::pair<int, int>> {};

}  // namespace

TEST_P(Recurrences, FirstDerivative) {
  const auto [la, lb] = GetParam();
  TestPair s(la, lb);
  ECPIntegral engine(std::max(la, lb), 2, 1);
  std::array<TwoIndex<double>, 3> QA;
  engine.left_shell_derivative(s.U, s.A, s.B, QA);
  for (int p = 0; p < 3; p++) {
    const auto ref = fd(s.A.localCenter[p], [&]() {
      TwoIndex<double> v;
      engine.compute_shell_pair(s.U, s.A, s.B, v);
      return v;
    });
    expect_close(QA[p], ref, 1e-8, "d/dA_" + std::to_string(p));
  }
}

TEST_P(Recurrences, SecondDerivatives) {
  const auto [la, lb] = GetParam();
  TestPair s(la, lb);
  ECPIntegral engine(std::max(la, lb), 2, 2);
  std::array<TwoIndex<double>, 6> QAA;
  std::array<TwoIndex<double>, 9> QAB;
  engine.left_shell_second_derivative(s.U, s.A, s.B, QAA);
  engine.mixed_second_derivative(s.U, s.A, s.B, QAB);
  const int sym[3][3] = {{0, 1, 2}, {1, 3, 4}, {2, 4, 5}};
  for (int q = 0; q < 3; q++) {
    // d/dA_q and d/dB_q of the first derivatives d/dA_p
    auto first = [&]() {
      std::array<TwoIndex<double>, 3> Q;
      engine.left_shell_derivative(s.U, s.A, s.B, Q);
      return Q;
    };
    for (int p = 0; p < 3; p++) {
      const auto dAA = fd(s.A.localCenter[q], [&]() { return first()[p]; });
      expect_close(QAA[sym[p][q]], dAA, 1e-7, "d2/dA_p dA_q");
      const auto dAB = fd(s.B.localCenter[q], [&]() { return first()[p]; });
      expect_close(QAB[3 * p + q], dAB, 1e-7, "d2/dA_p dB_q");
    }
  }
}

INSTANTIATE_TEST_SUITE_P(ShellPairs, Recurrences,
                         ::testing::Values(std::make_pair(0, 0), std::make_pair(1, 0),
                                           std::make_pair(0, 2), std::make_pair(2, 1),
                                           std::make_pair(1, 3), std::make_pair(6, 2),
                                           std::make_pair(4, 6)));

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
