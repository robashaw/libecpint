#include "mathutil.hpp"
#include "multiarr.hpp"
#include "gtest/gtest.h"
#include <cmath>
#include <cstdlib>
#include <iostream>

using namespace libecpint;

TEST(MathUtil, Defines) {
  EXPECT_EQ(MAX_FAC, 100);
  EXPECT_EQ(MAX_DFAC, 200);
}

TEST(MathUtil, Gamma) {
  EXPECT_EQ(GAMMA[0], 1.7724538509055);
  EXPECT_EQ(GAMMA[29], 8.7178291200e10);
}

TEST(MathUtil, Factorials) {
  bool check = false;
#ifdef FAC_INIT
  check = true;
#endif
  EXPECT_FALSE(check);
}

TEST(RealSphericalHarmonics, Leq0) {
  // result should always be 1/sqrt(4*pi)
  double x = 2.0 * rand() / RAND_MAX - 1.0;  // cos(theta) in [-1, 1]
  double phi = M_PI * rand() / RAND_MAX;     // phi in [0, pi]
  TwoIndex<double> result = realSphericalHarmonics(0, x, phi);
  EXPECT_EQ(result.data.size(), 1);
  EXPECT_DOUBLE_EQ(result(0, 0), 1.0 / std::sqrt(4.0 * M_PI));
}

TEST(RealSphericalHarmonics, Lgt0) {
  // try with L=1 first, theta=pi/2 and phi=pi/4
  TwoIndex<double> result = realSphericalHarmonics(1, 0.0, M_PI / 4.0);
  EXPECT_EQ(result.data.size(), 6);
  EXPECT_NEAR(result(0, 0), 0.282095, 1e-6);
  EXPECT_DOUBLE_EQ(result(1, 1), 0.0);

  // then L=3, theta=PI/6 and phi = 3PI/5
  result = realSphericalHarmonics(3, 0.5 * std::sqrt(3), 0.6 * M_PI);
  EXPECT_EQ(result.data.size(), 28);
  EXPECT_NEAR(result(0, 0), 0.282095, 1e-6);
  EXPECT_NEAR(result(1, 1), 0.423142, 1e-6);
  EXPECT_NEAR(result(2, 2), 0.394239, 1e-6);
}

TEST(FrobeniusNorm, ZeroMat) {
  TwoIndex<double> mat;
  EXPECT_DOUBLE_EQ(frobenius_norm(mat), 0.0);

  mat.assign(5, 5, 0.0);
  EXPECT_DOUBLE_EQ(frobenius_norm(mat), 0.0);
}

TEST(FrobeniusNorm, NonZeroMat) {
  TwoIndex<double> mat(4, 4);
  for (int i = 0; i < 4; i++) mat(i, i) = 0.5;
  EXPECT_DOUBLE_EQ(frobenius_norm(mat), 1.0);

  mat(0, 2) = 4.1;
  mat(1, 3) = 0.05;
  mat(3, 0) = 0.842;
  EXPECT_NEAR(frobenius_norm(mat), 4.3036570495, 1e-10);
}

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
