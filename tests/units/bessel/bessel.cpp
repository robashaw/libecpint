#include "bessel.hpp"

#include <cmath>
#include <iostream>

#include "../../support/radial_reference.hpp"
#include "gtest/gtest.h"
#include "mathutil.hpp"
#include "radial.hpp"

using namespace libecpint;

class BesselTest : public testing::Test {
 protected:
  BesselFunction bessie;

  void SetUp() override {
    initFactorials();
    bessie.init(2, 1600, 200, 1e-15);
  }
};

TEST_F(BesselTest, CalculateBigL) {
  std::vector<double> values(4, 0.0);
  // check it catches maxL > lMax
  bessie.calculate(SMALL, 3, values);
}

TEST_F(BesselTest, CalculateSmallZ) {
  // check negative
  std::vector<double> values(3, 0.0);
  bessie.calculate(-1.0, 2, values);
  EXPECT_DOUBLE_EQ(values[0], 1.0);

  values.assign(3, 0.0);
  bessie.calculate(SMALL / 2, 2, values);
  EXPECT_NEAR(values[0], 1.0, 1e-7);
  EXPECT_NEAR(values[1], 1.66667e-8, 1e-12);

  double v0 = bessie.calculate(SMALL / 2, 0);
  double v1 = bessie.calculate(SMALL / 2, 1);
  double v2 = bessie.calculate(SMALL / 2, 2);
  EXPECT_NEAR(values[0], v0, 1e-12);
  EXPECT_NEAR(values[1], v1, 1e-12);
  EXPECT_NEAR(values[2], v2, 1e-12);
}

TEST_F(BesselTest, SingleVsMultiSmallZ) {
  // Regression test: the single-L small-z branch must follow K_L(z) ~ (1-z) z^L / (2L+1)!!,
  // i.e. agree (in relative terms) with the multi-L routine and the analytic value, even though
  // the magnitudes are tiny. A previous bug used (z/(2L+1))^L, wrong for L >= 2.
  double z = 1e-8;
  std::vector<double> v(3, 0.0);
  bessie.calculate(z, 2, v);

  double analytic2 = (1.0 - z) * z * z / 15.0;  // (2*2+1)!! = 15
  EXPECT_NEAR(v[2] / analytic2, 1.0, 1e-6);
  EXPECT_NEAR(bessie.calculate(z, 2) / analytic2, 1.0, 1e-6);
  EXPECT_NEAR(bessie.calculate(z, 2) / v[2], 1.0, 1e-9);

  double analytic1 = (1.0 - z) * z / 3.0;  // (2*1+1)!! = 3
  EXPECT_NEAR(bessie.calculate(z, 1) / analytic1, 1.0, 1e-6);
}

TEST_F(BesselTest, CalculateBigZ) {
  std::vector<double> values(2, 0.0);
  bessie.calculate(17.0, 1, values);
  EXPECT_NEAR(values[0], 0.0294118, 1e-7);
  EXPECT_NEAR(values[1], 0.0276817, 1e-7);

  double v0 = bessie.calculate(17.0, 0);
  double v1 = bessie.calculate(17.0, 1);
  EXPECT_NEAR(values[0], v0, 1e-12);
  EXPECT_NEAR(values[1], v1, 1e-12);
}

TEST_F(BesselTest, CalculateMidZ) {
  std::vector<double> values(2, 0.0);
  bessie.calculate(5.0, 1, values);
  EXPECT_NEAR(values[0], 0.0999955, 1e-7);
  EXPECT_NEAR(values[1], 0.0800054, 1e-7);

  double v0 = bessie.calculate(5.0, 0);
  double v1 = bessie.calculate(5.0, 1);
  EXPECT_NEAR(values[0], v0, 1e-12);
  EXPECT_NEAR(values[1], v1, 1e-12);
}

TEST_F(BesselTest, UpperBoundSmallZ) {
  std::vector<double> values(3, 0.0);
  bessie.calculate(SMALL / 2, 2, values);
  double ub0 = bessie.upper_bound(SMALL / 2, 0);
  double ub1 = bessie.upper_bound(SMALL / 2, 1);
  double ub2 = bessie.upper_bound(SMALL / 2, 2);
  EXPECT_TRUE(ub0 >= values[0]);
  EXPECT_TRUE(ub1 >= values[1]);
  EXPECT_TRUE(ub2 >= values[2]);
}

TEST_F(BesselTest, UpperBoundBigZ) {
  std::vector<double> values(3, 0.0);
  bessie.calculate(17.0, 2, values);
  double ub0 = bessie.upper_bound(17.0, 0);
  double ub1 = bessie.upper_bound(17.0, 1);
  double ub2 = bessie.upper_bound(17.0, 2);
  EXPECT_TRUE(ub0 >= values[0]);
  EXPECT_TRUE(ub1 >= values[1]);
  EXPECT_TRUE(ub2 >= values[2]);
}

TEST_F(BesselTest, UpperBoundMidZ) {
  std::vector<double> values(3, 0.0);
  bessie.calculate(5.054, 2, values);
  double ub0 = bessie.upper_bound(5.054, 0);
  double ub1 = bessie.upper_bound(5.054, 1);
  double ub2 = bessie.upper_bound(5.054, 2);
  EXPECT_TRUE(ub0 >= values[0]);
  EXPECT_TRUE(ub1 >= values[1]);
  EXPECT_TRUE(ub2 >= values[2]);
}

TEST(BesselAccuracy, AllOrdersAgainstSeries) {
  // K_l(z) <= 1, from both forms of calculate (all orders up to lMax, and a single order) as the
  // radial integrals use them, against the reference (power series, exact to ~1e-16), over the
  // table (z < 16), the large-z forms beyond it, and the change between them; absolute error, as
  // the integrals sum K_l times Gaussian factors. Orders up to 24 = 2 LIBECPINT_MAX_L + 8, the
  // most the radial integrals are initialised for
  initFactorials();
  const int lMax = 24;
  BesselFunction bessie(lMax, 1600, 200, libecpint::BESSEL_ACCURACY);
  std::vector<double> values(lMax + 1);
  for (int l = 0; l <= lMax; l++) {
    double worst = 0.0, worst_z = 0.0, worst1 = 0.0, worst1_z = 0.0;
    for (double z = 1e-3; z < 40.0; z += 0.00731) {
      bessie.calculate(z, lMax, values);
      const double ref = radref::K(l, z);
      const double err = std::abs(values[l] - ref), err1 = std::abs(bessie.calculate(z, l) - ref);
      if (!(err <= worst)) {
        worst = err;
        worst_z = z;
      }
      if (!(err1 <= worst1)) {
        worst1 = err1;
        worst1_z = z;
      }
    }
    EXPECT_LE(worst, 2e-15) << "all orders, l=" << l << ": worst at z=" << worst_z;
    EXPECT_LE(worst1, 2e-15) << "single order, l=" << l << ": worst at z=" << worst1_z;
  }
}

TEST(BesselAccuracy, IndependentOfMaxL) {
  // Each K_l is bitwise the same whatever maximum order is asked for, as the shared grid's cache
  // reuses values computed for a larger one (with the downward recurrence started from the order
  // asked for, they differed in the last bit, and results depended on the order of evaluation)
  initFactorials();
  const int lMax = 24;
  BesselFunction bessie(lMax, 1600, 200, libecpint::BESSEL_ACCURACY);
  std::vector<double> all(lMax + 1), some(lMax + 1);
  int differ = 0;
  for (double z = 1e-3; z < 60.0; z += 0.0371) {
    bessie.calculate(z, lMax, all);
    for (int maxL = 0; maxL < lMax; maxL++) {
      bessie.calculate(z, maxL, some);
      for (int l = 0; l <= maxL; l++) differ += some[l] != all[l];
    }
  }
  EXPECT_EQ(differ, 0);
}

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
