#include "gaussquad.hpp"
#include "gtest/gtest.h"
#include <cmath>
#include <iostream>

using namespace libecpint;

using ::testing::TestWithParam;
using ::testing::Values;

typedef GCQuadrature* CreateGrid();

GCQuadrature* CreateBigGrid() {
  GCQuadrature* gc = new GCQuadrature;
  gc->initGrid(1024, ONEPOINT);
  return gc;
}

GCQuadrature* CreateSmallGrid() {
  GCQuadrature* gc = new GCQuadrature;
  gc->initGrid(256, TWOPOINT);
  return gc;
}

double polynomial(const double r, const double* p, const int ix) {
  double result = 0;
  int n = int(p[0]);
  for (int i = 0; i <= n; i++) result += p[i + 1] * std::pow(r, i);
  return result;
}

double gaussian(const double r, const double* p, const int ix) {
  return p[1] * std::exp(-r * r * p[0]);
}

class QuadTest : public TestWithParam<CreateGrid*> {
 public:
  void SetUp() override { grid_ = (*GetParam())(); }

 protected:
  GCQuadrature* grid_;
};

TEST_P(QuadTest, CheckInit) { EXPECT_EQ(grid_->getN(), grid_->getX().size()); }

TEST_P(QuadTest, IntegratePoly) {
  std::function<double(double, const double*, int)> intgd = polynomial;
  double params[5] = {3, 1.0, -2.0, 3.0, -4.0};
  const auto integral_and_test =
      grid_->integrate(intgd, params, 1e-6, 0, grid_->getN() - 1);
  EXPECT_TRUE(integral_and_test.second);
  EXPECT_NEAR(integral_and_test.first, 4.0, 1e-6);
}

TEST_P(QuadTest, IntegrateGauss) {
  std::function<double(double, const double*, int)> intgd = gaussian;
  // integrate on [0, inf) instead of [-1, 1]
  grid_->transformZeroInf();
  double params[2] = {0.5, 2.4};
  const auto integral_and_test =
      grid_->integrate(intgd, params, 1e-6, 0, grid_->getN() - 1);
  EXPECT_TRUE(integral_and_test.second);
  EXPECT_NEAR(integral_and_test.first, 1.2 * std::sqrt(2.0 * M_PI), 1e-6);
}

TEST_P(QuadTest, TransformPoly) {
  std::function<double(double, const double*, int)> intgd = polynomial;
  // integrate on [10, 20] instead of [-1, 1]
  grid_->transformRMinMax(2.56, 14.375);
  double params[4] = {2, 0.2, 0.1, -0.003};
  const auto integral_and_test =
      grid_->integrate(intgd, params, 1e-6, 0, grid_->getN() - 1);
  EXPECT_TRUE(integral_and_test.second);
  EXPECT_NEAR(integral_and_test.first, 10.0, 1e-6);
}

INSTANTIATE_TEST_SUITE_P(GaussQuad, QuadTest,
                         Values(&CreateBigGrid, &CreateSmallGrid));

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
