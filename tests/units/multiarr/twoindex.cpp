#include "multiarr.hpp"
#include "gtest/gtest.h"

using namespace libecpint;

class TwoIndexTest : public testing::Test {
 protected:
  TwoIndex<double> matrix_;

  void SetUp() override { matrix_ = TwoIndex<double>(2, 3, 0.0); }
};

TEST_F(TwoIndexTest, Assign) {
  matrix_.assign(5, 4, 1.0);
  ASSERT_EQ(matrix_.dims[0], 5);
  ASSERT_EQ(matrix_.dims[1], 4);
  EXPECT_EQ(matrix_.data.size(), 20);
  EXPECT_DOUBLE_EQ(matrix_(2, 3), 1.0);
}

TEST_F(TwoIndexTest, Accessor) {
  EXPECT_DOUBLE_EQ(matrix_(1, 1), 0.0);
  matrix_(1, 1) = 5.4;
  EXPECT_DOUBLE_EQ(matrix_(1, 1), 5.4);
  EXPECT_DOUBLE_EQ(matrix_(0, 2), 0.0);
}

TEST_F(TwoIndexTest, Transpose) {
  matrix_(1, 0) = 0.9;
  matrix_(0, 1) = -0.9;
  TwoIndex<double> tmp = matrix_.transpose();
  EXPECT_DOUBLE_EQ(tmp(1, 0), -0.9);
  EXPECT_DOUBLE_EQ(tmp(0, 1), 0.9);
}

TEST_F(TwoIndexTest, Multiply) {
  matrix_(0, 0) = 0.9;
  matrix_(0, 2) = -8.2;
  matrix_.multiply(-0.5);
  EXPECT_DOUBLE_EQ(matrix_(1, 1), 0.0);
  EXPECT_DOUBLE_EQ(matrix_(0, 0), -0.45);
  EXPECT_DOUBLE_EQ(matrix_(0, 2), 4.1);
}

TEST_F(TwoIndexTest, CopyCtor) {
  TwoIndex<double> tmp;
  EXPECT_EQ(tmp.dims[0], 0);
  EXPECT_EQ(tmp.dims[1], 0);

  tmp = matrix_;
  EXPECT_EQ(tmp.dims[0], 2);
  EXPECT_EQ(tmp.dims[1], 3);

  tmp(1, 2) = 4.4;
  EXPECT_DOUBLE_EQ(tmp(1, 2), 4.4);
  EXPECT_DOUBLE_EQ(matrix_(1, 2), 0.0);
}

TEST_F(TwoIndexTest, Add) {
  TwoIndex<double> m2(2, 3);
  m2(0, 1) = 1.4;
  m2(1, 0) = -0.8;
  m2(1, 2) = 3.12;

  matrix_.add(m2);
  EXPECT_DOUBLE_EQ(matrix_(0, 1), 1.4);
  EXPECT_DOUBLE_EQ(matrix_(1, 0), -0.8);
  EXPECT_DOUBLE_EQ(matrix_(1, 2), 3.12);
  EXPECT_DOUBLE_EQ(matrix_(0, 0), 0);
}