#include "multiarr.hpp"
#include "gtest/gtest.h"
#include <algorithm>

using namespace libecpint;

class SevenIndexTest : public testing::Test {
 protected:
  SevenIndex<double> matrix_;

  void SetUp() override {
    matrix_ = SevenIndex<double>(2, 9, 3, 7, 5, 4, 1);
    std::fill(matrix_.data.begin(), matrix_.data.end(), 0.0);
  }
};

TEST_F(SevenIndexTest, Accessor) {
  EXPECT_DOUBLE_EQ(matrix_(1, 1, 1, 2, 4, 1, 0), 0.0);
  matrix_(1, 1, 1, 2, 4, 1, 0) = 5.4;
  EXPECT_DOUBLE_EQ(matrix_(1, 1, 1, 2, 4, 1, 0), 5.4);
  EXPECT_DOUBLE_EQ(matrix_(0, 3, 1, 4, 2, 0, 0), 0.0);
}

TEST_F(SevenIndexTest, CopyCtor) {
  SevenIndex<double> tmp;
  EXPECT_EQ(tmp.dims[0], 0);
  EXPECT_EQ(tmp.dims[1], 0);
  EXPECT_EQ(tmp.dims[2], 0);
  EXPECT_EQ(tmp.dims[3], 0);
  EXPECT_EQ(tmp.dims[4], 0);
  EXPECT_EQ(tmp.dims[5], 0);
  EXPECT_EQ(tmp.dims[6], 0);

  tmp = matrix_;
  EXPECT_EQ(tmp.dims[0], 2);
  EXPECT_EQ(tmp.dims[1], 9);
  EXPECT_EQ(tmp.dims[2], 3);
  EXPECT_EQ(tmp.dims[3], 7);
  EXPECT_EQ(tmp.dims[4], 5);
  EXPECT_EQ(tmp.dims[5], 4);
  EXPECT_EQ(tmp.dims[6], 1);

  tmp(0, 1, 1, 4, 3, 0, 0) = 4.9;
  EXPECT_DOUBLE_EQ(tmp(0, 1, 1, 4, 3, 0, 0), 4.9);
  EXPECT_DOUBLE_EQ(matrix_(0, 0, 2, 1, 4, 0, 0), 0.0);
}