#include "ecp.hpp"
#include "gtest/gtest.h"
#include <iostream>

using namespace libecpint;

class ECPTest : public testing::Test {
 protected:
  ECP u;

  void SetUp() override {
    u.center_ = {0.1, -0.1, 0.5};
    u.addPrimitive(2, 0, 0.6, 1.2, false);
    u.addPrimitive(2, 1, 1.4, -0.12, false);
    u.addPrimitive(2, 2, 0.02, 0.44, false);
    u.addPrimitive(2, 3, 8.7, 0.679, true);
  }
};

TEST(ECP, DefaultCtor) {
  ECP newU;
  EXPECT_EQ(newU.N, 0);
  EXPECT_EQ(newU.L, -1);
  EXPECT_DOUBLE_EQ(newU.center_[0], 0.0);
  EXPECT_DOUBLE_EQ(newU.center_[1], 0.0);
  EXPECT_DOUBLE_EQ(newU.center_[2], 0.0);
}

TEST(ECP, PosCtor) {
  double C[3] = {0.9, 0.4, 0.1};
  ECP newU(C);
  EXPECT_EQ(newU.N, 0);
  EXPECT_EQ(newU.L, -1);
  EXPECT_DOUBLE_EQ(newU.center_[0], 0.9);
  EXPECT_DOUBLE_EQ(newU.center_[1], 0.4);
  EXPECT_DOUBLE_EQ(newU.center_[2], 0.1);
}

TEST_F(ECPTest, CopyCtor) {
  ECP newU(u);
  EXPECT_EQ(newU.N, u.N);
  EXPECT_EQ(newU.L, u.L);
  EXPECT_EQ(newU.gaussians.size(), 4);
  newU.setPos(-0.4, 0.0, 1.8);
  EXPECT_DOUBLE_EQ(newU.center_[2], 1.8);
  EXPECT_DOUBLE_EQ(u.center_[2], 0.5);
}

TEST_F(ECPTest, AddPrim) {
  EXPECT_EQ(u.getN(), 4);
  u.addPrimitive(1, 4, 0.94, 1.111, true);
  EXPECT_EQ(u.getN(), 5);
  EXPECT_EQ(u.getL(), 4);
  EXPECT_EQ(u.getGaussian(4).n, -1);
}

TEST_F(ECPTest, SetPos) {
  EXPECT_DOUBLE_EQ(u.center_[0], 0.1);
  u.setPos(5.4, -1.8, 0.2);
  EXPECT_DOUBLE_EQ(u.center_[0], 5.4);
  EXPECT_DOUBLE_EQ(u.center_[1], -1.8);
  EXPECT_DOUBLE_EQ(u.center_[2], 0.2);
}

TEST_F(ECPTest, Type1) {
  EXPECT_FALSE(u.noType1());
  u.getGaussian(3).d = 0.0;
  EXPECT_TRUE(u.noType1());
}

TEST_F(ECPTest, Evaluate) {
  EXPECT_NEAR(u.evaluate(0.5, 1), -0.0845626, 1e-7);
  EXPECT_NEAR(u.evaluate(1.8, 3), 3.89024e-13, 1e-13);
  EXPECT_NEAR(u.evaluate(5.0, 4), 0.0, 1e-15);
}

TEST(ECP, Sort) {
  ECP newU;
  newU.addPrimitive(2, 4, 0.3, 1.8, false);
  newU.addPrimitive(2, 1, 0.7, 0.22, false);
  newU.addPrimitive(2, 3, 0.1, -0.99, false);
  newU.sort();
  EXPECT_EQ(newU.getGaussian(0).l, 1);
  EXPECT_EQ(newU.getGaussian(1).l, 3);
  EXPECT_EQ(newU.getGaussian(2).l, 4);
}