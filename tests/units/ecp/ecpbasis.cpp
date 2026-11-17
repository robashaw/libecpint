#include <iostream>
#include <stdexcept>

#include "ecp.hpp"
#include "gtest/gtest.h"

using namespace libecpint;

class ECPBasisTest : public testing::Test {
 protected:
  ECPBasis basis;

  void SetUp() override {
    // load in ECP for silicon
    basis.addECP_from_file(14, {0.0, -0.5, 0.5}, "testbasis.ecp");
  }
};

TEST_F(ECPBasisTest, NMaxL) {
  EXPECT_EQ(basis.getMaxL(), 2);
  EXPECT_EQ(basis.getN(), 1);
}

TEST_F(ECPBasisTest, ECPCore) {
  EXPECT_EQ(basis.getECPCore(14), 10);
  EXPECT_EQ(basis.getECPCore(33), 0);
}

TEST_F(ECPBasisTest, AddFromFile) {
  basis.addECP_from_file(33, {0.0, -0.5, 2.4}, "testbasis.ecp");
  EXPECT_EQ(basis.getMaxL(), 3);
  EXPECT_EQ(basis.getN(), 2);
  EXPECT_EQ(basis.getECPCore(33), 28);
}

// testbasis_bse.ecp is the Basis Set Exchange MOLPRO export of LANL2DZ for Si and As: one
// primitive per line, comments, an orbital basis block, and no spin-orbit count in the header.
TEST_F(ECPBasisTest, BSELayout) {
  basis.addECP_from_file(33, {0.0, -0.5, 2.4}, "testbasis.ecp");
  ECPBasis bse;
  bse.addECP_from_file(14, {0.0, -0.5, 0.5}, "testbasis_bse.ecp");
  bse.addECP_from_file(33, {0.0, -0.5, 2.4}, "testbasis_bse.ecp");

  EXPECT_EQ(bse.getECPCore(14), 10);
  EXPECT_EQ(bse.getECPCore(33), 28);
  ASSERT_EQ(bse.getN(), 2);
  EXPECT_EQ(bse.getECP(0).getN(), 16);
  for (int i = 0; i < 2; i++) {
    const ECP& ref = basis.getECP(i);
    const ECP& U = bse.getECP(i);
    ASSERT_EQ(U.getN(), ref.getN());
    EXPECT_EQ(U.getL(), ref.getL());
    for (int j = 0; j < ref.getN(); j++) {
      EXPECT_EQ(U.getGaussian(j).n, ref.getGaussian(j).n);
      EXPECT_EQ(U.getGaussian(j).l, ref.getGaussian(j).l);
      EXPECT_EQ(U.getGaussian(j).a, ref.getGaussian(j).a);
      EXPECT_EQ(U.getGaussian(j).d, ref.getGaussian(j).d);
    }
  }
}

TEST_F(ECPBasisTest, SuperheavyElements) {
  basis.addECP_from_file(118, {0.0, 0.0, 2.0}, "testbasis_og.ecp");
  EXPECT_EQ(basis.getN(), 2);
  EXPECT_EQ(basis.getECPCore(118), 92);
  EXPECT_EQ(basis.getECP(1).getL(), 4);
  EXPECT_EQ(basis.getECP(1).getN(), 43);
  EXPECT_THROW(basis.addECP_from_file(119, {0.0, 0.0, 0.0}, "testbasis_og.ecp"),
               std::runtime_error);
}

TEST_F(ECPBasisTest, MissingFileOrAtom) {
  EXPECT_THROW(basis.addECP_from_file(14, {0.0, 0.0, 0.0}, "nonexistent.ecp"), std::runtime_error);
  EXPECT_THROW(basis.addECP_from_file(1, {0.0, 0.0, 0.0}, "testbasis.ecp"), std::runtime_error);
  EXPECT_THROW(basis.addECP_from_file(14, {0.0, 0.0, 0.0}, "testbasis_highl.ecp"),
               std::runtime_error);
  EXPECT_EQ(basis.getN(), 1);
}
