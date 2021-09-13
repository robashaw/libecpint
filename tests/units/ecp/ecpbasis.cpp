#include "ecp.hpp"
#include "gtest/gtest.h"
#include <iostream>

using namespace libecpint;
#ifdef HAS_PUGIXML

class ECPBasisTest : public testing::Test {
 protected:
  ECPBasis basis;

  void SetUp() override {
    // load in ECP for silicon
    basis.addECP_from_file(14, {0.0, -0.5, 0.5}, "testbasis.xml");
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
  basis.addECP_from_file(33, {0.0, -0.5, 2.4}, "testbasis.xml");
  EXPECT_EQ(basis.getMaxL(), 3);
  EXPECT_EQ(basis.getN(), 2);
  EXPECT_EQ(basis.getECPCore(33), 28);
}

#endif
