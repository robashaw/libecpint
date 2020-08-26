#include "gtest/gtest.h"
#include "gshell.hpp"
#include <iostream>

using namespace libecpint;

class GShellTest : public testing::Test {
protected:
	std::vector<GaussianShell> gvec; 
	double C[3] = {0.0, 0.5, -0.5};
	
	void SetUp() override {
		GaussianShell g(C, 2);
		g.addPrim(0.5, 0.5);
		g.addPrim(1.4, -0.021);
		g.addPrim(10.1, 0.011);
		gvec.push_back(g);
	}
};

TEST_F(GShellTest, NPrimitive) {
	EXPECT_EQ(gvec[0].nprimitive(), 3);
	gvec[0].addPrim(0.01, 1.4);
	EXPECT_EQ(gvec[0].nprimitive(), 4);
}

TEST_F(GShellTest, NCartesian) {
	EXPECT_EQ(gvec[0].ncartesian(), 6);
	gvec[0].l = 4;
	EXPECT_EQ(gvec[0].ncartesian(), 15);
}

TEST_F(GShellTest, Center) {
	double* C = gvec[0].center();
	EXPECT_DOUBLE_EQ(C[0], 0.0);
	EXPECT_DOUBLE_EQ(C[1], 0.5);
	EXPECT_DOUBLE_EQ(C[2], -0.5);
}

TEST_F(GShellTest, ExpCoef) {
	EXPECT_DOUBLE_EQ(gvec[0].exp(0), 0.5);
	EXPECT_DOUBLE_EQ(gvec[0].coef(2), 0.011);
}

TEST_F(GShellTest, L) {
	EXPECT_EQ(gvec[0].am(), 2);
	gvec[0].l = 4;
	EXPECT_EQ(gvec[0].am(), 4);
}

TEST_F(GShellTest, Copy) {
	GaussianShell g2 = gvec[0].copy();
	EXPECT_EQ(g2.l, 2);
	g2.l = 4;
	EXPECT_EQ(gvec[0].l, 2);
	EXPECT_EQ(g2.nprimitive(), 3);
	double* C = g2.center();
	EXPECT_DOUBLE_EQ(C[1], 0.5);
}

int main(int argc, char** argv)
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

