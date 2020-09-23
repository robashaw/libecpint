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
		
		std::array<double, 3> Carr = {C[0], C[1], C[2]}; 
		GaussianShell g2(Carr, 4);
		gvec.push_back(g2);
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

TEST_F(GShellTest, Center1) {
	EXPECT_FALSE(gvec[0].local_ptr);
	double* C = gvec[0].center();
	EXPECT_DOUBLE_EQ(C[0], 0.0);
	EXPECT_DOUBLE_EQ(C[1], 0.5);
	EXPECT_DOUBLE_EQ(C[2], -0.5);
}

TEST_F(GShellTest, Center2) {
	EXPECT_TRUE(gvec[1].local_ptr);
	double* C = gvec[1].center();
	EXPECT_DOUBLE_EQ(C[0], 0.0);
	EXPECT_DOUBLE_EQ(C[1], 0.5);
	EXPECT_DOUBLE_EQ(C[2], -0.5);
	gvec[1].localCenter[0] = 0.9;
	EXPECT_DOUBLE_EQ(C[0], 0.9);
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

TEST_F(GShellTest, Copy1) {
	GaussianShell g2 = gvec[0].copy();
	EXPECT_EQ(g2.l, 2);
	EXPECT_FALSE(g2.local_ptr);
	g2.l = 4;
	EXPECT_EQ(gvec[0].l, 2);
	EXPECT_EQ(g2.nprimitive(), 3);
	double* C = g2.center();
	EXPECT_DOUBLE_EQ(C[1], 0.5);
}

TEST_F(GShellTest, Copy2) {
	GaussianShell g2 = gvec[1].copy();
	EXPECT_TRUE(g2.local_ptr);
	EXPECT_DOUBLE_EQ(g2.center()[0], gvec[1].center()[0]);
	g2.localCenter[1] = 0.9;
	EXPECT_DOUBLE_EQ(g2.center()[1] - gvec[1].center()[1], 0.4);
}

int main(int argc, char** argv)
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

