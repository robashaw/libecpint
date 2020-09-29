#include "gtest/gtest.h"
#include "multiarr.hpp"

using namespace libecpint;

class ThreeIndexTest : public testing::Test {
protected:
	ThreeIndex<double> matrix_;	
	
	void SetUp() override {
		matrix_ = ThreeIndex<double>(3, 2, 5);
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 2; j++)
				for (int k = 0; k < 5; k++)
					matrix_(i, j, k) = 2.5;
	}
};

TEST_F(ThreeIndexTest, Accessor) {
	EXPECT_DOUBLE_EQ(matrix_(1, 1, 1), 2.5);
	matrix_(1, 1, 1) = 5.4;
	EXPECT_DOUBLE_EQ(matrix_(1, 1, 1), 5.4);
	EXPECT_DOUBLE_EQ(matrix_(0, 0, 4), 2.5);
}

TEST_F(ThreeIndexTest, Fill) {
	EXPECT_DOUBLE_EQ(matrix_(0, 0, 0), 2.5);
	matrix_.fill(-3.91);
	EXPECT_DOUBLE_EQ(matrix_(0, 0, 0), -3.91);
}

TEST_F(ThreeIndexTest, CopyCtor) {
	ThreeIndex<double> tmp; 
	EXPECT_EQ(tmp.dims[0], 0);
	EXPECT_EQ(tmp.dims[1], 0);
	EXPECT_EQ(tmp.dims[2], 0);
	
	tmp = matrix_;
	EXPECT_EQ(tmp.dims[0], 3);
	EXPECT_EQ(tmp.dims[1], 2);
	EXPECT_EQ(tmp.dims[2], 5);
	
	tmp(2, 1, 1) = 4.4;
	EXPECT_DOUBLE_EQ(tmp(2, 1, 1), 4.4);
	EXPECT_DOUBLE_EQ(matrix_(0, 0, 3), 2.5);
}