#include "gtest/gtest.h"
#include "multiarr.hpp"
#include <algorithm>

using namespace libecpint;

class FiveIndexTest : public testing::Test {
protected:
	FiveIndex<double> matrix_;	
	
	void SetUp() override {
		matrix_ = FiveIndex<double>(2, 9, 3, 7, 5);
		std::fill(matrix_.data.begin(), matrix_.data.end(), -0.7);
	}
};

TEST_F(FiveIndexTest, Accessor) {
	EXPECT_DOUBLE_EQ(matrix_(1, 1, 1, 2, 4), -0.7);
	matrix_(1, 1, 1, 2, 4) = 5.4;
	EXPECT_DOUBLE_EQ(matrix_(1, 1, 1, 2, 4), 5.4);
	EXPECT_DOUBLE_EQ(matrix_(0, 3, 1, 4, 2), -0.7);
}

TEST_F(FiveIndexTest, CopyCtor) {
	FiveIndex<double> tmp; 
	EXPECT_EQ(tmp.dims[0], 0);
	EXPECT_EQ(tmp.dims[1], 0);
	EXPECT_EQ(tmp.dims[2], 0);
	EXPECT_EQ(tmp.dims[3], 0);
	EXPECT_EQ(tmp.dims[4], 0);
	
	tmp = matrix_;
	EXPECT_EQ(tmp.dims[0], 2);
	EXPECT_EQ(tmp.dims[1], 9);
	EXPECT_EQ(tmp.dims[2], 3);
	EXPECT_EQ(tmp.dims[3], 7);
	EXPECT_EQ(tmp.dims[4], 5);
	
	tmp(0, 1, 1, 4, 3) = 4.9;
	EXPECT_DOUBLE_EQ(tmp(0, 1, 1, 4, 3), 4.9);
	EXPECT_DOUBLE_EQ(matrix_(0, 0, 2, 1, 4), -0.7);
}