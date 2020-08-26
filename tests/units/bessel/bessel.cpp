#include "gtest/gtest.h"
#include "bessel.hpp"
#include "mathutil.hpp"
#include <iostream>

using namespace libecpint;

class BesselTest : public testing::Test {
protected:
	BesselFunction bessie;
	
	void SetUp() override {
		initFactorials();
		bessie.init(2, 1600, 200, 1e-15);
	}
};

TEST_F(BesselTest, CalculateBigL) {
	std::vector<double> values;
	// check it catches maxL > lMax
	bessie.calculate(SMALL, 3, values);
	EXPECT_EQ(values.size(), 3);
}

TEST_F(BesselTest, CalculateSmallZ) {
	// check negative
	std::vector<double> values;
	bessie.calculate(-1.0, 2, values);
	EXPECT_EQ(values.size(), 3);
	EXPECT_DOUBLE_EQ(values[0], 1.0);
	
	bessie.calculate(SMALL/2, 2, values);
	EXPECT_NEAR(values[0], 1.0, 1e-7);
	EXPECT_NEAR(values[1], 1.66667e-8, 1e-12);
}

TEST_F(BesselTest, CalculateBigZ) {
	std::vector<double> values;
	bessie.calculate(17.0, 1, values);
	EXPECT_EQ(values.size(), 2);
	EXPECT_NEAR(values[0], 0.0294118, 1e-7);
	EXPECT_NEAR(values[1], 0.0276817, 1e-7);
}

TEST_F(BesselTest, CalculateMidZ) {
	std::vector<double> values;
	bessie.calculate(5.0, 1, values);
	EXPECT_EQ(values.size(), 2);
	EXPECT_NEAR(values[0], 0.0999955, 1e-7);
	EXPECT_NEAR(values[1], 0.0800054, 1e-7);
}

int main(int argc, char** argv)
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
