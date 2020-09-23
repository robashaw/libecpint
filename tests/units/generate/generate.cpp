#include "gtest/gtest.h"
#include "generate.hpp"
#include <iostream>

using namespace libecpint;

class SumTermTest : public testing::Test {
protected:
	std::vector<SumTerm> terms;
	
	void SetUp() override {
		int lam = 1;
		int nb = 0;		
										
		for (int lam1 = 0; lam1 <= lam; lam1++) {
			int lam2start = (lam1) % 2; 
			for (int lam2 = lam2start; lam2 <= lam; lam2+=2) {
												
				for (int mu1 = -lam1; mu1 <= lam1; mu1++) {
					for (int mu2 = -lam2; mu2 <= lam2; mu2++) {
																																													
						for (int mu = -lam; mu <= lam; mu++) {
							SumTerm newTerm; 
							newTerm.SA = Pair(lam1, lam1+mu1); 
							newTerm.SB = Pair(lam2, lam2+mu2);
							newTerm.radial = Triple(0, lam1, lam2);
							newTerm.CA = Quintuple(0, 0, 0, 0, 0); 
							newTerm.CB = Quintuple(0, 0, 0, 0, 0); 
							newTerm.ang = 1.0;
							newTerm.mu = lam+mu; 
							newTerm.na = 0;
							newTerm.nb = 0;
																
							terms.push_back(newTerm); 

						} 
					}
																
				}
			}
		}							
	}					
};

TEST_F(SumTermTest, Operators) {
	EXPECT_TRUE(terms[0] < terms[1]);
	EXPECT_FALSE(terms[2] < terms[3]);
	EXPECT_FALSE(terms[3] < terms[6]);
	EXPECT_TRUE(terms[3] <= terms[6]);
	EXPECT_TRUE(terms[3] == terms[6]);
}

TEST_F(SumTermTest, Indices) {
	EXPECT_EQ(terms[5].ca_index(), 0);
	EXPECT_EQ(terms[5].cb_index(), 0);
	
	terms[5].CA = Quintuple(5, 4, 3, 2, 1);
	terms[5].CB = Quintuple(0, 4, 3, 2, 1);
	EXPECT_EQ(terms[5].ca_index(), 612);
	EXPECT_EQ(terms[5].cb_index(), 612);
}

TEST_F(SumTermTest, Compare) {
	Heptuple result = terms[9].compare(terms[4]);
	EXPECT_TRUE(result == Heptuple(0, 1, 1, 0, 1, 1, 1));
	
	terms[4].CA = Quintuple(5, 4, 3, 2, 1);
	terms[4].ang = 4.4;
	result = terms[9].compare(terms[4]);
 	EXPECT_TRUE(result == Heptuple(0, 1, 1, 0, 0, 0, 1));
}


int main(int argc, char** argv)
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

