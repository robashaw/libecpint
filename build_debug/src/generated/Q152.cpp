// Generated as part of Libecpint, Copyright 2017 Robert A Shaw
#include "qgen.hpp"
namespace libecpint {
namespace qgen {
void Q1_5_2(ECP& U, GaussianShell& shellA, GaussianShell& shellB, FiveIndex<double> &CA, FiveIndex<double> &CB, TwoIndex<double> &SA, TwoIndex<double> &SB, double Am, double Bm, RadialIntegral &radint, AngularIntegral& angint, ThreeIndex<double> &values) {

	std::vector<Triple> radial_triples_A = {
		{0, 2, 2},
		{1, 1, 2},
		{1, 2, 3},
		{2, 1, 1},
		{2, 1, 3},
		{2, 2, 2},
		{2, 2, 4},
		{2, 3, 3},
		{3, 1, 2},
		{3, 1, 4},
		{3, 2, 3},
		{3, 2, 5},
		{3, 3, 4},
		{4, 1, 1},
		{4, 1, 3},
		{4, 1, 5},
		{4, 2, 2},
		{4, 2, 4},
		{4, 2, 6},
		{4, 3, 3},
		{4, 3, 5},
		{5, 1, 2},
		{5, 1, 4},
		{5, 1, 6},
		{5, 2, 3},
		{5, 2, 5},
		{5, 2, 7},
		{5, 3, 4},
		{5, 3, 6},
		{6, 1, 1},
		{6, 1, 3},
		{6, 1, 5},
		{6, 1, 7},
		{6, 3, 3},
		{6, 3, 5},
		{6, 3, 7}	};

	ThreeIndex<double> radials(9, 4, 8);
	radint.type2(radial_triples_A, 8, 2, U, shellA, shellB, Am, Bm, radials);

	std::vector<Triple> radial_triples_B = {
		{1, 1, 2},
		{1, 2, 3},
		{2, 0, 2},
		{2, 1, 3},
		{3, 0, 1},
		{3, 1, 2},
		{3, 0, 3},
		{3, 2, 3},
		{4, 0, 2},
		{4, 1, 3},
		{5, 0, 1},
		{5, 1, 2},
		{5, 0, 3},
		{5, 2, 3},
		{6, 1, 3}	};

	ThreeIndex<double> radials_B(9, 8, 4);
	radint.type2(radial_triples_B, 8, 2, U, shellB, shellA, Bm, Am, radials_B);
	for (Triple& t : radial_triples_B) radials(std::get<0>(t), std::get<2>(t), std::get<1>(t)) = radials_B(std::get<0>(t), std::get<1>(t), std::get<2>(t));

	rolled_up(2, 1, 5, radials, CA, CB, SA, SB, angint, values);
}
}
}
