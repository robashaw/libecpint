// Generated as part of Libecpint, Copyright 2017 Robert A Shaw
#include "qgen.hpp"
namespace libecpint {
namespace qgen {
void Q2_4_4(ECP& U, GaussianShell& shellA, GaussianShell& shellB, FiveIndex<double> &CA, FiveIndex<double> &CB, TwoIndex<double> &SA, TwoIndex<double> &SB, double Am, double Bm, RadialIntegral &radint, AngularIntegral& angint, ThreeIndex<double> &values) {

	std::vector<Triple> radial_triples_A = {
		{0, 4, 4},
		{1, 3, 4},
		{1, 4, 5},
		{2, 2, 4},
		{2, 3, 3},
		{2, 3, 5},
		{2, 4, 4},
		{2, 4, 6},
		{2, 5, 5},
		{3, 2, 3},
		{3, 2, 5},
		{3, 3, 4},
		{3, 3, 6},
		{3, 4, 5},
		{3, 4, 7},
		{3, 5, 6},
		{4, 2, 2},
		{4, 2, 4},
		{4, 2, 6},
		{4, 3, 3},
		{4, 3, 5},
		{4, 3, 7},
		{4, 4, 4},
		{4, 4, 6},
		{4, 4, 8},
		{4, 5, 5},
		{4, 5, 7},
		{4, 6, 6},
		{5, 2, 3},
		{5, 2, 5},
		{5, 2, 7},
		{5, 3, 4},
		{5, 3, 6},
		{5, 3, 8},
		{5, 4, 5},
		{5, 4, 7},
		{5, 5, 6},
		{5, 5, 8},
		{5, 6, 7},
		{6, 2, 2},
		{6, 2, 4},
		{6, 2, 6},
		{6, 2, 8},
		{6, 4, 4},
		{6, 4, 6},
		{6, 4, 8},
		{6, 6, 6},
		{6, 6, 8}	};

	ThreeIndex<double> radials(11, 7, 9);
	radint.type2(radial_triples_A, 11, 4, U, shellA, shellB, Am, Bm, radials);

	std::vector<Triple> radial_triples_B = {
		{1, 3, 4},
		{1, 4, 5},
		{2, 2, 4},
		{2, 3, 5},
		{2, 4, 6},
		{3, 2, 3},
		{3, 1, 4},
		{3, 3, 4},
		{3, 2, 5},
		{3, 4, 5},
		{3, 3, 6},
		{3, 5, 6},
		{4, 1, 3},
		{4, 0, 4},
		{4, 2, 4},
		{4, 1, 5},
		{4, 3, 5},
		{4, 2, 6},
		{4, 4, 6},
		{5, 1, 2},
		{5, 0, 3},
		{5, 2, 3},
		{5, 1, 4},
		{5, 3, 4},
		{5, 0, 5},
		{5, 2, 5},
		{5, 4, 5},
		{5, 1, 6},
		{5, 3, 6},
		{5, 5, 6},
		{6, 0, 2},
		{6, 0, 4},
		{6, 2, 4},
		{6, 0, 6},
		{6, 2, 6},
		{6, 4, 6}	};

	ThreeIndex<double> radials_B(11, 9, 7);
	radint.type2(radial_triples_B, 11, 4, U, shellB, shellA, Bm, Am, radials_B);
	for (Triple& t : radial_triples_B) radials(std::get<0>(t), std::get<2>(t), std::get<1>(t)) = radials_B(std::get<0>(t), std::get<1>(t), std::get<2>(t));

	rolled_up(4, 2, 4, radials, CA, CB, SA, SB, angint, values);
}
}
}
