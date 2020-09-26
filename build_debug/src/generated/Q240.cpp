// Generated as part of Libecpint, Copyright 2017 Robert A Shaw
#include "qgen.hpp"
namespace libecpint {
namespace qgen {
void Q2_4_0(ECP& U, GaussianShell& shellA, GaussianShell& shellB, FiveIndex<double> &CA, FiveIndex<double> &CB, TwoIndex<double> &SA, TwoIndex<double> &SB, double Am, double Bm, RadialIntegral &radint, AngularIntegral& angint, ThreeIndex<double> &values) {

	std::vector<Triple> radial_triples_A = {
		{0, 0, 0},
		{1, 0, 1},
		{2, 0, 0},
		{2, 0, 2},
		{2, 1, 1},
		{3, 0, 1},
		{3, 0, 3},
		{3, 1, 2},
		{4, 0, 0},
		{4, 0, 2},
		{4, 0, 4},
		{4, 1, 1},
		{4, 1, 3},
		{4, 2, 2},
		{5, 0, 1},
		{5, 0, 3},
		{5, 1, 2},
		{5, 1, 4},
		{5, 2, 3},
		{6, 0, 0},
		{6, 0, 2},
		{6, 0, 4},
		{6, 2, 2},
		{6, 2, 4}	};

	ThreeIndex<double> radials(7, 3, 5);
	radint.type2(radial_triples_A, 7, 0, U, shellA, shellB, Am, Bm, radials);

	std::vector<Triple> radial_triples_B = {
		{1, 0, 1},
		{2, 0, 2},
		{3, 0, 1},
		{3, 1, 2},
		{4, 0, 2},
		{5, 0, 1},
		{5, 1, 2},
		{6, 0, 2}	};

	ThreeIndex<double> radials_B(7, 5, 3);
	radint.type2(radial_triples_B, 7, 0, U, shellB, shellA, Bm, Am, radials_B);
	for (Triple& t : radial_triples_B) radials(std::get<0>(t), std::get<2>(t), std::get<1>(t)) = radials_B(std::get<0>(t), std::get<1>(t), std::get<2>(t));

	rolled_up(0, 2, 4, radials, CA, CB, SA, SB, angint, values);
}
}
}
