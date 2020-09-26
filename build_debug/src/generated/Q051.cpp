// Generated as part of Libecpint, Copyright 2017 Robert A Shaw
#include "qgen.hpp"
namespace libecpint {
namespace qgen {
void Q0_5_1(ECP& U, GaussianShell& shellA, GaussianShell& shellB, FiveIndex<double> &CA, FiveIndex<double> &CB, TwoIndex<double> &SA, TwoIndex<double> &SB, double Am, double Bm, RadialIntegral &radint, AngularIntegral& angint, ThreeIndex<double> &values) {

	std::vector<Triple> radial_triples_A = {
		{0, 1, 1},
		{1, 1, 2},
		{2, 1, 1},
		{2, 1, 3},
		{3, 1, 2},
		{3, 1, 4},
		{4, 1, 1},
		{4, 1, 3},
		{4, 1, 5},
		{5, 1, 2},
		{5, 1, 4},
		{5, 1, 6}	};

	ThreeIndex<double> radials(7, 2, 7);
	radint.type2(radial_triples_A, 5, 1, U, shellA, shellB, Am, Bm, radials);

	std::vector<Triple> radial_triples_B = {
		{1, 0, 1},
		{3, 0, 1},
		{5, 0, 1}	};

	ThreeIndex<double> radials_B(7, 7, 2);
	radint.type2(radial_triples_B, 5, 1, U, shellB, shellA, Bm, Am, radials_B);
	for (Triple& t : radial_triples_B) radials(std::get<0>(t), std::get<2>(t), std::get<1>(t)) = radials_B(std::get<0>(t), std::get<1>(t), std::get<2>(t));

	rolled_up(1, 0, 5, radials, CA, CB, SA, SB, angint, values);
}
}
}
