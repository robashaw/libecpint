// Generated as part of Libecpint, Copyright 2017 Robert A Shaw
#include "qgen.hpp"
namespace libecpint {
namespace qgen {
void Q0_4_3(ECP& U, GaussianShell& shellA, GaussianShell& shellB, FiveIndex<double> &CA, FiveIndex<double> &CB, TwoIndex<double> &SA, TwoIndex<double> &SB, double Am, double Bm, RadialIntegral &radint, AngularIntegral& angint, ThreeIndex<double> &values) {

	std::vector<Triple> radial_triples_A = {
		{0, 3, 3},
		{1, 3, 4},
		{2, 3, 3},
		{2, 3, 5},
		{3, 3, 4},
		{3, 3, 6},
		{4, 3, 3},
		{4, 3, 5},
		{4, 3, 7}	};

	ThreeIndex<double> radials(8, 4, 8);
	radint.type2(radial_triples_A, 6, 3, U, shellA, shellB, Am, Bm, radials);

	std::vector<Triple> radial_triples_B = {
		{1, 2, 3},
		{2, 1, 3},
		{3, 0, 3},
		{3, 2, 3},
		{4, 1, 3}	};

	ThreeIndex<double> radials_B(8, 8, 4);
	radint.type2(radial_triples_B, 6, 3, U, shellB, shellA, Bm, Am, radials_B);
	for (Triple& t : radial_triples_B) radials(std::get<0>(t), std::get<2>(t), std::get<1>(t)) = radials_B(std::get<0>(t), std::get<1>(t), std::get<2>(t));

	rolled_up(3, 0, 4, radials, CA, CB, SA, SB, angint, values);
}
}
}
