// Generated as part of Libecpint, Copyright 2017 Robert A Shaw
#include "qgen.hpp"
namespace libecpint {
namespace qgen {
void Q0_5_4(ECP& U, GaussianShell& shellA, GaussianShell& shellB, FiveIndex<double> &CA, FiveIndex<double> &CB, TwoIndex<double> &SA, TwoIndex<double> &SB, double Am, double Bm, RadialIntegral &radint, AngularIntegral& angint, ThreeIndex<double> &values) {

	std::vector<Triple> radial_triples_A = {
		{0, 4, 4},
		{1, 4, 5},
		{2, 4, 4},
		{2, 4, 6},
		{3, 4, 5},
		{3, 4, 7},
		{4, 4, 4},
		{4, 4, 6},
		{4, 4, 8},
		{5, 4, 5},
		{5, 4, 7},
		{5, 4, 9}	};

	ThreeIndex<double> radials(10, 5, 10);
	radint.type2(radial_triples_A, 8, 4, U, shellA, shellB, Am, Bm, radials);

	std::vector<Triple> radial_triples_B = {
		{1, 3, 4},
		{2, 2, 4},
		{3, 1, 4},
		{3, 3, 4},
		{4, 0, 4},
		{4, 2, 4},
		{5, 1, 4},
		{5, 3, 4}	};

	ThreeIndex<double> radials_B(10, 10, 5);
	radint.type2(radial_triples_B, 8, 4, U, shellB, shellA, Bm, Am, radials_B);
	for (Triple& t : radial_triples_B) radials(std::get<0>(t), std::get<2>(t), std::get<1>(t)) = radials_B(std::get<0>(t), std::get<1>(t), std::get<2>(t));

	rolled_up(4, 0, 5, radials, CA, CB, SA, SB, angint, values);
}
}
}
