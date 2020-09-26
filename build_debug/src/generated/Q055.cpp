// Generated as part of Libecpint, Copyright 2017 Robert A Shaw
#include "qgen.hpp"
namespace libecpint {
namespace qgen {
void Q0_5_5(ECP& U, GaussianShell& shellA, GaussianShell& shellB, FiveIndex<double> &CA, FiveIndex<double> &CB, TwoIndex<double> &SA, TwoIndex<double> &SB, double Am, double Bm, RadialIntegral &radint, AngularIntegral& angint, ThreeIndex<double> &values) {

	std::vector<Triple> radial_triples_A = {
		{0, 5, 5},
		{1, 5, 6},
		{2, 5, 5},
		{2, 5, 7},
		{3, 5, 6},
		{3, 5, 8},
		{4, 5, 5},
		{4, 5, 7},
		{4, 5, 9},
		{5, 5, 6},
		{5, 5, 8},
		{5, 5, 10}	};

	ThreeIndex<double> radials(11, 6, 11);
	radint.type2(radial_triples_A, 9, 5, U, shellA, shellB, Am, Bm, radials);

	std::vector<Triple> radial_triples_B = {
		{1, 4, 5},
		{2, 3, 5},
		{3, 2, 5},
		{3, 4, 5},
		{4, 1, 5},
		{4, 3, 5},
		{5, 0, 5},
		{5, 2, 5},
		{5, 4, 5}	};

	ThreeIndex<double> radials_B(11, 11, 6);
	radint.type2(radial_triples_B, 9, 5, U, shellB, shellA, Bm, Am, radials_B);
	for (Triple& t : radial_triples_B) radials(std::get<0>(t), std::get<2>(t), std::get<1>(t)) = radials_B(std::get<0>(t), std::get<1>(t), std::get<2>(t));

	rolled_up(5, 0, 5, radials, CA, CB, SA, SB, angint, values);
}
}
}
