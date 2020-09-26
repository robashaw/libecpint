// Generated as part of Libecpint, Copyright 2017 Robert A Shaw
#include "qgen.hpp"
namespace libecpint {
namespace qgen {
void Q0_1_0(ECP& U, GaussianShell& shellA, GaussianShell& shellB, FiveIndex<double> &CA, FiveIndex<double> &CB, TwoIndex<double> &SA, TwoIndex<double> &SB, double Am, double Bm, RadialIntegral &radint, AngularIntegral& angint, ThreeIndex<double> &values) {

	std::vector<Triple> radial_triples_A = {
		{0, 0, 0},
		{1, 0, 1}	};

	ThreeIndex<double> radials(2, 1, 2);
	radint.type2(radial_triples_A, 0, 0, U, shellA, shellB, Am, Bm, radials);

	std::vector<Triple> radial_triples_B = {
	};

	ThreeIndex<double> radials_B(2, 2, 1);
	radint.type2(radial_triples_B, 0, 0, U, shellB, shellA, Bm, Am, radials_B);
	for (Triple& t : radial_triples_B) radials(std::get<0>(t), std::get<2>(t), std::get<1>(t)) = radials_B(std::get<0>(t), std::get<1>(t), std::get<2>(t));

	values(0, 0, 0) += 157.914 * CA(0, 0, 0, 0, 0) * CB(0, 0, 0, 0, 0) * radials(0, 0, 0) * SA(0, 0) * SB(0, 0);
	values(0, 0, 0) += 91.1715 * CA(0, 0, 0, 0, 0) * CB(0, 0, 1, 0, 0) * radials(1, 0, 1) * SA(0, 0) * SB(1, 2);
	values(0, 1, 0) += 157.914 * CA(0, 0, 0, 0, 0) * CB(0, 1, 0, 0, 0) * radials(0, 0, 0) * SA(0, 0) * SB(0, 0);
	values(0, 1, 0) += 91.1715 * CA(0, 0, 0, 0, 0) * CB(0, 1, 0, 1, 0) * radials(1, 0, 1) * SA(0, 0) * SB(1, 0);
	values(0, 2, 0) += 157.914 * CA(0, 0, 0, 0, 0) * CB(0, 2, 0, 0, 0) * radials(0, 0, 0) * SA(0, 0) * SB(0, 0);
	values(0, 2, 0) += 91.1715 * CA(0, 0, 0, 0, 0) * CB(0, 2, 0, 0, 1) * radials(1, 0, 1) * SA(0, 0) * SB(1, 1);
}
}
}
