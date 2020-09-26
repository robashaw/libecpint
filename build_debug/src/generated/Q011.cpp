// Generated as part of Libecpint, Copyright 2017 Robert A Shaw
#include "qgen.hpp"
namespace libecpint {
namespace qgen {
void Q0_1_1(ECP& U, GaussianShell& shellA, GaussianShell& shellB, FiveIndex<double> &CA, FiveIndex<double> &CB, TwoIndex<double> &SA, TwoIndex<double> &SB, double Am, double Bm, RadialIntegral &radint, AngularIntegral& angint, ThreeIndex<double> &values) {

	std::vector<Triple> radial_triples_A = {
		{0, 1, 1},
		{1, 1, 2}	};

	ThreeIndex<double> radials(3, 2, 3);
	radint.type2(radial_triples_A, 1, 1, U, shellA, shellB, Am, Bm, radials);

	std::vector<Triple> radial_triples_B = {
		{1, 0, 1}	};

	ThreeIndex<double> radials_B(3, 3, 2);
	radint.type2(radial_triples_B, 1, 1, U, shellB, shellA, Bm, Am, radials_B);
	for (Triple& t : radial_triples_B) radials(std::get<0>(t), std::get<2>(t), std::get<1>(t)) = radials_B(std::get<0>(t), std::get<1>(t), std::get<2>(t));

	values(0, 0, 0) += 157.914 * CA(0, 0, 0, 0, 0) * CB(0, 0, 0, 0, 0) * radials(0, 1, 1) * SA(1, 0) * SB(1, 0);
	values(0, 0, 1) += 157.914 * CA(0, 0, 0, 0, 0) * CB(0, 0, 0, 0, 0) * radials(0, 1, 1) * SA(1, 1) * SB(1, 1);
	values(0, 0, 2) += 157.914 * CA(0, 0, 0, 0, 0) * CB(0, 0, 0, 0, 0) * radials(0, 1, 1) * SA(1, 2) * SB(1, 2);
	values(0, 0, 2) += 91.1715 * CA(0, 0, 0, 0, 0) * CB(0, 0, 1, 0, 0) * radials(1, 1, 0) * SA(1, 2) * SB(0, 0);
	values(0, 0, 0) += 70.6211 * CA(0, 0, 0, 0, 0) * CB(0, 0, 1, 0, 0) * radials(1, 1, 2) * SA(1, 0) * SB(2, 0);
	values(0, 0, 1) += 70.6211 * CA(0, 0, 0, 0, 0) * CB(0, 0, 1, 0, 0) * radials(1, 1, 2) * SA(1, 1) * SB(2, 3);
	values(0, 0, 2) += -40.7731 * CA(0, 0, 0, 0, 0) * CB(0, 0, 1, 0, 0) * radials(1, 1, 2) * SA(1, 2) * SB(2, 2);
	values(0, 0, 2) += 70.6211 * CA(0, 0, 0, 0, 0) * CB(0, 0, 1, 0, 0) * radials(1, 1, 2) * SA(1, 2) * SB(2, 4);
	values(0, 1, 0) += 157.914 * CA(0, 0, 0, 0, 0) * CB(0, 1, 0, 0, 0) * radials(0, 1, 1) * SA(1, 0) * SB(1, 0);
	values(0, 1, 1) += 157.914 * CA(0, 0, 0, 0, 0) * CB(0, 1, 0, 0, 0) * radials(0, 1, 1) * SA(1, 1) * SB(1, 1);
	values(0, 1, 2) += 157.914 * CA(0, 0, 0, 0, 0) * CB(0, 1, 0, 0, 0) * radials(0, 1, 1) * SA(1, 2) * SB(1, 2);
	values(0, 1, 0) += 91.1715 * CA(0, 0, 0, 0, 0) * CB(0, 1, 0, 1, 0) * radials(1, 1, 0) * SA(1, 0) * SB(0, 0);
	values(0, 1, 0) += -40.7731 * CA(0, 0, 0, 0, 0) * CB(0, 1, 0, 1, 0) * radials(1, 1, 2) * SA(1, 0) * SB(2, 2);
	values(0, 1, 0) += -70.6211 * CA(0, 0, 0, 0, 0) * CB(0, 1, 0, 1, 0) * radials(1, 1, 2) * SA(1, 0) * SB(2, 4);
	values(0, 1, 1) += 70.6211 * CA(0, 0, 0, 0, 0) * CB(0, 1, 0, 1, 0) * radials(1, 1, 2) * SA(1, 1) * SB(2, 1);
	values(0, 1, 2) += 70.6211 * CA(0, 0, 0, 0, 0) * CB(0, 1, 0, 1, 0) * radials(1, 1, 2) * SA(1, 2) * SB(2, 0);
	values(0, 2, 0) += 157.914 * CA(0, 0, 0, 0, 0) * CB(0, 2, 0, 0, 0) * radials(0, 1, 1) * SA(1, 0) * SB(1, 0);
	values(0, 2, 1) += 157.914 * CA(0, 0, 0, 0, 0) * CB(0, 2, 0, 0, 0) * radials(0, 1, 1) * SA(1, 1) * SB(1, 1);
	values(0, 2, 2) += 157.914 * CA(0, 0, 0, 0, 0) * CB(0, 2, 0, 0, 0) * radials(0, 1, 1) * SA(1, 2) * SB(1, 2);
	values(0, 2, 1) += 91.1715 * CA(0, 0, 0, 0, 0) * CB(0, 2, 0, 0, 1) * radials(1, 1, 0) * SA(1, 1) * SB(0, 0);
	values(0, 2, 0) += 70.6211 * CA(0, 0, 0, 0, 0) * CB(0, 2, 0, 0, 1) * radials(1, 1, 2) * SA(1, 0) * SB(2, 1);
	values(0, 2, 1) += 81.5463 * CA(0, 0, 0, 0, 0) * CB(0, 2, 0, 0, 1) * radials(1, 1, 2) * SA(1, 1) * SB(2, 2);
	values(0, 2, 2) += 70.6211 * CA(0, 0, 0, 0, 0) * CB(0, 2, 0, 0, 1) * radials(1, 1, 2) * SA(1, 2) * SB(2, 3);
}
}
}
