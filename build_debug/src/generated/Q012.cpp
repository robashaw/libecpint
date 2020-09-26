// Generated as part of Libecpint, Copyright 2017 Robert A Shaw
#include "qgen.hpp"
namespace libecpint {
namespace qgen {
void Q0_1_2(ECP& U, GaussianShell& shellA, GaussianShell& shellB, FiveIndex<double> &CA, FiveIndex<double> &CB, TwoIndex<double> &SA, TwoIndex<double> &SB, double Am, double Bm, RadialIntegral &radint, AngularIntegral& angint, ThreeIndex<double> &values) {

	std::vector<Triple> radial_triples_A = {
		{0, 2, 2},
		{1, 2, 3}	};

	ThreeIndex<double> radials(4, 3, 4);
	radint.type2(radial_triples_A, 2, 2, U, shellA, shellB, Am, Bm, radials);

	std::vector<Triple> radial_triples_B = {
		{1, 1, 2}	};

	ThreeIndex<double> radials_B(4, 4, 3);
	radint.type2(radial_triples_B, 2, 2, U, shellB, shellA, Bm, Am, radials_B);
	for (Triple& t : radial_triples_B) radials(std::get<0>(t), std::get<2>(t), std::get<1>(t)) = radials_B(std::get<0>(t), std::get<1>(t), std::get<2>(t));

	values(0, 0, 0) += 157.914 * CA(0, 0, 0, 0, 0) * CB(0, 0, 0, 0, 0) * radials(0, 2, 2) * SA(2, 0) * SB(2, 0);
	values(0, 0, 1) += 157.914 * CA(0, 0, 0, 0, 0) * CB(0, 0, 0, 0, 0) * radials(0, 2, 2) * SA(2, 1) * SB(2, 1);
	values(0, 0, 2) += 157.914 * CA(0, 0, 0, 0, 0) * CB(0, 0, 0, 0, 0) * radials(0, 2, 2) * SA(2, 2) * SB(2, 2);
	values(0, 0, 3) += 157.914 * CA(0, 0, 0, 0, 0) * CB(0, 0, 0, 0, 0) * radials(0, 2, 2) * SA(2, 3) * SB(2, 3);
	values(0, 0, 4) += 157.914 * CA(0, 0, 0, 0, 0) * CB(0, 0, 0, 0, 0) * radials(0, 2, 2) * SA(2, 4) * SB(2, 4);
	values(0, 0, 0) += 70.6211 * CA(0, 0, 0, 0, 0) * CB(0, 0, 1, 0, 0) * radials(1, 2, 1) * SA(2, 0) * SB(1, 0);
	values(0, 0, 2) += -40.7731 * CA(0, 0, 0, 0, 0) * CB(0, 0, 1, 0, 0) * radials(1, 2, 1) * SA(2, 2) * SB(1, 2);
	values(0, 0, 3) += 70.6211 * CA(0, 0, 0, 0, 0) * CB(0, 0, 1, 0, 0) * radials(1, 2, 1) * SA(2, 3) * SB(1, 1);
	values(0, 0, 4) += 70.6211 * CA(0, 0, 0, 0, 0) * CB(0, 0, 1, 0, 0) * radials(1, 2, 1) * SA(2, 4) * SB(1, 2);
	values(0, 0, 0) += 73.0998 * CA(0, 0, 0, 0, 0) * CB(0, 0, 1, 0, 0) * radials(1, 2, 3) * SA(2, 0) * SB(3, 0);
	values(0, 0, 0) += -18.8743 * CA(0, 0, 0, 0, 0) * CB(0, 0, 1, 0, 0) * radials(1, 2, 3) * SA(2, 0) * SB(3, 2);
	values(0, 0, 1) += 59.6858 * CA(0, 0, 0, 0, 0) * CB(0, 0, 1, 0, 0) * radials(1, 2, 3) * SA(2, 1) * SB(3, 1);
	values(0, 0, 2) += 65.3825 * CA(0, 0, 0, 0, 0) * CB(0, 0, 1, 0, 0) * radials(1, 2, 3) * SA(2, 2) * SB(3, 4);
	values(0, 0, 3) += -46.2324 * CA(0, 0, 0, 0, 0) * CB(0, 0, 1, 0, 0) * radials(1, 2, 3) * SA(2, 3) * SB(3, 3);
	values(0, 0, 3) += 59.6858 * CA(0, 0, 0, 0, 0) * CB(0, 0, 1, 0, 0) * radials(1, 2, 3) * SA(2, 3) * SB(3, 5);
	values(0, 0, 4) += -18.8743 * CA(0, 0, 0, 0, 0) * CB(0, 0, 1, 0, 0) * radials(1, 2, 3) * SA(2, 4) * SB(3, 4);
	values(0, 0, 4) += 73.0998 * CA(0, 0, 0, 0, 0) * CB(0, 0, 1, 0, 0) * radials(1, 2, 3) * SA(2, 4) * SB(3, 6);
	values(0, 1, 0) += 157.914 * CA(0, 0, 0, 0, 0) * CB(0, 1, 0, 0, 0) * radials(0, 2, 2) * SA(2, 0) * SB(2, 0);
	values(0, 1, 1) += 157.914 * CA(0, 0, 0, 0, 0) * CB(0, 1, 0, 0, 0) * radials(0, 2, 2) * SA(2, 1) * SB(2, 1);
	values(0, 1, 2) += 157.914 * CA(0, 0, 0, 0, 0) * CB(0, 1, 0, 0, 0) * radials(0, 2, 2) * SA(2, 2) * SB(2, 2);
	values(0, 1, 3) += 157.914 * CA(0, 0, 0, 0, 0) * CB(0, 1, 0, 0, 0) * radials(0, 2, 2) * SA(2, 3) * SB(2, 3);
	values(0, 1, 4) += 157.914 * CA(0, 0, 0, 0, 0) * CB(0, 1, 0, 0, 0) * radials(0, 2, 2) * SA(2, 4) * SB(2, 4);
	values(0, 1, 0) += 70.6211 * CA(0, 0, 0, 0, 0) * CB(0, 1, 0, 1, 0) * radials(1, 2, 1) * SA(2, 0) * SB(1, 2);
	values(0, 1, 1) += 70.6211 * CA(0, 0, 0, 0, 0) * CB(0, 1, 0, 1, 0) * radials(1, 2, 1) * SA(2, 1) * SB(1, 1);
	values(0, 1, 2) += -40.7731 * CA(0, 0, 0, 0, 0) * CB(0, 1, 0, 1, 0) * radials(1, 2, 1) * SA(2, 2) * SB(1, 0);
	values(0, 1, 4) += -70.6211 * CA(0, 0, 0, 0, 0) * CB(0, 1, 0, 1, 0) * radials(1, 2, 1) * SA(2, 4) * SB(1, 0);
	values(0, 1, 0) += -18.8743 * CA(0, 0, 0, 0, 0) * CB(0, 1, 0, 1, 0) * radials(1, 2, 3) * SA(2, 0) * SB(3, 4);
	values(0, 1, 0) += -73.0998 * CA(0, 0, 0, 0, 0) * CB(0, 1, 0, 1, 0) * radials(1, 2, 3) * SA(2, 0) * SB(3, 6);
	values(0, 1, 1) += -46.2324 * CA(0, 0, 0, 0, 0) * CB(0, 1, 0, 1, 0) * radials(1, 2, 3) * SA(2, 1) * SB(3, 3);
	values(0, 1, 1) += -59.6858 * CA(0, 0, 0, 0, 0) * CB(0, 1, 0, 1, 0) * radials(1, 2, 3) * SA(2, 1) * SB(3, 5);
	values(0, 1, 2) += 65.3825 * CA(0, 0, 0, 0, 0) * CB(0, 1, 0, 1, 0) * radials(1, 2, 3) * SA(2, 2) * SB(3, 2);
	values(0, 1, 3) += 59.6858 * CA(0, 0, 0, 0, 0) * CB(0, 1, 0, 1, 0) * radials(1, 2, 3) * SA(2, 3) * SB(3, 1);
	values(0, 1, 4) += 73.0998 * CA(0, 0, 0, 0, 0) * CB(0, 1, 0, 1, 0) * radials(1, 2, 3) * SA(2, 4) * SB(3, 0);
	values(0, 1, 4) += 18.8743 * CA(0, 0, 0, 0, 0) * CB(0, 1, 0, 1, 0) * radials(1, 2, 3) * SA(2, 4) * SB(3, 2);
	values(0, 2, 0) += 157.914 * CA(0, 0, 0, 0, 0) * CB(0, 2, 0, 0, 0) * radials(0, 2, 2) * SA(2, 0) * SB(2, 0);
	values(0, 2, 1) += 157.914 * CA(0, 0, 0, 0, 0) * CB(0, 2, 0, 0, 0) * radials(0, 2, 2) * SA(2, 1) * SB(2, 1);
	values(0, 2, 2) += 157.914 * CA(0, 0, 0, 0, 0) * CB(0, 2, 0, 0, 0) * radials(0, 2, 2) * SA(2, 2) * SB(2, 2);
	values(0, 2, 3) += 157.914 * CA(0, 0, 0, 0, 0) * CB(0, 2, 0, 0, 0) * radials(0, 2, 2) * SA(2, 3) * SB(2, 3);
	values(0, 2, 4) += 157.914 * CA(0, 0, 0, 0, 0) * CB(0, 2, 0, 0, 0) * radials(0, 2, 2) * SA(2, 4) * SB(2, 4);
	values(0, 2, 1) += 70.6211 * CA(0, 0, 0, 0, 0) * CB(0, 2, 0, 0, 1) * radials(1, 2, 1) * SA(2, 1) * SB(1, 0);
	values(0, 2, 2) += 81.5463 * CA(0, 0, 0, 0, 0) * CB(0, 2, 0, 0, 1) * radials(1, 2, 1) * SA(2, 2) * SB(1, 1);
	values(0, 2, 3) += 70.6211 * CA(0, 0, 0, 0, 0) * CB(0, 2, 0, 0, 1) * radials(1, 2, 1) * SA(2, 3) * SB(1, 2);
	values(0, 2, 0) += 59.6858 * CA(0, 0, 0, 0, 0) * CB(0, 2, 0, 0, 1) * radials(1, 2, 3) * SA(2, 0) * SB(3, 1);
	values(0, 2, 1) += 75.4972 * CA(0, 0, 0, 0, 0) * CB(0, 2, 0, 0, 1) * radials(1, 2, 3) * SA(2, 1) * SB(3, 2);
	values(0, 2, 2) += 80.0768 * CA(0, 0, 0, 0, 0) * CB(0, 2, 0, 0, 1) * radials(1, 2, 3) * SA(2, 2) * SB(3, 3);
	values(0, 2, 3) += 75.4972 * CA(0, 0, 0, 0, 0) * CB(0, 2, 0, 0, 1) * radials(1, 2, 3) * SA(2, 3) * SB(3, 4);
	values(0, 2, 4) += 59.6858 * CA(0, 0, 0, 0, 0) * CB(0, 2, 0, 0, 1) * radials(1, 2, 3) * SA(2, 4) * SB(3, 5);
}
}
}
