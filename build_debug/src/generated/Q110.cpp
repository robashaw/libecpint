// Generated as part of Libecpint, Copyright 2017 Robert A Shaw
#include "qgen.hpp"
namespace libecpint {
namespace qgen {
void Q1_1_0(ECP& U, GaussianShell& shellA, GaussianShell& shellB, FiveIndex<double> &CA, FiveIndex<double> &CB, TwoIndex<double> &SA, TwoIndex<double> &SB, double Am, double Bm, RadialIntegral &radint, AngularIntegral& angint, ThreeIndex<double> &values) {

	std::vector<Triple> radial_triples_A = {
		{0, 0, 0},
		{1, 0, 1},
		{2, 1, 1}	};

	ThreeIndex<double> radials(3, 2, 2);
	radint.type2(radial_triples_A, 2, 0, U, shellA, shellB, Am, Bm, radials);

	std::vector<Triple> radial_triples_B = {
		{1, 0, 1}	};

	ThreeIndex<double> radials_B(3, 2, 2);
	radint.type2(radial_triples_B, 2, 0, U, shellB, shellA, Bm, Am, radials_B);
	for (Triple& t : radial_triples_B) radials(std::get<0>(t), std::get<2>(t), std::get<1>(t)) = radials_B(std::get<0>(t), std::get<1>(t), std::get<2>(t));

	values(0, 0, 0) += 157.914 * CA(0, 0, 0, 0, 0) * CB(0, 0, 0, 0, 0) * radials(0, 0, 0) * SA(0, 0) * SB(0, 0);
	values(0, 0, 0) += 91.1715 * CA(0, 0, 0, 0, 0) * CB(0, 0, 1, 0, 0) * radials(1, 0, 1) * SA(0, 0) * SB(1, 2);
	values(0, 0, 0) += 91.1715 * CA(0, 0, 1, 0, 0) * CB(0, 0, 0, 0, 0) * radials(1, 1, 0) * SA(1, 2) * SB(0, 0);
	values(0, 0, 0) += 52.6379 * CA(0, 0, 1, 0, 0) * CB(0, 0, 1, 0, 0) * radials(2, 1, 1) * SA(1, 2) * SB(1, 2);
	values(0, 1, 0) += 157.914 * CA(0, 0, 0, 0, 0) * CB(0, 1, 0, 0, 0) * radials(0, 0, 0) * SA(0, 0) * SB(0, 0);
	values(0, 1, 0) += 91.1715 * CA(0, 0, 0, 0, 0) * CB(0, 1, 0, 1, 0) * radials(1, 0, 1) * SA(0, 0) * SB(1, 0);
	values(0, 1, 0) += 91.1715 * CA(0, 0, 1, 0, 0) * CB(0, 1, 0, 0, 0) * radials(1, 1, 0) * SA(1, 2) * SB(0, 0);
	values(0, 1, 0) += 52.6379 * CA(0, 0, 1, 0, 0) * CB(0, 1, 0, 1, 0) * radials(2, 1, 1) * SA(1, 2) * SB(1, 0);
	values(0, 2, 0) += 157.914 * CA(0, 0, 0, 0, 0) * CB(0, 2, 0, 0, 0) * radials(0, 0, 0) * SA(0, 0) * SB(0, 0);
	values(0, 2, 0) += 91.1715 * CA(0, 0, 0, 0, 0) * CB(0, 2, 0, 0, 1) * radials(1, 0, 1) * SA(0, 0) * SB(1, 1);
	values(0, 2, 0) += 91.1715 * CA(0, 0, 1, 0, 0) * CB(0, 2, 0, 0, 0) * radials(1, 1, 0) * SA(1, 2) * SB(0, 0);
	values(0, 2, 0) += 52.6379 * CA(0, 0, 1, 0, 0) * CB(0, 2, 0, 0, 1) * radials(2, 1, 1) * SA(1, 2) * SB(1, 1);
	values(1, 0, 0) += 157.914 * CA(0, 1, 0, 0, 0) * CB(0, 0, 0, 0, 0) * radials(0, 0, 0) * SA(0, 0) * SB(0, 0);
	values(1, 0, 0) += 91.1715 * CA(0, 1, 0, 0, 0) * CB(0, 0, 1, 0, 0) * radials(1, 0, 1) * SA(0, 0) * SB(1, 2);
	values(1, 0, 0) += 91.1715 * CA(0, 1, 0, 1, 0) * CB(0, 0, 0, 0, 0) * radials(1, 1, 0) * SA(1, 0) * SB(0, 0);
	values(1, 0, 0) += 52.6379 * CA(0, 1, 0, 1, 0) * CB(0, 0, 1, 0, 0) * radials(2, 1, 1) * SA(1, 0) * SB(1, 2);
	values(1, 1, 0) += 157.914 * CA(0, 1, 0, 0, 0) * CB(0, 1, 0, 0, 0) * radials(0, 0, 0) * SA(0, 0) * SB(0, 0);
	values(1, 1, 0) += 91.1715 * CA(0, 1, 0, 0, 0) * CB(0, 1, 0, 1, 0) * radials(1, 0, 1) * SA(0, 0) * SB(1, 0);
	values(1, 1, 0) += 91.1715 * CA(0, 1, 0, 1, 0) * CB(0, 1, 0, 0, 0) * radials(1, 1, 0) * SA(1, 0) * SB(0, 0);
	values(1, 1, 0) += 52.6379 * CA(0, 1, 0, 1, 0) * CB(0, 1, 0, 1, 0) * radials(2, 1, 1) * SA(1, 0) * SB(1, 0);
	values(1, 2, 0) += 157.914 * CA(0, 1, 0, 0, 0) * CB(0, 2, 0, 0, 0) * radials(0, 0, 0) * SA(0, 0) * SB(0, 0);
	values(1, 2, 0) += 91.1715 * CA(0, 1, 0, 0, 0) * CB(0, 2, 0, 0, 1) * radials(1, 0, 1) * SA(0, 0) * SB(1, 1);
	values(1, 2, 0) += 91.1715 * CA(0, 1, 0, 1, 0) * CB(0, 2, 0, 0, 0) * radials(1, 1, 0) * SA(1, 0) * SB(0, 0);
	values(1, 2, 0) += 52.6379 * CA(0, 1, 0, 1, 0) * CB(0, 2, 0, 0, 1) * radials(2, 1, 1) * SA(1, 0) * SB(1, 1);
	values(2, 0, 0) += 157.914 * CA(0, 2, 0, 0, 0) * CB(0, 0, 0, 0, 0) * radials(0, 0, 0) * SA(0, 0) * SB(0, 0);
	values(2, 0, 0) += 91.1715 * CA(0, 2, 0, 0, 0) * CB(0, 0, 1, 0, 0) * radials(1, 0, 1) * SA(0, 0) * SB(1, 2);
	values(2, 0, 0) += 91.1715 * CA(0, 2, 0, 0, 1) * CB(0, 0, 0, 0, 0) * radials(1, 1, 0) * SA(1, 1) * SB(0, 0);
	values(2, 0, 0) += 52.6379 * CA(0, 2, 0, 0, 1) * CB(0, 0, 1, 0, 0) * radials(2, 1, 1) * SA(1, 1) * SB(1, 2);
	values(2, 1, 0) += 157.914 * CA(0, 2, 0, 0, 0) * CB(0, 1, 0, 0, 0) * radials(0, 0, 0) * SA(0, 0) * SB(0, 0);
	values(2, 1, 0) += 91.1715 * CA(0, 2, 0, 0, 0) * CB(0, 1, 0, 1, 0) * radials(1, 0, 1) * SA(0, 0) * SB(1, 0);
	values(2, 1, 0) += 91.1715 * CA(0, 2, 0, 0, 1) * CB(0, 1, 0, 0, 0) * radials(1, 1, 0) * SA(1, 1) * SB(0, 0);
	values(2, 1, 0) += 52.6379 * CA(0, 2, 0, 0, 1) * CB(0, 1, 0, 1, 0) * radials(2, 1, 1) * SA(1, 1) * SB(1, 0);
	values(2, 2, 0) += 157.914 * CA(0, 2, 0, 0, 0) * CB(0, 2, 0, 0, 0) * radials(0, 0, 0) * SA(0, 0) * SB(0, 0);
	values(2, 2, 0) += 91.1715 * CA(0, 2, 0, 0, 0) * CB(0, 2, 0, 0, 1) * radials(1, 0, 1) * SA(0, 0) * SB(1, 1);
	values(2, 2, 0) += 91.1715 * CA(0, 2, 0, 0, 1) * CB(0, 2, 0, 0, 0) * radials(1, 1, 0) * SA(1, 1) * SB(0, 0);
	values(2, 2, 0) += 52.6379 * CA(0, 2, 0, 0, 1) * CB(0, 2, 0, 0, 1) * radials(2, 1, 1) * SA(1, 1) * SB(1, 1);
}
}
}
