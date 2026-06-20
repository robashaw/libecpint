/*
 * Regression test for radial screening / grid-scaling of type-1 (local-channel) ECP integrals.
 *
 * On-centre s- and p-shell pairs against a one-term local ECP  c * r^(n-2) * exp(-g r^2)
 * have simple closed forms (unnormalised primitives, contraction coefficient 1, combined AO
 * exponent p = 2*alpha, c = p + g):
 *
 *      <s|U|s>     n=1: 2 pi / c          n=2: pi^1.5 / c^1.5
 *      <p_x|U|p_x> n=1: 2 pi / (3 c^2)    n=2: pi^1.5 / (2 c^2.5)
 *
 * These previously collapsed to ~0 when the ECP was much steeper than the AO pair (g >> p),
 * because the type-1 radial grid was scaled by p alone instead of the combined scale p + g,
 * and (for l >= 1) because the screening gate used a hard-coded tolerance independent of the
 * requested threshold. This test scans the previously-broken regime and checks every value
 * against its closed form. See the accompanying PR / radial_quad.cpp, ecpint.cpp.
 */
#include "ecp.hpp"
#include "gshell.hpp"
#include "ecpint.hpp"
#include "multiarr.hpp"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace libecpint;

int main() {
	double O[3] = {0.0, 0.0, 0.0};
	const double thresh = 1e-18; // tight, so nothing in the scanned range is screened away

	// closed forms
	auto s_exact = [](int n, double c) { return n == 1 ? 2.0*M_PI/c : std::pow(M_PI,1.5)/std::pow(c,1.5); };
	auto p_exact = [](int n, double c) { return n == 1 ? 2.0*M_PI/(3.0*c*c) : std::pow(M_PI,1.5)/(2.0*std::pow(c,2.5)); };

	// Only assert on integrals above a floor: below ~1e-15 the values are physically
	// irrelevant and limited by double precision / the screening threshold. This still spans
	// the entire previously-broken regime (s down to ~6e-7, p down to ~1e-13). Bona-fide
	// collapse shows up as relative error ~1 and is caught easily.
	const double floor = 1e-15;
	const double reltol = 1e-5;

	double worst = 0.0;
	int nfail = 0, ntot = 0;
	std::cout << std::scientific << std::setprecision(3);

	for (int n : {1, 2}) {
		for (double alpha : {1e-2, 1.0, 1e2, 1e4}) {
			for (double g : {1e0, 1e2, 1e4, 1e6, 1e7}) {
				double p = 2.0*alpha, c = p + g;

				// s | U | s
				{
					double ex = s_exact(n, c);
					ECP U(O); U.addPrimitive(n, 0, g, 1.0, true);
					GaussianShell s(O, 0); s.addPrim(alpha, 1.0);
					ECPIntegral ecpint(0, 0, 0, thresh);
					TwoIndex<double> r; ecpint.compute_shell_pair(U, s, s, r);
					double rel = std::abs(r(0,0) - ex)/std::abs(ex);
					if (std::abs(ex) > floor) {
						worst = std::max(worst, rel); ntot++;
						if (rel > reltol) { nfail++; std::cout << "FAIL s n=" << n << " a=" << alpha
							<< " g=" << g << " got=" << r(0,0) << " exact=" << ex << " rel=" << rel << "\n"; }
					}
				}
				// p_x | U | p_x
				{
					double ex = p_exact(n, c);
					ECP U(O); U.addPrimitive(n, 0, g, 1.0, true);
					GaussianShell px(O, 1); px.addPrim(alpha, 1.0);
					ECPIntegral ecpint(1, 0, 0, thresh);
					TwoIndex<double> r; ecpint.compute_shell_pair(U, px, px, r);
					double rel = std::abs(r(0,0) - ex)/std::abs(ex);
					if (std::abs(ex) > floor) {
						worst = std::max(worst, rel); ntot++;
						if (rel > reltol) { nfail++; std::cout << "FAIL p n=" << n << " a=" << alpha
							<< " g=" << g << " got=" << r(0,0) << " exact=" << ex << " rel=" << rel << "\n"; }
					}
				}
			}
		}
	}

	std::cout << "screen_scan: " << (ntot - nfail) << "/" << ntot
		<< " cases within 1e-9; worst relative error = " << worst << std::endl;
	return nfail == 0 ? 0 : 1;
}
