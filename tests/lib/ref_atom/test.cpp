/*
 * On-centre (single-atom) ECP integral test against an independent closed-form reference.
 *
 * ref_atom.output is produced by generate_reference.py, which evaluates the same integrals
 * analytically with sympy (exact angular monomial integrals over the sphere and half-integer
 * Gamma radial integrals) - a completely different method from libecpint's Bessel/recursion
 * scheme. Three scenarios are covered (keep this file in sync with the generator):
 *   A: moderate exponents, mixed ECP channels (general sanity);
 *   B: a STEEP local channel (g up to 1e5) across s/p/d shells and diffuse AOs - the
 *      regression for the type-1 grid-scale collapse (g >> p) and the l>0 screening-gate
 *      collapse at large p+g;
 *   C: a diffuse AO pair vs a standard local primitive (g/p ~ 250) - the "ordinary"
 *      augmented-basis diffuse-diffuse regime that screening silently zeroed out.
 *
 * Because the reference is exact, every non-zero reference value is checked to a relative
 * tolerance regardless of magnitude (down to ~1e-13), and symmetry-zero references are checked
 * absolutely. A tight integral threshold is used so nothing in range is screened away.
 */
#include "ecp.hpp"
#include "gshell.hpp"
#include "ecpint.hpp"
#include "multiarr.hpp"
#include "testutil.hpp"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <array>

using namespace libecpint;

namespace {
	double O[3] = {0.0, 0.0, 0.0};

	GaussianShell makeShell(int l, const std::vector<std::array<double,2>>& prims) {
		GaussianShell s(O, l);
		for (auto& pc : prims) s.addPrim(pc[0], pc[1]);
		return s;
	}
}

int main() {
	const double thresh = 1e-18; // tight: do not screen away anything in the reference range
	ECPIntegral ecpint(2, 2, 0, thresh);

	// ---- Scenario A: moderate, mixed channels ----
	ECP UA(O);
	UA.addPrimitive(2, 0, 1.5,  2.0); UA.addPrimitive(3, 0, 0.7,  0.4);
	UA.addPrimitive(2, 1, 1.2, -1.0);
	UA.addPrimitive(2, 2, 0.9,  0.5); UA.addPrimitive(4, 2, 0.6,  0.25, true);
	GaussianShell sA = makeShell(0, {{1.1,1.0}});
	GaussianShell pA = makeShell(1, {{0.8,1.0},{2.5,0.7}});
	GaussianShell dA = makeShell(2, {{0.7,1.0}});

	// ---- Scenario B: steep local channel ----
	ECP UB(O);
	UB.addPrimitive(2, 0, 1000.0, 1.0); UB.addPrimitive(2, 0, 100000.0, 0.5, true);
	GaussianShell sB = makeShell(0, {{0.05,1.0}});
	GaussianShell pB = makeShell(1, {{0.05,1.0}});
	GaussianShell dB = makeShell(2, {{0.5,1.0}});

	// ---- Scenario C: diffuse AO pair, standard local primitive ----
	ECP UC(O);
	UC.addPrimitive(2, 0, 10.0, 1.0, true);
	GaussianShell sC = makeShell(0, {{0.02,1.0}});
	GaussianShell pC = makeShell(1, {{0.02,1.0}});

	struct Pair { ECP* U; GaussianShell* A; GaussianShell* B; };
	std::vector<Pair> pairs = {
		{&UA,&sA,&sA},{&UA,&pA,&sA},{&UA,&pA,&pA},{&UA,&dA,&sA},{&UA,&dA,&pA},{&UA,&dA,&dA},
		{&UB,&sB,&sB},{&UB,&pB,&sB},{&UB,&pB,&pB},{&UB,&dB,&sB},{&UB,&dB,&pB},{&UB,&dB,&dB},
		{&UC,&sC,&sC},{&UC,&pC,&sC},{&UC,&pC,&pC},
	};

	std::vector<double> got;
	TwoIndex<double> r;
	for (auto& p : pairs) {
		ecpint.compute_shell_pair(*p.U, *p.A, *p.B, r);
		for (double v : r.data) got.push_back(v);
	}

#ifdef WRITE_NEW_BENCHMARK
	for (double v : got) std::cout << std::setprecision(15) << v << std::endl;
#endif

	// Read exact reference and compare element-by-element.
	std::ifstream in("ref_atom.output");
	if (!in.is_open()) { std::cerr << "cannot open ref_atom.output" << std::endl; return 1; }
	std::vector<double> ref; std::string line;
	while (std::getline(in, line)) if (!line.empty()) ref.push_back(std::stod(line));

	if (ref.size() != got.size()) {
		std::cerr << "size mismatch: ref " << ref.size() << " vs got " << got.size() << std::endl;
		return 1;
	}

	const double reltol = 1e-7;   // relative tolerance for genuinely non-zero integrals
	const double abstol = 1e-11;  // absolute tolerance for symmetry-zero integrals
	double worst_rel = 0.0, worst_abs = 0.0;
	int nfail = 0;
	for (size_t i = 0; i < ref.size(); i++) {
		if (std::abs(ref[i]) > 1e-30) {
			double rel = std::abs(got[i] - ref[i]) / std::abs(ref[i]);
			worst_rel = std::max(worst_rel, rel);
			if (rel > reltol) { nfail++; std::cout << "FAIL line " << i << " ref=" << std::setprecision(8)
				<< ref[i] << " got=" << got[i] << " rel=" << rel << "\n"; }
		} else {
			worst_abs = std::max(worst_abs, std::abs(got[i]));
			if (std::abs(got[i]) > abstol) { nfail++; std::cout << "FAIL line " << i
				<< " expected 0 got=" << got[i] << "\n"; }
		}
	}
	std::cout << "ref_atom: " << (ref.size() - nfail) << "/" << ref.size()
		<< " ok; worst rel(nonzero)=" << worst_rel << " worst abs(zero)=" << worst_abs << std::endl;
	return nfail == 0 ? 0 : 1;
}
