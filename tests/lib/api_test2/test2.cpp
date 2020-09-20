#include "api.hpp"
#include "testutil.hpp"
#include <iostream>
#include <string>

double g_coords[51]  = {0.0, 0.0, 0.0, 
					  0.0, 0.0, 0.0,
					  0.0, 0.0, 0.0,
					  0.0, 0.0, 0.0,
					  0.0, 0.0, 0.0,
					  0.0, 0.0, 0.0,
					  0.0, 0.0, 0.0,
					  0.0, 0.0, 0.0,
					  0.0, 0.0, 0.0,
					  0.0, 0.0, 0.0,
					  0.0, 0.0, 0.0,
					  0.0, 0.0, 0.0,
					  0.0, 0.0, 2.36,
					  0.0, 0.0, 2.36,
					  0.0, 0.0, 2.36,
					  0.0, 0.0, 2.36,
					  0.0, 0.0, 2.36};
				  	
double g_exps[52]    = {2.808600E+03,4.211800E+02,5.034570E+01,1.791330E+01,3.805310E+00,1.749680E+00,4.485550E-01,1.644980E-01,
					  2.808600E+03,4.211800E+02,5.034570E+01,1.791330E+01,3.805310E+00,1.749680E+00,4.485550E-01,1.644980E-01,
				  	  4.485550E-01,1.644980E-01,5.020000E-02,
				  	  1.057520E+02,2.763680E+01,6.596560E+00,2.785220E+00,1.078120E+00,3.935370E-01,1.274690E-01,
				  	  1.057520E+02,2.763680E+01,6.596560E+00,2.785220E+00,1.078120E+00,3.935370E-01,1.274690E-01, 
				      1.274690E-01,3.940000E-02,
				  	  1.438650E+02,4.611630E+01,1.736940E+01,6.951070E+00,2.756070E+00,1.011780E+00,4.291000E-01, 
				  	  4.291000E-01,1.548000E-01,
					  1.301000E+01,1.962000E+00,4.446000E-01,1.220000E-01,
					  1.220000E-01,2.974000E-02,
					  7.270000E-01,1.410000E-01};
					  
double g_coefs[52]   = {1.606000E-03,8.393000E-03,6.957800E-02,-3.899080E-01,6.944970E-01,4.913540E-01,2.263700E-02,-3.723000E-03, 
					  -6.350000E-04,-3.492000E-03,-2.519500E-02,1.501130E-01,-3.662260E-01,-3.834220E-01,7.144680E-01,5.352530E-01, 
				  	  1.0, 1.0, 1.0, 
				  	  5.341000E-03,-8.308400E-02,4.477660E-01,5.506170E-01,1.235000E-01,-3.771000E-03,2.278000E-03, 
				      -1.308000E-03,2.292100E-02,-1.450290E-01,-2.090370E-01,9.373000E-02,6.050210E-01,4.571230E-01, 
				      1.0, 1.0, 
				      1.023700E-02,7.608300E-02,2.298070E-01,4.033470E-01,4.097280E-01,1.627900E-01,6.991000E-03, 
					  1.0, 1.0, 
					  1.968500E-02,1.379770E-01,4.781480E-01,5.012400E-01,
					  1.0, 1.0,
					  1.0, 1.0};
					  
int    g_ams[17]     = {0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 0, 0, 0, 1, 1};
int    g_lengths[17] = {8, 8, 1, 1, 1, 7, 7, 1, 1, 7, 1, 1, 4, 1, 1, 1, 1};

std::string get_deriv_id(int i, int natoms, int ncart) {
	std::string result = "";
	int nfuncs = ncart * ncart;
	int chunk = std::floor(i / nfuncs);
	int atom = std::floor(chunk / 3);
	int coord = chunk % 3; 
	int ix = i - chunk*nfuncs;
	int na = std::floor(ix / ncart);
	int nb = ix - na*ncart; 
	result = std::to_string(atom) + "." + std::to_string(coord) + "." + std::to_string(na) + "." + std::to_string(nb); 
	return result;
}

std::string get_hess_id(int i, int natoms, int ncart) {
	std::string result = "";
	int nfuncs = ncart * ncart; 
	int chunk = std::floor(i / nfuncs);
	int atom1 = (chunk < 15) ? 0 : 1;
	int atom2 = (chunk < 6) ? 0 : 1;
	int coord1, coord2;
	
	if (atom1 == atom2) {
		if (chunk % 15 < 3){
			coord1 = 0;
			coord2 = chunk % 15;
		} else if (chunk % 15 < 5) {
			coord1 = 1;
			coord2 = chunk % 15 - 2;
		} else {
			coord1 = 2;
			coord2 = 2;
		}
	} else {
		if (chunk - 6 < 3) {
			coord1 = 0;
			coord2 = chunk - 6;
		} else if (chunk - 6 < 6) {
			coord1 = 1;
			coord2 = chunk - 9;
		} else {
			coord1 = 2;
			coord2 = chunk - 12;
		}
	}
	int ix = i - chunk*nfuncs;
	int na = std::floor(ix / ncart);
	int nb = ix - na*ncart; 
	result = std::to_string(atom1) + "." + std::to_string(coord1) + ".";
	result += std::to_string(atom2) + "." + std::to_string(coord2) + ".";
	result += std::to_string(na) + "." + std::to_string(nb);
	return result;
}

void finite_diff(libecpint::ECPIntegrator& factory,  double h, int start1, int finish1, int start2, int finish2, int nshells,
	 			double* coords, double* ecp_center, int axis1, int axis2, bool shift_ecp_1, bool shift_ecp_2,
				 std::vector<double> &results) {
	std::vector<double> Ipp, Ipm, Imm, Imp;
	double rv;
	double o4h2 = 0.25 / (h * h);
	std::shared_ptr<std::vector<double>> ints;

	for (int i = start1; i < finish1; i++)
		coords[3*i+axis1] += h; 
	for (int i = start2; i < finish2; i++)
		coords[3*i+axis2] += h; 
	factory.update_gaussian_basis_coords(nshells, coords);
	if (shift_ecp_1) ecp_center[axis1] += h;
	if (shift_ecp_2) ecp_center[axis2] += h;
	factory.update_ecp_basis_coords(1, ecp_center);
	factory.compute_integrals();
	ints = factory.get_integrals();
	for (int i = 0; i < ints->size(); i++)
		Ipp.push_back((*ints)[i]);
	
	for (int i = start1; i < finish1; i++)
		coords[3*i+axis1] -= 2*h; 
	factory.update_gaussian_basis_coords(nshells, coords);
	if (shift_ecp_1) ecp_center[axis1] -= 2*h;
	factory.update_ecp_basis_coords(1, ecp_center);
	factory.compute_integrals();
	ints = factory.get_integrals();
	for (int i = 0; i < ints->size(); i++)
		Imp.push_back((*ints)[i]);
	
	for (int i = start2; i < finish2; i++)
		coords[3*i+axis2] -= 2*h; 
	factory.update_gaussian_basis_coords(nshells, coords);
	if (shift_ecp_2) ecp_center[axis2] -= 2*h;
	factory.update_ecp_basis_coords(1, ecp_center);
	factory.compute_integrals();
	ints = factory.get_integrals();
	for (int i = 0; i < ints->size(); i++)
		Imm.push_back((*ints)[i]);
	
	for (int i = start1; i < finish1; i++)
		coords[3*i+axis1] += 2*h; 
	factory.update_gaussian_basis_coords(nshells, coords);
	if (shift_ecp_1) ecp_center[axis1] += 2*h;
	factory.update_ecp_basis_coords(1, ecp_center);
	factory.compute_integrals();
	ints = factory.get_integrals();
	for (int i = 0; i < ints->size(); i++)
		Ipm.push_back((*ints)[i]);
	
	for (int i = start1; i < finish1; i++)
		coords[3*i+axis1] -= h; 
	for (int i = start2; i < finish2; i++)
		coords[3*i+axis2] += h; 
	factory.update_gaussian_basis_coords(nshells, coords);
	if (shift_ecp_1) ecp_center[axis1] -= h;
	if (shift_ecp_2) ecp_center[axis2] += h;
	factory.update_ecp_basis_coords(1, ecp_center);
	
	for (int i = 0; i < Ipp.size(); i++) {
		rv = Ipp[i] - Imp[i] - Ipm[i] + Imm[i];
		results.push_back(rv * o4h2);
	}
}

int main(int argc, char* argv[]) {
	using namespace libecpint; 
	
	double br_center[3] = {0.0, 0.0, 0.0};
	int charges[1] = {35};
	std::vector<std::string> names;
	names.push_back("ecp10mdf");
	std::string share_dir = argv[1];
	
	ECPIntegrator factory;
	factory.set_gaussian_basis(17, g_coords, g_exps, g_coefs, g_ams, g_lengths);
	factory.set_ecp_basis_from_library(1, br_center, charges, names, share_dir);
	factory.init(2);
	factory.compute_integrals();
	factory.compute_first_derivs();
	factory.compute_second_derivs();
	
	std::vector<double> flat_result, flat_finite;
	std::shared_ptr<std::vector<double>> ints = factory.get_integrals();
	for (int i = 0; i < ints->size(); i++) {
		double val = (*ints)[i];
		if (std::abs(val) > 1e-12) flat_result.push_back(val);
	}
	
	std::vector<std::shared_ptr<std::vector<double>>> first_derivs = factory.get_first_derivs();
	for (auto& s : first_derivs) {
		for (int i = 0; i < s->size(); i++) {
			double val = (*s)[i];
			if (std::abs(val) > 1e-12) flat_result.push_back(val);
		}
	}
	
	std::vector<std::shared_ptr<std::vector<double>>> second_derivs = factory.get_second_derivs();
	for (auto& s : second_derivs) {
		for (int i = 0; i < s->size(); i++) {
			double val = (*s)[i];
			if (std::abs(val) > 1e-12) flat_result.push_back(val);
		}
	}
	
	/*double h = 0.01;
	finite_diff(factory, h, 0, 12, 0, 12, 17, g_coords, br_center, 0, 0, true, true, flat_finite);
	finite_diff(factory, h, 0, 12, 0, 12, 17, g_coords, br_center, 0, 1, true, true, flat_finite);
	finite_diff(factory, h, 0, 12, 0, 12, 17, g_coords, br_center, 0, 2, true, true, flat_finite);
	finite_diff(factory, h, 0, 12, 0, 12, 17, g_coords, br_center, 1, 1, true, true, flat_finite);
	finite_diff(factory, h, 0, 12, 0, 12, 17, g_coords, br_center, 1, 2, true, true, flat_finite);
	finite_diff(factory, h, 0, 12, 0, 12, 17, g_coords, br_center, 2, 2, true, true, flat_finite);
	
	finite_diff(factory, h, 0, 12, 12, 17, 17, g_coords, br_center, 0, 0, true, false, flat_finite);
	finite_diff(factory, h, 0, 12, 12, 17, 17, g_coords, br_center, 0, 1, true, false, flat_finite);
	finite_diff(factory, h, 0, 12, 12, 17, 17, g_coords, br_center, 0, 2, true, false, flat_finite);
	finite_diff(factory, h, 0, 12, 12, 17, 17, g_coords, br_center, 1, 0, true, false, flat_finite);
	finite_diff(factory, h, 0, 12, 12, 17, 17, g_coords, br_center, 1, 1, true, false, flat_finite);
	finite_diff(factory, h, 0, 12, 12, 17, 17, g_coords, br_center, 1, 2, true, false, flat_finite);
	finite_diff(factory, h, 0, 12, 12, 17, 17, g_coords, br_center, 2, 0, true, false, flat_finite);
	finite_diff(factory, h, 0, 12, 12, 17, 17, g_coords, br_center, 2, 1, true, false, flat_finite);
	finite_diff(factory, h, 0, 12, 12, 17, 17, g_coords, br_center, 2, 2, true, false, flat_finite);
	
	finite_diff(factory, h, 12, 17, 12, 17, 17, g_coords, br_center, 0, 0, false, false, flat_finite);
	finite_diff(factory, h, 12, 17, 12, 17, 17, g_coords, br_center, 0, 1, false, false, flat_finite);
	finite_diff(factory, h, 12, 17, 12, 17, 17, g_coords, br_center, 0, 2, false, false, flat_finite);
	finite_diff(factory, h, 12, 17, 12, 17, 17, g_coords, br_center, 1, 1, false, false, flat_finite);
	finite_diff(factory, h, 12, 17, 12, 17, 17, g_coords, br_center, 1, 2, false, false, flat_finite);
	finite_diff(factory, h, 12, 17, 12, 17, 17, g_coords, br_center, 2, 2, false, false, flat_finite);

	for (int i = 0; i < flat_finite.size(); i++) {
		double diff = flat_result[i]-flat_finite[i];
		if (std::abs(diff) > 1e-2)
			std::cout << std::setw(10) <<  i << " " << std::setw(15) << get_hess_id(i, factory.natoms, factory.ncart) 
				      << " " <<  std::setw(12) << flat_result[i] <<  " " << std::setw(12) << flat_finite[i] 
					  << " " << std::setw(12) << diff << std::endl;
	}*/
		
	return check_file<double>("api_test2.output", flat_result);
}
