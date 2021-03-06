/*! EXAMPLE
	shows calculation of ECP integrals and first derivatives for HI
	in a cc-pVDZ(-PP) basis, using ECP28MDF 

	USAGE
	./example [LIBECPINT_SHARE_DIR_PATH]
 */

#include "libecpint.hpp"

// if you want to try out the Eigen matrix build
#ifdef _WITH_EIGEN
	#include <Eigen/Dense>
	#include <Eigen/Core>
#endif

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <array>
#include <algorithm>
#include <iterator>

using DVec = std::vector<double>;
using IVec = std::vector<int>;
using Coord = std::array<double, 3>;

void read_basis_file(std::string, DVec&, DVec&, DVec&, IVec&, IVec&, Coord&);

int main(int argc, char* argv[]) {
	using namespace libecpint; 
	
	// this code expects the path to the libecpint/share directory
	// to be passed as an argument to the executable
	std::string share_dir = argv[1]; 
	
	// Roughly equilibrium position HI molecule
	// NOTE: coordinate values must be in BOHR
	Coord H_pos = {0.0, 0.0, 0.0};
	Coord I_pos = {0.0, 0.0, 3.0};
	
	// Let's read in the GTO basis from file
	// store in standard vectors
	DVec g_exps, g_coeffs; // length = total number of primitives in basis
	IVec g_ams, g_lens; // length = total number of shells in basis
	DVec g_coords; // length = 3*total number of shells
	
	// hydrogen cc-pVDZ
	read_basis_file("hydrogen.bas", g_exps, g_coeffs, g_coords, g_ams, g_lens, H_pos);
	// iodine cc-pVDZ-PP
	read_basis_file("iodine.bas", g_exps, g_coeffs,  g_coords, g_ams, g_lens, I_pos);
	// check sizes, should read 'Basis read: 44, 12, 36'
	std::cout << "Basis read: " << g_exps.size() << ", " << g_ams.size() << ", " << g_coords.size() << std::endl; 
	
	
	// Now to perform a calculation using the high-level API we do the following:
	ECPIntegrator factory; // object used to build ECP integral matrices
	
	// Set the orbital basis sets 
	factory.set_gaussian_basis(g_ams.size(), g_coords.data(), g_exps.data(), g_coeffs.data(), g_ams.data(), g_lens.data());
	// put an ECP on the iodine (charge = 53, name="ecp28mdf")
	std::vector<std::string> names = {"ecp28mdf"};
	int charges[1] = {53};
	factory.set_ecp_basis_from_library(1, I_pos.data(), charges, names, share_dir);
	
	// initialise - we want integrals and derivatives, so deriv_order=1
	factory.init(1);
	
	// and compute both the integrals and derivatives
	std::cout << "Computing integrals..." << std::endl;
	factory.compute_integrals();
	std::cout << "Computing first derivs..." << std::endl;
	factory.compute_first_derivs();
	std::cout << "Done." << std::endl; 
	
	// we can now access the integrals
	// as an example, let's build them into an
	// eigen matrix
#ifdef _WITH_EIGEN
	 // grab the integrals
	 std::shared_ptr<DVec> ints = factory.get_integrals();
	 
	 // map into a square matrix
	 // note that we use ROW-MAJOR ORDERING
	 // 12 basis functions
	 Eigen::MatrixXd ecpints = Eigen::Map<Eigen::Matrix<double, 12, 12, Eigen::RowMajor>>(ints->data());
	 std::cout << "ECP Integrals:" << std::endl;
	 std::cout << ecpints << std::endl; 
	
	 // we could do the same with derivatives
	 // this returns pointers in order [Hx, Hy, Hz, Ix, Iy, Iz]
	 std::vector<std::shared_ptr<DVec>> derivs = factory.get_first_derivs();
	 // e.g. the Iodine y-derivative
	 Eigen::MatrixXd iodine_y_derivs = Eigen::Map<Eigen::Matrix<double, 12, 12, Eigen::RowMajor>>(derivs[4]->data());
	 std::cout << "y-derivs on iodine:" << std::endl;
	 std::cout << iodine_y_derivs << std::endl; 
#endif
	
	return 0;
}

/*! Reads in the Gaussian orbital basis from a file. 

	filename - the path/name of the basis file
	exps - vector to place exponents in
	coeffs - vector to place coefficients in
	coords - the vector of xyz coords for each shell
	ams - vector to place angular momenta in
	lens - vector to place shell lengths in
	atom - the position of the atom being read in  
 */
void read_basis_file(std::string filename, DVec& exps, DVec& coeffs, DVec& coords, IVec& ams, IVec& lens, Coord& atom) { 
	std::ifstream input_file(filename);
	if (input_file.is_open()) {
		// We expect the file to have the format
		// L; x, c; x, c; x, c; x, c;
		// with a line for each shell
		std::string line, token; 
		while (!input_file.eof()) {
			std::getline(input_file, line);
			
			// split the line by semicolon
			size_t sc_pos = line.find(';');
			int len = -1;
			int am;
			while (sc_pos != std::string::npos) {
				token = line.substr(0, sc_pos);
				line.erase(0, sc_pos+1); 
				
				// First bit is L
				if (len == -1) {
					am = std::stoi(token);
					len++; 
				} else {
					// subsequent bits are x,c
					size_t c_pos = token.find(',');
					if (c_pos != std::string::npos) {
						exps.push_back(std::stod(token.substr(0, c_pos)));
						coeffs.push_back(std::stod(token.substr(c_pos+1, token.length())));
						len++; 
					} 
				}
				sc_pos = line.find(';'); 
			}
			if (len > 0) {
				// non-empty shell found
				ams.push_back(am);
				lens.push_back(len);
				std::copy(atom.begin(), atom.end(), std::back_inserter(coords));
			}
		}
	}
}
