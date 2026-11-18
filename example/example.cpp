/*! EXAMPLE
        shows the calculation of ECP integrals and their gradient for HI
        in a cc-pVDZ(-PP) basis, using the ECP28MDF pseudopotential from
        the built-in library

        USAGE
        ./example [LIBECPINT_SHARE_DIR_PATH]
 */

#include "libecpint.hpp"

// if you want to try out the Eigen matrix build
#ifdef _WITH_EIGEN
#include <Eigen/Core>
#include <Eigen/Dense>
#endif

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>

using DVec = std::vector<double>;
using IVec = std::vector<int>;
using Coord = std::array<double, 3>;

void read_basis_file(std::string, DVec&, DVec&, DVec&, IVec&, IVec&, Coord&);

int main(int argc, char* argv[]) {
  using namespace libecpint;

  // this code expects the path to the libecpint/share directory
  // to be passed as an argument to the executable
  if (argc < 2 || argv[1] == nullptr) {
    std::cerr << "Usage: " << argv[0] << " [LIBECPINT_SHARE_DIR_PATH]" << std::endl;
    return 1;
  }
  std::string share_dir = argv[1];

  // Roughly equilibrium position HI molecule
  // NOTE: coordinate values must be in BOHR
  Coord H_pos = {0.0, 0.0, 0.0};
  Coord I_pos = {0.0, 0.0, 3.0};

  // Let's read in the GTO basis from file
  // store in standard vectors
  DVec g_exps, g_coeffs;  // length = total number of primitives in basis
  IVec g_ams, g_lens;     // length = total number of shells in basis
  DVec g_coords;          // length = 3*total number of shells

  // hydrogen cc-pVDZ
  read_basis_file("hydrogen.bas", g_exps, g_coeffs, g_coords, g_ams, g_lens, H_pos);
  // iodine cc-pVDZ-PP
  read_basis_file("iodine.bas", g_exps, g_coeffs, g_coords, g_ams, g_lens, I_pos);
  // check sizes, should read 'Basis read: 44, 12, 36'
  std::cout << "Basis read: " << g_exps.size() << ", " << g_ams.size() << ", " << g_coords.size()
            << std::endl;

  // Now to perform a calculation using the high-level API we do the following:
  ECPIntegrator factory;  // object used to build ECP integral matrices

  // Set the orbital basis sets
  factory.set_gaussian_basis(g_ams.size(), g_coords.data(), g_exps.data(), g_coeffs.data(),
                             g_ams.data(), g_lens.data());
  // put an ECP on the iodine (charge = 53, name="ecp28mdf")
  std::vector<std::string> names = {"ecp28mdf"};
  int charges[1] = {53};
  factory.set_ecp_basis_from_library(1, I_pos.data(), charges, names, share_dir);

  // initialise - we want integrals and first derivatives, so deriv_order=1
  factory.init(1);
  const int ncart = factory.ncart;  // number of Cartesian basis functions (30 here)
  std::cout << "Cartesian functions: " << ncart << ", atoms: " << factory.natoms << std::endl;

  // compute the integrals, stored row-major as an ncart x ncart matrix
  std::cout << "Computing integrals..." << std::endl;
  factory.compute_integrals();
  std::shared_ptr<DVec> ints = factory.get_integrals();
  std::cout << "V(0, 0) = " << (*ints)[0] << ", V(0, 1) = " << (*ints)[1] << std::endl;

  // The gradient: the derivatives contracted with a density matrix D (ncart x ncart,
  // row-major), without storing the derivative matrices. In a real program D is the
  // density of the calculation; here we make one up.
  DVec D(ncart * ncart);
  for (int i = 0; i < ncart; i++)
    for (int j = 0; j < ncart; j++) D[i * ncart + j] = 1.0 / (1.0 + i + j);
  std::cout << "Computing gradient..." << std::endl;
  DVec gradient = factory.compute_gradient(D.data());  // [Hx, Hy, Hz, Ix, Iy, Iz]
  std::cout << "Gradient: H " << gradient[0] << " " << gradient[1] << " " << gradient[2] << ", I "
            << gradient[3] << " " << gradient[4] << " " << gradient[5] << std::endl;

  // The integrals depend only on the relative positions of the atoms, so the
  // gradient sums to zero over the atoms: a quick check of the result
  double drift = 0.0;
  for (int k = 0; k < 3; k++) drift = std::max(drift, std::abs(gradient[k] + gradient[3 + k]));
  if (drift > 1e-10) {
    std::cerr << "Gradient does not sum to zero: " << drift << std::endl;
    return 1;
  }
  std::cout << "Done." << std::endl;

  // If the derivative matrices themselves are needed, factory.compute_first_derivs()
  // stores all 3*natoms of them (see the documentation for the memory needed)

  // we can now access the integrals
  // as an example, let's build them into an
  // eigen matrix
#ifdef _WITH_EIGEN
  // map into a square matrix
  // note that we use ROW-MAJOR ORDERING
  using RowMatrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
  Eigen::MatrixXd ecpints = Eigen::Map<RowMatrix>(ints->data(), ncart, ncart);
  std::cout << "ECP Integrals:" << std::endl;
  std::cout << ecpints << std::endl;
#endif

  return 0;
}

/*! Reads in the Gaussian orbital basis from a file.
        libecpint uses the coefficients exactly as given, so they are normalised here:
        each primitive, then the contraction, for the x^l component of the shell.

        filename - the path/name of the basis file
        exps - vector to place exponents in
        coeffs - vector to place (normalised) coefficients in
        coords - the vector of xyz coords for each shell
        ams - vector to place angular momenta in
        lens - vector to place shell lengths in
        atom - the position of the atom being read in
 */
void read_basis_file(std::string filename, DVec& exps, DVec& coeffs, DVec& coords, IVec& ams,
                     IVec& lens, Coord& atom) {
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
        line.erase(0, sc_pos + 1);

        // First bit is L
        if (len == -1) {
          am = std::stoi(token);
          len++;
        } else {
          // subsequent bits are x,c
          size_t c_pos = token.find(',');
          if (c_pos != std::string::npos) {
            exps.push_back(std::stod(token.substr(0, c_pos)));
            coeffs.push_back(std::stod(token.substr(c_pos + 1, token.length())));
            len++;
          }
        }
        sc_pos = line.find(';');
      }
      if (len > 0) {
        // non-empty shell found: normalise its coefficients
        const double pi = std::acos(-1.0);
        const size_t first = exps.size() - len;
        double dfac = 1.0;  // (2l - 1)!!
        for (int k = 2 * am - 1; k > 1; k -= 2) dfac *= k;
        for (size_t i = first; i < exps.size(); i++)
          coeffs[i] *= std::pow(2.0 * exps[i] / pi, 0.75) * std::pow(4.0 * exps[i], 0.5 * am) /
                       std::sqrt(dfac);
        double S = 0.0;
        for (size_t i = first; i < exps.size(); i++)
          for (size_t j = first; j < exps.size(); j++) {
            const double p = exps[i] + exps[j];
            S += coeffs[i] * coeffs[j] * std::pow(pi / p, 1.5) * dfac / std::pow(2.0 * p, am);
          }
        for (size_t i = first; i < exps.size(); i++) coeffs[i] /= std::sqrt(S);

        ams.push_back(am);
        lens.push_back(len);
        std::copy(atom.begin(), atom.end(), std::back_inserter(coords));
      }
    }
  }
}
