/*
 * Test systems: a molecule, a Gaussian basis and ECPs, read from the files in tests/data and listed
 * in the manifest tests/data/systems.txt. See tests/README.md for the formats.
 */
#ifndef LIBECPINT_TESTS_SYSTEM_HPP
#define LIBECPINT_TESTS_SYSTEM_HPP

#include <array>
#include <map>
#include <string>
#include <vector>

#include "api.hpp"

namespace testsupport {

/// Directory of the test data (tests/data) and of the library's ECPs (share/libecpint)
std::string data_dir();
std::string share_dir();

struct Atom {
  std::string symbol;
  int Z;
  std::array<double, 3> xyz;  // bohr
};

struct Shell {
  int atom;  // index into System::atoms
  int l;
  std::vector<double> exps, coefs;  // coefficients include the normalisation (see read_gbs)
  std::vector<double> input_coefs;  // coefficients as read, for System::hash
};

struct EcpPrimitive {
  int atom;
  int n, l;  // d r^(n-2) exp(-a r^2) in channel l
  double a, d;
};

/// One line of the manifest
struct SystemSpec {
  std::string name, xyz, basis;
  std::vector<std::pair<std::string, std::string>> ecps;  // (element, source), see load_system
  int deriv = 0;                                          // highest derivative order tested
  bool live = false;  // whether to compare with the reference engine in the test run
};

struct System {
  SystemSpec spec;
  std::vector<Atom> atoms;
  std::vector<Shell> shells;
  std::vector<EcpPrimitive> ecp_prims;  // grouped by ECP, in atom order

  /// Sets the basis and ECPs of f and initialises it for derivatives up to order deriv
  void setup(libecpint::ECPIntegrator& f, int deriv, double tol = 1e-12,
             double radial_tol = 1e-15) const;

  /// Hash (hexadecimal) of everything the integrals depend on: geometry, basis and ECPs. It uses
  /// only values as read (not the normalised coefficients, whose last bits depend on the platform's
  /// pow), and the ECP primitives in a canonical order (ECP::sort leaves the order within a channel
  /// to the standard library), so that it is the same on every platform.
  std::string hash() const;

  /// The atom (index into atoms) of each atom id assigned by ECPIntegrator::init
  std::vector<int> integrator_atoms() const;

  /// A description of the atom with ECPIntegrator atom id id, e.g. "O3"
  std::string atom_label(int id) const;

  /// A description of basis function k (in the order of the integrals), e.g. "O3 shell 5 p(y)"
  std::string function_label(int k) const;
};

/// Reads the manifest (one system per line as key=value pairs; '#' starts a comment)
std::vector<SystemSpec> read_manifest(const std::string& path);

/// Loads a system, with paths relative to data_dir(). Each ECP source is "lib/NAME" (the library
/// file share/libecpint/raw/NAME.ecp) or "data/NAME" (tests/data/ecp/NAME.ecp), and gives the ECP
/// for every atom of its element.
System load_system(const SystemSpec& spec);

/// Loads the system called name from the manifest tests/data/systems.txt
System load_system(const std::string& name);

/// Reads an xyz file. Coordinates are in angstrom unless the comment line contains "bohr".
std::vector<Atom> read_xyz(const std::string& path);

/// Reads a Gaussian94-format basis file. Blocks are keyed by element symbol (for every atom of
/// that element) or by atom number (1-based, for that atom only, taking precedence). Coefficients
/// are multiplied by the normalisation of their primitives, and each contraction is then
/// normalised (both for the x^l component), as quantum chemistry codes pass them to libecpint.
struct GbsShell {
  int l;
  std::vector<double> exps, coefs;
  std::vector<double> input_coefs;  // coefficients as read, before normalisation
};
std::map<std::string, std::vector<GbsShell>> read_gbs(const std::string& path);

/// Atomic number of an element symbol (case-insensitive), 0 if unknown
int atomic_number(const std::string& symbol);

}  // namespace testsupport

#endif
