/*
 * Stored reference values, tests/data/reference/<system>.ref, written by regenerate_references from
 * the reference engine. They hold the integrals in full, and the first and second derivatives and
 * the magnetic field derivatives contracted with fixed test densities (dense derivative matrices
 * would make the files too large, and are compared in full by the live tests instead), each with
 * its conditioning scale.
 */
#ifndef LIBECPINT_TESTS_REFERENCE_IO_HPP
#define LIBECPINT_TESTS_REFERENCE_IO_HPP

#include <string>
#include <vector>

#include "reference.hpp"
#include "system.hpp"

namespace testsupport {

struct ReferenceData {
  std::string system, hash, generator;
  int deriv = 0, ncart = 0, natoms = 0;
  std::vector<double> ints;      // lower triangle, row-major (i >= j)
  std::vector<double> gradient;  // per test density, 3 natoms contractions (compute_gradient)
  std::vector<double> hessian;   // per test density, (3 natoms)^2 contractions (compute_hessian)
  std::vector<double> bfield;    // per test density, 3 contractions (compute_bfield_derivs)
  // Conditioning scales of the values above (see reference_integrals and CONDITIONING_TOL), the
  // derivative ones contracted with the absolute values of the densities
  std::vector<double> ints_scale, gradient_scale, hessian_scale, bfield_scale;
};

/// The non-symmetric test densities (ncart x ncart, row-major) used for the contractions
std::vector<std::vector<double>> test_densities(int ncart);

/// The data to store for a system, from reference integrals with derivatives up to deriv and their
/// conditioning scales
ReferenceData make_reference_data(const System& sys, const Integrals& ref, const Integrals& scale,
                                  int deriv);

/// The same data from the library, using compute_gradient, compute_hessian and
/// compute_bfield_derivs for the derivatives
ReferenceData library_reference_data(const System& sys, libecpint::ECPIntegrator& f, int deriv);

/// tests/data/reference/<name>.ref
std::string reference_path(const std::string& name);

void write_reference(const std::string& path, const ReferenceData& data);

/// Reads a reference file; throws if it is missing or malformed
ReferenceData read_reference(const std::string& path);

}  // namespace testsupport

#endif
