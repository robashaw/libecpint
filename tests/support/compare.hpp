/*
 * Comparison of integrals and derivatives against references, reporting the worst element in terms
 * of the system (atoms, shells, Cartesian components) rather than raw indices.
 */
#ifndef LIBECPINT_TESTS_COMPARE_HPP
#define LIBECPINT_TESTS_COMPARE_HPP

#include <gtest/gtest.h>

#include <string>
#include <vector>

#include "multiarr.hpp"
#include "reference.hpp"
#include "system.hpp"

namespace testsupport {

/**
 * Tolerance for a set of values (e.g. all first derivatives): each value must satisfy
 *   |value - reference| <= rel * max |reference| + CONDITIONING_TOL * scale,
 * the largest reference value of the set setting the scale of the first term. Rounding noise
 * between the library and the reference engine is ~1e-14 of that scale; errors that the tests
 * exist to catch have been >= 1e-13.
 */
constexpr double DEFAULT_REL_TOL = 1e-13;

/**
 * The second term allows for the precision limit of the library's type 2 method: it expands each
 * shell about the ECP centre, and for a primitive with exponent a at distance d from it the terms
 * of that expansion exceed the result by up to ~(1 + d sqrt(2a))^l (l including the angular
 * momentum added by derivatives), so small relative errors of the radial integrals (rounding,
 * ~5e-15 for the primitives above SG_TIGHT) are amplified by that much. scale is a value's
 * conditioning scale (see reference_integrals). CONDITIONING_TOL is ~3x the largest
 * (error - rel * max |reference|) / scale over every element of the live systems with the library
 * unscreened: 3.2e-15, for an s primitive with exponent 2266 in the second derivatives of la2o3
 * (testpair needs 5e-17, the others nothing). Removing the limit needs a formulation without that
 * expansion (noted for later).
 */
constexpr double CONDITIONING_TOL = 1e-14;

/// The worst element of a comparison (the furthest beyond its tolerance), and whether all passed
struct Comparison {
  bool pass = true;
  double max_ref = 0.0, tol = 0.0;  // tol: rel * max_ref, the tolerance without conditioning
  double worst_err = 0.0, worst_tol = 0.0;
  size_t worst = 0, nbad = 0;
  std::string describe(const std::string& what, const std::string& worst_label) const;
};

/// Compares values elementwise (see DEFAULT_REL_TOL), with conditioning scales if given (see
/// CONDITIONING_TOL)
Comparison compare_values(const std::vector<double>& value, const std::vector<double>& ref,
                          double rel = DEFAULT_REL_TOL, const std::vector<double>* scale = nullptr);

/// Labels of derivative indices: 3 * atom id + coordinate (first), packed order (second), field
/// component (bfield)
std::string first_label(const System& sys, int index);
std::string second_label(const System& sys, int index, int natoms);
std::string bfield_label(int k);

/**
 * Compares integrals and derivatives (up to deriv; the magnetic field derivatives with the first)
 * as gtest assertions, each set (ints, first, second, bfield) with its own scale and, if given,
 * conditioning scales, describing the worst element of any failing set.
 */
::testing::AssertionResult compare_integrals(const System& sys, const Integrals& value,
                                             const Integrals& ref, int deriv,
                                             double rel = DEFAULT_REL_TOL,
                                             const Integrals* scale = nullptr);

}  // namespace testsupport

#endif
