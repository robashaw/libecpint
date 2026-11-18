/*
 * Reference integrals in the layout of ECPIntegrator, from the independent reference engine
 * (integral_reference.hpp) for each ECP and shell pair, placed by the library's own assembly code
 * (assembly.hpp) with no screening. The assembly is covered by its own unit tests, so a difference
 * between these and the library's integrals is a difference in the integrals themselves.
 */
#ifndef LIBECPINT_TESTS_REFERENCE_HPP
#define LIBECPINT_TESTS_REFERENCE_HPP

#include <vector>

#include "api.hpp"
#include "integral_reference.hpp"
#include "multiarr.hpp"

namespace testsupport {

/// Integrals and their derivatives, laid out as ECPIntegrator's integrals, first_derivs,
/// second_derivs and bfield_derivs (derivative vectors are empty if not computed; the magnetic
/// field derivatives are computed with the first derivatives, as they need the same engine)
struct Integrals {
  libecpint::TwoIndex<double> ints;
  std::vector<libecpint::TwoIndex<double>> first, second, bfield;
};

/// The reference engine's input for the shells and ECPs of an ECPIntegrator (on which init has been
/// called, so that the atom ids are set)
std::vector<intref::Shell> reference_shells(const libecpint::ECPIntegrator& f);
std::vector<intref::Ecp> reference_ecps(const libecpint::ECPIntegrator& f);

/**
 * Reference integrals and derivatives up to order deriv for the system of f (initialised), and, if
 * scale is given, their conditioning scales in the same layout (see CONDITIONING_TOL in
 * compare.hpp): each value's sum over the (ECP, shell pair) blocks contributing to it of
 *   max |<A^(o_A)|U|B^(o_B)>|
 * for the derivative blocks of order o_A, o_B on the shells before assembly, where s^(o) is shell s
 * raised to l_s + o with coefficients |c| (2a)^o (1 + d_s sqrt(2a))^(l_s + o), d_s its distance
 * from the ECP. The library expands the shells about the ECP centre, which cancels by up to the
 * last factor for each primitive, and the derivative recurrence and the contraction cancel too, so
 * rounding gives errors up to ~1e-16 times this scale. The magnetic field derivatives are built
 * from the shells and the shells raised by one (ECPIntegral::compute_shell_pair_derivative_Bfield),
 * so a block's scale for them is
 *   |(A - B) x C| max |<A|U|B>| + |A - C| max |<A^(1)|U|B>| + |B - C| max |<A|U|B^(1)>|,
 * with s^(1) as above but without the factor 2a, which only differentiation brings.
 */
Integrals reference_integrals(const libecpint::ECPIntegrator& f, int deriv,
                              Integrals* scale = nullptr);

/// Shell s raised to l + o with coefficients |c| (2a)^o (1 + d sqrt(2a))^(l + o), d its distance
/// from U (without the factor (2a)^o if !diff): the size of the terms of derivatives of order o in
/// the library (see reference_integrals)
libecpint::GaussianShell conditioning_shell(const libecpint::GaussianShell& s,
                                            const libecpint::ECP& U, int o, bool diff = true);

/// The conditioning scale of the magnetic field derivative block (A | B) over U (see
/// reference_integrals), with the integrals from engine (constructed with deriv >= 1)
double bfield_conditioning_scale(const libecpint::ECPIntegral& engine, const libecpint::ECP& U,
                                 const libecpint::GaussianShell& A,
                                 const libecpint::GaussianShell& B);

/// The library's integrals and derivatives up to order deriv, computed with f's settings
Integrals library_integrals(libecpint::ECPIntegrator& f, int deriv);

}  // namespace testsupport

#endif
