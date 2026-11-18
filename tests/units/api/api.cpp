/*
 * Unit tests of the parts of the ECPIntegrator API that the suite (which goes through
 * System::setup, set_gaussian_basis and set_ecp_basis) does not use: ECPs from the built-in
 * library, moving the basis and ECPs, and the getters. Each is checked for identical results to
 * the suite's path, so their accuracy is the suite's.
 */
#include "api.hpp"

#include <gtest/gtest.h>

#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "system.hpp"

using namespace testsupport;

namespace {

// The Gaussian basis of sys in set_gaussian_basis's arrays, with the shells displaced by shift
// (one xyz triple per atom)
struct Basis {
  std::vector<double> coords, exps, coefs;
  std::vector<int> ams, lens;
};
Basis basis_arrays(const System& sys, const std::vector<double>& shift = {}) {
  Basis b;
  for (const Shell& sh : sys.shells) {
    for (int k = 0; k < 3; k++)
      b.coords.push_back(sys.atoms[sh.atom].xyz[k] +
                         (shift.empty() ? 0.0 : shift[3 * sh.atom + k]));
    b.exps.insert(b.exps.end(), sh.exps.begin(), sh.exps.end());
    b.coefs.insert(b.coefs.end(), sh.coefs.begin(), sh.coefs.end());
    b.ams.push_back(sh.l);
    b.lens.push_back(sh.exps.size());
  }
  return b;
}

// ag3_mdf, whose ECPs (ecp28mdf, from the library) are on every atom
const System& ag3() {
  static const System sys = load_system("ag3_mdf");
  return sys;
}

}  // namespace

TEST(Api, EcpBasisFromLibrary) {
  const System& sys = ag3();
  libecpint::ECPIntegrator f, g;
  sys.setup(f, 1);
  const Basis b = basis_arrays(sys);
  g.set_gaussian_basis(sys.shells.size(), b.coords.data(), b.exps.data(), b.coefs.data(),
                       b.ams.data(), b.lens.data());
  std::vector<double> centres;
  std::vector<int> charges;
  std::vector<std::string> names;
  for (const Atom& a : sys.atoms) {
    centres.insert(centres.end(), a.xyz.begin(), a.xyz.end());
    charges.push_back(a.Z);
    names.push_back("ecp28mdf");
  }
  g.set_ecp_basis_from_library(sys.atoms.size(), centres.data(), charges.data(), names,
                               share_dir());
  g.init(1);
  f.compute_integrals();
  g.compute_integrals();
  EXPECT_EQ(f.integrals.data, g.integrals.data);
  f.compute_first_derivs();
  g.compute_first_derivs();
  ASSERT_EQ(f.first_derivs.size(), g.first_derivs.size());
  for (size_t i = 0; i < f.first_derivs.size(); i++)
    EXPECT_EQ(f.first_derivs[i].data, g.first_derivs[i].data) << "derivative " << i;
}

TEST(Api, UpdateCoordinates) {
  // Moving the shells and ECPs of an initialised integrator gives the integrals of a fresh one at
  // the new geometry
  const System& sys = ag3();
  const std::vector<double> shift = {0.01, -0.02, 0.03, 0.0, 0.015, -0.01, -0.005, 0.0, 0.02};
  ASSERT_EQ(shift.size(), 3 * sys.atoms.size());
  System moved = sys;
  for (size_t a = 0; a < sys.atoms.size(); a++)
    for (int k = 0; k < 3; k++) moved.atoms[a].xyz[k] += shift[3 * a + k];

  libecpint::ECPIntegrator f, fresh;
  sys.setup(f, 1);
  f.compute_integrals();  // at the old geometry first
  const Basis b = basis_arrays(sys, shift);
  f.update_gaussian_basis_coords(sys.shells.size(), b.coords.data());
  std::vector<double> centres;
  for (const Atom& a : moved.atoms) centres.insert(centres.end(), a.xyz.begin(), a.xyz.end());
  f.update_ecp_basis_coords(sys.atoms.size(), centres.data());
  moved.setup(fresh, 1);
  f.compute_integrals();
  fresh.compute_integrals();
  EXPECT_EQ(f.integrals.data, fresh.integrals.data);
  f.compute_first_derivs();
  fresh.compute_first_derivs();
  for (size_t i = 0; i < f.first_derivs.size(); i++)
    EXPECT_EQ(f.first_derivs[i].data, fresh.first_derivs[i].data) << "derivative " << i;
}

TEST(Api, Getters) {
  const System sys = load_system("testpair");
  libecpint::ECPIntegrator f;
  sys.setup(f, 2);
  f.compute_integrals();
  f.compute_first_derivs();
  f.compute_second_derivs();
  EXPECT_EQ(*f.get_integrals(), f.integrals.data);
  const auto first = f.get_first_derivs();
  ASSERT_EQ(first.size(), f.first_derivs.size());
  for (size_t i = 0; i < first.size(); i++) EXPECT_EQ(*first[i], f.first_derivs[i].data);
  const auto second = f.get_second_derivs();
  ASSERT_EQ(second.size(), f.second_derivs.size());
  for (size_t i = 0; i < second.size(); i++) EXPECT_EQ(*second[i], f.second_derivs[i].data);
}

TEST(Api, RepeatedDerivatives) {
  // Computing the derivatives again replaces the previous results rather than adding to them
  const System sys = load_system("testpair");
  libecpint::ECPIntegrator f;
  sys.setup(f, 2);
  f.compute_first_derivs();
  f.compute_second_derivs();
  const auto first = f.first_derivs;
  const auto second = f.second_derivs;
  f.compute_first_derivs();
  f.compute_second_derivs();
  ASSERT_EQ(f.first_derivs.size(), first.size());
  for (size_t i = 0; i < first.size(); i++) EXPECT_EQ(f.first_derivs[i].data, first[i].data);
  ASSERT_EQ(f.second_derivs.size(), second.size());
  for (size_t i = 0; i < second.size(); i++) EXPECT_EQ(f.second_derivs[i].data, second[i].data);
}

TEST(Api, AngularMomentumLimit) {
  // Basis angular momentum plus derivative order, and ECP angular momentum, up to
  // LIBECPINT_MAX_L (8) are supported; beyond, the library throws rather than computing wrong
  // values (release builds once did, the check being an assert)
  double coords[3] = {0.0, 0.0, 1.0}, exps[1] = {1.0}, coefs[1] = {1.0};
  int ams[1] = {7}, lens[1] = {1};
  double ecoords[3] = {0.0, 0.0, 0.0}, eexps[1] = {2.0}, ecoefs[1] = {1.0};
  int eams[1] = {0}, ens[1] = {2}, elens[1] = {1};
  auto integrator = [&]() {
    auto f = std::make_unique<libecpint::ECPIntegrator>();
    f->set_gaussian_basis(1, coords, exps, coefs, ams, lens);
    f->set_ecp_basis(1, ecoords, eexps, ecoefs, eams, ens, elens);
    return f;
  };
  EXPECT_NO_THROW(integrator()->init(1));                      // k shell, first derivatives
  EXPECT_THROW(integrator()->init(2), std::invalid_argument);  // k shell, second derivatives
  eams[0] = LIBECPINT_MAX_L + 1;
  libecpint::ECPIntegrator g;
  EXPECT_THROW(g.set_ecp_basis(1, ecoords, eexps, ecoefs, eams, ens, elens), std::invalid_argument);
}

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
