/*
 *      Copyright (c) 2020 Robert Shaw
 *		This file is a part of Libecpint.
 *
 *      Permission is hereby granted, free of charge, to any person obtaining
 *      a copy of this software and associated documentation files (the
 *      "Software"), to deal in the Software without restriction, including
 *      without limitation the rights to use, copy, modify, merge, publish,
 *      distribute, sublicense, and/or sell copies of the Software, and to
 *      permit persons to whom the Software is furnished to do so, subject to
 *      the following conditions:
 *
 *      The above copyright notice and this permission notice shall be
 *      included in all copies or substantial portions of the Software.
 *
 *      THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 *      EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 *      MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 *      NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
 *      LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 *      OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 *      WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#include "api.hpp"

#include <cassert>
#include <cmath>
#include <iostream>

#include "mathutil.hpp"

namespace libecpint {

void ECPIntegrator::set_gaussian_basis(const int nshells, const double* coords,
                                       const double* exponents, const double* coefs, const int* ams,
                                       const int* shell_lengths) {
  int ctr = 0;
  min_alpha = 100.0;
  for (int i = 0; i < nshells; i++) {
    ncart += (ams[i] + 1) * (ams[i] + 2) / 2;
    std::array<double, 3> cvec = {coords[3 * i], coords[3 * i + 1], coords[3 * i + 2]};
    GaussianShell newShell(cvec, ams[i]);
    if (ams[i] > maxLB) maxLB = ams[i];
    for (int n = 0; n < shell_lengths[i]; n++) {
      newShell.addPrim(exponents[ctr], coefs[ctr]);
      ctr++;
    }
    min_alpha = newShell.min_exp < min_alpha ? newShell.min_exp : min_alpha;
    shells.push_back(newShell);
  }
  basis_is_set = true;
}

void ECPIntegrator::set_ecp_basis(const int necps, const double* coords, const double* exponents,
                                  const double* coefs, const int* ams, const int* ns,
                                  const int* shell_lengths) {
  int ctr = 0;
  for (int i = 0; i < necps; i++) {
    ECP newU(&coords[3 * i]);
    for (int n = 0; n < shell_lengths[i]; n++) {
      newU.addPrimitive(ns[ctr], ams[ctr], exponents[ctr], coefs[ctr]);
      ctr++;
    }
    newU.sort();
    ecps.addECP(newU, 0);
  }
  ecp_is_set = true;
}

void ECPIntegrator::set_ecp_basis_from_library(const int necps, const double* coords,
                                               const int* charges,
                                               const std::vector<std::string>& names,
                                               const std::string& share_dir) {
  for (int i = 0; i < necps; i++) {
    std::array<double, 3> center = {coords[3 * i], coords[3 * i + 1], coords[3 * i + 2]};
    int q = charges[i];
    std::string filename = share_dir + "/raw/" + names[i] + ".ecp";
    ecps.addECP_from_file(q, center, filename);
  }
  ecp_is_set = true;
}

void ECPIntegrator::update_gaussian_basis_coords(const int nshells, const double* coords) {
  assert(nshells == shells.size());

  for (int i = 0; i < nshells; i++) {
    shells[i].localCenter[0] = coords[3 * i];
    shells[i].localCenter[1] = coords[3 * i + 1];
    shells[i].localCenter[2] = coords[3 * i + 2];
  }
}

void ECPIntegrator::update_ecp_basis_coords(const int necps, const double* coords) {
  assert(necps == ecps.getN());

  for (int i = 0; i < necps; i++)
    ecps.getECP(i).setPos(coords[3 * i], coords[3 * i + 1], coords[3 * i + 2]);
}

void ECPIntegrator::init(const int deriv_, const double tol, const double radial_tol,
                         const unsigned smallGrid, const unsigned bigGrid) {
  assert(ecp_is_set);
  assert(basis_is_set);
  deriv = std::max(0, std::min(2, deriv_));

  // Set the screening tolerance for shells and shell pairs
  tolerance = tol;
  shell_tolerance = screening::shell_tolerance(tolerance, maxLB, min_alpha);

  ecpint =
      std::make_shared<ECPIntegral>(maxLB, ecps.getMaxL(), deriv, radial_tol, smallGrid, bigGrid);

  // Determine the internal atom ids
  natoms = 0;
  std::vector<std::array<double, 3>> centers;
  for (auto& s : shells) {
    int i = 0;
    bool found = false;
    while (!found && (i < centers.size())) {
      double diff = std::abs(centers[i][0] - s.centerVec[0]);
      diff += std::abs(centers[i][1] - s.centerVec[1]);
      diff += std::abs(centers[i][2] - s.centerVec[2]);
      if (diff < 1e-4) {
        s.atom_id = i;
        found = true;
      }
      i++;
    }
    if (!found) {
      s.atom_id = natoms;
      natoms++;
      centers.push_back({s.centerVec[0], s.centerVec[1], s.centerVec[2]});
    }
  }

  for (int n = 0; n < ecps.getN(); n++) {
    ECP& U = ecps.getECP(n);
    int i = 0;
    bool found = false;
    while (!found && (i < centers.size())) {
      double diff = std::abs(centers[i][0] - U.center_[0]);
      diff += std::abs(centers[i][1] - U.center_[1]);
      diff += std::abs(centers[i][2] - U.center_[2]);
      if (diff < 1e-4) {
        U.atom_id = i;
        found = true;
      }
      i++;
    }
    if (!found) {
      U.atom_id = natoms;
      natoms++;
      centers.push_back({U.center_[0], U.center_[1], U.center_[2]});
    }
  }
}

ShellPairScreen ECPIntegrator::screen(const int d) const {
  return {screening::shell_ecp_bounds(shells, ecps, d), shell_tolerance, tolerance};
}

void ECPIntegrator::compute_integrals() {
  integrals.assign(ncart, ncart, 0.0);
  integrals_loop(
      shells, ecps, screen(0),
      [&](const ECP& U, const GaussianShell& A, const GaussianShell& B, TwoIndex<double>& v) {
        ecpint->compute_shell_pair(U, A, B, v);
      },
      [&](const int n1, const int n2, const TwoIndex<double>& block) {
        for (int i = n1; i < n1 + block.dims[0]; i++) {
          for (int j = n2; j < n2 + block.dims[1]; j++) {
            integrals(i, j) = block(i - n1, j - n2);
            integrals(j, i) = integrals(i, j);
          }
        }
      });
}

template <typename Add>
void ECPIntegrator::first_derivs_loop(Add add) {
  assert(deriv > 0);
  libecpint::first_derivs_loop(
      shells, ecps, screen(1),
      [&](const ECP& U, const GaussianShell& A, const GaussianShell& B,
          std::array<TwoIndex<double>, 9>& v) {
        ecpint->compute_shell_pair_derivative(U, A, B, v);
      },
      add);
}

template <typename Add>
void ECPIntegrator::second_derivs_loop(Add add) {
  assert(deriv > 1);
  libecpint::second_derivs_loop(
      shells, ecps, natoms, screen(2),
      [&](const ECP& U, const GaussianShell& A, const GaussianShell& B,
          std::array<TwoIndex<double>, 45>& v) {
        ecpint->compute_shell_pair_second_derivative(U, A, B, v);
      },
      add);
}

void ECPIntegrator::compute_first_derivs() {
  for (int n = 0; n < 3 * natoms; n++) first_derivs.push_back(TwoIndex<double>(ncart, ncart, 0.0));
  first_derivs_loop([&](const int i, const int k, const int l, const double v, const bool mirror) {
    TwoIndex<double>& M = first_derivs[i];
    M(k, l) += v;
    if (mirror) M(l, k) = M(k, l);
  });
}

// Only the block (k in shell A, l in shell B) of each shell pair with B <= A is calculated; the
// matrices are symmetric, so for B < A the block also stands for its transpose
std::vector<double> ECPIntegrator::compute_gradient(const double* density) {
  std::vector<double> gradient(3 * natoms, 0.0);
  first_derivs_loop([&](const int i, const int k, const int l, const double v, const bool mirror) {
    gradient[i] +=
        v * (mirror ? density[k * ncart + l] + density[l * ncart + k] : density[k * ncart + l]);
  });
  return gradient;
}

void ECPIntegrator::compute_second_derivs() {
  const int nhess = (3 * natoms * (3 * natoms + 1)) / 2;
  for (int n = 0; n < nhess; n++) second_derivs.push_back(TwoIndex<double>(ncart, ncart, 0.0));
  second_derivs_loop([&](const int i, const int k, const int l, const double v, const bool mirror) {
    TwoIndex<double>& M = second_derivs[i];
    M(k, l) += v;
    if (mirror) M(l, k) = M(k, l);
  });
}

void ECPIntegrator::compute_bfield_derivs() {
  assert(deriv > 0);
  bfield_derivs.assign(3, TwoIndex<double>(ncart, ncart, 0.0));
  bfield_derivs_loop(
      shells, ecps,
      ShellPairScreen{screening::shell_ecp_bfield_bounds(shells, ecps), shell_tolerance, tolerance},
      [&](const ECP& U, const GaussianShell& A, const GaussianShell& B,
          std::array<TwoIndex<double>, 3>& v) {
        ecpint->compute_shell_pair_derivative_Bfield(U, A, B, v);
      },
      [&](const int n1, const int n2, const std::array<TwoIndex<double>, 3>& block) {
        for (int k = 0; k < 3; k++)
          for (int i = n1; i < n1 + block[k].dims[0]; i++)
            for (int j = n2; j < n2 + block[k].dims[1]; j++) {
              bfield_derivs[k](i, j) = block[k](i - n1, j - n2);
              if (n1 != n2) bfield_derivs[k](j, i) = -bfield_derivs[k](i, j);
            }
      });
}

std::vector<double> ECPIntegrator::compute_hessian(const double* density) {
  // Contracted in the packed order of second_derivs, then unpacked
  std::vector<double> packed((3 * natoms * (3 * natoms + 1)) / 2, 0.0);
  second_derivs_loop([&](const int i, const int k, const int l, const double v, const bool mirror) {
    packed[i] +=
        v * (mirror ? density[k * ncart + l] + density[l * ncart + k] : density[k * ncart + l]);
  });
  return unpack_hessian(packed, natoms);
}

}  // namespace libecpint
