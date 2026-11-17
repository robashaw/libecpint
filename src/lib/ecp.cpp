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

#include "ecp.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include "mathutil.hpp"

namespace libecpint {

// GaussianECP constructor and copy constructor
GaussianECP::GaussianECP() : n(0), l(0), a(0), d(0) {}
GaussianECP::GaussianECP(const int _n, const int _l, const double _a, const double _d)
    : n(_n - 2), l(_l), a(_a), d(_d) {}
GaussianECP::GaussianECP(const GaussianECP& other)
    : n(other.n), l(other.l), a(other.a), d(other.d) {}

// class ECP

ECP::ECP() : N(0), L(-1) {
  center_[0] = center_[1] = center_[2] = 0.0;
  min_exp = 1000.0;
  for (int i = 0; i < LIBECPINT_MAX_L + 1; i++) {
    min_exp_l[i] = 1000.0;
    l_starts[i] = 0;
  }
  l_starts[LIBECPINT_MAX_L + 1] = 0;
}

ECP::ECP(const double* _center) : N(0), L(-1) {
  center_[0] = _center[0];
  center_[1] = _center[1];
  center_[2] = _center[2];
  min_exp = 1000.0;
  for (int i = 0; i < LIBECPINT_MAX_L + 1; i++) {
    min_exp_l[i] = 1000.0;
    l_starts[i] = 0;
  }
  l_starts[LIBECPINT_MAX_L + 1] = 0;
}

ECP::ECP(const ECP& other) {
  gaussians = other.gaussians;
  N = other.N;
  L = other.L;
  atom_id = other.atom_id;
  min_exp = other.min_exp;
  for (int i = 0; i < LIBECPINT_MAX_L + 1; i++) {
    min_exp_l[i] = other.min_exp_l[i];
    l_starts[i] = other.l_starts[i];
  }
  l_starts[LIBECPINT_MAX_L + 1] = other.l_starts[LIBECPINT_MAX_L + 1];
  center_ = other.center_;
}

void ECP::addPrimitive(const int n, const int l, const double a, const double d,
                       const bool needSort) {
  if (l < 0 || l > LIBECPINT_MAX_L)
    throw std::invalid_argument("libecpint: ECP angular momentum " + std::to_string(l) +
                                " outside 0.." + std::to_string(LIBECPINT_MAX_L));
  if (n < 0 || n > 4)
    throw std::invalid_argument("libecpint: ECP radial power n = " + std::to_string(n) +
                                " outside 0..4 (r^(n-2))");
  GaussianECP newEcp(n, l, a, d);
  gaussians.push_back(newEcp);
  N++;
  L = l > L ? l : L;
  min_exp = a < min_exp ? a : min_exp;
  min_exp_l[l] = a < min_exp_l[l] ? a : min_exp_l[l];
  for (int lx = l + 1; lx < LIBECPINT_MAX_L + 2; lx++) l_starts[lx] += 1;
  if (needSort) sort();
}

void ECP::sort() {
  std::stable_sort(gaussians.begin(), gaussians.end(),
                   [&](const GaussianECP& g1, const GaussianECP& g2) { return (g1.l < g2.l); });
}

bool ECP::noType1() const {
  bool zero = true;
  for (auto& g : gaussians)
    if (g.l == L && fabs(g.d) > 1e-12) zero = false;
  return zero;
}

// Evaluate U_l(r), assuming that gaussians sorted by angular momentum
double ECP::evaluate(const double r, const int l) const {
  double value = 0.0;
  double r2 = r * r;
  int p;
  for (int i = l_starts[l]; i < l_starts[l + 1]; i++) {
    p = gaussians[i].n > -1 ? gaussians[i].n : MAX_POW - gaussians[i].n;
    value += FAST_POW[p](r) * gaussians[i].d * exp(-gaussians[i].a * r2);
  }
  return value;
}

void ECP::setPos(const double x, const double y, const double z) {
  center_[0] = x;
  center_[1] = y;
  center_[2] = z;
}

ECPBasis::ECPBasis() : N(0), maxL(-1) {}

void ECPBasis::addECP(const ECP& U, const int atom) {
  basis.push_back(U);
  atomList.push_back(atom);
  N++;
  maxL = U.getL() > maxL ? U.getL() : maxL;
}

ECP& ECPBasis::getECP(const int i) { return basis[i]; }
const ECP& ECPBasis::getECP(const int i) const { return basis[i]; }

int ECPBasis::getECPCore(const int q) const {
  int core = 0;
  auto it = core_electrons.find(q);
  if (it != core_electrons.end()) core = it->second;
  return core;
}

static std::vector<std::string> split_commas(const std::string& s) {
  std::vector<std::string> tokens;
  std::stringstream ss(s);
  std::string token;
  while (std::getline(ss, token, ',')) tokens.push_back(token);
  return tokens;
}

void ECPBasis::addECP_from_file(const int q, const std::array<double, 3>& coords,
                                const std::string& filename) {
  const int n_atoms = sizeof(atom_names) / sizeof(atom_names[0]);
  if (q > n_atoms)
    throw std::runtime_error("libecpint: no element name for atomic number " + std::to_string(q));
  std::ifstream file(filename);
  if (!file) throw std::runtime_error("libecpint: cannot open ECP file " + filename);

  // Molpro format: '!' starts a comment, and statements end at ';' or a newline.
  // Statements are stored lower case with all whitespace removed.
  std::vector<std::string> statements;
  std::string line, statement;
  while (std::getline(file, line)) {
    std::stringstream ss(line.substr(0, line.find('!')));
    while (std::getline(ss, statement, ';')) {
      statement.erase(std::remove_if(statement.begin(), statement.end(),
                                     [](unsigned char c) { return std::isspace(c); }),
                      statement.end());
      std::transform(statement.begin(), statement.end(), statement.begin(),
                     [](unsigned char c) { return std::tolower(c); });
      if (!statement.empty()) statements.push_back(statement);
    }
  }

  const std::string atom_name = q < 1 ? "x" : atom_names[q - 1];
  for (size_t i = 0; i < statements.size(); i++) {
    // Header: ecp,name,ncore,lmax[,nso]. The lmax+1 blocks that follow are in the order
    // lmax, 0, 1, ..., lmax-1, each a primitive count k then k primitives "n,exponent,coef".
    // Any nso spin-orbit blocks come after these, and are not read.
    const std::vector<std::string> header = split_commas(statements[i]);
    if (header.size() < 4 || header[0] != "ecp" || header[1] != atom_name) continue;

    ECP newECP;
    newECP.center_ = coords;
    const int maxl = std::stoi(header[3]);
    if (maxl > LIBECPINT_MAX_L)
      throw std::runtime_error("libecpint: ECP for " + atom_name + " in " + filename +
                               " has lmax " + header[3] + " > LIBECPINT_MAX_L");
    for (int block = 0; block <= maxl; block++) {
      const int l = block == 0 ? maxl : block - 1;
      const int nprims = std::stoi(statements.at(++i));
      for (int p = 0; p < nprims; p++) {
        const std::vector<std::string> nxc = split_commas(statements.at(++i));
        newECP.addPrimitive(std::stoi(nxc.at(0)), l, std::stod(nxc.at(1)), std::stod(nxc.at(2)));
      }
    }

    auto it = core_electrons.find(q);
    if (it == core_electrons.end()) core_electrons[q] = std::stoi(header[2]);

    newECP.sort();
    addECP(newECP, 0);
    return;
  }
  throw std::runtime_error("libecpint: no ECP for " + atom_name + " in " + filename);
}
}  // namespace libecpint
