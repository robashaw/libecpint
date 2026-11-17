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

#include "shared_grid.hpp"

#include <algorithm>
#include <cmath>

#include "config.hpp"
#include "gaussquad.hpp"

namespace libecpint {
namespace detail {

std::shared_ptr<const SharedGrid> shared_grid(const ECP& U) {
  thread_local LruCache<std::shared_ptr<const SharedGrid>> grids(1u << 20);
  thread_local double next_id = 0.0;
  thread_local std::vector<double> key;
  key.assign(1, static_cast<double>(U.getL()));
  for (const auto& g : U.gaussians) key.insert(key.end(), {double(g.n), double(g.l), g.a, g.d});
  if (const auto* found = grids.find(key)) return *found;

  // Each ECP primitive is negligible beyond R_u, where |d| r^(k + n) exp(-u r^2) < 1e-25 for all
  // k = N + 2 in [2, 20]. The grid extends to the largest R_u, and the panels are graded: a panel
  // starting at r resolves z = 2 SG_TIGHT + (largest u with R_u > r), non-increasing in r.
  std::vector<std::pair<double, double>> extent;  // (R_u, u)
  for (const auto& g : U.gaussians) {
    if (g.d == 0.0) continue;
    const double kmin = 2.0 + g.n, kmax = 20.0 + g.n;
    auto log_size = [&](const double r) {
      return std::log(std::abs(g.d)) + (r > 1.0 ? kmax : kmin) * std::log(r) - g.a * r * r;
    };
    double r = std::sqrt(std::max(kmin, 0.0) / (2.0 * g.a));  // decreasing beyond here
    while (log_size(r + 0.01) > std::log(1e-25)) r += 0.01;
    extent.push_back({r + 0.01, g.a});
  }
  auto zres = [&](const double r) {
    double u = 0.0;
    for (const auto& e : extent)
      if (e.first > r) u = std::max(u, e.second);
    return 2.0 * SG_TIGHT + u;
  };
  double R = 0.0;
  for (const auto& e : extent) R = std::max(R, e.first);

  auto grid = std::make_shared<SharedGrid>();
  const auto& gl = gauss_legendre<SG_PANEL_POINTS>();
  std::vector<double> w;
  for (double r0 = 0.0; r0 < R;) {
    const double h = SG_PANEL_WIDTH / std::sqrt(zres(r0));
    for (int g = 0; g < SG_PANEL_POINTS; g++) {
      grid->r.push_back(r0 + 0.5 * h * (1.0 + gl.x[g]));
      w.push_back(0.5 * h * gl.w[g]);
    }
    r0 += h;
  }
  const size_t n = grid->r.size();
  grid->Uw.resize(U.getL() + 1, std::vector<double>(n));
  for (int lam = 0; lam <= U.getL(); lam++)
    for (size_t i = 0; i < n; i++) grid->Uw[lam][i] = w[i] * U.evaluate(grid->r[i], lam);
  grid->id = next_id++;
  grids.insert(key, grid, n * (U.getL() + 2) * sizeof(double));
  return grid;
}

std::shared_ptr<const SharedFactor> shared_factor(const SharedGrid& grid,
                                                  const GaussianShell& shell, const double A,
                                                  const int lmax, const BesselFunction& bessie) {
  thread_local LruCache<std::shared_ptr<const SharedFactor>> factors(
      size_t(LIBECPINT_RADIAL_CACHE_MB) << 20);
  thread_local std::vector<double> key;
  key.assign({grid.id, A});
  for (int i = 0; i < shell.nprimitive(); i++)
    if (shell.exp(i) <= SG_TIGHT) key.insert(key.end(), {shell.exp(i), shell.coef(i)});
  if (const auto* found = factors.find(key))
    if ((*found)->lmax >= lmax)
      return *found;  // values for l <= lmax do not depend on (*found)->lmax

  auto factor = std::make_shared<SharedFactor>();
  const size_t n = grid.r.size();
  factor->lmax = lmax;
  factor->F.assign((lmax + 1) * n, 0.0);
  std::vector<double> bes(lmax + 1);
  for (size_t k = 2; k < key.size(); k += 2) {
    const double a = key[k], c = key[k + 1], x = 2.0 * a * A;
    for (size_t i = 0; i < n; i++) {
      const double d = grid.r[i] - A, e = a * d * d;
      if (e > 708.0) continue;  // exp(-e) underflows
      const double f = c * std::exp(-e);
      bessie.calculate(x * grid.r[i], lmax, bes);
      for (int l = 0; l <= lmax; l++) factor->F[l * n + i] += f * bes[l];
    }
  }
  factors.insert(key, factor, factor->F.size() * sizeof(double));
  return factor;
}

std::shared_ptr<const SharedWeights> shared_weights(const SharedGrid& grid, const int lam,
                                                    const int Nmax) {
  thread_local LruCache<std::shared_ptr<const SharedWeights>> cache(1u << 20);
  thread_local std::vector<double> key;
  key.assign({grid.id, static_cast<double>(lam)});
  if (const auto* found = cache.find(key))
    if ((*found)->Nmax >= Nmax) return *found;  // values for N <= Nmax do not depend on its Nmax

  auto weights = std::make_shared<SharedWeights>();
  const size_t n = grid.r.size();
  weights->Nmax = Nmax;
  weights->W.resize((Nmax + 1) * n);
  weights->S.resize((Nmax + 1) * (n + 1));
  for (size_t i = 0; i < n; i++) {
    double w = grid.Uw[lam][i] * grid.r[i] * grid.r[i];
    for (int N = 0; N <= Nmax; N++, w *= grid.r[i]) weights->W[N * n + i] = w;
  }
  for (int N = 0; N <= Nmax; N++) {
    double* S = &weights->S[N * (n + 1)];
    S[n] = 0.0;
    for (size_t i = n; i-- > 0;) S[i] = std::max(S[i + 1], std::abs(weights->W[N * n + i]));
  }
  cache.insert(key, weights, (weights->W.size() + weights->S.size()) * sizeof(double));
  return weights;
}

}  // namespace detail
}  // namespace libecpint
