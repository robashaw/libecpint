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

#ifndef SHARED_GRID_HEAD
#define SHARED_GRID_HEAD

#include <cstddef>
#include <cstdint>
#include <cstring>
#include <list>
#include <memory>
#include <unordered_map>
#include <utility>
#include <vector>

#include "bessel.hpp"
#include "ecp.hpp"
#include "gshell.hpp"

namespace libecpint {

/**
 * Shared radial grid for the type 2 radial integrals (internal to RadialIntegral::type2).
 *
 * With shells A and B at distances A and B from the ECP centre (either may be 0, where only
 * K_0 = 1 is non-zero),
 *   sum_triples coef G = int r^(N+2) U_lam(r) F_A,l1(r) F_B,l2(r) dr,
 *   F_A,l(r) = sum_a c_a exp(-a (r - A)^2) K_l(2 a A r).
 * F_A depends only on the shell, its distance from the ECP and the grid, so with one grid per ECP
 * type it is computed once and reused for every partner shell and ECP angular momentum, instead
 * of evaluating Bessel functions per primitive triple and shell pair. The grid is a composite
 * Gauss-Legendre rule on panels resolving r^k exp(-z (r - p)^2) to ~1e-14 relative for
 * z <= 2 * SG_TIGHT + (largest ECP exponent not negligible there), any p and k <= 40, so it
 * integrates every primitive triple whose basis exponents are <= SG_TIGHT. (The powers of r narrow
 * the integrand: at p = 0 its peak is that of a Gaussian with exponent 2z, which panels sized for
 * exp(-z r^2) alone resolve only to ~5e-12.) Primitives with larger exponents are
 * integrated per primitive triple on their own window. Grids and factors are kept in small
 * per-thread LRU caches keyed by content, so results depend neither on the cache state nor on the
 * threads.
 */
namespace detail {

constexpr int SG_PANEL_POINTS = 48;
constexpr double SG_PANEL_WIDTH = 11.0;  ///< in units of 1/sqrt(z); 13 gives 3e-14 at k = 38
constexpr double SG_TIGHT = 30.0;        ///< largest basis exponent integrated on the shared grid

/// Least-recently-used cache from a key of doubles to V, limited to a total size in bytes
template <typename V>
class LruCache {
 public:
  explicit LruCache(const size_t budget) : budget_(budget) {}

  const V* find(const std::vector<double>& key) {
    // The two most recent items first (e.g. the factors of both shells of the current pair),
    // without hashing
    auto mru = items_.begin();
    for (int i = 0; i < 2 && mru != items_.end(); i++, ++mru)
      if (mru->key == key) {
        items_.splice(items_.begin(), items_, mru);
        return &mru->value;
      }
    const auto it = index_.find(hash(key));
    if (it == index_.end() || it->second->key != key) return nullptr;
    items_.splice(items_.begin(), items_, it->second);
    return &it->second->value;
  }

  // key is copied: callers reuse their key buffers so that lookups do not allocate
  void insert(const std::vector<double>& key, V value, const size_t bytes) {
    const uint64_t h = hash(key);
    const auto it = index_.find(h);
    if (it != index_.end()) {  // same key (being replaced), or a hash collision
      bytes_ -= it->second->bytes;
      items_.erase(it->second);
      index_.erase(it);
    }
    items_.push_front({key, std::move(value), bytes, h});
    index_[h] = items_.begin();
    bytes_ += bytes;
    while (bytes_ > budget_ && items_.size() > 1) {
      bytes_ -= items_.back().bytes;
      index_.erase(items_.back().hash);
      items_.pop_back();
    }
  }

  /// @return the number of items held
  size_t size() const { return items_.size(); }

 private:
  struct Item {
    std::vector<double> key;
    V value;
    size_t bytes;
    uint64_t hash;
  };
  static uint64_t hash(const std::vector<double>& key) {  // FNV-1a on 64-bit words, then mixed
    uint64_t h = 14695981039346656037ull;
    for (const double x : key) {
      uint64_t bits;
      std::memcpy(&bits, &x, sizeof bits);
      h = (h ^ bits) * 1099511628211ull;
    }
    h ^= h >> 32;  // the low bits pick the hash table bucket
    return h * 0x9e3779b97f4a7c15ull;
  }
  std::list<Item> items_;
  std::unordered_map<uint64_t, typename std::list<Item>::iterator> index_;
  size_t budget_, bytes_ = 0;
};

/// A composite Gauss-Legendre grid for one ECP
struct SharedGrid {
  std::vector<double> r;                ///< nodes
  std::vector<std::vector<double>> Uw;  ///< for each lam, weight * U_lam(r)
  double id;                            ///< distinct for each grid built by this thread
};

/**
 * @return the shared grid for U (from the per-thread cache, or built). The grid extends to the
 * largest R_u beyond which each ECP primitive is negligible (|d| r^(k + n) exp(-u r^2) < 1e-25 for
 * all k = N + 2 in [2, 20]), and its panels are graded: a panel starting at r resolves
 * z = 2 SG_TIGHT + (largest u with R_u > r), non-increasing in r.
 */
std::shared_ptr<const SharedGrid> shared_grid(const ECP& U);

/// F_A,l(r) on a shared grid, F[l * n + i] for l <= lmax, from the primitives with a <= SG_TIGHT
struct SharedFactor {
  int lmax;
  std::vector<double> F;
};

/// @return F_A,l on grid for l <= lmax (from the per-thread cache, or built)
std::shared_ptr<const SharedFactor> shared_factor(const SharedGrid& grid,
                                                  const GaussianShell& shell, double A, int lmax,
                                                  const BesselFunction& bessie);

/**
 * W_N(r) = weight * U_lam(r) * r^(N+2) on a shared grid, W[N * n + i] for N <= Nmax, and the suffix
 * maxima S[N * (n + 1) + i] = max_{j >= i} |W_N(r_j)| (0 at i = n), non-increasing in i, so that
 * the nodes beyond which the integrand is negligible are found by a binary search. These depend
 * only on the grid and lam, so are kept in a per-thread cache rather than rebuilt for every shell
 * pair.
 */
struct SharedWeights {
  int Nmax;
  std::vector<double> W, S;
};

/// @return W_N and its suffix maxima on grid for channel lam and N <= Nmax (cached, or built)
std::shared_ptr<const SharedWeights> shared_weights(const SharedGrid& grid, int lam, int Nmax);

}  // namespace detail
}  // namespace libecpint

#endif
