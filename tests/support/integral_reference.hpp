/*
 * Independent reference values for ECP integrals over contracted Cartesian Gaussian shells, and
 * their first and second derivatives with respect to the shell centres, for regression tests.
 * Deliberately shares no code with the library, and uses no expansion of the Gaussians about the
 * ECP centre (no binomial/Bessel expansions, which can cancel), and no screening beyond dropping
 * terms below ~1e-40:
 *
 *   <f|U|g> = int f U_L g d^3r + sum_{l<L} int r^2 U_l(r) sum_m P_f,lm(r) P_g,lm(r) dr,
 *   P_f,lm(r) = int f(r w) Y_lm(w) dw,
 *
 * with real orthonormal spherical harmonics Y_lm about the ECP centre, U_l(r) = sum d r^(n-2)
 * exp(-a r^2) over the ECP primitives of angular momentum l, and L the largest of these (the local
 * part). Basis functions and their derivatives with respect to their centres are evaluated
 * pointwise. Angular integrals use a grid aligned with the direction of the Gaussian's centre (for
 * P_f) or of the Gaussian product centre (local part), so the integrand is exp(k cos t) times a
 * polynomial: uniform nodes in the azimuth integrate the polynomial exactly, and panels in
 * 1 - cos t growing geometrically from 1/k resolve any concentration k. Radial integrals use
 * composite Gauss-Legendre panels no wider than 1.5/sqrt(exponent) within 12/sqrt(exponent) of
 * every Gaussian feature.
 *
 * The engine gives the integrals over each ECP and shell pair, and their derivatives with respect
 * to the two shell centres. Optionally it also gives the derivatives with respect to a uniform
 * magnetic field over gauge-including orbitals (as
 * ECPIntegral::compute_shell_pair_derivative_Bfield),
 *
 *   <f| [(A - C) x r]_k U |g> - <f| U [(B - C) x r']_k |g>,
 *
 * by multiplying each basis function pointwise by its factor (no raised shells). Derivatives with
 * respect to the ECP centre (translational invariance) and the placement in molecular matrices are
 * bookkeeping, not accuracy, and are left to the library's assembly code, which is unit-tested
 * separately (see reference.hpp).
 */
#ifndef LIBECPINT_TESTS_INTEGRAL_REFERENCE_HPP
#define LIBECPINT_TESTS_INTEGRAL_REFERENCE_HPP

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <thread>
#include <utility>
#include <vector>

namespace intref {

struct Shell {
  std::array<double, 3> centre;
  int l;
  std::vector<double> exps, coefs;  // coefficients used as given (no normalisation)
  int atom;
};

struct EcpPrimitive {
  int n, l;  // d r^(n-2) exp(-a r^2) in channel l
  double a, d;
};

struct Ecp {
  std::array<double, 3> centre;
  std::vector<EcpPrimitive> prims;
  int atom;
};

namespace detail {
using V3 = std::array<double, 3>;

// Gaussian factors exp(-x) with x above this (< 1e-52) are dropped: even multiplied by the largest
// polynomial and derivative factors arising here (~1e12), such terms are below 1e-40
constexpr double NEGLIGIBLE = 120.0;

inline void gauss_legendre(const int n, std::vector<double>& x, std::vector<double>& w) {
  x.assign(n, 0.0);
  w.assign(n, 0.0);
  for (int i = 0; i < n; i++) {
    double t = std::cos(M_PI * (i + 0.75) / (n + 0.5)), dp = 0.0;
    for (int it = 0; it < 100; it++) {
      double p0 = 1.0, p1 = t;
      for (int k = 2; k <= n; k++) {
        const double p2 = ((2 * k - 1) * t * p1 - (k - 1) * p0) / k;
        p0 = p1;
        p1 = p2;
      }
      dp = n * (t * p1 - p0) / (t * t - 1.0);
      const double dt = p1 / dp;
      t -= dt;
      if (std::abs(dt) < 1e-16) break;
    }
    x[i] = t;
    w[i] = 2.0 / ((1.0 - t * t) * dp * dp);
  }
}

struct GL16 {
  std::vector<double> x, w;
  GL16() { gauss_legendre(16, x, w); }
};
inline const GL16& gl16() {
  static const GL16 g;
  return g;
}

// Cartesian components (i, j, k) of angular momentum l, x-power descending then y-power descending
inline std::vector<std::array<int, 3>> cartesians(const int l) {
  std::vector<std::array<int, 3>> c;
  for (int i = l; i >= 0; i--)
    for (int j = l - i; j >= 0; j--) c.push_back({i, j, l - i - j});
  return c;
}

// Derivative multi-indices of order <= dmax: none; x, y, z; xx, xy, xz, yy, yz, zz
inline std::vector<std::array<int, 3>> variants(const int dmax) {
  std::vector<std::array<int, 3>> v{{0, 0, 0}};
  if (dmax >= 1) v.insert(v.end(), {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}});
  if (dmax >= 2)
    v.insert(v.end(), {{2, 0, 0}, {1, 1, 0}, {1, 0, 1}, {0, 2, 0}, {0, 1, 1}, {0, 0, 2}});
  return v;
}
inline int order(const std::array<int, 3>& v) { return v[0] + v[1] + v[2]; }
inline int var1(const int p) { return 1 + p; }
inline int var2(int p, int q) {
  if (p > q) std::swap(p, q);
  static const int ix[3][3] = {{0, 1, 2}, {1, 3, 4}, {2, 4, 5}};
  return 4 + ix[p][q];
}

// (-d/dA)^d [t^i exp(-a t^2)] / exp(-a t^2), with t = x - A, for d <= 2
inline double factor(const int i, const int d, const double a, const double t) {
  auto pw = [t](const int k) {
    if (k < 0) return 0.0;
    double v = 1.0;
    for (int m = 0; m < k; m++) v *= t;
    return v;
  };
  if (d == 0) return pw(i);
  if (d == 1) return 2.0 * a * pw(i + 1) - i * pw(i - 1);
  return 4.0 * a * a * pw(i + 2) - 2.0 * a * (2 * i + 1) * pw(i) + double(i) * (i - 1) * pw(i - 2);
}

// c * polynomial part of each variant of each Cartesian component of a primitive with exponent a,
// at displacement t from its centre: out[v * ncart + c]
inline void primitive_poly(const double a, const double coef,
                           const std::vector<std::array<int, 3>>& carts,
                           const std::vector<std::array<int, 3>>& vars, const V3& t,
                           std::vector<double>& out) {
  const int l = carts.empty() ? 0 : carts[0][0];
  std::vector<std::array<std::array<double, 3>, 3>> f(l + 1);  // [power][dim][derivative]
  for (int x = 0; x < 3; x++)
    for (int i = 0; i <= l; i++)
      for (int d = 0; d < 3; d++) f[i][x][d] = factor(i, d, a, t[x]);
  const size_t nc = carts.size();
  out.resize(vars.size() * nc);
  for (size_t v = 0; v < vars.size(); v++)
    for (size_t c = 0; c < nc; c++)
      out[v * nc + c] = coef * f[carts[c][0]][0][vars[v][0]] * f[carts[c][1]][1][vars[v][1]] *
                        f[carts[c][2]][2][vars[v][2]];
}

// Real orthonormal spherical harmonics Y[l*l + l + m] for l <= lmax at the unit vector u
inline void harmonics(const int lmax, const V3& u, std::vector<double>& Y) {
  Y.assign((lmax + 1) * (lmax + 1), 0.0);
  if (lmax < 0) return;
  const double z = u[2];
  // Pb[l][m] = P_l^m(z) / (1 - z^2)^(m/2) up to sign, by the standard recurrences
  std::vector<std::vector<double>> Pb(lmax + 1, std::vector<double>(lmax + 1, 0.0));
  double dfac = 1.0;
  for (int m = 0; m <= lmax; m++) {
    if (m > 0) dfac *= 2 * m - 1;
    Pb[m][m] = dfac;
    if (m + 1 <= lmax) Pb[m + 1][m] = (2 * m + 1) * z * dfac;
    for (int l = m + 2; l <= lmax; l++)
      Pb[l][m] = ((2 * l - 1) * z * Pb[l - 1][m] - (l + m - 1) * Pb[l - 2][m]) / (l - m);
  }
  // (x + i y)^m = (1 - z^2)^(m/2) (cos m phi + i sin m phi)
  std::vector<double> cm(lmax + 1), sm(lmax + 1);
  cm[0] = 1.0;
  sm[0] = 0.0;
  for (int m = 1; m <= lmax; m++) {
    cm[m] = cm[m - 1] * u[0] - sm[m - 1] * u[1];
    sm[m] = sm[m - 1] * u[0] + cm[m - 1] * u[1];
  }
  for (int l = 0; l <= lmax; l++) {
    double ratio = 1.0;  // (l - m)! / (l + m)!
    for (int m = 0; m <= l; m++) {
      if (m > 0) ratio /= double(l + m) * (l - m + 1);
      const double N = std::sqrt((2 * l + 1) / (4.0 * M_PI) * ratio);
      if (m == 0)
        Y[l * l + l] = N * Pb[l][0];
      else {
        Y[l * l + l + m] = std::sqrt(2.0) * N * Pb[l][m] * cm[m];
        Y[l * l + l - m] = std::sqrt(2.0) * N * Pb[l][m] * sm[m];
      }
    }
  }
}

// Nodes and weights on the unit sphere for integrands exp(kappa w.axis) * (polynomial of degree
// <= degree in the components of w). axis must be a unit vector.
inline void sphere_rule(const V3& axis, const double kappa, const int degree, std::vector<V3>& dir,
                        std::vector<double>& w) {
  dir.clear();
  w.clear();
  // orthonormal frame (e1, e2, axis)
  V3 e1 = std::abs(axis[0]) < 0.9 ? V3{1.0, 0.0, 0.0} : V3{0.0, 1.0, 0.0};
  const double proj = e1[0] * axis[0] + e1[1] * axis[1] + e1[2] * axis[2];
  for (int i = 0; i < 3; i++) e1[i] -= proj * axis[i];
  const double n1 = std::sqrt(e1[0] * e1[0] + e1[1] * e1[1] + e1[2] * e1[2]);
  for (int i = 0; i < 3; i++) e1[i] /= n1;
  const V3 e2{axis[1] * e1[2] - axis[2] * e1[1], axis[2] * e1[0] - axis[0] * e1[2],
              axis[0] * e1[1] - axis[1] * e1[0]};
  // panels in s = 1 - cos t on [0, 2]
  std::vector<double> edges{0.0};
  if (2.0 * kappa <= 1.0)
    edges.push_back(2.0);
  else
    for (double s = 1.0 / kappa; edges.back() < 2.0; s *= 2.0) edges.push_back(std::min(s, 2.0));
  const int nphi = degree + 2;
  const auto& g = gl16();
  for (size_t p = 0; p + 1 < edges.size(); p++) {
    const double h = 0.5 * (edges[p + 1] - edges[p]), mid = 0.5 * (edges[p + 1] + edges[p]);
    for (int i = 0; i < 16; i++) {
      const double x = 1.0 - (mid + h * g.x[i]), st = std::sqrt(std::max(0.0, 1.0 - x * x));
      for (int k = 0; k < nphi; k++) {
        const double phi = 2.0 * M_PI * k / nphi, c = st * std::cos(phi), s = st * std::sin(phi);
        dir.push_back({c * e1[0] + s * e2[0] + x * axis[0], c * e1[1] + s * e2[1] + x * axis[1],
                       c * e1[2] + s * e2[2] + x * axis[2]});
        w.push_back(h * g.w[i] * 2.0 * M_PI / nphi);
      }
    }
  }
}

// Composite 16-point Gauss-Legendre rule on [0, R]: panels no wider than hmax, and no wider than
// c/sqrt(e) within 12/sqrt(e) of the centre of each feature (centre, exponent e)
inline void radial_rule(const std::vector<std::pair<double, double>>& features, const double R,
                        std::vector<double>& r, std::vector<double>& w, const double hmax = 0.25,
                        const double c = 1.5) {
  r.clear();
  w.clear();
  const auto& g = gl16();
  for (double r0 = 0.0; r0 < R;) {
    double h = hmax;
    for (bool changed = true; changed;) {
      changed = false;
      for (const auto& f : features) {
        const double half = 12.0 / std::sqrt(f.second), hf = c / std::sqrt(f.second);
        if (hf < h && r0 < f.first + half && r0 + h > f.first - half) {
          h = hf;
          changed = true;
        }
      }
    }
    for (int i = 0; i < 16; i++) {
      r.push_back(r0 + 0.5 * h * (1.0 + g.x[i]));
      w.push_back(0.5 * h * g.w[i]);
    }
    r0 += h;
  }
}

// Radius beyond which every primitive satisfies |d| r^(n-2+k) exp(-a r^2) < 1e-35 for k in [2, 24]
inline double extent(const std::vector<EcpPrimitive>& prims) {
  double R = 0.0;
  for (const auto& p : prims) {
    if (p.d == 0.0) continue;
    for (double r = 0.005; r < 60.0; r += 0.005) {
      const int k = (r > 1.0 ? 24 : 2) + p.n - 2;
      if (std::log(std::abs(p.d)) + k * std::log(r) - p.a * r * r > std::log(1e-35))
        R = std::max(R, r);
    }
  }
  return R + 0.01;
}

inline double U(const std::vector<EcpPrimitive>& prims, const int l, const double r) {
  double v = 0.0;
  for (const auto& p : prims)
    if (p.l == l) v += p.d * std::pow(r, p.n - 2) * std::exp(-p.a * r * r);
  return v;
}

// Compensated (Neumaier) summation: sums of millions of terms of both signs otherwise lose
// ~1e-12 relative to rounding
struct Sum {
  double s = 0.0, c = 0.0;
  void add(const double x) {
    const double t = s + x;
    c += std::abs(s) >= std::abs(x) ? (s - t) + x : (x - t) + s;
    s = t;
  }
  double value() const { return s + c; }
};

inline V3 sub(const V3& a, const V3& b) { return {a[0] - b[0], a[1] - b[1], a[2] - b[2]}; }
inline double norm(const V3& a) { return std::sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]); }
inline V3 unit_or_z(const V3& a) {
  const double n = norm(a);
  return n < 1e-12 ? V3{0.0, 0.0, 1.0} : V3{a[0] / n, a[1] / n, a[2] / n};
}

// Run f(i) for i in [0, n) on all hardware threads
template <typename F>
void parallel_for(const int n, F f) {
  const int nt = std::max(1u, std::thread::hardware_concurrency());
  std::atomic<int> next{0};
  std::vector<std::thread> pool;
  for (int t = 0; t < nt; t++)
    pool.emplace_back([&]() {
      for (int i; (i = next++) < n;) f(i);
    });
  for (auto& th : pool) th.join();
}

}  // namespace detail

/**
 * Integrals over one ECP and one shell pair, and their derivatives with respect to the two shell
 * centres up to order deriv, for every ECP and every shell pair (sa >= sb), computed in parallel on
 * construction. Derivatives with respect to the ECP centre, and the placement of the blocks in
 * molecular matrices, are left to the caller (the library's assembly code in the tests).
 */
class Engine {
 public:
  /// With bfield, also the magnetic field derivatives (get_bfield)
  Engine(const std::vector<Shell>& shells, const std::vector<Ecp>& ecps, const int deriv,
         const bool bfield = false)
      : shells_(shells),
        ecps_(ecps),
        deriv_(deriv),
        vars_(detail::variants(deriv)),
        nb_(bfield ? 3 : 0) {
    const int nv = vars_.size(), nt = nv + nb_;
    for (int a = 0; a < nv; a++)
      for (int b = 0; b < nv; b++)
        if (detail::order(vars_[a]) + detail::order(vars_[b]) <= deriv) vpairs_.push_back({a, b});
    // the field factors, variants nv + k, on either side of an undifferentiated function
    for (int k = 0; k < nb_; k++) vpairs_.insert(vpairs_.end(), {{nv + k, 0}, {0, nv + k}});
    vpair_index_.assign(nt * nt, -1);
    for (size_t k = 0; k < vpairs_.size(); k++)
      vpair_index_[vpairs_[k].first * nt + vpairs_[k].second] = k;
    for (const Shell& sh : shells_) carts_.push_back(detail::cartesians(sh.l));
    for (size_t e = 0; e < ecps_.size(); e++) blocks_.push_back(compute_ecp(e));
  }

  int deriv() const { return deriv_; }

  /// Index of the derivative variant with respect to coordinate p (first) or p and q (second) of
  /// one centre, for get(); 0 is no derivative
  static int var1(const int p) { return detail::var1(p); }
  static int var2(const int p, const int q) { return detail::var2(p, q); }

  /// d^va/dA d^vb/dB <ca|U_e|cb> for shells sa >= sb, Cartesian components ca, cb (in library
  /// order) and derivative variants va (on A) and vb (on B) of total order <= deriv
  double get(const int e, const int sa, const int sb, const int va, const int vb, const int ca,
             const int cb) const {
    const int nca = carts_[sa].size(), ncb = carts_[sb].size();
    const int k = vpair_index_[va * (vars_.size() + nb_) + vb];
    return blocks_[e][sa * (sa + 1) / 2 + sb][(k * nca + ca) * ncb + cb];
  }

  /// Component k of the magnetic field derivative of <ca|U_e|cb>, for shells sa >= sb (the engine
  /// must be constructed with bfield)
  double get_bfield(const int e, const int sa, const int sb, const int k, const int ca,
                    const int cb) const {
    const int nv = vars_.size();
    return get(e, sa, sb, nv + k, 0, ca, cb) - get(e, sa, sb, 0, nv + k, ca, cb);
  }

 private:
  using V3 = detail::V3;

  // Blocks T[vpair][ca][cb] for every shell pair (sa >= sb, in the order of pair_index) over ECP e
  std::vector<std::vector<double>> compute_ecp(const size_t e) const {
    using namespace detail;
    const Ecp& ecp = ecps_[e];
    const auto& vars = vars_;
    const auto& vpairs = vpairs_;
    const auto& carts = carts_;
    const auto& shells = shells_;
    const int deriv = deriv_;
    const int nv = vars.size(), nb = nb_, nt = nv + nb;
    const int nshell = shells.size();
    // polynomial degree added to the basis functions by the derivatives or field factors
    const int extra = std::max(deriv, nb > 0 ? 1 : 0);

    int L = 0;
    for (const auto& p : ecp.prims) L = std::max(L, p.l);
    std::vector<EcpPrimitive> local, semilocal;
    for (const auto& p : ecp.prims) (p.l == L ? local : semilocal).push_back(p);
    bool has_local = false;
    for (const auto& p : local) has_local = has_local || p.d != 0.0;
    std::vector<V3> centre(nshell);  // shell centres relative to the ECP
    for (int s = 0; s < nshell; s++) centre[s] = sub(shells[s].centre, ecp.centre);
    // Field factor of shell s at pt (relative to the ECP): [(A - C) x r]_k for the absolute
    // position r = pt + C (so exactly 0 for a shell on the ECP); appended to val as variants nv + k
    auto add_field = [&](const int s, const V3& pt, const int nc, double* val) {
      if (nb == 0) return;
      const V3& D = centre[s];
      const V3 r{pt[0] + ecp.centre[0], pt[1] + ecp.centre[1], pt[2] + ecp.centre[2]};
      const double m[3] = {D[1] * r[2] - D[2] * r[1], D[2] * r[0] - D[0] * r[2],
                           D[0] * r[1] - D[1] * r[0]};
      for (int k = 0; k < 3; k++)
        for (int c = 0; c < nc; c++) val[(nv + k) * nc + c] = m[k] * val[c];
    };
    const int nlm = L * L;  // harmonics l < L

    // Semi-local part: projections of every shell on a radial grid shared by all shells
    std::vector<double> rs, ws;
    std::vector<std::vector<double>> Wl;         // w r^2 U_l(r) for each l < L
    std::vector<std::vector<double>> P(nshell);  // [node][v][c][lm], v < nt
    if (L > 0) {
      std::vector<std::pair<double, double>> features;
      for (int s = 0; s < nshell; s++)
        for (const double a : shells[s].exps) features.push_back({norm(centre[s]), a});
      for (const auto& p : semilocal) features.push_back({0.0, p.a});
      radial_rule(features, extent(semilocal), rs, ws);
      Wl.assign(L, std::vector<double>(rs.size()));
      for (int l = 0; l < L; l++)
        for (size_t i = 0; i < rs.size(); i++)
          Wl[l][i] = ws[i] * rs[i] * rs[i] * U(semilocal, l, rs[i]);
      parallel_for(nshell, [&](const int s) {
        const Shell& sh = shells[s];
        const int nc = carts[s].size();
        const double A = norm(centre[s]);
        const V3 axis = unit_or_z(centre[s]);
        P[s].assign(rs.size() * nt * nc * nlm, 0.0);
        std::vector<Sum> acc(nt * nc * nlm);
        std::vector<V3> dir;
        std::vector<double> dw, Y, poly, val(nt * nc);
        for (size_t i = 0; i < rs.size(); i++) {
          const double r = rs[i];
          double kappa = 0.0;
          for (const double a : sh.exps)
            if (a * (r - A) * (r - A) < 100.0) kappa = std::max(kappa, 2.0 * a * A * r);
          sphere_rule(axis, kappa, sh.l + extra + L - 1, dir, dw);
          std::fill(acc.begin(), acc.end(), Sum());
          for (size_t k = 0; k < dir.size(); k++) {
            const V3 pt{r * dir[k][0], r * dir[k][1], r * dir[k][2]};
            const V3 t = sub(pt, centre[s]);
            const double t2 = t[0] * t[0] + t[1] * t[1] + t[2] * t[2];
            std::fill(val.begin(), val.end(), 0.0);
            for (size_t m = 0; m < sh.exps.size(); m++) {
              if (sh.exps[m] * t2 > NEGLIGIBLE) continue;
              const double e = std::exp(-sh.exps[m] * t2);
              primitive_poly(sh.exps[m], sh.coefs[m] * e, carts[s], vars, t, poly);
              for (int x = 0; x < nv * nc; x++) val[x] += poly[x];
            }
            add_field(s, pt, nc, val.data());
            harmonics(L - 1, dir[k], Y);
            for (int x = 0; x < nt * nc; x++) {
              const double f = dw[k] * val[x];
              for (int lm = 0; lm < nlm; lm++) acc[x * nlm + lm].add(f * Y[lm]);
            }
          }
          double* Pi = &P[s][i * nt * nc * nlm];
          for (int x = 0; x < nt * nc * nlm; x++) Pi[x] = acc[x].value();
        }
      });
    }

    // Shell pairs (sa >= sb): blocks T[vpair][ca][cb]
    std::vector<std::pair<int, int>> pairs;
    for (int sa = 0; sa < nshell; sa++)
      for (int sb = 0; sb <= sa; sb++) pairs.push_back({sa, sb});
    std::vector<std::vector<double>> blocks(pairs.size());
    parallel_for(pairs.size(), [&](const int ip) {
      const int sa = pairs[ip].first, sb = pairs[ip].second;
      const int nca = carts[sa].size(), ncb = carts[sb].size();
      std::vector<double>& T = blocks[ip];
      std::vector<Sum> acc(vpairs.size() * nca * ncb);
      auto at = [&](const int k, const int ca, const int cb) -> Sum& {
        return acc[(k * nca + ca) * ncb + cb];
      };
      auto store = [&]() {
        T.resize(acc.size());
        for (size_t x = 0; x < acc.size(); x++) T[x] = acc[x].value();
      };
      // semi-local
      for (size_t i = 0; i < rs.size(); i++) {
        const double* Pa = &P[sa][i * nt * nca * nlm];
        const double* Pb = &P[sb][i * nt * ncb * nlm];
        for (size_t k = 0; k < vpairs.size(); k++)
          for (int ca = 0; ca < nca; ca++)
            for (int cb = 0; cb < ncb; cb++) {
              const double* pa = &Pa[(vpairs[k].first * nca + ca) * nlm];
              const double* pb = &Pb[(vpairs[k].second * ncb + cb) * nlm];
              double sum = 0.0;
              for (int l = 0; l < L; l++) {
                double sl = 0.0;
                for (int lm = l * l; lm < (l + 1) * (l + 1); lm++) sl += pa[lm] * pb[lm];
                sum += Wl[l][i] * sl;
              }
              at(k, ca, cb).add(sum);
            }
      }
      // local, per primitive pair, about the Gaussian product centre
      if (!has_local) {
        store();
        return;
      }
      const Shell &shA = shells[sa], &shB = shells[sb];
      const V3 &A = centre[sa], &B = centre[sb];
      const double RL = extent(local);
      std::vector<double> rl, wl, pa, pb;
      std::vector<V3> dir;
      std::vector<double> dw;
      for (size_t ia = 0; ia < shA.exps.size(); ia++)
        for (size_t ib = 0; ib < shB.exps.size(); ib++) {
          const double a = shA.exps[ia], b = shB.exps[ib], p = a + b;
          const V3 AB = sub(A, B);
          const double AB2 = AB[0] * AB[0] + AB[1] * AB[1] + AB[2] * AB[2];
          if (a * b / p * AB2 > NEGLIGIBLE) continue;
          const V3 Pc{(a * A[0] + b * B[0]) / p, (a * A[1] + b * B[1]) / p,
                      (a * A[2] + b * B[2]) / p};
          const double Pm = norm(Pc);
          const V3 axis = unit_or_z(Pc);
          std::vector<std::pair<double, double>> features{{Pm, p}};
          for (const auto& u : local) features.push_back({0.0, u.a});
          radial_rule(features, RL, rl, wl);
          for (size_t i = 0; i < rl.size(); i++) {
            const double r = rl[i];
            // |r w - A|^2 a + |r w - B|^2 b >= p (r - |P|)^2 + ab/p |A - B|^2
            if (p * (r - Pm) * (r - Pm) + a * b / p * AB2 > NEGLIGIBLE) continue;
            const double Wr = wl[i] * r * r * U(local, L, r);
            if (Wr == 0.0) continue;
            sphere_rule(axis, 2.0 * p * Pm * r, shA.l + shB.l + extra, dir, dw);
            for (size_t k = 0; k < dir.size(); k++) {
              const V3 pt{r * dir[k][0], r * dir[k][1], r * dir[k][2]};
              const V3 ta = sub(pt, A), tb = sub(pt, B);
              const double arg = a * (ta[0] * ta[0] + ta[1] * ta[1] + ta[2] * ta[2]) +
                                 b * (tb[0] * tb[0] + tb[1] * tb[1] + tb[2] * tb[2]);
              if (arg > NEGLIGIBLE) continue;
              const double e = std::exp(-arg);
              primitive_poly(a, shA.coefs[ia], carts[sa], vars, ta, pa);
              primitive_poly(b, shB.coefs[ib], carts[sb], vars, tb, pb);
              pa.resize(nt * nca);
              pb.resize(nt * ncb);
              add_field(sa, pt, nca, pa.data());
              add_field(sb, pt, ncb, pb.data());
              const double f = Wr * dw[k] * e;
              for (size_t q = 0; q < vpairs.size(); q++)
                for (int ca = 0; ca < nca; ca++) {
                  const double fa = f * pa[vpairs[q].first * nca + ca];
                  const double* pbv = &pb[vpairs[q].second * ncb];
                  for (int cb = 0; cb < ncb; cb++) at(q, ca, cb).add(fa * pbv[cb]);
                }
            }
          }
        }
      store();
    });

    return blocks;
  }

  std::vector<Shell> shells_;
  std::vector<Ecp> ecps_;
  int deriv_;
  std::vector<std::array<int, 3>> vars_;
  int nb_;  // number of field-factor variants (3 with bfield, else 0), after vars_
  std::vector<std::pair<int, int>> vpairs_;
  std::vector<int> vpair_index_;
  std::vector<std::vector<std::array<int, 3>>> carts_;
  std::vector<std::vector<std::vector<double>>> blocks_;  // [ecp][pair][vpair][ca][cb]
};

}  // namespace intref

#endif
