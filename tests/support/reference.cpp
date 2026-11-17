#include "reference.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <map>

#include "assembly.hpp"

namespace testsupport {

using libecpint::ECP;
using libecpint::GaussianShell;
using libecpint::TwoIndex;

std::vector<intref::Shell> reference_shells(const libecpint::ECPIntegrator& f) {
  std::vector<intref::Shell> shells;
  for (const GaussianShell& s : f.shells) {
    intref::Shell r;
    r.centre = {s.center()[0], s.center()[1], s.center()[2]};
    r.l = s.am();
    for (int i = 0; i < s.nprimitive(); i++) {
      r.exps.push_back(s.exp(i));
      r.coefs.push_back(s.coef(i));
    }
    r.atom = s.atom_id;
    shells.push_back(r);
  }
  return shells;
}

std::vector<intref::Ecp> reference_ecps(const libecpint::ECPIntegrator& f) {
  std::vector<intref::Ecp> ecps;
  for (int i = 0; i < f.ecps.getN(); i++) {
    const ECP& U = f.ecps.getECP(i);
    intref::Ecp r;
    r.centre = {U.center()[0], U.center()[1], U.center()[2]};
    for (const auto& g : U.gaussians) r.prims.push_back({g.n + 2, g.l, g.a, g.d});  // n holds n - 2
    r.atom = U.atom_id;
    ecps.push_back(r);
  }
  return ecps;
}

namespace {
// Integrals and derivatives assembled by the library's loops from the blocks
// q(e, sa, sb, va, vb, ca, cb) (see intref::Engine::get), adding absolute values if absolute
template <typename Q>
Integrals assemble(const libecpint::ECPIntegrator& f, const int deriv, Q block_of,
                   const bool absolute) {
  const auto screen = libecpint::ShellPairScreen::none(f.shells.size(), f.ecps.getN());
  // The assembly loops pass references into f.shells and f.ecps
  auto shell_index = [&](const GaussianShell& s) { return int(&s - f.shells.data()); };
  auto ecp_index = [&](const ECP& U) { return int(&U - &f.ecps.getECP(0)); };
  using V = intref::Engine;

  Integrals res;
  const int ncart = f.ncart, natoms = f.natoms;
  res.ints.assign(ncart, ncart, 0.0);
  libecpint::integrals_loop(
      f.shells, f.ecps, screen,
      [&](const ECP& U, const GaussianShell& A, const GaussianShell& B, TwoIndex<double>& v) {
        const int e = ecp_index(U), sa = shell_index(A), sb = shell_index(B);
        v.assign(A.ncartesian(), B.ncartesian(), 0.0);
        for (int ca = 0; ca < A.ncartesian(); ca++)
          for (int cb = 0; cb < B.ncartesian(); cb++) v(ca, cb) = block_of(e, sa, sb, 0, 0, ca, cb);
      },
      [&](const int n1, const int n2, const TwoIndex<double>& block) {
        for (int i = n1; i < n1 + block.dims[0]; i++)
          for (int j = n2; j < n2 + block.dims[1]; j++) {
            res.ints(i, j) = absolute ? std::abs(block(i - n1, j - n2)) : block(i - n1, j - n2);
            res.ints(j, i) = res.ints(i, j);
          }
      });

  // As ECPIntegrator::compute_first_derivs / compute_second_derivs
  auto add_to = [absolute](std::vector<TwoIndex<double>>& mats) {
    return [&mats, absolute](const int i, const int k, const int l, const double v,
                             const bool mirror) {
      TwoIndex<double>& M = mats[i];
      M(k, l) += absolute ? std::abs(v) : v;
      if (mirror) M(l, k) = M(k, l);
    };
  };

  if (deriv >= 1) {
    res.first.assign(3 * natoms, TwoIndex<double>(ncart, ncart, 0.0));
    libecpint::first_derivs_loop(
        f.shells, f.ecps, screen,
        [&](const ECP& U, const GaussianShell& A, const GaussianShell& B,
            std::array<TwoIndex<double>, 9>& out) {
          const int e = ecp_index(U), sa = shell_index(A), sb = shell_index(B);
          const int na = A.ncartesian(), nb = B.ncartesian();
          std::array<TwoIndex<double>, 3> QA, QB;
          for (int p = 0; p < 3; p++) {
            QA[p].assign(na, nb, 0.0);
            QB[p].assign(na, nb, 0.0);
            for (int ca = 0; ca < na; ca++)
              for (int cb = 0; cb < nb; cb++) {
                QA[p](ca, cb) = block_of(e, sa, sb, V::var1(p), 0, ca, cb);
                QB[p](ca, cb) = block_of(e, sa, sb, 0, V::var1(p), ca, cb);
              }
          }
          libecpint::pack_first_derivatives(QA, QB, libecpint::on_ecp(A, U),
                                            libecpint::on_ecp(B, U), na, nb, out);
        },
        add_to(res.first));
  }

  if (deriv >= 2) {
    res.second.assign(3 * natoms * (3 * natoms + 1) / 2, TwoIndex<double>(ncart, ncart, 0.0));
    const int sym[6][2] = {{0, 0}, {0, 1}, {0, 2}, {1, 1}, {1, 2}, {2, 2}};
    libecpint::second_derivs_loop(
        f.shells, f.ecps, natoms, screen,
        [&](const ECP& U, const GaussianShell& A, const GaussianShell& B,
            std::array<TwoIndex<double>, 45>& out) {
          const int e = ecp_index(U), sa = shell_index(A), sb = shell_index(B);
          const int na = A.ncartesian(), nb = B.ncartesian();
          std::array<TwoIndex<double>, 6> QAA, QBB;
          std::array<TwoIndex<double>, 9> QAB;
          for (int k = 0; k < 6; k++) {
            const int vk = V::var2(sym[k][0], sym[k][1]);
            QAA[k].assign(na, nb, 0.0);
            QBB[k].assign(na, nb, 0.0);
            for (int ca = 0; ca < na; ca++)
              for (int cb = 0; cb < nb; cb++) {
                QAA[k](ca, cb) = block_of(e, sa, sb, vk, 0, ca, cb);
                QBB[k](ca, cb) = block_of(e, sa, sb, 0, vk, ca, cb);
              }
          }
          for (int p = 0; p < 3; p++)
            for (int q = 0; q < 3; q++) {
              TwoIndex<double>& M = QAB[3 * p + q];
              M.assign(na, nb, 0.0);
              for (int ca = 0; ca < na; ca++)
                for (int cb = 0; cb < nb; cb++)
                  M(ca, cb) = block_of(e, sa, sb, V::var1(p), V::var1(q), ca, cb);
            }
          libecpint::pack_second_derivatives(QAA, QAB, QBB, libecpint::on_ecp(A, U),
                                             libecpint::on_ecp(B, U), na, nb, out);
        },
        add_to(res.second));
  }
  return res;
}

}  // namespace

Integrals reference_integrals(const libecpint::ECPIntegrator& f, const int deriv,
                              Integrals* scale) {
  const intref::Engine engine(reference_shells(f), reference_ecps(f), deriv);
  auto get = [&](const int e, const int sa, const int sb, const int va, const int vb, const int ca,
                 const int cb) { return engine.get(e, sa, sb, va, vb, ca, cb); };
  if (scale) {
    // The size of the terms of a derivative of order o on shell s about ECP e: the shell raised to
    // l + o (the leading term of the derivative recurrence) with coefficients
    //   |c| (2a)^o (1 + d sqrt(2a))^(l + o),
    // d the distance of s from the ECP, the last factor the cancellation of the expansion about
    // the ECP centre for that primitive, and nothing cancelling between primitives (computed by the
    // library: only a magnitude)
    auto raised = [&](const int e, const int s, const int o) {
      const GaussianShell& sh = f.shells[s];
      const ECP& U = f.ecps.getECP(e);
      double d2 = 0.0;
      for (int i = 0; i < 3; i++) d2 += std::pow(sh.center()[i] - U.center()[i], 2);
      GaussianShell r(sh.center(), sh.am() + o);
      for (int i = 0; i < sh.nprimitive(); i++) {
        const double a = sh.exp(i);
        r.addPrim(a, std::abs(sh.coef(i)) * std::pow(2.0 * a, o) *
                         std::pow(1.0 + std::sqrt(2.0 * a * d2), sh.am() + o));
      }
      return r;
    };
    std::map<std::array<int, 5>, double> terms;  // (e, sa, sb, oA, oB) -> max |block|
    auto term_scale = [&](const int e, const int sa, const int sb, const int oA, const int oB) {
      const std::array<int, 5> key{e, sa, sb, oA, oB};
      auto it = terms.find(key);
      if (it != terms.end()) return it->second;
      TwoIndex<double> v;
      f.ecpint->compute_shell_pair(f.ecps.getECP(e), raised(e, sa, oA), raised(e, sb, oB), v);
      double m = 0.0;
      for (const double x : v.data) m = std::max(m, std::abs(x));
      return terms[key] = m;
    };
    auto order = [](const int v) { return v == 0 ? 0 : (v <= 3 ? 1 : 2); };  // variants(2)
    *scale = assemble(
        f, deriv,
        [&](const int e, const int sa, const int sb, const int va, const int vb, const int,
            const int) { return term_scale(e, sa, sb, order(va), order(vb)); },
        true);
  }
  return assemble(f, deriv, get, false);
}

Integrals library_integrals(libecpint::ECPIntegrator& f, const int deriv) {
  Integrals res;
  f.compute_integrals();
  res.ints = f.integrals;
  if (deriv >= 1) {
    f.first_derivs.clear();
    f.compute_first_derivs();
    res.first = f.first_derivs;
  }
  if (deriv >= 2) {
    f.second_derivs.clear();
    f.compute_second_derivs();
    res.second = f.second_derivs;
  }
  return res;
}

}  // namespace testsupport
