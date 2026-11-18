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

#ifndef ASSEMBLY_HEAD
#define ASSEMBLY_HEAD

#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

#include "ecp.hpp"
#include "gshell.hpp"
#include "multiarr.hpp"
#include "screening.hpp"

/// returns the index for the ij-th second derivative for a system with N atoms
#define H_START(i, j, N) (9 * j + 3 * (3 * N - 1) * i - (9 * i * (i + 1)) / 2 - 3)

namespace libecpint {

/**
 * Bookkeeping that turns per-(ECP, shell pair) results into whole-molecule integrals and
 * derivatives: packing the shell-centre derivatives of one ECP block into derivatives with respect
 * to the three centres (using translational invariance), and the loops over shell pairs and ECPs
 * that place them. None of this affects accuracy, so it is shared by ECPIntegrator and the test
 * reference, and tested on its own.
 */

/// @return true if shell is on the centre of U (to 1e-6 in each coordinate, summed)
inline bool on_ecp(const GaussianShell& shell, const ECP& U) {
  const double* A = shell.center();
  const double* C = U.center();
  return std::abs(A[0] - C[0]) + std::abs(A[1] - C[1]) + std::abs(A[2] - C[2]) <= 1e-6;
}

/**
 * Packs the first derivatives of one ECP block <A|U_C|B> with respect to the shell centres into
 * derivatives with respect to all three centres, results = [Ax, Ay, Az, Bx, By, Bz, Cx, Cy, Cz],
 * each ncartA x ncartB. d/dC = -(d/dA + d/dB) by translational invariance. A shell on the ECP
 * moves with it, so its derivative is not needed (and is not read): if A is on the ECP,
 * d/dA = -d/dB and d/dC = 0 (the atom's total derivative is then d/dA + d/dC), and likewise for B.
 *
 * @param QA - d/dA_q <A|U|B> for q = x, y, z (ncartA x ncartB); read only if !onA
 * @param QB - d/dB_q <A|U|B> for q = x, y, z (ncartA x ncartB); read only if !onB
 * @param onA, onB - whether shell A, B is on the ECP centre (on_ecp)
 */
void pack_first_derivatives(const std::array<TwoIndex<double>, 3>& QA,
                            const std::array<TwoIndex<double>, 3>& QB, bool onA, bool onB,
                            int ncartA, int ncartB, std::array<TwoIndex<double>, 9>& results);

/**
 * Packs the second derivatives of one ECP block into results = [AA, AB, AC, BB, BC, CC], with
 * components [xx, xy, xz, yy, yz, zz] for AA, BB and CC, and [xx, xy, xz, yx, yy, yz, zx, zy, zz]
 * (first index on the first centre) for AB, AC and BC: 6 + 9 + 9 + 6 + 9 + 6 = 45 matrices, each
 * ncartA x ncartB. The C blocks follow from translational invariance. As for the first derivatives,
 * a shell on the ECP moves with it: if A is on the ECP, AA = BB = QBB, AB = -QBB and the C blocks
 * are 0, and likewise for B.
 *
 * @param QAA - d2/dA_p dA_q <A|U|B>, [xx, xy, xz, yy, yz, zz]; read only if !onA
 * @param QAB - d2/dA_p dB_q <A|U|B>, [xx, xy, xz, yx, yy, yz, zx, zy, zz]; read if !onA && !onB
 * @param QBB - d2/dB_p dB_q <A|U|B>, [xx, xy, xz, yy, yz, zz]; read only if !onB
 * @param onA, onB - whether shell A, B is on the ECP centre (on_ecp)
 */
void pack_second_derivatives(const std::array<TwoIndex<double>, 6>& QAA,
                             const std::array<TwoIndex<double>, 9>& QAB,
                             const std::array<TwoIndex<double>, 6>& QBB, bool onA, bool onB,
                             int ncartA, int ncartB, std::array<TwoIndex<double>, 45>& results);

/**
 * Screening of (shell, ECP) and (shell pair, ECP) combinations in the assembly loops: the shell
 * s is skipped against ECP i if bounds[s * necp + i] <= shell_tolerance, and the pair (s1, s2) if
 * bounds[s1 * necp + i] * bounds[s2 * necp + i] <= tolerance.
 */
struct ShellPairScreen {
  std::vector<double> bounds;
  double shell_tolerance, tolerance;

  /// @return a screen that skips nothing, for nshells shells and necp ECPs
  static ShellPairScreen none(const size_t nshells, const int necp) {
    return {std::vector<double>(nshells * necp, 1.0), -1.0, -1.0};
  }

  /// @return the ECPs (index, bound) not screened out for shell s
  std::vector<std::pair<int, double>> ecps_for(const int s, const int necp) const {
    return screening::screen_shell_ecps(bounds, s, necp, shell_tolerance);
  }
};

/**
 * Loops over the shell pairs (s2 <= s1) and ECPs for the integrals. For each pair with at least
 * one ECP not screened out, the blocks kernel(U, shellA, shellB, block) are summed over the ECPs
 * and passed to store(n1, n2, sum), n1 and n2 being the offsets of the shells' functions. Pairs
 * screened out against every ECP are not passed (most pairs, in a large system).
 */
template <typename Kernel, typename Store>
void integrals_loop(const std::vector<GaussianShell>& shells, const ECPBasis& ecps,
                    const ShellPairScreen& screen, Kernel&& kernel, Store&& store) {
  TwoIndex<double> tempValues, shellPairInts;
  const int nshells = shells.size(), necp = ecps.getN();

  int n1 = 0;
  for (auto s1 = 0; s1 < nshells; ++s1) {
    const GaussianShell& shellA = shells[s1];
    int ncartA = shellA.ncartesian();
    auto ns = screen.ecps_for(s1, necp);

    if (ns.size() > 0) {
      int n2 = 0;
      for (auto s2 = 0; s2 <= s1; ++s2) {
        const GaussianShell& shellB = shells[s2];
        int ncartB = shellB.ncartesian();

        bool computed = false;
        for (auto& [i, sbA] : ns) {
          const ECP& U = ecps.getECP(i);
          if (sbA * screen.bounds[s2 * necp + i] <= screen.tolerance) continue;

          if (!computed) shellPairInts.assign(ncartA, ncartB, 0.0);
          computed = true;
          kernel(U, shellA, shellB, tempValues);
          shellPairInts.add(tempValues);
        }
        if (computed) store(n1, n2, shellPairInts);
        n2 += ncartB;
      }
    }
    n1 += ncartA;
  }
}

/**
 * As integrals_loop, for the magnetic field derivatives: kernel(U, shellA, shellB, block) gives
 * the [x, y, z] components of one ECP block (as from compute_shell_pair_derivative_Bfield), and
 * store(n1, n2, sum) receives their sum over the ECPs. The matrices are antisymmetric, so for
 * s2 < s1 the block also stands for minus its transpose.
 */
template <typename Kernel, typename Store>
void bfield_derivs_loop(const std::vector<GaussianShell>& shells, const ECPBasis& ecps,
                        const ShellPairScreen& screen, Kernel&& kernel, Store&& store) {
  std::array<TwoIndex<double>, 3> tempValues, shellPairInts;
  const int nshells = shells.size(), necp = ecps.getN();

  int n1 = 0;
  for (auto s1 = 0; s1 < nshells; ++s1) {
    const GaussianShell& shellA = shells[s1];
    int ncartA = shellA.ncartesian();
    auto ns = screen.ecps_for(s1, necp);

    int n2 = 0;
    for (auto s2 = 0; s2 <= s1; ++s2) {
      const GaussianShell& shellB = shells[s2];
      int ncartB = shellB.ncartesian();

      bool computed = false;
      for (auto& [i, sbA] : ns) {
        const ECP& U = ecps.getECP(i);
        if (sbA * screen.bounds[s2 * necp + i] <= screen.tolerance) continue;

        if (!computed)
          for (auto& m : shellPairInts) m.assign(ncartA, ncartB, 0.0);
        computed = true;
        kernel(U, shellA, shellB, tempValues);
        for (int k = 0; k < 3; k++) shellPairInts[k].add(tempValues[k]);
      }
      if (computed) store(n1, n2, shellPairInts);
      n2 += ncartB;
    }
    n1 += ncartA;
  }
}

/**
 * Loops over the shell pairs (s2 <= s1) and ECPs for the first derivatives, with
 * kernel(U, shellA, shellB, results) giving the packed derivatives of one ECP block (as from
 * pack_first_derivatives). Each contribution goes to add(index, k, l, value, mirror): value is the
 * contribution to element (k, l) of derivative matrix index (3 * atom + coordinate), and mirror is
 * true if it also stands for (l, k) (the matrices are symmetric and only one triangle of shell
 * pairs is computed). Atoms are the atom_id of the shells and ECPs.
 */
template <typename Kernel, typename Add>
void first_derivs_loop(const std::vector<GaussianShell>& shells, const ECPBasis& ecps,
                       const ShellPairScreen& screen, Kernel&& kernel, Add&& add) {
  std::array<TwoIndex<double>, 9> tempValues;
  const int nshells = shells.size(), necp = ecps.getN();

  int n1 = 0;
  int Aix, Bix, Cix;
  for (auto s1 = 0; s1 < nshells; ++s1) {
    const GaussianShell& shellA = shells[s1];
    int ncartA = shellA.ncartesian();
    Aix = shellA.atom_id;

    auto ns = screen.ecps_for(s1, necp);

    int n2 = 0;
    for (auto s2 = 0; s2 <= s1; ++s2) {
      const GaussianShell& shellB = shells[s2];
      int ncartB = shellB.ncartesian();
      Bix = shellB.atom_id;

      for (auto& [i, sbA] : ns) {
        const ECP& U = ecps.getECP(i);
        Cix = U.atom_id;
        if (sbA * screen.bounds[s2 * necp + i] <= screen.tolerance) continue;

        kernel(U, shellA, shellB, tempValues);

        // work out where to put them
        for (int n = 0; n < 3; n++) {
          for (int k = n1; k < n1 + ncartA; k++) {
            for (int l = n2; l < n2 + ncartB; l++) {
              add(3 * Aix + n, k, l, tempValues[n](k - n1, l - n2), s2 < s1);
              add(3 * Bix + n, k, l, tempValues[n + 3](k - n1, l - n2), s2 < s1);
              add(3 * Cix + n, k, l, tempValues[n + 6](k - n1, l - n2), s2 < s1);
            }
          }
        }
      }

      n2 += ncartB;
    }

    n1 += ncartA;
  }
}

/**
 * As first_derivs_loop, for the second derivatives: kernel gives the packed second derivatives of
 * one ECP block (as from pack_second_derivatives), and index is in the packed order of
 * ECPIntegrator::second_derivs (see H_START) for natoms atoms.
 */
template <typename Kernel, typename Add>
void second_derivs_loop(const std::vector<GaussianShell>& shells, const ECPBasis& ecps,
                        const int natoms, const ShellPairScreen& screen, Kernel&& kernel,
                        Add&& add) {
  std::array<TwoIndex<double>, 45> tempValues;
  const int nshells = shells.size(), necp = ecps.getN();

  int n1 = 0;
  int Aix, Bix, Cix;
  int saa, sab, sac, sbb, sbc, scc;
  int ixes[6] = {0, 1, 2, 4, 5, 8};
  int back_ixes[6] = {0, 3, 6, 4, 7, 8};
  int jxes[9] = {0, 3, 6, 1, 4, 7, 2, 5, 8};
  for (auto s1 = 0; s1 < nshells; ++s1) {
    const GaussianShell& shellA = shells[s1];
    int ncartA = shellA.ncartesian();
    Aix = shellA.atom_id;
    auto ns = screen.ecps_for(s1, necp);

    int n2 = 0;
    for (auto s2 = 0; s2 <= s1; ++s2) {
      const GaussianShell& shellB = shells[s2];
      int ncartB = shellB.ncartesian();
      Bix = shellB.atom_id;

      saa = H_START(Aix, Aix, natoms) + 3;
      sbb = H_START(Bix, Bix, natoms) + 3;
      sab = H_START(std::min(Aix, Bix), std::max(Aix, Bix), natoms);
      sab = Aix == Bix ? sab + 3 : sab;

      for (auto& [i, sbA] : ns) {
        const ECP& U = ecps.getECP(i);
        Cix = U.atom_id;
        if (sbA * screen.bounds[s2 * necp + i] <= screen.tolerance) continue;

        kernel(U, shellA, shellB, tempValues);

        // work out where to put them
        scc = H_START(Cix, Cix, natoms) + 3;
        sac = H_START(std::min(Aix, Cix), std::max(Aix, Cix), natoms);
        sac = Aix == Cix ? sac + 3 : sac;
        sbc = H_START(std::min(Bix, Cix), std::max(Bix, Cix), natoms);
        sbc = Bix == Cix ? sbc + 3 : sbc;

        if ((Aix == Cix) || (Bix == Cix)) {
          if (Bix != Aix) {
            // two distinct atoms
            // only need to worry about AA, AB, and BB blocks
            for (int n = 0; n < 6; n++) {
              for (int k = n1; k < n1 + ncartA; k++) {
                for (int l = n2; l < n2 + ncartB; l++) {
                  add(saa + n, k, l, tempValues[n](k - n1, l - n2), s1 != s2);
                  add(sbb + n, k, l, tempValues[n + 24](k - n1, l - n2), s1 != s2);
                }
              }
            }

            for (int n = 0; n < 9; n++) {
              for (int k = n1; k < n1 + ncartA; k++) {
                for (int l = n2; l < n2 + ncartB; l++) {
                  if (Aix > Bix) {
                    add(sab + n, k, l, tempValues[jxes[n] + 6](k - n1, l - n2), s1 != s2);
                  } else {
                    add(sab + n, k, l, tempValues[n + 6](k - n1, l - n2), s1 != s2);
                  }
                }
              }
            }
          }  // else everything is zero
        } else if (Aix == Bix) {
          // two distinct atoms, need to worry about everything
          for (int n = 0; n < 6; n++) {
            for (int k = n1; k < n1 + ncartA; k++) {
              for (int l = n2; l < n2 + ncartB; l++) {
                add(saa + n, k, l, tempValues[n](k - n1, l - n2), s1 != s2);            // aa
                add(saa + n, k, l, tempValues[n + 24](k - n1, l - n2), s1 != s2);       // bb = aa
                add(scc + n, k, l, tempValues[n + 39](k - n1, l - n2), s1 != s2);       // cc
                add(saa + n, k, l, tempValues[ixes[n] + 6](k - n1, l - n2), s1 != s2);  // ab = aa
                add(saa + n, k, l, tempValues[back_ixes[n] + 6](k - n1, l - n2),
                    s1 != s2);  // ba = aa
              }
            }
          }

          for (int n = 0; n < 9; n++) {
            for (int k = n1; k < n1 + ncartA; k++) {
              for (int l = n2; l < n2 + ncartB; l++) {
                if (Aix > Cix) {
                  add(sac + n, k, l, tempValues[jxes[n] + 15](k - n1, l - n2), s1 != s2);
                  add(sac + n, k, l, tempValues[jxes[n] + 30](k - n1, l - n2),
                      s1 != s2);  // bc = ac
                } else {
                  add(sac + n, k, l, tempValues[n + 15](k - n1, l - n2), s1 != s2);
                  add(sac + n, k, l, tempValues[n + 30](k - n1, l - n2), s1 != s2);  // bc = ac
                }
              }
            }
          }
        } else {
          for (int n = 0; n < 6; n++) {
            for (int k = n1; k < n1 + ncartA; k++) {
              for (int l = n2; l < n2 + ncartB; l++) {
                add(saa + n, k, l, tempValues[n](k - n1, l - n2), s1 != s2);
                add(sbb + n, k, l, tempValues[n + 24](k - n1, l - n2), s1 != s2);
                add(scc + n, k, l, tempValues[n + 39](k - n1, l - n2), s1 != s2);
              }
            }
          }

          for (int n = 0; n < 9; n++) {
            for (int k = n1; k < n1 + ncartA; k++) {
              for (int l = n2; l < n2 + ncartB; l++) {
                if (Aix > Bix) {
                  add(sab + n, k, l, tempValues[jxes[n] + 6](k - n1, l - n2), s1 != s2);
                } else {
                  add(sab + n, k, l, tempValues[n + 6](k - n1, l - n2), s1 != s2);
                }

                if (Aix > Cix) {
                  add(sac + n, k, l, tempValues[jxes[n] + 15](k - n1, l - n2), s1 != s2);
                } else {
                  add(sac + n, k, l, tempValues[n + 15](k - n1, l - n2), s1 != s2);
                }

                if (Bix > Cix) {
                  add(sbc + n, k, l, tempValues[jxes[n] + 30](k - n1, l - n2), s1 != s2);
                } else {
                  add(sbc + n, k, l, tempValues[n + 30](k - n1, l - n2), s1 != s2);
                }
              }
            }
          }
        }
      }

      n2 += ncartB;
    }

    n1 += ncartA;
  }
}

/**
 * Unpacks second derivatives (or their contractions) from the packed order of
 * ECPIntegrator::second_derivs into a full symmetric 3 natoms x 3 natoms matrix, row-major.
 */
std::vector<double> unpack_hessian(const std::vector<double>& packed, int natoms);

}  // namespace libecpint

#endif
