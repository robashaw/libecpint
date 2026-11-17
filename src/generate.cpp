/*
*      Copyright (c) 2020 Robert Shaw
       This file is a part of Libecpint.
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

#include "generate.hpp"

#include <algorithm>
#include <iomanip>
#include <limits>
#include <set>

/**
 * \file generate.cpp
 * \brief ECP integral code generator
 */

/**
 * For a given ECP integral, Q(LA, LB, lam), generates the unrolled integral code.
 * This determines which radial integrals are necessary, based on angular integral screening, and
 * unrolls the angular integration.
 *
 * @param LA - the shellA angular momentum
 * @param LB - the shellB angular momentum
 * @param lam - the ECP angular momentum
 * @param angInts - the angular integrals
 */
void generate_lists(int LA, int LB, int lam, libecpint::AngularIntegral& angInts) {
  using namespace libecpint;

  // Create the code file
  std::string ofname =
      "generated/Q" + std::to_string(LA) + std::to_string(LB) + std::to_string(lam) + ".cpp";
  std::ofstream outfile(ofname);
  // Unrolled angular coefficients are written as literals; the default 6 significant digits
  // (e.g. 16pi^2 -> 157.914) made the unrolled s/p integrals wrong at ~1e-6 relative.
  outfile << std::setprecision(std::numeric_limits<double>::max_digits10);

  if (!outfile.is_open())
    std::cerr << "Problems writing to file!" << std::endl;
  else {
#ifdef DEBUG
    std::cout << "Generating Q(" << LA << ", " << LB << ", " << lam << ")... " << std::flush;
#endif
    // Top matter
    outfile << "// Generated as part of Libecpint, Copyright 2017 Robert A Shaw" << std::endl;
    outfile << "#include \"qgen.hpp\"" << std::endl;
    outfile << "namespace libecpint {" << std::endl << "namespace qgen {" << std::endl;
    outfile << "void Q" << LA << "_" << LB << "_" << lam
            << "(const ECP& U, const GaussianShell& shellA, const GaussianShell& shellB, "
            << "const FiveIndex<double> &CA, const FiveIndex<double> &CB, const TwoIndex<double> "
               "&SA, const TwoIndex<double> &SB, const double Am, const double Bm, "
            << "const RadialIntegral &radint, const AngularIntegral& angint, const "
               "RadialIntegral::Parameters& parameters, ThreeIndex<double> &values) {"
            << std::endl
            << std::endl;

    // Store the terms and radials
    std::vector<SumTerm> terms;
    std::set<Triple> radial_triples;

    int z1, z2, w_l;
    int w_l1;
    int w_l2;
    double ang;
    double prefac = 16.0 * M_PI * M_PI;
    int w_size = 2 * lam + 1;
    std::vector<double> w1_contr(w_size * (lam + LA + 1));
    std::vector<double> w2_contr(w_size * (lam + LB + 1));

    // w_contr[l * w_size + lam + mu] = sum_mu1 Omega(x, y, z, lam, mu, l, mu1), for l <= lam + lmax
    auto contract = [&](const int x, const int y, const int z, const int lmax,
                        std::vector<double>& w_contr) {
      for (int l1 = 0; l1 <= lam + lmax; l1++) {
        w_l = l1 * w_size + lam;
        for (int mu = -lam; mu <= lam; mu++) {
          const double* om = angInts.getOmegaRow(x, y, z, lam, mu) + l1 * l1 + l1;
          w_contr[w_l + mu] = 0.0;
          for (int mu1 = -l1; mu1 <= l1; mu1++) w_contr[w_l + mu] += om[mu1];
        }
      }
    };

    // Loop over cartesian shell functions in alpha order, e.g. {xx xy, xz, yy, yz, zz} for l=2
    int na = 0;  // Rows are shellA
    for (int x1 = LA; x1 >= 0; x1--) {
      for (int r1 = LA - x1; r1 >= 0; r1--) {
        z1 = LA - x1 - r1;

        int nb = 0;  // Cols are shellB
        for (int x2 = LB; x2 >= 0; x2--) {
          for (int y2 = LB - x2; y2 >= 0; y2--) {
            z2 = LB - x2 - y2;

            // Begin full ECP integral expansion
            for (int alpha_x = 0; alpha_x <= x1; alpha_x++) {
              for (int alpha_y = 0; alpha_y <= r1; alpha_y++) {
                for (int alpha_z = 0; alpha_z <= z1; alpha_z++) {
                  int alpha = alpha_x + alpha_y + alpha_z;

                  for (int beta_x = 0; beta_x <= x2; beta_x++) {
                    for (int beta_y = 0; beta_y <= y2; beta_y++) {
                      for (int beta_z = 0; beta_z <= z2; beta_z++) {
                        int beta = beta_x + beta_y + beta_z;
                        int N = alpha + beta;

                        contract(alpha_x, alpha_y, alpha_z, alpha, w1_contr);
                        contract(beta_x, beta_y, beta_z, beta, w2_contr);

                        for (int lam1 = 0; lam1 <= lam + alpha; lam1++) {
                          w_l1 = lam1 * w_size + lam;
                          int lam2start = (lam1 + N) % 2;
                          for (int lam2 = lam2start; lam2 <= lam + beta; lam2 += 2) {
                            w_l2 = lam2 * w_size + lam;

                            ang = 0.0;
                            for (int mu = -lam; mu <= lam; mu++)
                              ang += prefac * w1_contr[w_l1 + mu] * w2_contr[w_l2 + mu];

                            if (fabs(ang) > 1e-15) {
                              radial_triples.insert(Triple{N, lam1, lam2});
                              for (int mu = -lam; mu <= lam; mu++) {
                                for (int mu1 = -lam1; mu1 <= lam1; mu1++) {
                                  for (int mu2 = -lam2; mu2 <= lam2; mu2++) {
                                    SumTerm newTerm;
                                    newTerm.SA = Pair(lam1, lam1 + mu1);
                                    newTerm.SB = Pair(lam2, lam2 + mu2);
                                    newTerm.radial = Triple(N, lam1, lam2);
                                    newTerm.CA = Quintuple(0, na, alpha_x, alpha_y, alpha_z);
                                    newTerm.CB = Quintuple(0, nb, beta_x, beta_y, beta_z);
                                    newTerm.ang =
                                        prefac * angInts.getIntegral(alpha_x, alpha_y, alpha_z, lam,
                                                                     mu, lam1, mu1);
                                    newTerm.ang *= angInts.getIntegral(beta_x, beta_y, beta_z, lam,
                                                                       mu, lam2, mu2);
                                    newTerm.mu = lam + mu;
                                    newTerm.na = na;
                                    newTerm.nb = nb;

                                    // exact zeros would still be evaluated (0 * x is not
                                    // folded without fast-math), so are not emitted
                                    if (newTerm.ang != 0.0) terms.push_back(newTerm);
                                  }
                                }
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }

            nb++;
          }
        }

        na++;
      }
    }

    // Determine the maximum number of base integrals needed across the set of all radial integrals
    int nbase = 0;
    if (!radial_triples.empty()) {
      const Triple& tmax = *radial_triples.crbegin();
      nbase = std::get<0>(tmax) + std::get<1>(tmax) - 1;
      nbase = nbase < 0 ? 0 : nbase;
    }

    // Sort the radials into two lists, depending on whether l1 <= l2 (radial_A), or l2 > l1
    // (radial_B) swapping the order of l1/l2 in the latter case
    std::vector<Triple> radial_A, radial_B;
    for (const Triple& t : radial_triples) {
      if (std::get<1>(t) <= std::get<2>(t))
        radial_A.push_back(t);
      else
        radial_B.push_back(Triple{std::get<0>(t), std::get<2>(t), std::get<1>(t)});
    }

    // Compute the correctly ordered radials first
    outfile << "\tstd::vector<Triple> radial_triples_A = {" << std::endl;
    bool first = true;
    for (Triple& t : radial_A) {
      if (!first)
        outfile << "," << std::endl;
      else
        first = false;
      outfile << "\t\tTriple{" + std::to_string(std::get<0>(t)) + ", " +
                     std::to_string(std::get<1>(t)) + ", " + std::to_string(std::get<2>(t)) + "}";
    }
    outfile << "\t};" << std::endl << std::endl;

    outfile << "\tThreeIndex<double> radials(" << lam + LA + LB + 1 << ", " << lam + LA + 1 << ", "
            << lam + LB + 1 << ");" << std::endl;
    outfile << "\tradint.type2(radial_triples_A, " << nbase << ", " << lam
            << ", U, shellA, shellB, Am, Bm, radials);" << std::endl
            << std::endl;

    // Now compute the reverse-ordered radials
    outfile << "\tstd::vector<Triple> radial_triples_B = {" << std::endl;
    first = true;
    for (Triple& t : radial_B) {
      if (!first)
        outfile << "," << std::endl;
      else
        first = false;
      outfile << "\t\tTriple{" + std::to_string(std::get<0>(t)) + ", " +
                     std::to_string(std::get<1>(t)) + ", " + std::to_string(std::get<2>(t)) + "}";
    }
    outfile << "\t};" << std::endl << std::endl;

    outfile << "\tThreeIndex<double> radials_B(" << lam + LA + LB + 1 << ", " << lam + LB + 1
            << ", " << lam + LA + 1 << ");" << std::endl;
    outfile << "\tradint.type2(radial_triples_B, " << nbase << ", " << lam
            << ", U, shellB, shellA, Bm, Am, radials_B);" << std::endl;
    // These need to be compressed into the radials array, with l1/l2 reversed back
    outfile << "\tfor (Triple& t : radial_triples_B) radials(std::get<0>(t), std::get<2>(t), "
               "std::get<1>(t)) = radials_B(std::get<0>(t), std::get<1>(t), std::get<2>(t));"
            << std::endl
            << std::endl;

    // Print out the unrolled angular integral code
    for (auto& term : terms) outfile << "\t" << term << std::endl;
    outfile << "}" << std::endl << "}" << std::endl << "}" << std::endl;
#ifdef DEBUG
    std::cout << "done." << std::endl;
#endif
    outfile.close();
  }
}

int main(int argc, char* argv[]) {
  // Factorial singletons will not have been initialised
  libecpint::initFactorials();
  // Only the small cases are unrolled: LA <= LB <= maxUnrol, LA + LB + lam <= 3 maxUnrol (and
  // src/makelist.py lists the files). The library integrates all others with its rolled-up code.
  const int maxUnrol = libecpint::maxUnrol, maxLam = 3 * maxUnrol;
  libecpint::AngularIntegral angInts(maxUnrol, maxLam);
  angInts.compute();

  // Generate the qgen.hpp header file
  std::string header_name;
  if (argc > 1) {
    header_name = argv[1];
    header_name += "qgen.hpp";
  } else {
    header_name = "generated/qgen";
  }

  std::ifstream qgen_part("generated/qgen.part");
  std::ofstream qgen_head(header_name);
  if (!qgen_part.is_open() || !qgen_head.is_open())
    std::cerr << "Problem creating qgen header file!" << std::endl;
  else {
    std::string line;
    while (!qgen_part.eof()) {
      std::getline(qgen_part, line);
      qgen_head << line << std::endl;
    }
    qgen_part.close();

    // Loop over the unrolled (l1, l2, lam) integrals with l1 <= l2, generating the code and adding
    // the function to the header file.
    for (int j = 0; j <= maxUnrol; j++) {
      for (int i = 0; i <= j; i++) {
        for (int k = 0; i + j + k <= maxLam; k++) {
          generate_lists(i, j, k, angInts);
          qgen_head
              << "\tvoid Q" << i << "_" << j << "_" << k << "("
              << "const ECP&, const GaussianShell&, const GaussianShell&, const "
                 "FiveIndex<double>&, const FiveIndex<double>&, "
              << "const TwoIndex<double>&, const TwoIndex<double>&, double, double, const "
                 "RadialIntegral&, "
              << "const AngularIntegral&, const RadialIntegral::Parameters&, ThreeIndex<double>&);"
              << std::endl;
        }
      }
    }
    qgen_head << std::endl << "}" << std::endl << "}" << std::endl;
    qgen_head << "#endif" << std::endl;
    qgen_head.close();

    // Now generate the kernel lookup in ecpint_gen.cpp
    std::ifstream ecpgen_part("generated/ecpint_gen.part");
    std::ofstream ecpgen_head("generated/ecpint_gen.cpp");

    if (!ecpgen_part.is_open() || !ecpgen_head.is_open())
      std::cerr << "Problem reading/writing ecpgen file!" << std::endl;
    else {
      while (!ecpgen_part.eof()) {
        std::getline(ecpgen_part, line);
        ecpgen_head << line << std::endl;
      }
      ecpgen_part.close();

      // qgen::unrolled(LA, LB, lam): the unrolled kernel, or nullptr
      for (int j = 0; j <= maxUnrol; j++)
        for (int i = 0; i <= j; i++)
          for (int k = 0; i + j + k <= maxLam; k++)
            ecpgen_head << "\tif (LA == " << i << " && LB == " << j << " && lam == " << k
                        << ") return Q" << i << "_" << j << "_" << k << ";" << std::endl;
      ecpgen_head << "\treturn nullptr;" << std::endl << "}" << std::endl;
      ecpgen_head << "}  // namespace qgen" << std::endl
                  << "}  // namespace libecpint" << std::endl;
      ecpgen_head.close();
    }
  }
  return 0;
}
