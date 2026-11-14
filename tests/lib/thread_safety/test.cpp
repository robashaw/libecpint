/*
 * Thread-safety test: one ECPIntegral shared by several threads.
 *
 * compute_shell_pair is const and documented as thread safe, so callers may share a single
 * ECPIntegral between threads. Scratch buffers held as mutable members of RadialIntegral were
 * written during the type 1 and on-centre type 2 radial integrals, so concurrent calls raced on
 * them. This computes shell pairs covering type 1, type 2 with a shell on the ECP centre, and type
 * 2 with both shells off the centre, serially and then repeatedly from several threads sharing the
 * engine, and requires identical results. (A data race need not change results in a given run;
 * build with -fsanitize=thread to detect it reliably.)
 */
#include <array>
#include <atomic>
#include <cstdlib>
#include <iostream>
#include <thread>
#include <vector>

#include "ecp.hpp"
#include "ecpint.hpp"
#include "gshell.hpp"
#include "multiarr.hpp"

using namespace libecpint;

namespace {
// La2O3 fragment (bohr): La (ECP centre), La, O, O
const double X[4][3] = {{0.0, 0.944863, -0.3286479},
                        {0.0, -0.944863, -0.3286479},
                        {0.0, 1.88973, 1.5610781},
                        {0.0, 0.0, 1.5610781}};

// La def2 ECP (46 core electrons): n, l, exponent, coefficient
const double LA_ECP[10][4] = {{2, 3, 4.0286, -36.010016}, {2, 0, 3.3099, 91.932177},
                              {2, 0, 1.655, -3.788764},   {2, 0, 4.0286, 36.010016},
                              {2, 1, 2.8368, 63.759486},  {2, 1, 1.4184, -0.647958},
                              {2, 1, 4.0286, 36.010016},  {2, 2, 2.0213, 36.116173},
                              {2, 2, 1.0107, 0.219114},   {2, 2, 4.0286, 36.010016}};

struct Shell {
  int centre, l;
  std::vector<double> exps, coefs;
};

const Shell SHELLS[] = {
    {0, 2, {1.5, 0.5}, {1.0, 0.6}},             // La d, on the ECP centre
    {0, 0, {2.0, 0.3}, {0.8, 0.4}},             // La s, on the ECP centre
    {1, 3, {0.25683}, {0.069128541039011268}},  // La f, other La
    {2, 1, {17.721504317, 3.863550544, 1.0480920883}, {3.24601, 2.57350, 1.12084}},     // O p
    {3, 2, {1.2}, {2.2645224825069943}},                                                // O d
    {2, 0, {2266.1767785, 340.87010191, 21.47964494}, {-1.26159, -2.27485, -3.32999}},  // O 1s
};
constexpr int NSHELL = sizeof(SHELLS) / sizeof(SHELLS[0]);

GaussianShell make(const Shell& s) {
  GaussianShell g(std::array<double, 3>{X[s.centre][0], X[s.centre][1], X[s.centre][2]}, s.l);
  for (size_t i = 0; i < s.exps.size(); i++) g.addPrim(s.exps[i], s.coefs[i]);
  return g;
}
}  // namespace

int main(int argc, char** argv) {
  ECP U(X[0]);
  for (int k = 0; k < 10; k++)
    U.addPrimitive((int)LA_ECP[k][0], (int)LA_ECP[k][1], LA_ECP[k][2], LA_ECP[k][3], k == 9);
  std::vector<GaussianShell> shells;
  for (const Shell& s : SHELLS) shells.push_back(make(s));

  const ECPIntegral eng(3, 3, 0);
  std::vector<std::vector<double>> ref;
  for (int a = 0; a < NSHELL; a++)
    for (int b = 0; b <= a; b++) {
      TwoIndex<double> v;
      eng.compute_shell_pair(U, shells[a], shells[b], v);
      ref.push_back(v.data);
    }

  // optional argument: number of repetitions (fewer under ThreadSanitizer, which is slow)
  const int nthreads = 8, nrep = argc > 1 ? std::atoi(argv[1]) : 20;
  std::atomic<int> nbad{0};
  std::vector<std::thread> pool;
  for (int t = 0; t < nthreads; t++)
    pool.emplace_back([&, t]() {
      for (int rep = 0; rep < nrep; rep++) {
        int ix = 0;
        // each thread walks the pairs in a different order, so different code paths overlap
        for (int a0 = 0; a0 < NSHELL; a0++) {
          const int a = (a0 + t) % NSHELL;
          for (int b = 0; b <= a; b++) {
            TwoIndex<double> v;
            eng.compute_shell_pair(U, shells[a], shells[b], v);
            ix = a * (a + 1) / 2 + b;
            if (v.data != ref[ix]) nbad++;
          }
        }
      }
    });
  for (auto& th : pool) th.join();

  std::cout << nthreads << " threads x " << nrep << " repetitions x " << ref.size()
            << " shell pairs sharing one ECPIntegral: " << nbad << " results differ from serial"
            << std::endl;
  if (nbad) {
    std::cout << "FAIL" << std::endl;
    return 1;
  }
  std::cout << "PASS" << std::endl;
  return 0;
}
