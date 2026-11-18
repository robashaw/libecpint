/*
 * Integral checks over every system of the manifest tests/data/systems.txt (see tests/README.md):
 *
 *   Stored     the library, with its default settings, against the stored reference values
 *              (tests/data/reference), through compute_integrals, compute_gradient,
 *              compute_hessian and compute_bfield_derivs
 *   (Stored and Live allow for the precision limit of the library's method: CONDITIONING_TOL)
 *   Live       the library against the reference engine computed now, every element of the
 *              integrals and of the dense first, second and magnetic field derivatives (slow:
 *              ctest -L reference)
 *   Screening  the library with its default screening against the library with none: a
 *              difference is a screening error, whatever the accuracy of the integrals
 *   Threads    one ECPIntegral shared by several threads gives results identical to serial ones,
 *              for every per-block kernel (integrals, first, second and magnetic field
 *              derivatives)
 */
#include <gtest/gtest.h>

#include <array>
#include <cstdio>
#include <functional>
#include <string>
#include <thread>
#include <tuple>
#include <vector>

#include "compare.hpp"
#include "reference.hpp"
#include "reference_io.hpp"
#include "system.hpp"

using namespace testsupport;

namespace {

std::vector<std::string> system_names() {
  std::vector<std::string> names;
  for (const SystemSpec& spec : read_manifest(data_dir() + "/systems.txt"))
    names.push_back(spec.name);
  return names;
}

std::string param_name(const ::testing::TestParamInfo<std::string>& info) { return info.param; }

class SuiteTest : public ::testing::TestWithParam<std::string> {
 protected:
  void SetUp() override { sys_ = load_system(GetParam()); }
  System sys_;
};

// Labels of the stored values: lower triangle of the integrals, contractions per density
std::string ints_label(const System& sys, size_t k) {
  size_t i = 0;
  while ((i + 1) * (i + 2) / 2 <= k) i++;
  const size_t j = k - i * (i + 1) / 2;
  return "V [" + sys.function_label(i) + ", " + sys.function_label(j) + "]";
}

}  // namespace

using Stored = SuiteTest;
TEST_P(Stored, Values) {
  const std::string path = reference_path(sys_.spec.name);
  ReferenceData ref;
  try {
    ref = read_reference(path);
  } catch (const std::exception& e) {
    FAIL() << e.what() << "\nrun: cmake --build <build> --target regenerate_references";
  }
  ASSERT_EQ(ref.hash, sys_.hash())
      << "references stale (the system's inputs changed since " << path
      << " was written)\nrun: cmake --build <build> --target regenerate_references";
  const int deriv = sys_.spec.deriv;
  ASSERT_EQ(ref.deriv, deriv) << "references stale: regenerate_references";
  ASSERT_EQ(ref.ints_scale.size(), ref.ints.size())
      << "references predate conditioning scales: regenerate_references";

  libecpint::ECPIntegrator f;
  sys_.setup(f, deriv);
  const ReferenceData lib = library_reference_data(sys_, f, deriv);
  const int n3 = 3 * ref.natoms;

  const Comparison ci = compare_values(lib.ints, ref.ints, DEFAULT_REL_TOL, &ref.ints_scale);
  EXPECT_TRUE(ci.pass) << ci.describe("integrals", ints_label(sys_, ci.worst));
  if (deriv >= 1) {
    const Comparison cg =
        compare_values(lib.gradient, ref.gradient, DEFAULT_REL_TOL, &ref.gradient_scale);
    EXPECT_TRUE(cg.pass) << cg.describe(
        "contracted first derivatives",
        first_label(sys_, cg.worst % n3) + " (density " + std::to_string(cg.worst / n3) + ")");
    ASSERT_EQ(ref.bfield.size(), lib.bfield.size())
        << "references predate magnetic field derivatives: regenerate_references";
    const Comparison cb =
        compare_values(lib.bfield, ref.bfield, DEFAULT_REL_TOL, &ref.bfield_scale);
    EXPECT_TRUE(cb.pass) << cb.describe(
        "contracted magnetic field derivatives",
        bfield_label(cb.worst % 3) + " (density " + std::to_string(cb.worst / 3) + ")");
  }
  if (deriv >= 2) {
    const Comparison ch =
        compare_values(lib.hessian, ref.hessian, DEFAULT_REL_TOL, &ref.hessian_scale);
    const int k = ch.worst % (n3 * n3), row = k / n3, col = k % n3;
    EXPECT_TRUE(ch.pass) << ch.describe(
        "contracted second derivatives",
        "d2/d" + sys_.atom_label(row / 3) + "_" + "xyz"[row % 3] + " d" + sys_.atom_label(col / 3) +
            "_" + "xyz"[col % 3] + " (density " + std::to_string(ch.worst / (n3 * n3)) + ")");
  }
}

using Live = SuiteTest;
TEST_P(Live, Dense) {
  const int deriv = sys_.spec.deriv;
  libecpint::ECPIntegrator f;
  sys_.setup(f, deriv);
  Integrals scale;
  const Integrals ref = reference_integrals(f, deriv, &scale);
  const Integrals lib = library_integrals(f, deriv);
  EXPECT_TRUE(compare_integrals(sys_, lib, ref, deriv, DEFAULT_REL_TOL, &scale));
}

using Screening = SuiteTest;
TEST_P(Screening, OnOff) {
  const int deriv = sys_.spec.deriv;
  libecpint::ECPIntegrator on, off;
  sys_.setup(on, deriv);
  // A negative shell-pair tolerance skips no shell pairs, and a zero radial tolerance no channels
  // or primitive triples
  sys_.setup(off, deriv, -1.0, 0.0);
  EXPECT_TRUE(
      compare_integrals(sys_, library_integrals(on, deriv), library_integrals(off, deriv), deriv));
}

using Threads = SuiteTest;
TEST_P(Threads, Deterministic) {
  // Every (ECP, shell pair) block of each per-block kernel (integrals; first and second
  // derivatives, up to the system's order; magnetic field derivatives), computed serially and then
  // concurrently from several threads sharing one ECPIntegral, must be identical. The derivative
  // kernels reach the shifted shells (l - 2 to l + 2, exponent-scaled coefficients) and the
  // radial integrals and caches for them, which the integrals alone do not.
  const int deriv = sys_.spec.deriv;
  libecpint::ECPIntegrator f;
  sys_.setup(f, deriv);
  const libecpint::ECPIntegral& engine = *f.ecpint;
  struct Task {
    int e, a, b;
  };
  std::vector<Task> tasks;
  for (int e = 0; e < f.ecps.getN(); e++)
    for (int a = 0; a < int(f.shells.size()); a++)
      for (int b = 0; b <= a; b++) tasks.push_back({e, a, b});

  // Each kernel writes its block(s) for a task, flattened, into out
  using libecpint::TwoIndex;
  using Kernel = std::function<void(const Task&, std::vector<double>&)>;
  auto flatten = [](const auto& blocks, std::vector<double>& out) {
    out.clear();
    for (const TwoIndex<double>& m : blocks) out.insert(out.end(), m.data.begin(), m.data.end());
  };
  auto args = [&](const Task& t) {
    return std::tie(f.ecps.getECP(t.e), f.shells[t.a], f.shells[t.b]);
  };
  std::vector<std::pair<const char*, Kernel>> kernels;
  kernels.push_back({"integrals", [&](const Task& t, std::vector<double>& out) {
                       TwoIndex<double> v;
                       const auto [U, A, B] = args(t);
                       engine.compute_shell_pair(U, A, B, v);
                       out = v.data;
                     }});
  if (deriv >= 1) {
    kernels.push_back({"first derivatives", [&](const Task& t, std::vector<double>& out) {
                         std::array<TwoIndex<double>, 9> v;
                         const auto [U, A, B] = args(t);
                         engine.compute_shell_pair_derivative(U, A, B, v);
                         flatten(v, out);
                       }});
    kernels.push_back({"magnetic field derivatives", [&](const Task& t, std::vector<double>& out) {
                         std::array<TwoIndex<double>, 3> v;
                         const auto [U, A, B] = args(t);
                         engine.compute_shell_pair_derivative_Bfield(U, A, B, v);
                         flatten(v, out);
                       }});
  }
  if (deriv >= 2)
    kernels.push_back({"second derivatives", [&](const Task& t, std::vector<double>& out) {
                         std::array<TwoIndex<double>, 45> v;
                         const auto [U, A, B] = args(t);
                         engine.compute_shell_pair_second_derivative(U, A, B, v);
                         flatten(v, out);
                       }});

  const int nthreads = 8;
  for (const auto& [name, run] : kernels) {
    std::vector<std::vector<double>> serial(tasks.size());
    for (size_t i = 0; i < tasks.size(); i++) run(tasks[i], serial[i]);

    std::vector<int> differ(nthreads, 0);
    std::vector<std::thread> pool;
    for (int t = 0; t < nthreads; t++)
      pool.emplace_back([&, t]() {
        std::vector<double> out;
        for (int rep = 0; rep < 2; rep++)  // each thread in a different order
          for (size_t k = 0; k < tasks.size(); k++) {
            const size_t i = (k * (2 * t + 1) + rep * 7) % tasks.size();
            run(tasks[i], out);
            differ[t] += out != serial[i];
          }
      });
    for (auto& th : pool) th.join();
    int total = 0;
    for (const int d : differ) total += d;
    EXPECT_EQ(total, 0) << name << ": of " << nthreads * 2 * tasks.size()
                        << " concurrent shell pairs";
  }
}

INSTANTIATE_TEST_SUITE_P(Systems, Stored, ::testing::ValuesIn(system_names()), param_name);
INSTANTIATE_TEST_SUITE_P(Systems, Live, ::testing::ValuesIn(system_names()), param_name);
INSTANTIATE_TEST_SUITE_P(Systems, Screening, ::testing::ValuesIn(system_names()), param_name);
INSTANTIATE_TEST_SUITE_P(Systems, Threads, ::testing::ValuesIn(system_names()), param_name);

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
