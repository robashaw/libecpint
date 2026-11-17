#include "compare.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>

#include "assembly.hpp"

namespace testsupport {

namespace {
const char* const XYZ = "xyz";

// Flattens matrices, recording where each starts
std::vector<double> flatten(const std::vector<libecpint::TwoIndex<double>>& mats) {
  std::vector<double> v;
  for (const auto& m : mats) v.insert(v.end(), m.data.begin(), m.data.end());
  return v;
}
}  // namespace

std::string Comparison::describe(const std::string& what, const std::string& worst_label) const {
  char buf[320];
  std::snprintf(buf, sizeof buf,
                "%s: %zu values beyond tolerance (%.2e = rel * max |reference| %.3e, plus any "
                "conditioning allowance); worst error %.3e (tolerance %.2e) at ",
                what.c_str(), nbad, tol, max_ref, worst_err, worst_tol);
  return buf + worst_label;
}

Comparison compare_values(const std::vector<double>& value, const std::vector<double>& ref,
                          const double rel, const std::vector<double>* scale) {
  Comparison c;
  if (value.size() != ref.size() || (scale && scale->size() != ref.size())) {
    c.pass = false;
    c.nbad = std::max(value.size(), ref.size());
    return c;
  }
  for (const double r : ref) c.max_ref = std::max(c.max_ref, std::abs(r));
  c.tol = rel * c.max_ref;
  double worst_ratio = -1.0;
  for (size_t i = 0; i < ref.size(); i++) {
    const double e = std::abs(value[i] - ref[i]);
    const double t = c.tol + (scale ? CONDITIONING_TOL * (*scale)[i] : 0.0);
    if (!(e <= t)) c.nbad++;  // NaN counts as bad
    const double ratio = t > 0.0 ? e / t : (e > 0.0 ? HUGE_VAL : 0.0);
    if (!(ratio <= worst_ratio)) {
      worst_ratio = ratio;
      c.worst = i;
      c.worst_err = e;
      c.worst_tol = t;
    }
  }
  c.pass = c.nbad == 0;
  return c;
}

std::string first_label(const System& sys, const int index) {
  return "d/d" + sys.atom_label(index / 3) + "_" + XYZ[index % 3];
}

std::string second_label(const System& sys, const int index, const int natoms) {
  const char* const sym[6] = {"xx", "xy", "xz", "yy", "yz", "zz"};
  for (int a = 0; a < natoms; a++)
    for (int b = a; b < natoms; b++) {
      const int start = a == b ? H_START(a, a, natoms) + 3 : H_START(a, b, natoms);
      const int n = a == b ? 6 : 9;
      if (index < start || index >= start + n) continue;
      const int k = index - start;
      const std::string pq =
          a == b ? sym[k] : std::string(1, XYZ[k / 3]) + std::string(1, XYZ[k % 3]);
      return "d2/d" + sys.atom_label(a) + "_" + pq[0] + " d" + sys.atom_label(b) + "_" + pq[1];
    }
  return "?";
}

::testing::AssertionResult compare_integrals(const System& sys, const Integrals& value,
                                             const Integrals& ref, const int deriv,
                                             const double rel, const Integrals* scale) {
  std::string failures;
  const int ncart = ref.ints.dims[0];
  const int natoms = ref.first.size() / 3;
  auto element = [&](const size_t flat, const std::string& matrix) {
    const int k = (flat / ncart) % ncart, l = flat % ncart;
    return matrix + " [" + sys.function_label(k) + ", " + sys.function_label(l) + "]";
  };
  auto check = [&](const char* what, const std::vector<double>& v, const std::vector<double>& r,
                   const std::vector<double>* sc, auto label) {
    const Comparison c = compare_values(v, r, rel, sc);
    if (c.pass) return;
    if (v.size() != r.size())
      failures += std::string(what) + ": sizes differ\n";
    else
      failures += c.describe(what, label(c.worst)) + "\n";
  };
  check("integrals", value.ints.data, ref.ints.data, scale ? &scale->ints.data : nullptr,
        [&](const size_t i) { return element(i, "V"); });
  const size_t block = size_t(ncart) * ncart;
  if (deriv >= 1) {
    const std::vector<double> sc = scale ? flatten(scale->first) : std::vector<double>();
    check("first derivatives", flatten(value.first), flatten(ref.first), scale ? &sc : nullptr,
          [&](const size_t i) { return element(i, first_label(sys, i / block) + " V"); });
  }
  if (deriv >= 2) {
    const std::vector<double> sc = scale ? flatten(scale->second) : std::vector<double>();
    check("second derivatives", flatten(value.second), flatten(ref.second), scale ? &sc : nullptr,
          [&](const size_t i) { return element(i, second_label(sys, i / block, natoms) + " V"); });
  }
  if (failures.empty()) return ::testing::AssertionSuccess();
  return ::testing::AssertionFailure() << failures;
}

}  // namespace testsupport
