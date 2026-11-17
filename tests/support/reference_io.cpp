#include "reference_io.hpp"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <stdexcept>

#include "assembly.hpp"

namespace testsupport {

namespace {
double contract(const std::vector<double>& D, const libecpint::TwoIndex<double>& M) {
  double sum = 0.0;
  for (size_t k = 0; k < M.data.size(); k++) sum += D[k] * M.data[k];
  return sum;
}

std::vector<double> lower_triangle(const libecpint::TwoIndex<double>& M) {
  std::vector<double> v;
  for (int i = 0; i < M.dims[0]; i++)
    for (int j = 0; j <= i; j++) v.push_back(M(i, j));
  return v;
}
}  // namespace

std::vector<std::vector<double>> test_densities(const int ncart) {
  std::vector<std::vector<double>> D(2, std::vector<double>(ncart * ncart));
  for (int k = 0; k < ncart * ncart; k++) {
    D[0][k] = std::sin(1.0 + 0.37 * k);
    D[1][k] = std::cos(0.3 + 0.61 * k);
  }
  return D;
}

ReferenceData make_reference_data(const System& sys, const Integrals& ref, const Integrals& scale,
                                  const int deriv) {
  ReferenceData data;
  data.system = sys.spec.name;
  data.hash = sys.hash();
  data.generator = LIBECPINT_TEST_GENERATOR;
  data.deriv = deriv;
  data.ncart = ref.ints.dims[0];
  data.natoms = sys.integrator_atoms().size();
  data.ints = lower_triangle(ref.ints);
  data.ints_scale = lower_triangle(scale.ints);
  // The derivatives contracted with each density, and their scales with its absolute values
  auto contract_all = [&](const Integrals& in, const std::vector<double>& D,
                          std::vector<double>& gradient, std::vector<double>& hessian) {
    if (deriv >= 1)
      for (const auto& M : in.first) gradient.push_back(contract(D, M));
    if (deriv >= 2) {
      std::vector<double> packed;
      for (const auto& M : in.second) packed.push_back(contract(D, M));
      const auto h = libecpint::unpack_hessian(packed, data.natoms);
      hessian.insert(hessian.end(), h.begin(), h.end());
    }
  };
  for (const auto& D : test_densities(data.ncart)) {
    contract_all(ref, D, data.gradient, data.hessian);
    std::vector<double> absD(D.size());
    for (size_t k = 0; k < D.size(); k++) absD[k] = std::abs(D[k]);
    contract_all(scale, absD, data.gradient_scale, data.hessian_scale);
  }
  return data;
}

ReferenceData library_reference_data(const System& sys, libecpint::ECPIntegrator& f,
                                     const int deriv) {
  ReferenceData data;
  data.system = sys.spec.name;
  data.hash = sys.hash();
  data.deriv = deriv;
  data.ncart = f.ncart;
  data.natoms = f.natoms;
  f.compute_integrals();
  data.ints = lower_triangle(f.integrals);
  for (const auto& D : test_densities(data.ncart)) {
    if (deriv >= 1) {
      const auto g = f.compute_gradient(D.data());
      data.gradient.insert(data.gradient.end(), g.begin(), g.end());
    }
    if (deriv >= 2) {
      const auto h = f.compute_hessian(D.data());
      data.hessian.insert(data.hessian.end(), h.begin(), h.end());
    }
  }
  return data;
}

std::string reference_path(const std::string& name) {
  return data_dir() + "/reference/" + name + ".ref";
}

void write_reference(const std::string& path, const ReferenceData& data) {
  std::FILE* fp = std::fopen(path.c_str(), "w");
  if (!fp) throw std::runtime_error("cannot write " + path);
  std::fprintf(fp, "# libecpint reference values from the independent reference engine.\n");
  std::fprintf(fp, "# Regenerate with: cmake --build <build> --target regenerate_references\n");
  std::fprintf(fp, "system %s\nhash %s\ngenerator %s\nderiv %d\nncart %d\nnatoms %d\n",
               data.system.c_str(), data.hash.c_str(), data.generator.c_str(), data.deriv,
               data.ncart, data.natoms);
  auto section = [&](const char* name, const std::vector<double>& v) {
    std::fprintf(fp, "%s %zu\n", name, v.size());
    for (const double x : v) std::fprintf(fp, "%.17g\n", x);
  };
  section("integrals", data.ints);
  section("gradient", data.gradient);
  section("hessian", data.hessian);
  section("integrals_scale", data.ints_scale);
  section("gradient_scale", data.gradient_scale);
  section("hessian_scale", data.hessian_scale);
  std::fclose(fp);
}

ReferenceData read_reference(const std::string& path) {
  std::ifstream file(path);
  if (!file) throw std::runtime_error("cannot open " + path);
  ReferenceData data;
  std::string line, key;
  while (std::getline(file, line)) {
    if (line.empty() || line[0] == '#') continue;
    std::istringstream ss(line);
    ss >> key;
    if (key == "system")
      ss >> data.system;
    else if (key == "hash")
      ss >> data.hash;
    else if (key == "generator")
      ss >> data.generator;
    else if (key == "deriv")
      ss >> data.deriv;
    else if (key == "ncart")
      ss >> data.ncart;
    else if (key == "natoms")
      ss >> data.natoms;
    else if (key == "integrals" || key == "gradient" || key == "hessian" ||
             key == "integrals_scale" || key == "gradient_scale" || key == "hessian_scale") {
      size_t n;
      ss >> n;
      std::vector<double>& v = key == "integrals"         ? data.ints
                               : key == "gradient"        ? data.gradient
                               : key == "hessian"         ? data.hessian
                               : key == "integrals_scale" ? data.ints_scale
                               : key == "gradient_scale"  ? data.gradient_scale
                                                          : data.hessian_scale;
      v.resize(n);
      for (size_t i = 0; i < n; i++)
        if (!(file >> v[i])) throw std::runtime_error(path + ": truncated " + key);
      if (n > 0) std::getline(file, line);  // rest of the last value's line
    } else
      throw std::runtime_error(path + ": unknown entry " + key);
  }
  return data;
}

}  // namespace testsupport
