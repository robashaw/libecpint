/*
 * Writes the stored reference values tests/data/reference/<system>.ref from the independent
 * reference engine, for every system of the manifest or for those named on the command line.
 * Run through the regenerate_references target: cmake --build <build> --target
 * regenerate_references. It uses all hardware threads, and takes minutes.
 */
#include <chrono>
#include <cstdio>
#include <exception>
#include <string>
#include <vector>

#include "reference.hpp"
#include "reference_io.hpp"
#include "system.hpp"

using namespace testsupport;

int main(int argc, char** argv) {
  std::vector<std::string> names(argv + 1, argv + argc);
  if (names.empty())
    for (const SystemSpec& spec : read_manifest(data_dir() + "/systems.txt"))
      names.push_back(spec.name);

  int failed = 0;
  for (const std::string& name : names) {
    try {
      const System sys = load_system(name);
      const auto start = std::chrono::steady_clock::now();
      libecpint::ECPIntegrator f;
      sys.setup(f, sys.spec.deriv);
      Integrals scale;
      const Integrals ref = reference_integrals(f, sys.spec.deriv, &scale);
      const std::string path = reference_path(name);
      write_reference(path, make_reference_data(sys, ref, scale, sys.spec.deriv));
      const double secs =
          std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
      std::printf("%-14s %4d functions, derivatives to order %d: %7.1f s -> %s\n", name.c_str(),
                  f.ncart, sys.spec.deriv, secs, path.c_str());
    } catch (const std::exception& e) {
      std::printf("%-14s FAILED: %s\n", name.c_str(), e.what());
      failed++;
    }
  }
  return failed > 0;
}
