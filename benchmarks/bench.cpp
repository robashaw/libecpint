/*
 * Timing benchmark of libecpint (see benchmarks/README.md). Three groups:
 *
 *   kernels  ECPIntegral::compute_shell_pair for each angular class (LA <= LB <= 4) against the
 *            def2 ECP of Au (s, p, d projectors and an f local channel), with the shells off the
 *            ECP centre (the type 2 code that LIBECPINT_MAX_UNROL changes), and, as controls that
 *            it does not change, one or both on it. Warm: the same pair is called repeatedly, so
 *            the per-thread radial caches hit, and the times weigh the angular part more than in
 *            whole calculations.
 *   systems  ECPIntegrator init, compute_integrals, compute_gradient and compute_hessian (if the
 *            manifest tests second derivatives) for test systems from tests/data/systems.txt, and
 *            for an iodoalkane chain C16 (as in scaling). Cold: every run is on a fresh thread,
 *            with empty per-thread caches, as for the first calculation at a geometry; each sample
 *            is the mean of as many runs as take 50 ms.
 *   scaling  compute_integrals and compute_gradient for iodoalkane chains of doubling length
 *            (def2-SVP, ecp28mdf for I), cold, with local exponents of the growth with size.
 *
 * Each measurement is repeated (after one discarded run), and reported as the median, minimum and
 * spread (interquartile range / median) of the repetitions, with a checksum of the results (a
 * weighted sum of every value), so that builds that compute different values can be detected.
 *
 * usage: ecpint_bench [kernels] [systems] [scaling] [options]   (default: kernels systems)
 *   --reps N         repetitions (default: 11 for kernels, 5 for systems, 3 for scaling)
 *   --csv FILE       also write the results to FILE as CSV
 *   --tag NAME       label for the CSV rows, e.g. the build (default "default")
 *   --max-carbons N  largest chain for scaling (default 64)
 *   --wasted         scaling: also count the screened-in shell-pair calls of compute_integrals
 *                    that return only values below the tolerance (slack in the screening)
 */
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <memory>
#include <string>
#include <thread>
#include <tuple>
#include <utility>
#include <vector>

#include "api.hpp"
#include "config.hpp"
#include "screening.hpp"
#include "system.hpp"

namespace {

using Clock = std::chrono::steady_clock;
double seconds_since(const Clock::time_point t0) {
  return std::chrono::duration<double>(Clock::now() - t0).count();
}

/// sum_k sin(1 + 0.37 k) v_k: changes with any value
double checksum(const std::vector<double>& v) {
  double s = 0.0;
  for (size_t k = 0; k < v.size(); k++) s += std::sin(1.0 + 0.37 * k) * v[k];
  return s;
}

struct Measurement {
  std::string group, name, unit;
  std::vector<double> samples;
  double checksum = 0.0;

  double quantile(const double q) const {
    std::vector<double> s = samples;
    std::sort(s.begin(), s.end());
    const double x = q * (s.size() - 1);
    const size_t i = static_cast<size_t>(x);
    return i + 1 < s.size() ? s[i] + (x - i) * (s[i + 1] - s[i]) : s[i];
  }
  double median() const { return quantile(0.5); }
  double min() const { return *std::min_element(samples.begin(), samples.end()); }
  double spread() const { return (quantile(0.75) - quantile(0.25)) / median(); }
};

struct Options {
  bool kernels = false, systems = false, scaling = false, wasted = false;
  int reps = 0, max_carbons = 64;
  std::string csv, tag = "default";
};

std::vector<Measurement> results;

void report(const Measurement& m) {
  std::printf("  %-34s %12.4g %12.4g %7.1f%%   %-7s %+.15e\n", m.name.c_str(), m.median(), m.min(),
              100.0 * m.spread(), m.unit.c_str(), m.checksum);
  std::fflush(stdout);
  results.push_back(m);
}

void header(const char* title) {
  std::printf("\n%s\n  %-34s %12s %12s %8s   %-7s %s\n", title, "", "median", "min", "spread",
              "unit", "checksum");
}

/// Runs f on a fresh thread (empty per-thread caches) and returns its result
template <class F>
auto on_fresh_thread(F f) -> decltype(f()) {
  decltype(f()) r;
  std::thread t([&]() { r = f(); });
  t.join();
  return r;
}

// ---------------------------------------------------------------------------------------------
// kernels

void kernels(const int reps) {
  using namespace libecpint;
  header("kernels: compute_shell_pair, Au def2 ECP (s, p, d projectors), warm, per call");
  ECPBasis basis;
  basis.addECP_from_file(79, {0.0, 0.0, 0.0}, testsupport::share_dir() + "/raw/def2-ecp.ecp");
  const ECP& U = basis.getECP(0);
  const int maxL = 4;
  const ECPIntegral engine(maxL, U.getL(), 0);
  // Shell centres (bohr) about 3.8 bohr from the ECP, as for a ligand atom, and the ECP centre
  const std::array<double, 3> origin = {0.0, 0.0, 0.0}, A = {0.9, -0.6, 3.6}, B = {-2.8, 1.5, 2.2};
  struct Case {
    const char* name;
    std::array<double, 3> a, b;
  };
  const Case cases[3] = {{"off", A, B}, {"A on ECP", origin, B}, {"both on ECP", origin, origin}};
  const char* L = "spdfghik";
  for (const Case& c : cases)
    for (int lb = 0; lb <= maxL; lb++)
      for (int la = 0; la <= lb; la++) {
        GaussianShell sa(c.a, la), sb(c.b, lb);  // two primitives each, valence-like
        sa.addPrim(2.4, 0.4);
        sa.addPrim(0.45, 0.7);
        sb.addPrim(1.9, 0.5);
        sb.addPrim(0.32, 0.6);
        TwoIndex<double> v;
        auto run = [&](const long n) {
          const auto t0 = Clock::now();
          for (long i = 0; i < n; i++) engine.compute_shell_pair(U, sa, sb, v);
          return seconds_since(t0);
        };
        run(1);  // fills the caches
        long n = 1;
        while (run(n) < 0.01) n *= 2;  // samples of at least 10 ms
        Measurement m{"kernels", std::string(c.name) + " " + L[la] + L[lb], "us/call"};
        for (int r = 0; r < reps; r++) m.samples.push_back(1e6 * run(n) / n);
        m.checksum = checksum(v.data);
        report(m);
      }
}

// ---------------------------------------------------------------------------------------------
// systems

std::vector<double> test_density(const int n) {
  std::vector<double> d(size_t(n) * n);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++) d[size_t(i) * n + j] = std::sin(1.0 + 0.37 * (i + j)) / (1 + i + j);
  return d;
}

// ---------------------------------------------------------------------------------------------
// scaling

struct ChainAtom {
  std::string symbol;
  std::array<double, 3> r;  // bohr
};

// Zigzag chain C_n H_(2n+2) in the xy plane (C-C 1.54 A, angle 112 degrees). Each carbon has two
// substituents at +-z, tilted away from the chain axis (tetrahedral), C-H 1.09 A and C-I 2.14 A;
// carbons 2, 6, 10, ... carry I at +z, so the density of ECP centres is constant. The end carbons
// have a third hydrogen along the chain axis.
std::vector<ChainAtom> make_chain(const int ncarbon) {
  const double bohr = 1.8897261246;
  const double cc = 1.54, half = 0.5 * 112.0 * M_PI / 180.0, tet = 0.5 * 109.47 * M_PI / 180.0;
  const double dx = cc * std::sin(half), dy = 0.5 * cc * std::cos(half);
  std::vector<ChainAtom> atoms;
  auto add = [&](const char* s, const double x, const double y, const double z) {
    atoms.push_back({s, {x * bohr, y * bohr, z * bohr}});
  };
  for (int i = 0; i < ncarbon; i++) {
    const double side = i % 2 == 0 ? 1.0 : -1.0, x = i * dx, y = side * dy;
    add("C", x, y, 0.0);
    for (const double sz : {1.0, -1.0}) {
      const bool iodine = sz > 0 && i % 4 == 2;
      const double len = iodine ? 2.14 : 1.09;
      add(iodine ? "I" : "H", x, y + side * len * std::cos(tet), sz * len * std::sin(tet));
    }
    if (i == 0) add("H", x - 1.09, y, 0.0);
    if (i == ncarbon - 1) add("H", x + 1.09, y, 0.0);
  }
  return atoms;
}

void setup_chain(libecpint::ECPIntegrator& f, const std::vector<ChainAtom>& atoms,
                 const int deriv) {
  static const auto basis = testsupport::read_gbs(testsupport::data_dir() + "/basis/def2-svp.gbs");
  std::vector<double> coords, exps, coefs, ecp_coords;
  std::vector<int> ams, lengths, charges;
  for (const auto& a : atoms) {
    for (const auto& s : basis.at(a.symbol)) {
      coords.insert(coords.end(), a.r.begin(), a.r.end());
      exps.insert(exps.end(), s.exps.begin(), s.exps.end());
      coefs.insert(coefs.end(), s.coefs.begin(), s.coefs.end());
      ams.push_back(s.l);
      lengths.push_back(s.exps.size());
    }
    if (a.symbol == "I") {
      ecp_coords.insert(ecp_coords.end(), a.r.begin(), a.r.end());
      charges.push_back(53);
    }
  }
  f.set_gaussian_basis(ams.size(), coords.data(), exps.data(), coefs.data(), ams.data(),
                       lengths.data());
  f.set_ecp_basis_from_library(charges.size(), ecp_coords.data(), charges.data(),
                               std::vector<std::string>(charges.size(), "ecp28mdf"),
                               testsupport::share_dir());
  f.init(deriv);
}

/// Screened-in shell-pair calls of compute_integrals, and how many return only values below the
/// tolerance
std::pair<long long, long long> count_calls(libecpint::ECPIntegrator& f) {
  long long calls = 0, wasted = 0;
  libecpint::TwoIndex<double> values;
  // the same screen as compute_integrals (see integrals_loop in assembly.hpp)
  const int necp = f.ecps.getN();
  const auto screen = libecpint::ShellPairScreen::from_bounds(
      libecpint::screening::shell_ecp_bounds(f.shells, f.ecps, 0), necp, f.tolerance);
  for (size_t s1 = 0; s1 < f.shells.size(); s1++)
    for (const auto& [i, sbA] : screen.ecps_for(s1, necp))
      for (size_t s2 = 0; s2 <= s1; s2++) {
        const libecpint::ECP& U = f.ecps.getECP(i);
        if (sbA * screen.bounds[s2 * necp + i] <= f.tolerance) continue;
        calls++;
        f.ecpint->compute_shell_pair(U, f.shells[s1], f.shells[s2], values);
        double vmax = 0.0;
        for (const double v : values.data) vmax = std::max(vmax, std::abs(v));
        if (vmax < f.tolerance) wasted++;
      }
  return {calls, wasted};
}

/// Sets up an integrator (basis, ECPs, init) for derivatives up to the given order
using Setup = std::function<void(libecpint::ECPIntegrator&, int)>;

/// One cold run, on a fresh thread: what = 0 times the setup (for derivatives to init_deriv);
/// otherwise the setup is untimed, then compute_integrals (1), compute_gradient (2) or
/// compute_hessian (3) is timed. Returns (seconds, checksum of the results).
std::pair<double, double> cold_run(const Setup& setup, const int what, const int init_deriv) {
  return on_fresh_thread([&]() {
    libecpint::ECPIntegrator f;
    auto t0 = Clock::now();
    setup(f, what == 0 ? init_deriv : what - 1);
    if (what == 0) return std::make_pair(seconds_since(t0), double(f.ncart));
    const std::vector<double> d = test_density(f.ncart);
    std::vector<double> out;
    t0 = Clock::now();
    if (what == 1) {
      f.compute_integrals();
      out = f.integrals.data;
    } else if (what == 2)
      out = f.compute_gradient(d.data());
    else
      out = f.compute_hessian(d.data());
    return std::make_pair(seconds_since(t0), checksum(out));
  });
}

/// Cold timings of `what` (see cold_run), after one discarded run; each sample is the mean of as
/// many cold runs as take at least 50 ms, so that short calculations are resolved
Measurement measure_cold(const std::string& group, const std::string& name, const Setup& setup,
                         const int what, const int init_deriv, const int reps) {
  static const char* names[4] = {"init", "integrals", "gradient", "hessian"};
  Measurement m{group, name + " " + names[what], "s"};
  const auto first = cold_run(setup, what, init_deriv);
  const int k = std::max(1, static_cast<int>(std::ceil(0.05 / first.first)));
  for (int r = 0; r < reps; r++) {
    double total = 0.0;
    for (int i = 0; i < k; i++) {
      const auto run = cold_run(setup, what, init_deriv);
      total += run.first;
      m.checksum = run.second;
    }
    m.samples.push_back(total / k);
  }
  return m;
}

void systems(const int reps) {
  header("systems: cold (fresh thread per run), per calculation");
  struct Workload {
    std::string name;
    int deriv;  // highest derivative order timed
    Setup setup;
  };
  std::vector<Workload> workloads;
  for (const char* name : {"ag3_mdf", "la2o3", "ch3i_svp", "au2_tzvpp"}) {
    auto sys = std::make_shared<const testsupport::System>(testsupport::load_system(name));
    workloads.push_back({name, sys->spec.deriv,
                         [sys](libecpint::ECPIntegrator& f, const int d) { sys->setup(f, d); }});
  }
  const auto chain = make_chain(16);
  workloads.push_back({"chain_c16", 2, [chain](libecpint::ECPIntegrator& f, const int d) {
                         setup_chain(f, chain, d);
                       }});
  for (const Workload& w : workloads)
    for (int what = 0; what <= 1 + w.deriv; what++)
      report(measure_cold("systems", w.name, w.setup, what, w.deriv, reps));
}

void scaling(const int reps, const int max_carbons, const bool wasted) {
  header("scaling: iodoalkane chains (I on every 4th carbon), def2-SVP + ecp28mdf, cold");
  struct Size {
    int natoms;
    double ints, grad;
    long long calls, wasted;
  };
  std::vector<Size> sizes;
  for (int n = 8; n <= max_carbons; n *= 2) {
    const auto atoms = make_chain(n);
    Size s{int(atoms.size()), 0.0, 0.0, 0, 0};
    const Setup setup = [&atoms](libecpint::ECPIntegrator& f, const int d) {
      setup_chain(f, atoms, d);
    };
    for (int what = 1; what <= 2; what++) {
      const Measurement m =
          measure_cold("scaling", "chain_c" + std::to_string(n), setup, what, 1, reps);
      report(m);
      (what == 1 ? s.ints : s.grad) = m.median();
    }
    if (wasted) {
      libecpint::ECPIntegrator f;
      setup_chain(f, atoms, 0);
      std::tie(s.calls, s.wasted) = count_calls(f);
      std::printf("  %-34s %lld calls, %.1f%% wasted\n", "", s.calls,
                  100.0 * s.wasted / std::max(1LL, s.calls));
    }
    sizes.push_back(s);
  }
  // log(y2 / y1) / log(N2 / N1) between successive sizes, N the number of atoms
  std::printf("\n  local exponents (y ~ N^p)\n  %13s %8s %8s%s\n", "atoms", "ints", "grad",
              wasted ? "    calls   needed" : "");
  for (size_t i = 1; i < sizes.size(); i++) {
    const Size &a = sizes[i - 1], &b = sizes[i];
    const double ln = std::log(double(b.natoms) / a.natoms);
    auto p = [&](const double ya, const double yb) { return std::log(yb / ya) / ln; };
    std::printf("  %6d-%-6d %8.2f %8.2f", a.natoms, b.natoms, p(a.ints, b.ints), p(a.grad, b.grad));
    if (wasted)
      std::printf(" %8.2f %8.2f", p(double(a.calls), double(b.calls)),
                  p(double(a.calls - a.wasted), double(b.calls - b.wasted)));
    std::printf("\n");
  }
}

void write_csv(const std::string& path, const std::string& tag) {
  FILE* f = std::fopen(path.c_str(), "w");
  if (!f) {
    std::fprintf(stderr, "cannot write %s\n", path.c_str());
    std::exit(1);
  }
  std::fprintf(f, "tag,group,name,unit,reps,median,min,spread,checksum\n");
  for (const Measurement& m : results)
    std::fprintf(f, "%s,%s,%s,%s,%zu,%.6e,%.6e,%.4f,%.15e\n", tag.c_str(), m.group.c_str(),
                 m.name.c_str(), m.unit.c_str(), m.samples.size(), m.median(), m.min(), m.spread(),
                 m.checksum);
  std::fclose(f);
}

}  // namespace

int main(int argc, char* argv[]) {
  Options o;
  for (int i = 1; i < argc; i++) {
    const std::string a = argv[i];
    auto value = [&]() {
      if (i + 1 >= argc) {
        std::fprintf(stderr, "%s needs a value\n", a.c_str());
        std::exit(1);
      }
      return std::string(argv[++i]);
    };
    if (a == "kernels")
      o.kernels = true;
    else if (a == "systems")
      o.systems = true;
    else if (a == "scaling")
      o.scaling = true;
    else if (a == "--wasted")
      o.wasted = true;
    else if (a == "--reps")
      o.reps = std::atoi(value().c_str());
    else if (a == "--max-carbons")
      o.max_carbons = std::atoi(value().c_str());
    else if (a == "--csv")
      o.csv = value();
    else if (a == "--tag")
      o.tag = value();
    else {
      std::fprintf(stderr,
                   "usage: %s [kernels] [systems] [scaling] [--reps N] [--csv FILE] [--tag NAME] "
                   "[--max-carbons N] [--wasted]\n",
                   argv[0]);
      return 1;
    }
  }
  if (!o.kernels && !o.systems && !o.scaling) o.kernels = o.systems = true;

  std::printf("libecpint benchmark, LIBECPINT_MAX_UNROL = %d\n", LIBECPINT_MAX_UNROL);
  if (o.kernels) kernels(o.reps > 0 ? o.reps : 11);
  if (o.systems) systems(o.reps > 0 ? o.reps : 5);
  if (o.scaling) scaling(o.reps > 0 ? o.reps : 3, o.max_carbons, o.wasted);
  if (!o.csv.empty()) write_csv(o.csv, o.tag);
  return 0;
}
