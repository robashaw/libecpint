#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "api.hpp"
#include "config.hpp"
#include "mathutil.hpp"

// Silver valence basis: 4s4p3d1f = 12 shells, 61 primitives total
// Shells: s(8) s(8) s(8) s(1) p(7) p(7) p(7) p(1) d(6) d(6) d(1) f(1)
// ECP28MDF: 28-electron ECP, maxL=4, channels l=0..3 plus local l=4
static const int SHELLS_PER_ATOM = 12;
static const int PRIMS_PER_ATOM = 61;
static const double TWO_C_TOLERANCE = 1e-12;  // screening tolerance for shell pairs

static const double ag_exps[61] = {
    1.800750E+02, 2.189870E+01, 1.386700E+01, 6.142630E+00, 1.438140E+00, 6.483820E-01,
    1.288820E-01, 4.573800E-02, 1.800750E+02, 2.189870E+01, 1.386700E+01, 6.142630E+00,
    1.438140E+00, 6.483820E-01, 1.288820E-01, 4.573800E-02, 1.800750E+02, 2.189870E+01,
    1.386700E+01, 6.142630E+00, 1.438140E+00, 6.483820E-01, 1.288820E-01, 4.573800E-02,
    4.573800E-02, 1.187510E+01, 8.002450E+00, 2.017660E+00, 9.542300E-01, 4.231180E-01,
    1.358850E-01, 4.540100E-02, 1.187510E+01, 8.002450E+00, 2.017660E+00, 9.542300E-01,
    4.231180E-01, 1.358850E-01, 4.540100E-02, 1.187510E+01, 8.002450E+00, 2.017660E+00,
    9.542300E-01, 4.231180E-01, 1.358850E-01, 4.540100E-02, 4.540100E-02, 2.643200E+01,
    1.103450E+01, 2.737870E+00, 1.195750E+00, 4.820420E-01, 1.729080E-01, 2.643200E+01,
    1.103450E+01, 2.737870E+00, 1.195750E+00, 4.820420E-01, 1.729080E-01, 1.729080E-01,
    1.329900E+00};

static const double ag_coefs[61] = {8.490000E-04,
                                    -6.545000E-02,
                                    2.977650E-01,
                                    -7.531210E-01,
                                    8.811750E-01,
                                    4.351760E-01,
                                    1.473800E-02,
                                    -1.399000E-03,
                                    -2.030000E-04,
                                    1.572300E-02,
                                    -7.922900E-02,
                                    2.226340E-01,
                                    -3.491960E-01,
                                    -2.559780E-01,
                                    5.486660E-01,
                                    6.181030E-01,
                                    -8.620000E-04,
                                    5.254600E-02,
                                    -2.081000E-01,
                                    5.249720E-01,
                                    -1.270173E+00,
                                    5.393730E-01,
                                    1.653391E+00,
                                    -1.583975E+00,
                                    1.0,
                                    1.162480E-01,
                                    -3.072860E-01,
                                    5.157360E-01,
                                    5.031040E-01,
                                    1.420950E-01,
                                    5.153000E-03,
                                    1.040000E-04,
                                    -2.828400E-02,
                                    7.834700E-02,
                                    -1.567410E-01,
                                    -1.886250E-01,
                                    6.356500E-02,
                                    5.817780E-01,
                                    4.999110E-01,
                                    -5.285700E-02,
                                    1.470390E-01,
                                    -3.129730E-01,
                                    -3.708280E-01,
                                    4.236290E-01,
                                    7.536550E-01,
                                    7.087000E-02,
                                    1.0,
                                    3.479000E-03,
                                    -1.384800E-02,
                                    2.545990E-01,
                                    4.498490E-01,
                                    3.757380E-01,
                                    1.458790E-01,
                                    -4.733000E-03,
                                    1.907000E-02,
                                    -4.332980E-01,
                                    -4.445680E-01,
                                    4.914420E-01,
                                    5.728660E-01,
                                    1.0,
                                    1.0};

static const int ag_shell_lengths[12] = {8, 8, 8, 1, 7, 7, 7, 1, 6, 6, 1, 1};
static const int ag_shell_ams[12] = {0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 3};

static std::vector<double> make_grid(int n, double spacing) {
  int side = 1;
  while (side * side * side < n) side++;

  std::vector<double> coords;
  coords.reserve(3 * n);
  int placed = 0;
  for (int ix = 0; ix < side && placed < n; ix++) {
    for (int iy = 0; iy < side && placed < n; iy++) {
      for (int iz = 0; iz < side && placed < n; iz++) {
        coords.push_back(ix * spacing);
        coords.push_back(iy * spacing);
        coords.push_back(iz * spacing);
        placed++;
      }
    }
  }
  return coords;
}

// Mirrors the screening logic in compute_integrals() to count how many
// compute_shell_pair calls actually happen after screening.
static long long count_ecp_calls(libecpint::ECPIntegrator& factory) {
  using namespace libecpint;
  initFactorials();

  int maxLB = factory.maxLB;
  double min_alpha = factory.min_alpha;
  auto& shells = factory.shells;
  auto& ecps = factory.ecps;
  int nshells = shells.size();

  double thresh =
      FAST_POW[maxLB + 3]((maxLB + 3.0) / min_alpha) * FAST_POW[3](M_PI / (2 * maxLB + 3.0));
  thresh /= FAST_POW[maxLB](2.0 * M_EULER);
  thresh = TWO_C_TOLERANCE / std::sqrt(thresh);

  double product_thresh = TWO_C_TOLERANCE;

  long long calls = 0;
  for (int s1 = 0; s1 < nshells; ++s1) {
    GaussianShell& shellA = shells[s1];
    std::vector<std::pair<int, double>> ns;

    for (int i = 0; i < ecps.getN(); i++) {
      ECP& U = ecps.getECP(i);
      double acx = shellA.center()[0] - U.center_[0];
      double acy = shellA.center()[1] - U.center_[1];
      double acz = shellA.center()[2] - U.center_[2];
      double A2 = acx * acx + acy * acy + acz * acz;
      double sb = shell_bound(shellA.l, shellA.min_exp, A2, U.min_exp);
      if (sb > thresh) ns.push_back({i, sb});
    }

    if (ns.size() > 0) {
      for (int s2 = 0; s2 <= s1; ++s2) {
        GaussianShell& shellB = shells[s2];
        for (auto& [i, sbA] : ns) {
          ECP& U = ecps.getECP(i);
          double bcx = shellB.center()[0] - U.center_[0];
          double bcy = shellB.center()[1] - U.center_[1];
          double bcz = shellB.center()[2] - U.center_[2];
          double B2 = bcx * bcx + bcy * bcy + bcz * bcz;
          double sbB = shell_bound(shellB.l, shellB.min_exp, B2, U.min_exp);
          if (sbA * sbB > product_thresh) calls++;
        }
      }
    }
  }
  return calls;
}

static void dump_screening_diagnostics(libecpint::ECPIntegrator& factory) {
  using namespace libecpint;
  initFactorials();

  int maxLB = factory.maxLB;
  double min_alpha = factory.min_alpha;
  auto& shells = factory.shells;
  auto& ecps = factory.ecps;

  double thresh =
      FAST_POW[maxLB + 3]((maxLB + 3.0) / min_alpha) * FAST_POW[3](M_PI / (2 * maxLB + 3.0));
  thresh /= FAST_POW[maxLB](2.0 * M_EULER);
  thresh = TWO_C_TOLERANCE / std::sqrt(thresh);

  std::cout << "\n=== Screening Diagnostics ===\n";
  std::cout << "thresh (independent) = " << std::scientific << thresh << "\n";
  std::cout << "thresh^2 (effective product) = " << thresh * thresh << "\n";
  std::cout << "TWO_C_TOLERANCE = " << TWO_C_TOLERANCE << "\n\n";

  // For the first ECP, show shell_bound vs distance for each shell type
  ECP& U0 = ecps.getECP(0);
  std::cout << "shell_bound values vs distance (ECP #0, eta_min=" << U0.min_exp << "):\n";
  std::cout << std::setw(5) << "L" << std::setw(10) << "alpha" << std::setw(10) << "Dist";
  std::cout << std::setw(15) << "sb" << std::setw(15) << "sb*sb_near" << std::setw(10) << "ind?"
            << std::setw(10) << "prod?" << "\n";
  std::cout << std::string(75, '-') << "\n";

  // Collect unique (l, min_exp) shell types
  struct ShType {
    int l;
    double alpha;
  };
  std::vector<ShType> types;
  for (auto& sh : shells) {
    bool found = false;
    for (auto& t : types)
      if (t.l == sh.l && std::abs(t.alpha - sh.min_exp) < 1e-10) {
        found = true;
        break;
      }
    if (!found) types.push_back({sh.l, sh.min_exp});
  }

  for (auto& t : types) {
    double sb_near = shell_bound(t.l, t.alpha, 0.0, U0.min_exp);
    for (double dist = 0; dist <= 25; dist += 2.5) {
      double A2 = dist * dist;
      double sb = shell_bound(t.l, t.alpha, A2, U0.min_exp);
      double product = sb_near * sb;
      std::cout << std::setw(5) << t.l << std::setw(10) << std::scientific << std::setprecision(2)
                << t.alpha << std::setw(10) << std::fixed << std::setprecision(1) << dist
                << std::setw(15) << std::scientific << std::setprecision(3) << sb << std::setw(15)
                << product << std::setw(10) << (sb > thresh ? "pass" : "SCREEN") << std::setw(10)
                << (product > TWO_C_TOLERANCE ? "pass" : "SCREEN") << "\n";
    }
    std::cout << "\n";
  }
}

struct BenchResult {
  int n_atoms;
  int n_shells;
  int n_cart;
  long long n_shell_pairs;
  long long n_ecp_calls;
  double t_init;
  double t_integrals;
};

static BenchResult run_bench(int n_atoms, double spacing, const std::string& share_dir) {
  auto coords = make_grid(n_atoms, spacing);

  int total_shells = SHELLS_PER_ATOM * n_atoms;
  int total_prims = PRIMS_PER_ATOM * n_atoms;

  std::vector<double> g_coords(3 * total_shells);
  std::vector<double> g_exps(total_prims);
  std::vector<double> g_coefs(total_prims);
  std::vector<int> g_lengths(total_shells);
  std::vector<int> g_ams(total_shells);
  std::vector<int> charges(n_atoms, 47);
  std::vector<std::string> names(n_atoms, "ecp28mdf");

  int pctr = 0;
  for (int i = 0; i < n_atoms; i++) {
    double cx = coords[3 * i];
    double cy = coords[3 * i + 1];
    double cz = coords[3 * i + 2];
    int lctr = 0;

    for (int j = 0; j < SHELLS_PER_ATOM; j++) {
      int si = SHELLS_PER_ATOM * i + j;
      g_coords[3 * si] = cx;
      g_coords[3 * si + 1] = cy;
      g_coords[3 * si + 2] = cz;
      g_lengths[si] = ag_shell_lengths[j];
      g_ams[si] = ag_shell_ams[j];

      for (int k = 0; k < ag_shell_lengths[j]; k++) {
        g_exps[pctr] = ag_exps[lctr];
        g_coefs[pctr] = ag_coefs[lctr];
        pctr++;
        lctr++;
      }
    }
  }

  BenchResult res;
  res.n_atoms = n_atoms;
  res.n_shells = total_shells;
  res.n_shell_pairs = (long long)total_shells * (total_shells + 1) / 2;

  auto t0 = std::chrono::steady_clock::now();
  libecpint::ECPIntegrator factory;
  factory.set_gaussian_basis(total_shells, g_coords.data(), g_exps.data(), g_coefs.data(),
                             g_ams.data(), g_lengths.data());
  factory.set_ecp_basis_from_library(n_atoms, coords.data(), charges.data(), names, share_dir);
  factory.init(0);
  auto t1 = std::chrono::steady_clock::now();
  res.t_init = std::chrono::duration<double>(t1 - t0).count();
  res.n_cart = factory.ncart;

  res.n_ecp_calls = count_ecp_calls(factory);

  t0 = std::chrono::steady_clock::now();
  factory.compute_integrals();
  t1 = std::chrono::steady_clock::now();
  res.t_integrals = std::chrono::duration<double>(t1 - t0).count();

  return res;
}

int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " <share_dir> [spacing]\n";
    return 1;
  }
  std::string share_dir = argv[1];
  double spacing = argc > 2 ? std::atof(argv[2]) : 5.0;

  std::vector<int> sizes = {2, 4, 8, 12, 16, 20, 25};

  std::cout << "Scaling test: Ag clusters, spacing=" << spacing << " angstrom"
            << ", LIBECPINT_MAX_L=" << LIBECPINT_MAX_L << "\n\n";

  std::cout << std::fixed << std::setprecision(4);
  std::cout << std::setw(6) << "Natom" << std::setw(8) << "Nshell" << std::setw(8) << "Ncart"
            << std::setw(10) << "Pairs" << std::setw(12) << "ECPcalls" << std::setw(12) << "Init(s)"
            << std::setw(12) << "Integ(s)" << std::setw(12) << "us/pair" << std::setw(12)
            << "us/call" << "\n";
  std::cout << std::string(92, '-') << "\n";

  // Run screening diagnostic on a small 2-atom system
  {
    auto diag_res = run_bench(2, spacing, share_dir);
    // Need to rebuild factory to call dump_screening_diagnostics
    auto coords2 = make_grid(2, spacing);
    int ts = SHELLS_PER_ATOM * 2;
    int tp = PRIMS_PER_ATOM * 2;
    std::vector<double> gc(3 * ts), ge(tp), gco(tp);
    std::vector<int> gl(ts), ga(ts), ch(2, 47);
    std::vector<std::string> nm(2, "ecp28mdf");
    int pc = 0;
    for (int i = 0; i < 2; i++) {
      int lc = 0;
      for (int j = 0; j < SHELLS_PER_ATOM; j++) {
        int si = SHELLS_PER_ATOM * i + j;
        gc[3 * si] = coords2[3 * i];
        gc[3 * si + 1] = coords2[3 * i + 1];
        gc[3 * si + 2] = coords2[3 * i + 2];
        gl[si] = ag_shell_lengths[j];
        ga[si] = ag_shell_ams[j];
        for (int k = 0; k < ag_shell_lengths[j]; k++) {
          ge[pc] = ag_exps[lc];
          gco[pc] = ag_coefs[lc];
          pc++;
          lc++;
        }
      }
    }
    libecpint::ECPIntegrator diag_factory;
    diag_factory.set_gaussian_basis(ts, gc.data(), ge.data(), gco.data(), ga.data(), gl.data());
    diag_factory.set_ecp_basis_from_library(2, coords2.data(), ch.data(), nm, share_dir);
    diag_factory.init(0);
    dump_screening_diagnostics(diag_factory);
  }

  std::vector<BenchResult> results;
  for (int n : sizes) {
    std::cerr << "Running N=" << n << "..." << std::flush;
    auto res = run_bench(n, spacing, share_dir);

    double us_per_pair = res.t_integrals * 1e6 / res.n_shell_pairs;
    double us_per_call = res.n_ecp_calls > 0 ? res.t_integrals * 1e6 / res.n_ecp_calls : 0;

    std::cout << std::setw(6) << res.n_atoms << std::setw(8) << res.n_shells << std::setw(8)
              << res.n_cart << std::setw(10) << res.n_shell_pairs << std::setw(12)
              << res.n_ecp_calls << std::setw(12) << res.t_init << std::setw(12) << res.t_integrals
              << std::setw(12) << us_per_pair << std::setw(12) << us_per_call << "\n";
    std::cerr << " done (" << std::setprecision(1) << res.t_integrals << "s)\n";
    results.push_back(res);
  }

  // Scaling fits (skip smallest point)
  std::cout << "\n--- Scaling Analysis ---\n";
  auto fit_slope = [&](auto get_y) {
    int n_fit = 0;
    double sx = 0, sy = 0, sxx = 0, sxy = 0;
    for (size_t i = 1; i < results.size(); i++) {
      double lx = std::log(results[i].n_atoms);
      double ly = std::log(get_y(results[i]));
      sx += lx;
      sy += ly;
      sxx += lx * lx;
      sxy += lx * ly;
      n_fit++;
    }
    return (n_fit * sxy - sx * sy) / (n_fit * sxx - sx * sx);
  };

  double slope_t = fit_slope([](const BenchResult& r) { return r.t_integrals; });
  double slope_c = fit_slope([](const BenchResult& r) { return (double)r.n_ecp_calls; });

  std::cout << "Integral wall time: N^" << std::setprecision(2) << slope_t
            << " (fit over N=" << results[1].n_atoms << ".." << results.back().n_atoms << ")\n";
  std::cout << "ECP call count:     N^" << std::setprecision(2) << slope_c << "\n";
  std::cout << "Avg cost per call:  N^" << std::setprecision(2) << slope_t - slope_c
            << " (should be ~0 if calls dominate)\n";

  // Screening analysis
  std::cout << "\n--- Screening Analysis ---\n";
  std::cout << std::setw(6) << "Natom" << std::setw(12) << "MaxCalls" << std::setw(12)
            << "ActualCalls" << std::setw(10) << "Screened" << std::setw(12) << "calls/pair"
            << "\n";
  for (auto& r : results) {
    long long max_calls = r.n_shell_pairs * r.n_atoms;
    double screen_pct = 100.0 * (1.0 - (double)r.n_ecp_calls / max_calls);
    std::cout << std::setw(6) << r.n_atoms << std::setw(12) << max_calls << std::setw(12)
              << r.n_ecp_calls << std::setw(9) << std::setprecision(1) << screen_pct << "%"
              << std::setw(12) << std::setprecision(2) << (double)r.n_ecp_calls / r.n_shell_pairs
              << "\n";
  }

  return 0;
}
