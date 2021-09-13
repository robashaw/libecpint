#include "api.hpp"
#include "config.hpp"
#include "testutil.hpp"
#include <chrono>
#include <cmath>
#include <ctime>
#include <iostream>
#include <string>

double ag_exps[61] = {
    1.800750E+02, 2.189870E+01, 1.386700E+01, 6.142630E+00, 1.438140E+00,
    6.483820E-01, 1.288820E-01, 4.573800E-02, 1.800750E+02, 2.189870E+01,
    1.386700E+01, 6.142630E+00, 1.438140E+00, 6.483820E-01, 1.288820E-01,
    4.573800E-02, 1.800750E+02, 2.189870E+01, 1.386700E+01, 6.142630E+00,
    1.438140E+00, 6.483820E-01, 1.288820E-01, 4.573800E-02, 4.573800E-02,
    1.187510E+01, 8.002450E+00, 2.017660E+00, 9.542300E-01, 4.231180E-01,
    1.358850E-01, 4.540100E-02, 1.187510E+01, 8.002450E+00, 2.017660E+00,
    9.542300E-01, 4.231180E-01, 1.358850E-01, 4.540100E-02, 1.187510E+01,
    8.002450E+00, 2.017660E+00, 9.542300E-01, 4.231180E-01, 1.358850E-01,
    4.540100E-02, 4.540100E-02, 2.643200E+01, 1.103450E+01, 2.737870E+00,
    1.195750E+00, 4.820420E-01, 1.729080E-01, 2.643200E+01, 1.103450E+01,
    2.737870E+00, 1.195750E+00, 4.820420E-01, 1.729080E-01, 1.729080E-01,
    1.329900E+00};

double ag_coefs[61] = {8.490000E-04,
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

double ag_2_centers[6] = {0.0, 0.0, -2.4, 0.0, 0.0, 2.4};

double ag_4_centers[12] = {0.0, 0.0, 4.5, 0.0, 0.0,  -4.5,
                           0.0, 2.4, 0.0, 0.0, -2.4, 0.0};

double ag_6_centers[18] = {0.0, 3.0,  0.0, 2.6,  -1.5, 0.0, 5.0,  2.9,  0.0,
                           0.0, -5.8, 0.0, -4.0, 2.9,  0.0, -2.6, -1.5, 0.0};

int ag_charges[6] = {47, 47, 47, 47, 47, 47};
int ag_lengths[12] = {8, 8, 8, 1, 7, 7, 7, 1, 6, 6, 1, 1};
int ag_ams[12] = {0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 3};

using namespace libecpint;

std::vector<double> run_bench(int n_centers, double* centers,
                              std::string& share_dir) {
  std::cout << "N: " << n_centers << std::endl;
  std::vector<std::string> names;
  std::chrono::steady_clock::time_point first_time, last_time;
  std::vector<double> timings;

  for (int i = 0; i < n_centers; i++) names.push_back("ecp28mdf");
  double g_coords[3 * 12 * n_centers];
  double g_coefs[61 * n_centers];
  double g_exps[61 * n_centers];
  int g_lengths[12 * n_centers];
  int g_ams[12 * n_centers];
  double x, y, z;
  int gctr = 0;
  int charges[n_centers];
  for (int i = 0; i < n_centers; i++) {
    charges[i] = 47;
    x = centers[3 * i];
    y = centers[3 * i + 1];
    z = centers[3 * i + 2];
    int lctr = 0;
    for (int j = 0; j < 12; j++) {
      g_coords[36 * i + 3 * j] = x;
      g_coords[36 * i + 3 * j + 1] = y;
      g_coords[36 * i + 3 * j + 2] = z;
      g_lengths[12 * i + j] = ag_lengths[j];
      g_ams[12 * i + j] = ag_ams[j];

      for (int k = 0; k < ag_lengths[j]; k++) {
        g_coefs[gctr] = ag_coefs[lctr];
        g_exps[gctr] = ag_exps[lctr];
        gctr++;
        lctr++;
      }
    }
  }
  std::cout << std::setw(20) << "Initialisation... ";
  first_time = std::chrono::steady_clock::now();
  ECPIntegrator factory;
  factory.set_gaussian_basis(12 * n_centers, g_coords, g_exps, g_coefs, g_ams,
                             g_lengths);
  factory.set_ecp_basis_from_library(n_centers, centers, charges, names,
                                     share_dir);

  int init_val = LIBECPINT_MAX_L > 4 ? 2 : (LIBECPINT_MAX_L > 3 ? 1 : 0);
  factory.init(init_val);
  last_time = std::chrono::steady_clock::now();
  timings.push_back(
      std::chrono::duration<double, std::deci>(last_time - first_time).count() /
      10.0);
  std::cout << std::setw(20) << "done. TIME TAKEN: " << std::setw(15)
            << timings[0] << " seconds" << std::endl;

  std::cout << std::setw(20) << "Integrals... ";
  first_time = std::chrono::steady_clock::now();
  factory.compute_integrals();
  last_time = std::chrono::steady_clock::now();
  timings.push_back(
      std::chrono::duration<double, std::deci>(last_time - first_time).count() /
      10.0);
  std::cout << std::setw(20) << "done. TIME TAKEN: " << std::setw(15)
            << timings[1] << " seconds" << std::endl;

  if (LIBECPINT_MAX_L > 3) {
    std::cout << std::setw(20) << "1st derivs... ";
    first_time = std::chrono::steady_clock::now();
    factory.compute_first_derivs();
    last_time = std::chrono::steady_clock::now();
    timings.push_back(
        std::chrono::duration<double, std::deci>(last_time - first_time)
            .count() /
        10.0);
    std::cout << std::setw(20) << "done. TIME TAKEN: " << std::setw(15)
              << timings[2] << " seconds" << std::endl;
    if (LIBECPINT_MAX_L > 4) {
      if (n_centers < 4) {
        std::cout << std::setw(20) << "2nd derivs... ";
        first_time = std::chrono::steady_clock::now();
        factory.compute_second_derivs();
        last_time = std::chrono::steady_clock::now();
        timings.push_back(
            std::chrono::duration<double, std::deci>(last_time - first_time)
                .count() /
            10.0);
        std::cout << std::setw(20) << "done. TIME TAKEN: " << std::setw(15)
                  << timings[3] << " seconds" << std::endl;
      }
    } else {
      std::cout << "Insufficient LIBECPINT_MAX_L for 2nd derivatives"
                << std::endl;
    }
  } else {
    std::cout << "Insufficient LIBECPINT_MAX_L for derivatives" << std::endl;
  }
  std::cout << std::endl;
  return timings;
}

int main(int argc, char* argv[]) {
  std::string share_dir = argv[1];

  std::vector<double> t2 = run_bench(2, ag_2_centers, share_dir);
  std::vector<double> t4 = run_bench(4, ag_4_centers, share_dir);
  std::vector<double> t6 = run_bench(6, ag_6_centers, share_dir);

  double nsum = 0.0, nsum2 = 0.0;
  double ysum_int = 0.0, ysum_deriv = 0.0;
  double ynsum_int = 0.0, ynsum_deriv = 0.0;

  double lnN = std::log(2);
  double lnY = std::log(t2[1]);
  nsum += lnN;
  nsum2 += lnN * lnN;
  ysum_int += lnY;
  ynsum_int += lnN * lnY;
  lnY = std::log(t2[2]);
  ysum_deriv += lnY;
  ynsum_deriv += lnN * lnY;

  lnN = std::log(4);
  lnY = std::log(t4[1]);
  nsum += lnN;
  nsum2 += lnN * lnN;
  ysum_int += lnY;
  ynsum_int += lnN * lnY;
  lnY = std::log(t4[2]);
  ysum_deriv += lnY;
  ynsum_deriv += lnN * lnY;

  lnN = std::log(6);
  lnY = std::log(t6[1]);
  nsum += lnN;
  nsum2 += lnN * lnN;
  ysum_int += lnY;
  ynsum_int += lnN * lnY;
  lnY = std::log(t6[2]);
  ysum_deriv += lnY;
  ynsum_deriv += lnN * lnY;

  nsum2 = nsum2 - (nsum * nsum / 3.0);
  ysum_int = ynsum_int - (nsum * ysum_int / 3.0);
  ysum_deriv = ynsum_deriv - (nsum * ysum_deriv / 3.0);
  std::cout << std::setw(30) << "Scaling of integrals: N**"
            << std::setprecision(3) << ysum_int / nsum2 << std::endl;
  std::cout << std::setw(30) << "Scaling of 1st derivs: N**"
            << std::setprecision(3) << ysum_deriv / nsum2 << std::endl;
}
