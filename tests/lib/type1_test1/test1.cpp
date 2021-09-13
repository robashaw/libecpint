#include "api.hpp"
#include "testutil.hpp"
#include <iostream>

double g_coords[24] = {
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 5.0, 0.0, 0.0, 5.0, 0.0, 0.0, 5.0, 0.0, 0.0, 5.0,
};

double g_exps[12] = {2.6130000, 0.5736000, 0.2014000, 7.8600000,
                     0.7387000, 0.2081000, 2.6130000, 0.5736000,
                     0.2014000, 7.8600000, 0.7387000, 0.2081000};

double g_coefs[12] = {-0.5110463, 1.2701236, 1.0, -0.0555167, 1.0115982, 1.0,
                      -0.5110463, 1.2701236, 1.0, -0.0555167, 1.0115982, 1.0};

int g_ams[8] = {0, 0, 1, 1, 0, 0, 1, 1};
int g_lengths[8] = {2, 1, 2, 1, 2, 1, 2, 1};

double u_coords[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 5.0};

double u_exps[32] = {
    152.8553033, 82.1424792,  83.7154800, 23.8557161,  4.3128823,   212.5573747,
    266.4884712, 139.9222477, 50.1097659, 14.5537276,  3.4828623,   711.5242175,
    144.6708689, 32.9246992,  9.9103877,  3.1328926,   152.8553033, 82.1424792,
    83.7154800,  23.8557161,  4.3128823,  212.5573747, 266.4884712, 139.9222477,
    50.1097659,  14.5537276,  3.4828623,  711.5242175, 144.6708689, 32.9246992,
    9.9103877,   3.1328926};

double u_coefs[32] = {3.0,         11.2465621,  250.7412812, 139.1606543,
                      41.2897981,  5.0,         5.2787358,   631.7135166,
                      305.1649809, 106.4807615, 17.0215765,  -10.0,
                      -99.0606669, -35.2711767, -11.8151947, -1.0453382,
                      3.0,         11.2465621,  250.7412812, 139.1606543,
                      41.2897981,  5.0,         5.2787358,   631.7135166,
                      305.1649809, 106.4807615, 17.0215765,  -10.0,
                      -99.0606669, -35.2711767, -11.8151947, -1.0453382};

int u_ams[32] = {0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2,
                 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2};
int u_ns[32] = {0, 1, 2, 2, 2, 0, 1, 2, 2, 2, 2, 1, 2, 2, 2, 2,
                0, 1, 2, 2, 2, 0, 1, 2, 2, 2, 2, 1, 2, 2, 2, 2};
int u_lengths[2] = {16, 16};

int main(int argc, char* argv[]) {
  using namespace libecpint;

  ECPIntegrator factory;
  factory.set_gaussian_basis(8, g_coords, g_exps, g_coefs, g_ams, g_lengths);
  factory.set_ecp_basis(2, u_coords, u_exps, u_coefs, u_ams, u_ns, u_lengths);
  factory.init();
  factory.compute_integrals();

  std::vector<double> flat_result;
  std::shared_ptr<std::vector<double>> ints = factory.get_integrals();
  for (int i = 0; i < ints->size(); i++) flat_result.push_back((*ints)[i]);

#ifdef WRITE_NEW_BENCHMARK
  for (auto& v : flat_result)
    std::cout << std::setprecision(15) << v << std::endl;
#endif

  return check_file<double>("type1_test1.output", flat_result);
}
