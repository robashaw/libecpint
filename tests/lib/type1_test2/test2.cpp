#include "api.hpp"
#include "testutil.hpp"
#include <iostream>

double g_coords[66] = {
    0.0,   0.0,   0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
    0.0,   0.0,   0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
    0.0,   0.0,   0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
    0.0,   0.0,   0.0, 0.0,   0.0,   4.545, 0.0,   0.0,   4.545, 0.0,   0.0,
    4.545, 0.0,   0.0, 4.545, 0.0,   0.0,   4.545, 0.0,   0.0,   4.545, 0.0,
    0.0,   4.545, 0.0, 0.0,   4.545, 0.0,   0.0,   4.545, 0.0,   0.0,   4.545};

double g_exps[42] = {
    20.115299000,  12.193477000,   6.0735294368,  1.3174451569,  0.58596768244,
    0.13875427354, 0.048876985527, 8.6096650000,  7.3353260000,  1.6575296365,
    0.78159310216, 0.32384840661,  0.0540000,     4.1439490000,  3.5682570000,
    1.2345757130,  0.48190232338,  0.16490636769, 0.7248200,     445.90489176,
    23.336842412,  19.583446104,   8.5112089186,  2.1996161861,  1.0668970454,
    11.720547572,  1.7625986449,   0.28697476318, 0.11227609607, 20.499027254,
    10.558754576,  1.5015895485,   0.64597173095, 3.0288656771,  0.34506591841,
    0.11134513813, 51.235354920,   15.616239482,  4.5266002139,  2.0529808766,
    0.87640281623, 0.30900000000};

double g_coefs[42] = {-0.15910719389,
                      0.79105526778,
                      1.0,
                      1.0,
                      1.0,
                      1.0,
                      1.0,
                      0.50053018599,
                      -0.72681584494,
                      0.57315511417,
                      0.49579068859,
                      1.0,
                      1.0,
                      -0.37099566643,
                      0.40197233762,
                      0.46001988624,
                      0.46152130957882,
                      1.0,
                      1.0,
                      0.20037290389E-02,
                      -0.14989397324,
                      0.36436474913,
                      -0.75538138199,
                      0.82816359356,
                      0.42161048110,
                      -0.16160234511E-01,
                      0.35271525193,
                      1.0,
                      1.0,
                      0.077159335656,
                      -0.38370372780,
                      0.83554174719,
                      0.16787259488,
                      1.0,
                      1.0,
                      1.0,
                      0.40345130346E-02,
                      -0.50768108881E-02,
                      0.29151065356,
                      0.51145783605,
                      0.31232025297,
                      1.0};

int g_ams[22] = {0, 0, 0, 0, 0, 0, 1, 1, 1, 2, 2,
                 3, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2};
int g_lengths[22] = {2, 1, 1, 1, 1, 1, 4, 1, 1, 4, 1,
                     1, 6, 2, 1, 1, 4, 1, 1, 1, 5, 1};

double u_coords[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 4.545};

double u_exps[30] = {4.78982,   2.39491,   13.2051,   6.60255,   4.78982,
                     2.39491,   10.45202,  5.22601,   4.78982,   2.39491,
                     7.8511,    3.92555,   4.78982,   2.39491,   1.0,
                     40.033376, 17.300576, 8.85172,   15.720141, 15.208222,
                     8.294186,  7.753949,  13.817751, 13.587805, 6.94763,
                     6.960099,  18.522950, 18.251035, 7.557901,  7.597404};

double u_coefs[30] = {
    30.4900889,   5.17107381,   426.8466792, 37.007082885, -30.4900889,
    -5.17107381,  261.19958038, 26.96249604, -30.4900889,  -5.17107381,
    124.79066561, 16.30072573,  -30.4900889, -5.17107381,  0.0,
    49.989649,    281.006556,   61.416739,   67.416239,    134.807696,
    14.566548,    28.968422,    35.538756,   53.339759,    9.716466,
    14.977500,    -20.176618,   -26.088077,  -0.220434,    -0.221646};

int u_ams[30] = {3, 3, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 4,
                 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3};
int u_ns[30] = {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2};
int u_lengths[2] = {14, 16};

int main(int argc, char* argv[]) {
  using namespace libecpint;

  ECPIntegrator factory;
  factory.set_gaussian_basis(22, g_coords, g_exps, g_coefs, g_ams, g_lengths);
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

  return check_file<double>("type1_test2.output", flat_result);
}
