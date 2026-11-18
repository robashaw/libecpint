/*
 * Unit tests of the test support code that is not about accuracy: the readers of the manifest,
 * xyz and Gaussian94 basis files, the system hash, and the round trip of reference files.
 */
#include <gtest/gtest.h>

#include <cmath>
#include <fstream>
#include <stdexcept>
#include <string>

#include "reference_io.hpp"
#include "system.hpp"

using namespace testsupport;

namespace {

std::string write_file(const std::string& name, const std::string& text) {
  const std::string path = ::testing::TempDir() + name;
  std::ofstream(path) << text;
  return path;
}

// int (x^l sum_i c_i exp(-a_i r^2))^2 d^3r by the trapezium rule (exact to rounding for these
// smooth, fast-decaying integrands), independently of the normalisation formula
double self_overlap(const GbsShell& sh) {
  auto integral = [](const int power, const double p) {
    const double h = 0.002 / std::sqrt(p), xmax = 40.0 / std::sqrt(p);
    double sum = 0.0;
    for (double x = -xmax; x <= xmax; x += h) sum += std::pow(x, power) * std::exp(-p * x * x);
    return sum * h;
  };
  double S = 0.0;
  for (size_t i = 0; i < sh.exps.size(); i++)
    for (size_t j = 0; j < sh.exps.size(); j++) {
      const double p = sh.exps[i] + sh.exps[j];
      S += sh.coefs[i] * sh.coefs[j] * integral(2 * sh.l, p) * std::pow(integral(0, p), 2);
    }
  return S;
}

}  // namespace

TEST(Manifest, ParsesKeysAndComments) {
  const std::string path = write_file("manifest.txt",
                                      "# comment line\n"
                                      "\n"
                                      "name=a xyz=x.xyz basis=b.gbs ecp=Ag:lib/ecp28mdf,Li:data/t "
                                      "deriv=2 live=yes  # trailing comment\n"
                                      "name=b xyz=y.xyz basis=c.gbs ecp=I:lib/def2-ecp\n");
  const auto specs = read_manifest(path);
  ASSERT_EQ(specs.size(), 2u);
  EXPECT_EQ(specs[0].name, "a");
  EXPECT_EQ(specs[0].xyz, "x.xyz");
  EXPECT_EQ(specs[0].basis, "b.gbs");
  ASSERT_EQ(specs[0].ecps.size(), 2u);
  EXPECT_EQ(specs[0].ecps[0], std::make_pair(std::string("Ag"), std::string("lib/ecp28mdf")));
  EXPECT_EQ(specs[0].ecps[1], std::make_pair(std::string("Li"), std::string("data/t")));
  EXPECT_EQ(specs[0].deriv, 2);
  EXPECT_TRUE(specs[0].live);
  EXPECT_EQ(specs[1].deriv, 0);
  EXPECT_FALSE(specs[1].live);
}

TEST(Manifest, RejectsUnknownKeys) {
  EXPECT_THROW(read_manifest(write_file("bad.txt", "name=a colour=red\n")), std::runtime_error);
}

TEST(Xyz, AngstromAndBohr) {
  const auto ang = read_xyz(write_file("a.xyz", "2\nangstrom\nAg 1.0 0.0 -2.0\nli 0 0.5 0\n"));
  ASSERT_EQ(ang.size(), 2u);
  EXPECT_EQ(ang[0].symbol, "Ag");
  EXPECT_EQ(ang[0].Z, 47);
  EXPECT_EQ(ang[1].Z, 3);
  EXPECT_NEAR(ang[0].xyz[0], 1.8897261246, 1e-9);
  EXPECT_NEAR(ang[0].xyz[2], -2.0 * 1.8897261246, 1e-9);
  const auto bohr = read_xyz(write_file("b.xyz", "1\nunits: Bohr\nAg 1.0 0.0 -2.0\n"));
  EXPECT_DOUBLE_EQ(bohr[0].xyz[0], 1.0);
  EXPECT_DOUBLE_EQ(bohr[0].xyz[2], -2.0);
  EXPECT_THROW(read_xyz(write_file("c.xyz", "2\n\nAg 0 0 0\n")), std::runtime_error);
  EXPECT_THROW(read_xyz(write_file("d.xyz", "1\n\nQq 0 0 0\n")), std::runtime_error);
}

TEST(AtomicNumber, Symbols) {
  EXPECT_EQ(atomic_number("H"), 1);
  EXPECT_EQ(atomic_number("ag"), 47);
  EXPECT_EQ(atomic_number("AU"), 79);
  EXPECT_EQ(atomic_number("Xx"), 0);
}

TEST(Gbs, KeysShellsAndNormalisation) {
  const auto basis = read_gbs(write_file("t.gbs",
                                         "! comment\n"
                                         "Ag     0\n"
                                         "S   2   1.00\n"
                                         "      2.5D+00   0.4\n"
                                         "      0.3       0.7\n"
                                         "F   1   1.00\n"
                                         "      0.8       1.0\n"
                                         "****\n"
                                         "2     0\n"
                                         "D   3   1.00\n"
                                         "      9.0       0.1\n"
                                         "      1.2d0     0.5\n"
                                         "      0.2       0.6\n"
                                         "****\n"));
  ASSERT_EQ(basis.size(), 2u);
  ASSERT_EQ(basis.count("Ag"), 1u);
  ASSERT_EQ(basis.count("2"), 1u);
  const auto& ag = basis.at("Ag");
  ASSERT_EQ(ag.size(), 2u);
  EXPECT_EQ(ag[0].l, 0);
  EXPECT_EQ(ag[1].l, 3);
  EXPECT_DOUBLE_EQ(ag[0].exps[0], 2.5);  // Fortran exponent
  const auto& d = basis.at("2")[0];
  EXPECT_EQ(d.l, 2);
  EXPECT_DOUBLE_EQ(d.exps[1], 1.2);
  // Contractions normalised for the x^l component, with the ratios of the coefficients set by
  // the primitive normalisation (2a/pi)^(3/4) (4a)^(l/2)
  for (const GbsShell* sh : {&ag[0], &ag[1], &d}) EXPECT_NEAR(self_overlap(*sh), 1.0, 1e-12);
  const double r = (0.4 / 0.7) * std::pow(2.5 / 0.3, 0.75);
  EXPECT_NEAR(ag[0].coefs[0] / ag[0].coefs[1], r, 1e-14 * r);
  EXPECT_THROW(read_gbs(write_file("u.gbs", "Ag 0\nSP 1 1.00\n 1.0 1.0 1.0\n****\n")),
               std::runtime_error);
}

TEST(SystemHash, ChangesWithInputs) {
  const System sys = load_system("testpair");
  System moved = sys, scaled = sys;
  EXPECT_EQ(sys.hash(), System(sys).hash());
  moved.atoms[0].xyz[1] += 1e-12;
  EXPECT_NE(moved.hash(), sys.hash());
  scaled.ecp_prims[0].d *= 1.0 + 1e-15;
  EXPECT_NE(scaled.hash(), sys.hash());
}

TEST(ReferenceFile, RoundTripIsExact) {
  ReferenceData data;
  data.system = "x";
  data.hash = "0123abcd";
  data.generator = "unit-test";
  data.deriv = 2;
  data.ncart = 2;
  data.natoms = 1;
  data.ints = {1.0 / 3.0, -2.5e-300, M_PI};
  data.gradient = {1e-17, -7.0, 0.1, 2.0, 3.0, 4.0};
  for (int k = 0; k < 18; k++) data.hessian.push_back(std::sin(1.0 + k) * std::pow(10.0, k - 9));
  data.ints_scale = {3.0, 1e-290, 7.25};
  data.gradient_scale = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0 / 7.0};
  for (int k = 0; k < 18; k++) data.hessian_scale.push_back(std::exp(0.7 * k));
  const std::string path = ::testing::TempDir() + "x.ref";
  write_reference(path, data);
  const ReferenceData back = read_reference(path);
  EXPECT_EQ(back.system, data.system);
  EXPECT_EQ(back.hash, data.hash);
  EXPECT_EQ(back.generator, data.generator);
  EXPECT_EQ(back.deriv, data.deriv);
  EXPECT_EQ(back.ncart, data.ncart);
  EXPECT_EQ(back.natoms, data.natoms);
  EXPECT_EQ(back.ints, data.ints);  // bitwise
  EXPECT_EQ(back.gradient, data.gradient);
  EXPECT_EQ(back.hessian, data.hessian);
  EXPECT_EQ(back.ints_scale, data.ints_scale);
  EXPECT_EQ(back.gradient_scale, data.gradient_scale);
  EXPECT_EQ(back.hessian_scale, data.hessian_scale);
  // empty sections (no second derivatives) followed by others
  data.hessian.clear();
  data.ints_scale.clear();
  write_reference(path, data);
  const ReferenceData back2 = read_reference(path);
  EXPECT_TRUE(back2.hessian.empty());
  EXPECT_TRUE(back2.ints_scale.empty());
  EXPECT_EQ(back2.gradient_scale, data.gradient_scale);
  EXPECT_EQ(back2.hessian_scale, data.hessian_scale);
}

TEST(ReferenceFile, MissingOrTruncatedThrows) {
  EXPECT_THROW(read_reference(::testing::TempDir() + "does-not-exist.ref"), std::runtime_error);
  EXPECT_THROW(read_reference(write_file("t.ref", "ncart 2\nintegrals 3\n1.0\n2.0\n")),
               std::runtime_error);
  EXPECT_THROW(read_reference(write_file("u.ref", "ncart 2\ncolour red\n")), std::runtime_error);
}

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
