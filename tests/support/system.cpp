#include "system.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <tuple>

namespace testsupport {

namespace {
const char* const ELEMENTS[] = {
    "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne", "Na", "Mg", "Al", "Si", "P",
    "S",  "Cl", "Ar", "K",  "Ca", "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y",  "Zr", "Nb", "Mo", "Tc", "Ru", "Rh",
    "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I",  "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
    "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W",  "Re",
    "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
    "Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db",
    "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"};

constexpr double BOHR_PER_ANGSTROM = 1.0 / 0.529177210903;

std::string lower(std::string s) {
  std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return std::tolower(c); });
  return s;
}

std::string capitalised(std::string s) {
  s = lower(s);
  if (!s.empty()) s[0] = std::toupper(static_cast<unsigned char>(s[0]));
  return s;
}

std::ifstream open(const std::string& path) {
  std::ifstream file(path);
  if (!file) throw std::runtime_error("cannot open " + path);
  return file;
}

double parse_double(std::string s) {  // accepts Fortran exponents (1.0D-01)
  std::replace(s.begin(), s.end(), 'D', 'E');
  std::replace(s.begin(), s.end(), 'd', 'e');
  return std::stod(s);
}

double double_factorial(const int n) {  // n!! for n >= -1
  double v = 1.0;
  for (int k = n; k > 1; k -= 2) v *= k;
  return v;
}

const char* const ANGULAR = "spdfghik";
}  // namespace

std::string data_dir() { return LIBECPINT_TEST_DATA_DIR; }
std::string share_dir() { return LIBECPINT_SHARE_DIR; }

int atomic_number(const std::string& symbol) {
  const std::string s = capitalised(symbol);
  for (int z = 0; z < int(sizeof(ELEMENTS) / sizeof(ELEMENTS[0])); z++)
    if (s == ELEMENTS[z]) return z + 1;
  return 0;
}

std::vector<Atom> read_xyz(const std::string& path) {
  std::ifstream file = open(path);
  std::string line;
  std::getline(file, line);
  const int n = std::stoi(line);
  std::getline(file, line);
  const double scale = lower(line).find("bohr") != std::string::npos ? 1.0 : BOHR_PER_ANGSTROM;
  std::vector<Atom> atoms;
  for (int i = 0; i < n; i++) {
    if (!std::getline(file, line)) throw std::runtime_error(path + ": too few atoms");
    std::istringstream ss(line);
    Atom a;
    ss >> a.symbol >> a.xyz[0] >> a.xyz[1] >> a.xyz[2];
    if (!ss) throw std::runtime_error(path + ": cannot read atom line '" + line + "'");
    a.symbol = capitalised(a.symbol);
    a.Z = atomic_number(a.symbol);
    if (a.Z == 0) throw std::runtime_error(path + ": unknown element " + a.symbol);
    for (double& x : a.xyz) x *= scale;
    atoms.push_back(a);
  }
  return atoms;
}

std::map<std::string, std::vector<GbsShell>> read_gbs(const std::string& path) {
  std::ifstream file = open(path);
  std::map<std::string, std::vector<GbsShell>> basis;
  std::string line, key;
  auto next = [&]() {  // next non-blank, non-comment line
    while (std::getline(file, line)) {
      const auto p = line.find_first_not_of(" \t\r");
      if (p != std::string::npos && line[p] != '!') return true;
    }
    return false;
  };
  while (next()) {
    std::istringstream ss(line);
    std::string first;
    ss >> first;
    if (first == "****") {
      key.clear();
      continue;
    }
    if (key.empty()) {  // centre line: symbol or atom number, then 0
      key = std::isdigit(static_cast<unsigned char>(first[0])) ? first : capitalised(first);
      basis[key];
      continue;
    }
    // shell line: type, number of primitives, scale factor
    const std::string type = lower(first);
    int nprim;
    ss >> nprim;
    const char* p = std::strchr(ANGULAR, type[0]);
    if (type.size() != 1 || p == nullptr || !ss)
      throw std::runtime_error(path + ": unsupported shell line '" + line + "'");
    GbsShell sh;
    sh.l = p - ANGULAR;
    for (int i = 0; i < nprim; i++) {
      if (!next()) throw std::runtime_error(path + ": truncated shell");
      std::istringstream ps(line);
      std::string e, c;
      ps >> e >> c;
      sh.exps.push_back(parse_double(e));
      sh.coefs.push_back(parse_double(c));
    }
    sh.input_coefs = sh.coefs;
    // Normalise the primitives, then the contraction, for the x^l component
    const int l = sh.l;
    const double dfac = double_factorial(2 * l - 1);
    for (int i = 0; i < nprim; i++) {
      const double a = sh.exps[i];
      sh.coefs[i] *= std::pow(2.0 * a / M_PI, 0.75) * std::pow(4.0 * a, 0.5 * l) / std::sqrt(dfac);
    }
    double S = 0.0;
    for (int i = 0; i < nprim; i++)
      for (int j = 0; j < nprim; j++) {
        const double p2 = sh.exps[i] + sh.exps[j];
        S += sh.coefs[i] * sh.coefs[j] * std::pow(M_PI / p2, 1.5) * dfac / std::pow(2.0 * p2, l);
      }
    for (double& c : sh.coefs) c /= std::sqrt(S);
    basis[key].push_back(sh);
  }
  return basis;
}

std::vector<SystemSpec> read_manifest(const std::string& path) {
  std::ifstream file = open(path);
  std::vector<SystemSpec> specs;
  std::string line;
  while (std::getline(file, line)) {
    line = line.substr(0, line.find('#'));
    std::istringstream ss(line);
    std::string token;
    SystemSpec spec;
    bool any = false;
    while (ss >> token) {
      any = true;
      const auto eq = token.find('=');
      if (eq == std::string::npos)
        throw std::runtime_error(path + ": expected key=value, got " + token);
      const std::string key = token.substr(0, eq), value = token.substr(eq + 1);
      if (key == "name")
        spec.name = value;
      else if (key == "xyz")
        spec.xyz = value;
      else if (key == "basis")
        spec.basis = value;
      else if (key == "ecp") {
        std::istringstream vs(value);
        std::string item;
        while (std::getline(vs, item, ',')) {
          const auto colon = item.find(':');
          if (colon == std::string::npos) throw std::runtime_error(path + ": bad ecp " + item);
          spec.ecps.push_back({capitalised(item.substr(0, colon)), item.substr(colon + 1)});
        }
      } else if (key == "deriv")
        spec.deriv = std::stoi(value);
      else if (key == "live")
        spec.live = value == "yes";
      else
        throw std::runtime_error(path + ": unknown key " + key);
    }
    if (any) specs.push_back(spec);
  }
  return specs;
}

System load_system(const SystemSpec& spec) {
  System sys;
  sys.spec = spec;
  sys.atoms = read_xyz(data_dir() + "/" + spec.xyz);
  const auto basis = read_gbs(data_dir() + "/" + spec.basis);
  for (size_t i = 0; i < sys.atoms.size(); i++) {
    auto it = basis.find(std::to_string(i + 1));
    if (it == basis.end()) it = basis.find(sys.atoms[i].symbol);
    if (it == basis.end()) continue;  // no basis functions (e.g. an ECP alone)
    for (const GbsShell& sh : it->second)
      sys.shells.push_back({int(i), sh.l, sh.exps, sh.coefs, sh.input_coefs});
  }
  for (size_t i = 0; i < sys.atoms.size(); i++) {
    const Atom& atom = sys.atoms[i];
    for (const auto& [element, source] : spec.ecps) {
      if (element != atom.symbol) continue;
      std::string path;
      if (source.rfind("lib/", 0) == 0)
        path = share_dir() + "/raw/" + source.substr(4) + ".ecp";
      else if (source.rfind("data/", 0) == 0)
        path = data_dir() + "/ecp/" + source.substr(5) + ".ecp";
      else
        throw std::runtime_error(spec.name + ": ECP source must start lib/ or data/: " + source);
      libecpint::ECPBasis ecps;
      ecps.addECP_from_file(atom.Z, atom.xyz, path);
      for (const auto& g : ecps.getECP(0).gaussians)
        sys.ecp_prims.push_back({int(i), g.n + 2, g.l, g.a, g.d});  // n holds n - 2
    }
  }
  return sys;
}

namespace {
SystemSpec manifest_spec(const std::string& name) {
  for (const SystemSpec& spec : read_manifest(data_dir() + "/systems.txt"))
    if (spec.name == name) return spec;
  throw std::runtime_error("no system " + name + " in the manifest");
}
}  // namespace

System load_system(const std::string& name) { return load_system(manifest_spec(name)); }

void System::setup(libecpint::ECPIntegrator& f, const int deriv, const double tol,
                   const double radial_tol) const {
  std::vector<double> coords, exps, coefs, ecoords, eexps, ecoefs;
  std::vector<int> ams, lens, eams, ens, elens;
  for (const Shell& sh : shells) {
    coords.insert(coords.end(), atoms[sh.atom].xyz.begin(), atoms[sh.atom].xyz.end());
    exps.insert(exps.end(), sh.exps.begin(), sh.exps.end());
    coefs.insert(coefs.end(), sh.coefs.begin(), sh.coefs.end());
    ams.push_back(sh.l);
    lens.push_back(sh.exps.size());
  }
  for (size_t k = 0; k < ecp_prims.size(); k++) {
    const EcpPrimitive& p = ecp_prims[k];
    if (k == 0 || p.atom != ecp_prims[k - 1].atom) {
      ecoords.insert(ecoords.end(), atoms[p.atom].xyz.begin(), atoms[p.atom].xyz.end());
      elens.push_back(0);
    }
    eexps.push_back(p.a);
    ecoefs.push_back(p.d);
    eams.push_back(p.l);
    ens.push_back(p.n);
    elens.back()++;
  }
  f.set_gaussian_basis(shells.size(), coords.data(), exps.data(), coefs.data(), ams.data(),
                       lens.data());
  f.set_ecp_basis(elens.size(), ecoords.data(), eexps.data(), ecoefs.data(), eams.data(),
                  ens.data(), elens.data());
  f.init(deriv, tol, radial_tol);
}

std::string System::hash() const {
  uint64_t h = 14695981039346656037ull;  // FNV-1a over a canonical text form
  auto mix = [&](const std::string& s) {
    for (const unsigned char c : s) h = (h ^ c) * 1099511628211ull;
  };
  char buf[64];
  auto num = [&](const double x) {
    std::snprintf(buf, sizeof buf, "%.17g ", x);
    mix(buf);
  };
  for (const Atom& a : atoms) {
    mix(a.symbol + " ");
    for (const double x : a.xyz) num(x);
  }
  for (const Shell& sh : shells) {
    num(sh.atom);
    num(sh.l);
    for (size_t i = 0; i < sh.exps.size(); i++) {
      num(sh.exps[i]);
      num(sh.input_coefs[i]);
    }
  }
  std::vector<EcpPrimitive> prims = ecp_prims;
  std::sort(prims.begin(), prims.end(), [](const EcpPrimitive& x, const EcpPrimitive& y) {
    return std::tie(x.atom, x.l, x.n, x.a, x.d) < std::tie(y.atom, y.l, y.n, y.a, y.d);
  });
  for (const EcpPrimitive& p : prims) {
    num(p.atom);
    num(p.n);
    num(p.l);
    num(p.a);
    num(p.d);
  }
  std::snprintf(buf, sizeof buf, "%016llx", static_cast<unsigned long long>(h));
  return buf;
}

std::vector<int> System::integrator_atoms() const {
  // ECPIntegrator::init numbers atoms by first appearance among the shell centres, then the ECPs
  std::vector<int> order;
  auto add = [&](const int atom) {
    if (std::find(order.begin(), order.end(), atom) == order.end()) order.push_back(atom);
  };
  for (const Shell& sh : shells) add(sh.atom);
  for (const EcpPrimitive& p : ecp_prims) add(p.atom);
  return order;
}

std::string System::atom_label(const int id) const {
  const int a = integrator_atoms().at(id);
  return atoms[a].symbol + std::to_string(a + 1);
}

std::string System::function_label(int k) const {
  for (size_t s = 0; s < shells.size(); s++) {
    const int l = shells[s].l, n = (l + 1) * (l + 2) / 2;
    if (k >= n) {
      k -= n;
      continue;
    }
    // Cartesian components in library order: x power descending, then y power descending
    int c = 0;
    for (int i = l; i >= 0; i--)
      for (int j = l - i; j >= 0; j--, c++)
        if (c == k) {
          std::string comp =
              std::string(i, 'x') + std::string(j, 'y') + std::string(l - i - j, 'z');
          const Atom& a = atoms[shells[s].atom];
          return a.symbol + std::to_string(shells[s].atom + 1) + " shell " + std::to_string(s) +
                 " " + ANGULAR[l] + (l > 0 ? "(" + comp + ")" : "");
        }
  }
  return "?";
}

}  // namespace testsupport
