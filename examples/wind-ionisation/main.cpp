// wind-ionisation: ion fractions of a stellar wind in photoionisation
// equilibrium at its boundary, for the WIND_<k>_TR<n> tracer values of a
// PION parameter file.
//
// PION sets the gas inside a wind-source boundary (radius WIND_<k>_radius)
// to fixed values, including the ion mass fractions WIND_<k>_TR<n>.  If the
// star's own radiation would keep the gas in a different ionisation state,
// the cells just outside the boundary change their ionisation state
// abruptly.  This program reads the parameter file, sets up one cell just
// outside the boundary of each wind source (the density and temperature
// that PION sets there, and the radiation of the star at the same position,
// unattenuated, as PION treats the gas inside the boundary as transparent),
// evolves the ion fractions with NEMO until they no longer change, and
// prints new WIND_<k>_TR<n> lines for the parameter file.
//
// The cell size matters for very optically thick ions (NEMO uses the
// photon-conserving rate (1 - exp(-dtau))), so by default it is the cell
// size of the finest grid level that contains the star.
//
// usage: wind-ionisation <parameter file> [options]
//   --wind <k>     only wind source k (default: all)
//   --ds <cm>      cell size (default: finest level that contains the star)
//   --T <K>        temperature (default: as PION sets it at the boundary)
//   --thermal      evolve the temperature too (photoionisation and thermal
//                  equilibrium) instead of keeping it fixed
//   --all-sources  also include the other radiation sources, unattenuated,
//                  at their distance from the wind boundary

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "nemo/nemo.h"  // everything in NEMO is in namespace NEMO

namespace {

const double msun      = 1.9891e33;    // g
const double year      = 3.1558150e7;  // s
const double gamma_eos = 5.0 / 3.0;

/// The parameter file as a map from key to value (first word after key).
std::map<std::string, std::string> read_parameters(const std::string &file)
{
  std::map<std::string, std::string> par;
  std::ifstream in(file);
  if (!in) {
    std::fprintf(stderr, "Cannot open parameter file %s\n", file.c_str());
    std::exit(1);
  }
  std::string line;
  while (std::getline(in, line)) {
    const size_t hash = line.find('#');
    if (hash != std::string::npos) line.erase(hash);
    std::istringstream ls(line);
    std::string key, value;
    if (ls >> key >> value) par[key] = value;
  }
  return par;
}

double get(
    const std::map<std::string, std::string> &par,
    const std::string &key,
    const double def)
{
  auto it = par.find(key);
  return (it == par.end()) ? def : std::atof(it->second.c_str());
}

double get_required(
    const std::map<std::string, std::string> &par, const std::string &key)
{
  auto it = par.find(key);
  if (it == par.end()) {
    std::fprintf(stderr, "Parameter %s not found\n", key.c_str());
    std::exit(1);
  }
  return std::atof(it->second.c_str());
}

std::string key_n(const std::string &a, int n, const std::string &b = "")
{
  return a + std::to_string(n) + b;
}

/// Element and charge of an ion tracer name ("C", "C2+"), or false.
bool parse_ion(const std::string &name, std::string &element, int &charge)
{
  size_t i = 0;
  while (i < name.size() && std::isalpha(static_cast<unsigned char>(name[i])))
    i++;
  if (i == 0 || i > 2 || !std::isupper(static_cast<unsigned char>(name[0])))
    return false;
  element = name.substr(0, i);
  if (i == name.size()) {
    charge = 0;
    return true;
  }
  if (name.back() != '+') return false;
  const std::string q = name.substr(i, name.size() - i - 1);
  if (q.empty() || !std::all_of(q.begin(), q.end(), ::isdigit)) return false;
  charge = std::atoi(q.c_str());
  return true;
}

/// Cell size of the finest nested-grid level that contains the position x.
double finest_cell_size(
    const std::map<std::string, std::string> &par, const double *x)
{
  const int ndim     = static_cast<int>(get_required(par, "ndim"));
  const int nlev     = static_cast<int>(get(par, "grid_nlevels", 1));
  const char *ax[3]  = {"X", "Y", "Z"};
  const char *axx[3] = {"XX", "YY", "ZZ"};
  double lo[3], hi[3];
  for (int i = 0; i < ndim; i++) {
    lo[i] = get_required(par, std::string(ax[i]) + "min");
    hi[i] = get_required(par, std::string(ax[i]) + "max");
  }
  const double dx0 = (hi[0] - lo[0]) / get_required(par, "NGridX");
  int finest       = 0;
  for (int l = 1; l < nlev; l++) {
    bool inside = true;
    for (int i = 0; i < ndim; i++) {
      // each refined level is half the size of the level above, centred on
      // NG_centre but kept inside the level above.
      if (get(par, std::string("NG_refine_") + axx[i], 1) == 0) continue;
      const double w = 0.5 * (hi[i] - lo[i]);
      double l0 = get(par, std::string("NG_centre_") + axx[i], 0.0) - 0.5 * w;
      l0        = std::max(l0, lo[i]);
      l0        = std::min(l0, hi[i] - w);
      lo[i]     = l0;
      hi[i]     = l0 + w;
      if (x[i] < lo[i] || x[i] > hi[i]) inside = false;
    }
    if (!inside) break;
    finest = l;
  }
  return dx0 / std::pow(2.0, finest);
}

}  // namespace

int main(int argc, char **argv)
{
  if (argc < 2) {
    std::fprintf(
        stderr,
        "usage: %s <parameter file> [--wind k] [--ds cm] [--T K] "
        "[--thermal] [--all-sources]\n",
        argv[0]);
    return 1;
  }
  int only_wind = -1;
  double ds_opt = 0.0, T_opt = 0.0;
  bool thermal = false, all_sources = false;
  for (int i = 2; i < argc; i++) {
    const std::string a(argv[i]);
    if (a == "--wind" && i + 1 < argc)
      only_wind = std::atoi(argv[++i]);
    else if (a == "--ds" && i + 1 < argc)
      ds_opt = std::atof(argv[++i]);
    else if (a == "--T" && i + 1 < argc)
      T_opt = std::atof(argv[++i]);
    else if (a == "--thermal")
      thermal = true;
    else if (a == "--all-sources")
      all_sources = true;
    else {
      std::fprintf(stderr, "Unknown option %s\n", a.c_str());
      return 1;
    }
  }
  // Only show errors: NEMO warns about intermediate states of the ODE
  // solver and about the temperature floor, which are harmless here.
  spdlog::set_level(spdlog::level::err);
  const auto par = read_parameters(argv[1]);
  const int ndim = static_cast<int>(get_required(par, "ndim"));

  // Tracers.  NEMO uses the element (X_El) and ion (El, El1+, ...) tracers;
  // other tracers are passive.
  const int ntr = static_cast<int>(get_required(par, "ntracer"));
  std::vector<std::string> tracers(ntr);
  for (int t = 0; t < ntr; t++) {
    char key[32];
    std::snprintf(key, sizeof(key), "Tracer%03d", t);
    auto it = par.find(key);
    if (it == par.end()) {
      std::fprintf(stderr, "Parameter %s not found\n", key);
      return 1;
    }
    tracers[t] = it->second;
  }

  // NEMO options, as PION sets them from the EP_ parameters.
  NEMO::nemo_options o{};
  o.nemo_reltol             = get(par, "EP_nemo_reltol", 1.0e-6);
  o.nemo_abstol             = get(par, "EP_nemo_abstol", 1.0e-12);
  o.nemo_minerg             = get(par, "EP_nemo_minerg", 1.0e-20);
  o.Min_Temperature         = get(par, "EP_Min_Temperature", 10.0);
  o.Max_Temperature         = get(par, "EP_Max_Temperature", 1.0e9);
  o.nv_prim                 = 2 + ntr;  // density, pressure, tracers
  o.n_scalar                = ntr;
  o.prim_dens_idx           = 0;
  o.prim_pres_idx           = 1;
  o.timestep_limiter        = 0;
  o.cooling                 = get(par, "EP_cooling", 1) != 0;
  o.update_energy           = thermal;
  o.collisional_ionization  = get(par, "EP_coll_ionisation", 1) != 0;
  o.radiative_recombination = get(par, "EP_rad_recombination", 1) != 0;
  o.charge_exchange =
      get(par, "EP_charge_exchange", 0) != 0 && o.collisional_ionization;
  o.cosmicray_ionization =
      get(par, "EP_cr_ionisation", 0) != 0 && o.collisional_ionization;
  o.onthespot_radiation = get(par, "EP_onthespot_radiation", 0) != 0;
  o.photoionization     = true;

  // Radiation sources that NEMO can use (stellar atmosphere models).
  std::vector<NEMO::nemo_star_properties> stars;
  std::vector<std::array<double, 3> > star_pos;
  const int nrt = static_cast<int>(get(par, "RT_Nsources", 0));
  for (int s = 0; s < nrt; s++) {
    const int effect =
        static_cast<int>(get_required(par, key_n("RT_effect___", s)));
    if (effect < RT_EFFECT_MP10_ATLAS_O3V || effect > RT_EFFECT_MP10_BLACKBODY)
      continue;
    NEMO::nemo_star_properties p;
    p.Luminosity = get_required(par, key_n("RT_strength_", s));
    p.Teff       = get_required(par, key_n("RT_Tstar____", s));
    p.Radius     = get_required(par, key_n("RT_Rstar____", s));
    p.SED_type   = effect;
    p.id         = static_cast<int>(stars.size());
    p.photo_emax = 0.0;
    p.Nbins      = 0;
    p.SED_info   = "";
    stars.push_back(p);
    std::array<double, 3> x = {0.0, 0.0, 0.0};
    for (int i = 0; i < ndim; i++)
      x[i] = get(par, key_n("RT_position_", s, "_" + std::to_string(i)), 0.0);
    star_pos.push_back(x);
  }
  if (stars.empty()) {
    std::fprintf(stderr, "No radiation sources that NEMO can use\n");
    return 1;
  }

  NEMO::nemo MP(o, tracers, stars, gamma_eos);
  MP.get_radiation_src_info(stars);
  NEMO::atomic_physics_database atomic;

  const int nwind = static_cast<int>(get(par, "WIND_NSRC", 0));
  for (int k = 0; k < nwind; k++) {
    if (only_wind >= 0 && k != only_wind) continue;
    const std::string w      = "WIND_" + std::to_string(k) + "_";
    std::array<double, 3> xw = {0.0, 0.0, 0.0};
    for (int i = 0; i < ndim; i++)
      xw[i] = get(par, w + "pos" + std::to_string(i), 0.0);
    const double R_wind = get_required(par, w + "radius");
    const double R_star = get_required(par, w + "Rstr");
    const double Mdot   = get_required(par, w + "mdot") * msun / year;
    const double vinf   = get_required(par, w + "vinf") * 1.0e5;
    const double Twind  = get_required(par, w + "temp");
    const bool accel    = get(par, w + "acceleration", 0) != 0;

    // the radiation source at the position of the star
    int own     = -1;
    double dmin = 1.0e99;
    for (size_t s = 0; s < stars.size(); s++) {
      double d2 = 0.0;
      for (int i = 0; i < 3; i++)
        d2 += (star_pos[s][i] - xw[i]) * (star_pos[s][i] - xw[i]);
      if (std::sqrt(d2) < dmin) {
        dmin = std::sqrt(d2);
        own  = static_cast<int>(s);
      }
    }
    if (own < 0 || dmin > 0.01 * R_wind) {
      std::printf(
          "# wind %d: no radiation source at the position of the star; "
          "skipped\n",
          k);
      continue;
    }

    // The cell just outside the boundary: radius r, size ds.
    const double ds =
        (ds_opt > 0.0) ? ds_opt : finest_cell_size(par, xw.data());
    const double r = R_wind + 0.5 * ds;

    // Density and pressure as PION sets them in the wind (3D geometry):
    // rho = Mdot / (4 pi r^2 v), and adiabatic expansion from the stellar
    // radius, where T = WIND_k_temp.
    double v = vinf;
    if (accel) {
      const double v0 = 4.0 * std::sqrt(atomic.kB() * Twind / atomic.m_p());
      v               = v0 + (vinf - v0) * (1.0 - R_star / r);
    }
    std::vector<double> state(2 + ntr, 0.0);
    state[0] = Mdot / (4.0 * M_PI * r * r * v);
    state[1] =
        atomic.kB() * Twind / atomic.m_p()
        * std::pow(4.0 * M_PI * R_star * R_star * v / Mdot, gamma_eos - 1.0)
        * std::pow(state[0], gamma_eos);
    for (int t = 0; t < ntr; t++)
      state[2 + t] = get(par, w + "TR" + std::to_string(t), 0.0);
    const std::vector<double> input = state;  // values in the parameter file

    // Element tracers (X_El), by element name.
    std::map<std::string, int> element_tracer;
    for (int t = 0; t < ntr; t++)
      if (tracers[t].rfind("X_", 0) == 0)
        element_tracer[tracers[t].substr(2)] = t;
    if (element_tracer.empty()) {
      std::fprintf(stderr, "wind %d: no element tracers (X_El)\n", k);
      return 1;
    }

    // The element mass fractions must add up to 1: if not, change He (or,
    // without He, the most abundant element) by the difference, and scale
    // the ions of that element with it, so their ion fractions stay the same.
    double sumX = 0.0;
    for (const auto &el : element_tracer)
      sumX += state[2 + el.second];
    std::string fix_el;
    if (std::fabs(sumX - 1.0) > 1.0e-12) {
      if (element_tracer.count("He"))
        fix_el = "He";
      else {
        double xmax = -1.0;
        for (const auto &el : element_tracer)
          if (state[2 + el.second] > xmax) {
            xmax   = state[2 + el.second];
            fix_el = el.first;
          }
      }
      const int te      = element_tracer[fix_el];
      const double Xold = state[2 + te];
      const double Xnew = Xold + (1.0 - sumX);
      if (Xnew <= 0.0) {
        std::fprintf(
            stderr,
            "wind %d: element mass fractions add up to %.12f, which cannot "
            "be corrected by changing X_%s = %.6e\n",
            k, sumX, fix_el.c_str(), Xold);
        return 1;
      }
      state[2 + te] = Xnew;
      for (int t = 0; t < ntr; t++) {
        std::string e;
        int q;
        if (parse_ion(tracers[t], e, q) && e == fix_el && Xold > 0.0)
          state[2 + t] *= Xnew / Xold;
      }
    }
    const std::vector<double> composition = state;  // corrected X_El

    // Temperature: PION sets at least the minimum temperature in the wind.
    double T = MP.Temperature(state, gamma_eos);
    if (T_opt > 0.0)
      T = T_opt;
    else
      T = std::max(T, o.Min_Temperature);
    MP.Set_Temp(state, T, gamma_eos);

    // Radiation field at the cell: the star, unattenuated (PION treats the
    // gas inside the wind boundary as transparent); optionally the other
    // sources, unattenuated, at their distance from the boundary.
    std::vector<NEMO::nemo_cell_rt_data> rt(stars.size());
    for (size_t s = 0; s < stars.size(); s++) {
      rt[s].star_id = static_cast<int>(s);
      rt[s].tau.assign(stars[s].Nbins, 0.0);
      rt[s].dtau.assign(stars[s].Nbins, 0.0);
      double d = r;
      if (static_cast<int>(s) != own) {
        double d2 = 0.0;
        for (int i = 0; i < 3; i++)
          d2 += (star_pos[s][i] - xw[i]) * (star_pos[s][i] - xw[i]);
        d = std::max(std::sqrt(d2) - R_wind, 0.5 * ds);
      }
      if (static_cast<int>(s) == own || all_sources) {
        rt[s].dS     = ds;
        rt[s].Vshell = 4.0 * M_PI * d * d * ds;
      }
      else {
        rt[s].dS     = 0.0;  // no flux from this source
        rt[s].Vshell = 1.0e200;
      }
    }

    // Evolve to equilibrium: steps of increasing length until the ion
    // fractions no longer change.
    std::vector<double> out(state.size());
    double dt = 1.0, change = 1.0, time = 0.0, dummy = 0.0;
    int nstep = 0;
    for (; nstep < 400 && (change > 1.0e-9 || nstep < 20); nstep++) {
      if (MP.update_single_cell(state, rt, out, dt, gamma_eos, &dummy)) {
        std::fprintf(stderr, "wind %d: NEMO update failed\n", k);
        return 1;
      }
      change = 0.0;
      for (size_t i = 2; i < state.size(); i++)
        change = std::max(change, std::fabs(out[i] - state[i]));
      if (thermal)
        change = std::max(
            change, std::fabs(
                        MP.Temperature(out, gamma_eos)
                            / MP.Temperature(state, gamma_eos)
                        - 1.0));
      state = out;
      if (!thermal) MP.Set_Temp(state, T, gamma_eos);
      time += dt;
      dt *= 1.5;
    }

    std::printf(
        "# Wind source %d at x = (%.4e, %.4e, %.4e) cm: boundary radius "
        "%.4e cm, cell size %.4e cm\n",
        k, xw[0], xw[1], xw[2], R_wind, ds);
    std::printf(
        "#   rho = %.4e g/cm^3, n_e = %.4e cm^-3, T = %.4e K (%s), "
        "radiation of source %d%s\n",
        state[0], MP.get_n_elec(state), MP.Temperature(state, gamma_eos),
        thermal ? "thermal equilibrium" : "fixed", own,
        all_sources ? " and the other sources" : "");
    std::printf(
        "#   %s after %d steps (%.3e s), largest change in the last step "
        "%.1e\n",
        change <= 1.0e-9 ? "converged" : "NOT CONVERGED", nstep, time, change);
    if (!fix_el.empty())
      std::printf(
          "#   element mass fractions added up to %.12f: X_%s changed from "
          "%.12e to %.12e\n",
          sumX, fix_el.c_str(), input[2 + element_tracer[fix_el]],
          state[2 + element_tracer[fix_el]]);
    else
      std::printf("#   element mass fractions add up to 1\n");

    // Ion fractions (of each element), input and equilibrium, by element.
    std::printf("#   %-6s %12s %12s\n", "ion", "input", "equilibrium");
    for (const auto &el : element_tracer) {
      const double X = input[2 + el.second];
      if (X <= 0.0) continue;
      const int Z = atomic.get_atomic_number(el.first);
      std::vector<double> f_in(Z + 1, 0.0), f_eq(Z + 1, 0.0);
      for (int t = 0; t < ntr; t++) {
        std::string e;
        int q;
        if (!parse_ion(tracers[t], e, q) || e != el.first) continue;
        f_in[q] = input[2 + t] / X;
        f_eq[q] = state[2 + t] / X;
      }
      double s_in = 0.0, s_eq = 0.0;
      for (int q = 0; q < Z; q++) {
        s_in += f_in[q];
        s_eq += f_eq[q];
      }
      f_in[Z] = std::max(0.0, 1.0 - s_in);  // bare ion (not a tracer)
      f_eq[Z] = std::max(0.0, 1.0 - s_eq);
      for (int q = 0; q <= Z; q++) {
        if (f_in[q] < 1.0e-4 && f_eq[q] < 1.0e-4) continue;
        const std::string name =
            (q == 0) ? el.first : el.first + std::to_string(q) + "+";
        std::printf("#   %-6s %12.4e %12.4e\n", name.c_str(), f_in[q], f_eq[q]);
      }
    }

    // New parameter-file lines for all tracers of this wind source: the
    // element mass fractions (corrected to add up to 1), the ions in
    // equilibrium, and the other tracers unchanged.  Ion mass fractions
    // below 1e-12 are set to 1e-12, and the most abundant ion is reduced so
    // that the ions of each element add up to at most (1 - 1e-11) times the
    // element.  The margin is larger than the rounding error of the printed
    // values (13 significant digits), so that the printed ions never add up
    // to more than the element.
    std::printf(
        "# New tracer values for the parameter file (all tracers of wind "
        "source %d):\n",
        k);
    std::vector<double> newtr(state.begin() + 2, state.end());
    for (const auto &el : element_tracer) {
      const double X   = composition[2 + el.second];
      newtr[el.second] = X;
      double sum       = 0.0;
      int imax         = -1;
      for (int t = 0; t < ntr; t++) {
        std::string e;
        int q;
        if (!parse_ion(tracers[t], e, q) || e != el.first) continue;
        newtr[t] = std::max(newtr[t], 1.0e-12);
        sum += newtr[t];
        if (imax < 0 || newtr[t] > newtr[imax]) imax = t;
      }
      const double Xmax = X * (1.0 - 1.0e-11);
      if (imax >= 0 && sum > Xmax) newtr[imax] -= (sum - Xmax);
    }
    for (int t = 0; t < ntr; t++) {
      std::string e;
      int q;
      const bool is_ion = parse_ion(tracers[t], e, q)
                          && element_tracer.find(e) != element_tracer.end();
      const bool is_element = tracers[t].rfind("X_", 0) == 0;
      // other (passive) tracers keep the value of the parameter file.
      const double value = (is_ion || is_element) ? newtr[t] : input[2 + t];
      std::printf(
          "%sTR%-3d %.12e   # %s (was %.6e)\n", w.c_str(), t, value,
          tracers[t].c_str(), input[2 + t]);
    }
    std::printf("\n");
  }
  return 0;
}
