/// isochoric-nei.cpp
///
/// Non-equilibrium isochoric cooling: a regression test.
///
/// The gas has the same composition and density as in the isochoric-cooling
/// test (H, He, C, N, O, Ne, Si, S, Fe; Asplund et al. 2009 mass fractions).
/// It starts in collisional ionisation equilibrium (CIE) at 1e9 K, and cools
/// radiatively at constant density while the ion fractions evolve with
/// collisional ionisation, radiative and dielectronic recombination, and
/// charge exchange.  Recombination cannot keep up with the cooling, so the
/// gas becomes over-ionised as it cools to ~1e4 K.
///
/// This is a regression test: the reference data are NEMO's own solution,
/// written by this program with
///     isochoric-nei --write-reference <file>
/// at a 100 times smaller relative tolerance than the test uses.  Write new
/// reference data (tests/data/isochoric_nei_reference.txt) when a change to
/// NEMO is meant to change the results, e.g. new atomic data.
///
/// The reference data give, at times t_i, the temperature and the fraction
/// of each element in each ion, and their time derivatives.  For each
/// interval, the test starts from the reference state at t_{i-1}, evolves it
/// to t_i, and compares the temperature and the ion fractions with the
/// reference.  A value passes if it agrees to a relative tolerance, or if
/// its error corresponds to a time error |y - y_ref| / |dy/dt| of at most
/// time_tolerance times the interval.  As in the isochoric-cooling test,
/// this keeps the test equally strict where the gas cools slowly and where
/// it cools and recombines very quickly.

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "nemo/nemo.h"
#include "nemo_git_version.h"
using namespace NEMO;

namespace {

const std::vector<std::string> elements = {"H",  "He", "C", "N", "O",
                                           "Ne", "Si", "S", "Fe"};
// Elemental mass fractions adopted from Asplund et al. (2009), with helium
// adjusted so the included elements sum to one.
const std::vector<double> element_mass_fractions = {
    0.7381,    0.2495736, 0.002368,  0.000693, 0.005739,
    0.0012579, 0.0006656, 0.0003096, 0.0012933};
const double density   = 3.847988437033092e-22;  // g/cm^3
const double T_start   = 1.0e9;                  // K
const double T_end     = 1.5e4;                  // K, last checkpoint
const int n_checkpoint = 40;
const double eos_gamma = 5.0 / 3.0;

const double test_reltol      = 1.0e-8;   // CVODE tolerance in the test
const double reference_reltol = 1.0e-10;  // CVODE tolerance for reference

// Test tolerances: relative tolerance on the temperature, relative and
// absolute tolerances on ion fractions (the absolute one is for trace ions,
// whose errors are set by CVODE's absolute tolerance), or time error as a
// fraction of the interval.
const double temperature_tolerance = 1.0e-6;
const double ion_tolerance         = 1.0e-6;
const double ion_abs_tolerance     = 1.0e-10;
const double time_tolerance        = 1.0e-5;
// Ion fractions below this are not compared.
const double min_ion_fraction = 1.0e-6;

/// State vector layout and tracer names.
struct problem {
  std::vector<std::string> tracers;
  std::vector<double> state;  ///< density, pressure, X_el..., ions...
  std::vector<int> ion_element;
  size_t first_ion = 0;  ///< index of the first ion in the state vector
  size_t n_ions    = 0;
};

/// Fully ionised gas at T_start (the pressure is set later).
problem setup_problem()
{
  atomic_physics_database atomic_data;
  problem p;
  p.state = {density, 0.0};
  for (size_t e = 0; e < elements.size(); ++e) {
    p.tracers.push_back("X_" + elements[e]);
    p.state.push_back(element_mass_fractions[e]);
  }
  p.first_ion = p.state.size();
  for (size_t e = 0; e < elements.size(); ++e) {
    const int Z = atomic_data.get_atomic_number(elements[e]);
    for (int q = 0; q < Z; ++q) {
      p.tracers.push_back(
          q == 0 ? elements[e] : elements[e] + std::to_string(q) + "+");
      p.state.push_back(0.0);  // only the top (bare) ion is present
      p.ion_element.push_back(e);
    }
  }
  p.n_ions = p.ion_element.size();
  // Pressure for a fully ionised gas at T_start; set_temperature() corrects
  // it for the actual ion fractions.
  double n = 0.0;
  for (size_t e = 0; e < elements.size(); ++e)
    n += density * element_mass_fractions[e]
         / atomic_data.get_atomic_mass(elements[e])
         * (1.0 + atomic_data.get_atomic_number(elements[e]));
  p.state[1] = n * atomic_data.kB() * T_start;
  return p;
}

nemo_options make_flags(const problem &p, double reltol, bool cooling)
{
  nemo_options flags{};
  flags.Min_Temperature         = 1.0e3;
  flags.Max_Temperature         = 1.0e9;
  flags.nemo_reltol             = reltol;
  flags.nemo_abstol             = 1.0e-12;
  flags.nemo_minerg             = 1.0e-20;
  flags.nv_prim                 = p.state.size();
  flags.n_scalar                = p.tracers.size();
  flags.prim_dens_idx           = 0;
  flags.prim_pres_idx           = 1;
  flags.timestep_limiter        = 0;
  flags.update_energy           = cooling;
  flags.cooling                 = cooling;
  flags.collisional_ionization  = 1;
  flags.photoionization         = 0;
  flags.charge_exchange         = 1;
  flags.cosmicray_ionization    = 0;
  flags.radiative_recombination = 1;
  flags.onthespot_radiation     = 0;
  return flags;
}

/// Set the pressure so that the temperature is T (the pressure is
/// proportional to T at fixed composition).
void set_temperature(nemo &MP, std::vector<double> &state, double T)
{
  state[1] *= T / MP.Temperature(state, eos_gamma);
}

/// Fraction of its element in each ion (not including the bare ion).
std::vector<double> ion_fractions(
    const problem &p, const std::vector<double> &state)
{
  std::vector<double> y(p.n_ions);
  for (size_t i = 0; i < p.n_ions; ++i)
    y[i] = state[p.first_ion + i] / element_mass_fractions[p.ion_element[i]];
  return y;
}

void set_ion_fractions(
    const problem &p, const std::vector<double> &y, std::vector<double> &state)
{
  for (size_t i = 0; i < p.n_ions; ++i)
    state[p.first_ion + i] = y[i] * element_mass_fractions[p.ion_element[i]];
}

/// Advance the state by dt; exits on failure.
void advance(nemo &MP, std::vector<double> &state, double dt)
{
  std::vector<nemo_cell_rt_data> cell_rt;
  std::vector<double> out(state.size());
  double dt_min = 0.0;
  const int err =
      MP.update_single_cell(state, cell_rt, out, dt, eos_gamma, &dt_min);
  if (err) {
    spdlog::error("update_single_cell failed with code {}", err);
    exit_nemo(1);
  }
  state.swap(out);
}

/// CIE state at T_start: evolve the ion fractions at fixed temperature
/// until they no longer change.
std::vector<double> cie_state(const problem &p)
{
  std::vector<struct nemo_star_properties> rs;
  nemo MP(make_flags(p, reference_reltol, false), p.tracers, rs, eos_gamma);
  std::vector<double> state = p.state;
  for (int k = 0; k < 10; ++k) {
    std::vector<double> old = state;
    advance(MP, state, 1.0e14);
    set_temperature(MP, state, T_start);
    double change = 0.0;
    for (size_t i = p.first_ion; i < state.size(); ++i)
      change = std::max(change, std::abs(state[i] - old[i]));
    if (k > 0 && change < 1.0e-14) return state;
  }
  spdlog::error("CIE state at {:.3e} K did not converge", T_start);
  exit_nemo(1);
  return state;
}

/// Reference data: one row per checkpoint.
struct reference_row {
  double t, T, dTdt;
  std::vector<double> y, dydt;
};

/// Write the reference data, with NEMO at reference_reltol.
int write_reference(const char *filename)
{
  problem p = setup_problem();
  std::vector<struct nemo_star_properties> rs;
  nemo MP(make_flags(p, reference_reltol, true), p.tracers, rs, eos_gamma);
  const std::vector<double> start = cie_state(p);

  // 1. Dense run (temperature changes of ~1% per step) to find the times
  //    at which the gas reaches the checkpoint temperatures.
  std::vector<double> dense_t = {0.0}, dense_logT = {std::log(T_start)};
  {
    std::vector<double> state = start;
    double t = 0.0, dt = 1.0e9;
    while (std::exp(dense_logT.back()) > 0.9 * T_end) {
      advance(MP, state, dt);
      t += dt;
      const double T      = MP.Temperature(state, eos_gamma);
      const double change = std::abs(std::log(T) - dense_logT.back()) + 1.0e-6;
      dense_t.push_back(t);
      dense_logT.push_back(std::log(T));
      dt *= std::min(2.0, std::max(0.5, 0.01 / change));
      if (dense_t.size() > 100000) {
        spdlog::error("dense run did not reach {:.3e} K", T_end);
        return 1;
      }
    }
  }
  std::vector<double> t_check = {0.0};
  for (int c = 1; c < n_checkpoint; ++c) {
    const double logT =
        std::log(T_start)
        + c * (std::log(T_end) - std::log(T_start)) / (n_checkpoint - 1);
    size_t k = 1;
    while (dense_logT[k] > logT)
      ++k;
    const double f =
        (logT - dense_logT[k - 1]) / (dense_logT[k] - dense_logT[k - 1]);
    t_check.push_back(dense_t[k - 1] + f * (dense_t[k] - dense_t[k - 1]));
  }

  // 2. Run from checkpoint to checkpoint, and get the time derivatives
  //    from a short extra step (a small fraction of the interval).
  std::vector<reference_row> rows;
  std::vector<double> state = start;
  for (int c = 0; c < n_checkpoint; ++c) {
    if (c > 0) advance(MP, state, t_check[c] - t_check[c - 1]);
    reference_row r;
    r.t = t_check[c];
    r.T = MP.Temperature(state, eos_gamma);
    r.y = ion_fractions(p, state);
    const double delta =
        1.0e-4 * (c > 0 ? t_check[c] - t_check[c - 1] : t_check[1]);
    std::vector<double> next = state;
    advance(MP, next, delta);
    r.dTdt = (MP.Temperature(next, eos_gamma) - r.T) / delta;
    const std::vector<double> y_next = ion_fractions(p, next);
    for (size_t i = 0; i < p.n_ions; ++i)
      r.dydt.push_back((y_next[i] - r.y[i]) / delta);
    rows.push_back(r);
  }

  FILE *f = std::fopen(filename, "w");
  if (!f) {
    spdlog::error("Can't open {} for writing", filename);
    return 1;
  }
  std::fprintf(
      f,
      "# Non-equilibrium isochoric cooling from CIE at %.3e K: reference\n"
      "# data for NEMO's isochoric-nei regression test, written by\n"
      "# 'isochoric-nei --write-reference' (NEMO git %s,\n"
      "# CVODE relative tolerance %.1e).\n#\n"
      "# Columns: t [s], T [K], dT/dt [K/s], then the fraction of its\n"
      "# element in each ion below (in this order), then the time\n"
      "# derivatives of these fractions [1/s].\n# ions:",
      T_start, git_version.c_str(), reference_reltol);
  for (size_t i = 0; i < p.n_ions; ++i)
    std::fprintf(f, " %s", p.tracers[elements.size() + i].c_str());
  std::fprintf(f, "\n");
  for (const auto &r : rows) {
    std::fprintf(f, "%.15e %.15e %.15e", r.t, r.T, r.dTdt);
    for (double v : r.y)
      std::fprintf(f, " %.15e", v);
    for (double v : r.dydt)
      std::fprintf(f, " %.15e", v);
    std::fprintf(f, "\n");
  }
  std::fclose(f);
  spdlog::info("Wrote {} checkpoints to {}", rows.size(), filename);
  return 0;
}

std::vector<reference_row> read_reference(size_t n_ions)
{
  std::ifstream file(ISOCHORIC_NEI_REFERENCE_FILE);
  if (!file.is_open()) {
    spdlog::error("Failed to open {}", ISOCHORIC_NEI_REFERENCE_FILE);
    exit_nemo(1);
  }
  std::vector<reference_row> rows;
  std::string line;
  while (std::getline(file, line)) {
    std::istringstream row(line);
    row >> std::ws;
    if (row.eof() || row.peek() == '#') continue;
    reference_row r;
    r.y.resize(n_ions);
    r.dydt.resize(n_ions);
    row >> r.t >> r.T >> r.dTdt;
    for (auto &v : r.y)
      row >> v;
    for (auto &v : r.dydt)
      row >> v;
    if (!row) {
      spdlog::error("Invalid reference row: {}", line);
      exit_nemo(1);
    }
    rows.push_back(r);
  }
  return rows;
}

/// True if value agrees with ref to rel_tol (plus abs_tol), or if the error
/// corresponds to a time error of at most time_tolerance * interval.
/// time_error is set to the time error, as a fraction of the interval, if
/// the value is outside rel_tol and abs_tol, and to zero otherwise.
bool agrees(
    double value,
    double ref,
    double dref_dt,
    double interval,
    double rel_tol,
    double abs_tol,
    double &time_error)
{
  const double err = std::abs(value - ref);
  if (std::isfinite(value) && err <= rel_tol * std::abs(ref) + abs_tol) {
    time_error = 0.0;
    return true;
  }
  time_error = (dref_dt != 0.0) ? err / std::abs(dref_dt) / interval : HUGE_VAL;
  return std::isfinite(value) && time_error <= time_tolerance;
}

int run_test()
{
  problem p = setup_problem();
  std::vector<struct nemo_star_properties> rs;
  nemo MP(make_flags(p, test_reltol, true), p.tracers, rs, eos_gamma);
  const std::vector<reference_row> ref = read_reference(p.n_ions);
  if (ref.size() < 2) {
    spdlog::error("Not enough reference data");
    return 1;
  }

  bool test_status = true;
  for (size_t c = 1; c < ref.size(); ++c) {
    // Start from the reference state at the start of the interval.
    std::vector<double> state = p.state;
    set_ion_fractions(p, ref[c - 1].y, state);
    set_temperature(MP, state, ref[c - 1].T);
    const double interval = ref[c].t - ref[c - 1].t;
    advance(MP, state, interval);

    // Temperature.
    double time_error = 0.0;
    const double T    = MP.Temperature(state, eos_gamma);
    bool passed       = agrees(
        T, ref[c].T, ref[c].dTdt, interval, temperature_tolerance, 0.0,
        time_error);
    double worst_time_error = time_error;
    std::string worst       = "T";

    // Ion fractions above min_ion_fraction.
    const std::vector<double> y = ion_fractions(p, state);
    for (size_t i = 0; i < p.n_ions; ++i) {
      if (ref[c].y[i] < min_ion_fraction) continue;
      const bool ion_passed = agrees(
          y[i], ref[c].y[i], ref[c].dydt[i], interval, ion_tolerance,
          ion_abs_tolerance, time_error);
      if (!ion_passed) {
        spdlog::error(
            "Interval {}: {} fraction {:.6e}, reference {:.6e}", c,
            p.tracers[elements.size() + i], y[i], ref[c].y[i]);
      }
      if (time_error > worst_time_error) {
        worst_time_error = time_error;
        worst            = p.tracers[elements.size() + i];
      }
      passed = passed && ion_passed;
    }

    test_status = test_status && passed;
    spdlog::info(
        "Interval {:2d}: T_ref = {:.4e} K, T error = {:.2e}, largest time "
        "error = {:.2e} of interval ({}): {}",
        c, ref[c].T, std::abs(T - ref[c].T) / ref[c].T, worst_time_error, worst,
        passed ? "PASSED" : "FAILED");
  }
  spdlog::info(
      "NEMO isochoric NEI cooling test: {}", test_status ? "PASSED" : "FAILED");
  return test_status ? 0 : 1;
}

}  // namespace

int main(int argc, char **argv)
{
  spdlog::set_level(spdlog::level::info);
  spdlog::flush_on(spdlog::level::err);
  spdlog::info("NEMO isochoric NEI cooling test: git {}", git_version);
  if (argc == 3 && std::strcmp(argv[1], "--write-reference") == 0)
    return write_reference(argv[2]);
  if (argc != 1) {
    spdlog::error("usage: {} [--write-reference <file>]", argv[0]);
    return 1;
  }
  return run_test();
}
