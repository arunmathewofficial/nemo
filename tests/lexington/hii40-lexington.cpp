/// hii40-lexington.cpp
///
/// Author: Arun Mathew
/// Created: 12 October 2026
///
/// Test NEMO's photoionisation and thermal evolution using the HII40
/// Lexington benchmark (Test 2 in Mathew et al. 2025, Section 3.4).
///
/// A central 40,000 K blackbody irradiates initially neutral, uniform gas
/// containing H, He, C, N, O and Ne on a one-dimensional spherical grid.
/// Includes radial attenuation, heating and cooling, without hydrodynamics.
///

#ifdef SPDLOG_FWD
#include <spdlog/fwd.h>
#endif
#include <fmt/ranges.h>
#include <spdlog/sinks/rotating_file_sink.h>
#include <spdlog/spdlog.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <vector>

#ifdef NEMO_OMP
#include <omp.h>
#endif

#include "nemo-datastructs.h"
#include "nemo/nemo.h"
#include "nemo_git_version.h"

#include "../common_functions.h"


std::pair<std::vector<double>, std::vector<std::string> > setup_initial_state()
{
  // Initialize the initial state vector with all ions zero (neutral state,
  // enforced by the conservation equations for each element)
  std::vector<double> initial_vector(59, 0.0);
  initial_vector[0] = 2.3589999999999998e-22;  // density
  initial_vector[1] = 1.0e-12;                 // pressure (~1e5 K)
  // Elemental abundances (from Asplund 2010/Solar System)
  initial_vector[2] = 7.0933e-01;  // H (72%)
  initial_vector[3] = 2.8373e-01;  // He (26%)
  initial_vector[4] = 1.8725e-03;  // C (0.25%)
  initial_vector[5] = 3.9726e-04;  // N (0.15%)
  initial_vector[6] = 3.7452e-03;  // O (0.15%)
  initial_vector[7] = 7.1647e-04;  // Ne (0.05%)
  initial_vector[8] = 2.0492e-04;  // S (0.015%)

  // Ion concentrations (all set to zero for neutral plasma)
  // H: 1 ion (H0) - index 9
  initial_vector[9] = 7.0933e-01;  // H0
  // He: 2 ions (He0, He+) - indices 10-11
  initial_vector[10] = 2.8373e-01;  // He0
  initial_vector[11] = 0.0;         // He1+
  // C: 6 ions (C0-C6+) - indices 12-17
  initial_vector[12] = 1.0000e-03;  // C0
  initial_vector[13] = 0.8725e-3;   // C1+
  initial_vector[14] = 0.0;         // C2+
  initial_vector[15] = 0.0;         // C3+
  initial_vector[16] = 0.0;         // C4+
  initial_vector[17] = 0.0;         // C5+
  // N: 7 ions (N0-N7+) - indices 18-24
  initial_vector[18] = 3.9726e-04;  // N0
  initial_vector[19] = 0.0;         // N1+
  initial_vector[20] = 0.0;         // N2+
  initial_vector[21] = 0.0;         // N3+
  initial_vector[22] = 0.0;         // N4+
  initial_vector[23] = 0.0;         // N5+
  initial_vector[24] = 0.0;         // N6+
  // O: 8 ions (O0-O8+) - indices 25-32
  initial_vector[25] = 3.7452e-03;  // O0
  initial_vector[26] = 0.0;         // O1+
  initial_vector[27] = 0.0;         // O2+
  initial_vector[28] = 0.0;         // O3+
  initial_vector[29] = 0.0;         // O4+
  initial_vector[30] = 0.0;         // O5+
  initial_vector[31] = 0.0;         // O6+
  initial_vector[32] = 0.0;         // O7+
  // Ne: 10 ions (Ne0-Ne9+) - indices 33-42
  initial_vector[33] = 7.1647e-04;  // Ne0
  initial_vector[34] = 0.0;         // Ne1+
  initial_vector[35] = 0.0;         // Ne2+
  initial_vector[36] = 0.0;         // Ne3+
  initial_vector[37] = 0.0;         // Ne4+
  initial_vector[38] = 0.0;         // Ne5+
  initial_vector[39] = 0.0;         // Ne6+
  initial_vector[40] = 0.0;         // Ne7+
  initial_vector[41] = 0.0;         // Ne8+
  initial_vector[42] = 0.0;         // Ne9+
  // S: 16 ions (S0-S15+) - indices 43-58
  initial_vector[43] = 2.0492e-04;  // S0
  initial_vector[44] = 0.0;         // S1+
  initial_vector[45] = 0.0;         // S2+
  initial_vector[46] = 0.0;         // S3+
  initial_vector[47] = 0.0;         // S4+
  initial_vector[48] = 0.0;         // S5+
  initial_vector[49] = 0.0;         // S6+
  initial_vector[50] = 0.0;         // S7+
  initial_vector[51] = 0.0;         // S8+
  initial_vector[52] = 0.0;         // S9+
  initial_vector[53] = 0.0;         // S10+
  initial_vector[54] = 0.0;         // S11+
  initial_vector[55] = 0.0;         // S12+
  initial_vector[56] = 0.0;         // S13+
  initial_vector[57] = 0.0;         // S14+
  initial_vector[58] = 0.0;         // S15+

  // Create tracer names vector
  std::vector<std::string> tracer_names(57);

  // Element names (X_ prefix)
  tracer_names[0] = "X_H";
  tracer_names[1] = "X_He";
  tracer_names[2] = "X_C";
  tracer_names[3] = "X_N";
  tracer_names[4] = "X_O";
  tracer_names[5] = "X_Ne";
  tracer_names[6] = "X_S";

  // Ion names (element name + ion charge)
  // H: 1 ion (H0)
  tracer_names[7] = "H";

  // He: 2 ions (He0, He+)
  tracer_names[8] = "He";
  tracer_names[9] = "He1+";

  // C: 6 ions (C0-C6+)
  tracer_names[10] = "C";
  tracer_names[11] = "C1+";
  tracer_names[12] = "C2+";
  tracer_names[13] = "C3+";
  tracer_names[14] = "C4+";
  tracer_names[15] = "C5+";

  // N: 7 ions (N0-N7+)
  tracer_names[16] = "N";
  tracer_names[17] = "N1+";
  tracer_names[18] = "N2+";
  tracer_names[19] = "N3+";
  tracer_names[20] = "N4+";
  tracer_names[21] = "N5+";
  tracer_names[22] = "N6+";

  // O: 8 ions (O0-O8+)
  tracer_names[23] = "O";
  tracer_names[24] = "O1+";
  tracer_names[25] = "O2+";
  tracer_names[26] = "O3+";
  tracer_names[27] = "O4+";
  tracer_names[28] = "O5+";
  tracer_names[29] = "O6+";
  tracer_names[30] = "O7+";

  // Ne: 10 ions (Ne0-Ne9+)
  tracer_names[31] = "Ne";
  tracer_names[32] = "Ne1+";
  tracer_names[33] = "Ne2+";
  tracer_names[34] = "Ne3+";
  tracer_names[35] = "Ne4+";
  tracer_names[36] = "Ne5+";
  tracer_names[37] = "Ne6+";
  tracer_names[38] = "Ne7+";
  tracer_names[39] = "Ne8+";
  tracer_names[40] = "Ne9+";

  // S: 16 ions (S0-S15+)
  tracer_names[41] = "S";
  tracer_names[42] = "S1+";
  tracer_names[43] = "S2+";
  tracer_names[44] = "S3+";
  tracer_names[45] = "S4+";
  tracer_names[46] = "S5+";
  tracer_names[47] = "S6+";
  tracer_names[48] = "S7+";
  tracer_names[49] = "S8+";
  tracer_names[50] = "S9+";
  tracer_names[51] = "S10+";
  tracer_names[52] = "S11+";
  tracer_names[53] = "S12+";
  tracer_names[54] = "S13+";
  tracer_names[55] = "S14+";
  tracer_names[56] = "S15+";

  return std::make_pair(initial_vector, tracer_names);
}



// *******************************************************************
// Load reference data file
std::vector<std::vector<double> > load_reference_state(const std::string &path)
{
  std::ifstream file(path);
  if (!file) {
    spdlog::error("Cannot open reference file: {}", path);
    return {};
  }

  std::vector<std::vector<double> > reference_state;
  std::string line;
  while (std::getline(file, line)) {
    line = line.substr(0, line.find('#'));
    std::istringstream row(line);
    std::vector<double> values;
    double value;
    while (row >> value) {
      values.push_back(value);
    }
    if (!row.eof()) {
      spdlog::error("Invalid reference data: {}", path);
      return {};
    }
    if (!values.empty()) reference_state.push_back(std::move(values));
  }
  if (file.bad() || reference_state.empty()) {
    spdlog::error("Cannot read reference data: {}", path);
    return {};
  }
  return reference_state;
}



// *******************************************************************
// Compare every cell and variable, including density and temperature.
bool compare_state(
    const std::vector<std::vector<double> > &result_state,
    const std::vector<std::vector<double> > &reference_state,
    double relative_tolerance = 1e-4,
    double absolute_tolerance = 1e-8)
{
  if (result_state.empty() || result_state.size() != reference_state.size()) {
    spdlog::error("States are empty or have different cell counts");
    return false;
  }

  bool test_status = true;
  for (size_t i = 0; i < result_state.size(); ++i) {
    if (result_state[i].size() != reference_state[i].size()) {
      spdlog::error("Different vector sizes in cell {}", i + 1);
      test_status = false;
      continue;
    }
    for (size_t j = 0; j < result_state[i].size(); ++j) {
      const double result    = result_state[i][j];
      const double reference = reference_state[i][j];
      const double diff      = std::abs(result - reference);
      const double absSum    = std::abs(result) + std::abs(reference);
      const double relDiff   = (absSum == 0.0) ? 0.0 : diff / absSum;
      if (diff > absolute_tolerance && relDiff > relative_tolerance) {
        spdlog::error(
            "Cell {}, element {}: result={:.6e}, reference={:.6e}, rel = {:.6e}, abs = {:.6e}",
            i + 1, j, result, reference, relDiff, diff);
        test_status = false;
      }
    }
  }
  return test_status;
}
// *******************************************************************


// *******************************************************************
int main(void)
{
#ifdef NDEBUG
  spdlog::set_level(spdlog::level::info);
  spdlog::flush_on(spdlog::level::err);
#else
  spdlog::set_level(spdlog::level::trace);
  spdlog::flush_on(spdlog::level::trace);
#endif

  spdlog::info("[NEMO::hii40-lexington] Initialising H II region simulation");
#ifdef NEMO_OMP
  spdlog::info(
      "NEMO::hii40-lexington test: git {} running, nth={}", git_version,
      omp_get_max_threads());
#else
  spdlog::info(
      "NEMO::hii40-lexington test: git {} running, nth=1", git_version);
#endif

  double Gamma                       = 5.0 / 3.0;
  auto init                          = setup_initial_state();
  std::vector<double> initial_vector = init.first;
  std::vector<string> tracers        = init.second;
  const int n_state_vars             = initial_vector.size();
  const int n_tracers                = n_state_vars - 2;

  // Radiation source(s)
  std::vector<struct nemo_star_properties> rs;
  struct nemo_star_properties s;
  s.SED_type   = RT_EFFECT_MP10_BLACKBODY;
  s.Luminosity = 3.07e39;  // erg/s
  s.Teff       = 4.0e4;    // K
  s.Radius     = 18.6517;  // solar radii
  s.SED_info   = "";       //? do i need to set this?
  s.id         = 0;
  rs.push_back(s);

  struct nemo_options flags;
  flags.Min_Temperature         = 1.0e3;
  flags.Max_Temperature         = 1.0e9;
  flags.nemo_reltol             = 1.0e-6;
  flags.nemo_abstol             = 1.0e-12;
  flags.nemo_minerg             = 1.0e-20;
  flags.timestep_limiter        = 1;
  flags.update_energy           = 1;
  flags.cooling                 = 1;
  flags.collisional_ionization  = 1;
  flags.photoionization         = 1;
  flags.charge_exchange         = 1;
  flags.cosmicray_ionization    = 1;
  flags.radiative_recombination = 1;
  flags.onthespot_radiation     = 1;

  // Initialise the NEMO class.
  static class nemo *MP;
#ifdef NEMO_OMP
  #pragma omp threadprivate(MP)

  #pragma omp parallel
  {
    MP = new class nemo(
        n_state_vars, n_tracers, 0, 1, tracers, flags, rs, Gamma);
  }
#else
  MP = new class nemo(n_state_vars, n_tracers, 0, 1, tracers, flags, rs, Gamma);
#endif
  MP->get_radiation_src_info(rs);


  // Set up the spherical shell grid
  const int n_cells          = 20;  // Number of spherical shells
  const double cm_per_parsec = 3.085677581e18;
  const double pi            = 3.141592653589793;
  // Simulation domain: radial boundaries
  const double r_min = 0.0 * cm_per_parsec;   // Inner radius: 0.01 pc
  const double r_max = 5.76 * cm_per_parsec;  // Outer radius: 6.0 pc
  const double dr    = (r_max - r_min) / n_cells;

  std::vector<double> r_inner(n_cells);   // Inner radii
  std::vector<double> r_outer(n_cells);   // Outer radii
  std::vector<double> r_center(n_cells);  // Radial midpoints
  std::vector<double> thickness(
      n_cells);  // Shell thicknesses //todo: remove this later
  std::vector<double> volume(n_cells);  // Shell volumes

  for (int i = 0; i < n_cells; ++i) {
    r_inner[i]   = r_min + i * dr;
    r_outer[i]   = r_min + (i + 1) * dr;
    r_center[i]  = 0.5 * (r_inner[i] + r_outer[i]);  // todo: remove this later
    thickness[i] = dr;                               // uniform thickness
    volume[i] =
        (4.0 * pi / 3.0) * (std::pow(r_outer[i], 3) - std::pow(r_inner[i], 3));
  }

  // Initial gas state of the simulation
  std::vector<std::vector<double> > initial_state(n_cells, initial_vector);

  // Simulation time settings
  const double start_time = 0.0;       // seconds
  const double end_time   = 3.156e13;  // seconds
  double time             = start_time;
  int step                = 0;
  std::vector<double> advanced_temperature(n_cells, 0.0);


  std::vector<std::vector<struct nemo_cell_rt_data> > cell_rt(n_cells);
  for (int i = 0; i < n_cells; ++i) {
    cell_rt[i].resize(1);
    struct nemo_cell_rt_data &rtdata = cell_rt[i][0];
    rtdata.star_id                   = rs[0].id;  // Radiation-source identifier
    rtdata.dS                        = thickness[i];  // Shell thickness, cm
    rtdata.Vshell                    = volume[i];     // Shell volume, cm3
    rtdata.tau.resize(rs[0].Nbins);                   // Number of energy bins
    rtdata.dtau.resize(rs[0].Nbins);                  // Number of energy bins
  }
  double x                    = 0.0;
  const auto simulation_start = std::chrono::steady_clock::now();

  while (time < end_time) {

    // ---------------------------------------------------------------
    // 1. Calculate attenuation from the current states of all shells.
    spherical_optical_depth(*MP, rs[0], initial_state, cell_rt);

    // ---------------------------------------------------------------
    // 2. Select the global timestep from the minimum shell timescale.
    const double remaining_time = end_time - time;
    double dt                   = remaining_time;
#ifdef NEMO_OMP
#pragma omp parallel for schedule(static) reduction(min: dt)
#endif
    for (int i = 0; i < n_cells; ++i) {
      const double cell_dt =
          MP->timescales(initial_state[i].data(), cell_rt[i], Gamma);
      dt = std::min(dt, cell_dt);
    }
    dt = std::min(std::max(dt, 0.1 * time), remaining_time);
    dt = std::min(dt, remaining_time);

    // ---------------------------------------------------------------
    // 3. Advance every shell through the same time interval.
#ifdef NEMO_OMP
    #pragma omp parallel
    {
      #pragma omp for
#endif
      for (int i = 0; i < n_cells; ++i) {
        std::vector<double> advanced_vector(initial_vector.size(), 0.0);

        // set_shell_radiation(i);

        // update all cells
        const int error = MP->update_single_cell(
            initial_state[i].data(),  // Current state of shell i
            cell_rt[i],               // Radiation data for this shell
            advanced_vector.data(),   // Advanced state
            dt,                       // Timestep, seconds
            Gamma,                    // Adiabatic index
            &x                        // Auxiliary output pointer
        );

        if (error != 0) {
          throw std::runtime_error(
              "[NEMO::hii40-lexington] NEMO update failed in shell "
              + std::to_string(i));
        }

        // Replace ith shell current state for the next timestep.
        initial_state[i] = advanced_vector;
      }
#ifdef NEMO_OMP
    }
#endif

    // ---------------------------------------------------------------
    // 4. Advance simulation time
    time += dt;
    ++step;

    const double runtime =
        std::chrono::duration<double>(
            std::chrono::steady_clock::now() - simulation_start)
            .count();
    spdlog::info(
        "New time: {:.6e}   dt: {:.6e}   steps: {:8d}   runtime: {:12.4e} s",
        time, dt, step, runtime);
  }


  // Calculate final shell temperatures for output.
#ifdef NEMO_OMP
  #pragma omp parallel for
#endif
  for (int i = 0; i < n_cells; ++i) {
    advanced_temperature[i] = MP->Temperature(initial_state[i].data(), Gamma);
  }

  // Assign Temperature to pressure variable for saving
  for (int i = 0; i < n_cells; ++i) {
    initial_state[i][1] = advanced_temperature[i];
    initial_state[i][0] = r_center[i];
  }

  for (int i = 0; i < n_cells; ++i) {
    spdlog::info(
        "optical_depth[{}] = {:8.3e}", i, fmt::join(cell_rt[i][0].tau, ", "));
  }

  // Load reference data
  std::vector<std::vector<double> > reference_state;
  reference_state = load_reference_state("data/hii40_lexington_reference.txt");

  // Compare the final state with reference state
  const bool test_status = compare_state(initial_state, reference_state);
  spdlog::info("Gas state comparison: {}", test_status ? "PASSED" : "FAILED");



  // Open file for writing final state
  std::ofstream outfile("output_states_hii40.txt");
  if (!outfile.is_open()) {
    spdlog::error("Failed to open output file");
    return 1;
  }
  // Write header line
  std::string header = "radius temperature";
  for (int i = 2; i < n_state_vars; i++)
    header += " " + tracers[i - 2];
  outfile << header << "\n";
  outfile.setf(std::ios::scientific);
  outfile.precision(6);  // 6 decimal places in scientific notation
  // save data in sequence (not in OMP loop)
  for (int i = 0; i < n_cells; i++) {
    for (int j = 0; j < n_state_vars; j++) {
      outfile << initial_state[i][j] << "  ";
    }
    outfile << "\n";
  }
  outfile.close();
  spdlog::info("Output saved to output_states_hii40.txt");

  // TODO: Remove the rest later
  // Save the final temperature at each shell's radial midpoint.
  const std::string profile_path = "temperature_profile.txt";
  std::ofstream file(profile_path);
  if (!file) {
    throw std::runtime_error("Cannot open " + profile_path);
  }

  file << "# Radius_pc  Temperature_K\n";
  file << std::scientific << std::setprecision(8);
  for (int i = 0; i < n_cells; ++i) {
    const double temperature = MP->Temperature(initial_state[i].data(), Gamma);
    file << r_center[i] / cm_per_parsec << "  " << temperature << '\n';
  }

  file.close();
  if (!file) {
    throw std::runtime_error("Failed to write " + profile_path);
  }
  spdlog::debug(
      "[NEMO::hii40-lexington] Saved temperature profile for {} cells",
      n_cells);



#ifdef NEMO_OMP
  #pragma omp parallel
  {
#endif
    delete MP;
#ifdef NEMO_OMP
  }
#endif

  return test_status ? 0 : 1;
}
// *******************************************************************
