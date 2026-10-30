/// isochoric-cooling.cpp
///
/// Isochoric radiative cooling test.
///
/// The gas initially has a temperature of 1e9 K and a mass density chosen
/// to give an electron number density of 100 cm^-3. The composition includes
/// H, He, C, N, O, Ne, Si, S, and Fe, with elemental mass fractions adopted
/// from Asplund et al. (2009). Each element is distributed equally among
/// its Z+1 charge states.
///
/// Radiative cooling updates the gas pressure while mass density and all
/// elemental and ion mass fractions remain fixed.
///
/// NEMO advances the gas to each elapsed time in the reference data, cooling
/// the gas using the same time intervals. At each step, the test calculates
/// the temperature and compares it with the reference temperature.

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "nemo/nemo.h"
#include "nemo_git_version.h"

bool compare_scalar(
    double result,
    double reference,
    double relative_tolerance = 1.0e-6,
    double absolute_tolerance = 1.0e-14)
{
  if (!std::isfinite(result) || !std::isfinite(reference)) return false;
  const double diff    = std::abs(result - reference);
  const double absSum  = std::abs(result) + std::abs(reference);
  const double relDiff = (absSum == 0.0) ? 0.0 : diff / absSum;
  return diff <= absolute_tolerance || relDiff <= relative_tolerance;
}


std::pair<std::vector<double>, std::vector<std::string> >
setup_isochoric_cooling_test()
{
  // Initially, every charge state has an equal fraction of its element.
  // Cool from T = 1e9 K at constant density, with n_e = 100 cm^-3.
  spdlog::info("Initialising NEMO isochoric cooling test");

  const std::vector<std::string> elements = {"H",  "He", "C", "N", "O",
                                             "Ne", "Si", "S", "Fe"};
  // Elemental mass fractions adopted from Asplund et al. (2009),
  // with helium adjusted so the included elements sum to one.
  std::vector<double> element_mass_fractions = {
      0.7381,     // H
      0.2495736,  // He
      0.002368,   // C
      0.000693,   // N
      0.005739,   // O
      0.0012579,  // Ne
      0.0006656,  // Si
      0.0003096,  // S
      0.0012933   // Fe
  };

  // Initial temperature and density
  const double initial_temperature = 1.0e9;  // K
  std::vector<double> initial_vector(2 + elements.size(), 0.0);
  // Density [g/cm^3] giving ne = 100 cm^-3 for the ion mass fractions below.
  initial_vector[0] = 3.847988437033092e-22;
  // tracer elements
  atomic_physics_database atomic_data;
  std::vector<std::string> tracer_names;
  for (size_t e = 0; e < elements.size(); ++e) {
    tracer_names.push_back("X_" + elements[e]);
    initial_vector[2 + e] = element_mass_fractions[e];
  }

  // Calculating total number density
  double total_number_density = 0.0;
  for (size_t e = 0; e < elements.size(); ++e) {

    const std::string &element = elements[e];
    const double element_mass  = atomic_data.get_atomic_mass(element);
    const int atomic_number    = atomic_data.get_atomic_number(element);
    const double elemental_density =
        initial_vector[0] * element_mass_fractions[e] / element_mass;
    const double ion_mass_fraction =
        element_mass_fractions[e] / (atomic_number + 1);

    // Equal charge-state fractions give Z/2 free electrons per atom or ion.
    total_number_density += elemental_density * (1.0 + atomic_number / 2.0);

    tracer_names.push_back(element);  // Neutral ion, charge 0.
    initial_vector.push_back(ion_mass_fraction);
    for (int charge = 1; charge < atomic_number; ++charge) {
      tracer_names.push_back(element + std::to_string(charge) + "+");
      initial_vector.push_back(ion_mass_fraction);
    }
  }

  // Ideal-gas pressure
  const double pressure =
      total_number_density * atomic_data.kB() * initial_temperature;
  initial_vector[1] = pressure;  // Pressure [erg/cm^3]

  return {initial_vector, tracer_names};
}


std::pair<std::vector<double>, std::vector<double> > load_reference_data()
{
  // Load NebulaPy reference times and temperatures.
  std::ifstream reference_file(ISOCHORIC_REFERENCE_FILE);
  if (!reference_file.is_open()) {
    spdlog::error("Failed to open {}", ISOCHORIC_REFERENCE_FILE);
    exit_nemo(1);
  }
  std::vector<double> reference_elapsed_times;
  std::vector<double> reference_temperatures;
  std::string line;
  while (std::getline(reference_file, line)) {
    std::istringstream row(line);
    row >> std::ws;
    if (row.eof() || row.peek() == '#') continue;
    double time, temperature;
    if (!(row >> time >> temperature)) {
      spdlog::error("Invalid reference row: {}", line);
      exit_nemo(1);
    }
    reference_elapsed_times.push_back(time);
    reference_temperatures.push_back(temperature);
  }
  if (reference_file.bad()) {
    spdlog::error("Error reading reference data");
    exit_nemo(1);
  }

  return {reference_elapsed_times, reference_temperatures};
}

int main(void)
{
  spdlog::set_level(spdlog::level::info);
  spdlog::flush_on(spdlog::level::err);
  spdlog::info(
      "NEMO isochoric cooling test: git {} running, nth=1", git_version);

  struct nemo_options flags {};
  flags.Min_Temperature         = 1.0e3;
  flags.Max_Temperature         = 1.0e9;
  flags.nemo_reltol             = 1.0e-8;
  flags.nemo_abstol             = 1.0e-12;
  flags.nemo_minerg             = 1.0e-20;
  flags.timestep_limiter        = 0;
  flags.update_energy           = 1;
  flags.cooling                 = 1;
  flags.collisional_ionization  = 0;
  flags.photoionization         = 0;
  flags.charge_exchange         = 0;
  flags.cosmicray_ionization    = 0;
  flags.radiative_recombination = 0;
  flags.onthespot_radiation     = 0;

  std::vector<struct nemo_star_properties> rs;
  std::vector<struct nemo_cell_rt_data> cell_rt;
  const double eos_gamma = 5.0 / 3.0;

  // Initialise the test state and load reference checkpoints.
  auto initial_state                     = setup_isochoric_cooling_test();
  std::vector<double> &initial_vector    = initial_state.first;
  std::vector<std::string> &tracer_names = initial_state.second;
  auto reference_data                    = load_reference_data();
  const std::vector<double> &reference_elapsed_times = reference_data.first;
  const std::vector<double> &reference_temperatures  = reference_data.second;
  const int N_state_vars                             = initial_vector.size();
  const int n_tracers                                = tracer_names.size();

  flags.nv_prim       = N_state_vars;
  flags.n_scalar      = n_tracers;
  flags.prim_dens_idx = 0;
  flags.prim_pres_idx = 1;
  class nemo MP(flags, tracer_names, rs, eos_gamma);

  std::vector<double> input_vector = initial_vector;
  std::vector<double> output_vector(N_state_vars);

  bool test_status = true;
  double x         = 0.0;

  for (size_t step = 0; step < reference_elapsed_times.size(); ++step) {

    if (step > 0) {
      // CVODE advances from the previous reference time to this exact time.
      const double timestep =
          reference_elapsed_times[step] - reference_elapsed_times[step - 1];
      const int err = MP.update_single_cell(
          input_vector, cell_rt, output_vector, timestep, eos_gamma, &x);
      if (err) {
        spdlog::error("Error updating cooling at step {}: code {}", step, err);
        return 1;
      }
      input_vector.swap(output_vector);
    }

    // Isochoric evolution requires constant mass density (fixed cell volume).
    // Compare the density ratio to detect relative changes in the small CGS
    // density.
    const bool density_passed =
        compare_scalar(input_vector[0] / initial_vector[0], 1.0, 1.0e-8, 0.0);
    if (!density_passed) {
      spdlog::error(
          "Step {}: density changed: initial={:.10e}, current={:.10e} g/cm^3",
          step, initial_vector[0], input_vector[0]);
    }

    // Verify that elemental and ion mass fractions remain unchanged during
    // cooling.
    bool mass_fractions_passed = true;
    for (int i = 2; i < N_state_vars; ++i) {
      if (!compare_scalar(input_vector[i], initial_vector[i], 1.0e-8, 0.0)) {
        spdlog::error(
            "Step {}: mass fraction {} changed", step, tracer_names[i - 2]);
        mass_fractions_passed = false;
      }
    }

    const bool constraints_passed = density_passed && mass_fractions_passed;
    test_status                   = test_status && constraints_passed;

    // Compare temperature with the reference only if both constraints pass.
    if (constraints_passed) {
      const double temperature = MP.Temperature(input_vector, eos_gamma);
      const double reference_temperature = reference_temperatures[step];
      const bool test_passed =
          compare_scalar(temperature, reference_temperature);
      test_status = test_status && test_passed;
      spdlog::info(
          "Time step: {:.6e} s | Constraints check: {} | Temperature check: {}",
          reference_elapsed_times[step],
          constraints_passed ? "PASSED" : "FAILED",
          test_passed ? "PASSED" : "FAILED");
    }
  }


  spdlog::info(
      "NEMO isochoric cooling test: {}", test_status ? "PASSED" : "FAILED");
  return test_status ? 0 : 1;
}
