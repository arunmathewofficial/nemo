/// cx-test.cpp
///
/// Single-cell CX equilibrium tests
///

#ifdef SPDLOG_FWD
#include <spdlog/fwd.h>
#endif
#include <spdlog/spdlog.h>

#include <cmath>
#include <limits>
#include <string>
#include <vector>

#include "nemo-datastructs.h"
#include "nemo/nemo.h"
#include "nemo_git_version.h"

struct test_configuration {
  std::vector<double> initial_vector;
  std::vector<std::string> tracer_names;
  double temperature;
  double finish_time;
  const double reference_value;
  std::vector<double> reference_vector;
};


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

// Compare every component and report all mismatches.
bool compare_vector(
    const std::vector<double> &result, const std::vector<double> &reference)
{
  if (result.size() != reference.size()) {
    spdlog::error(
        "Reference vector size {} differs from output size {}",
        reference.size(), result.size());
    return false;
  }
  bool passed = true;
  for (size_t i = 0; i < result.size(); ++i) {
    if (!compare_scalar(result[i], reference[i])) {
      spdlog::error(
          "Vector component {}: result={:.10e}, reference={:.10e}", i,
          result[i], reference[i]);
      passed = false;
    }
  }
  return passed;
}

// Positive ionic charge density in elementary charges per cm^3.
// Assumes the NEMO layout: density, pressure, then named tracers, with
// the fully stripped stage inferred from each elemental mass fraction.
double charge_density(
    const std::vector<double> &input_vector,
    const std::vector<std::string> &tracer_names)
{
  if (input_vector.size() != tracer_names.size() + 2) {
    spdlog::error(
        "Cannot calculate charge density: vector/tracer size mismatch");
    return std::numeric_limits<double>::quiet_NaN();
  }
  atomic_physics_database atomic_data;
  for (const std::string &name : tracer_names) {
    const std::string species =
        name.substr(0, 2) == "X_" ? name.substr(2) : name;
    if (!atomic_data.find_species(species)) {
      return std::numeric_limits<double>::quiet_NaN();
    }
  }

  double total_charge = 0.0;
  for (size_t e = 0; e < tracer_names.size(); ++e) {
    if (tracer_names[e].substr(0, 2) != "X_") continue;
    const std::string element = tracer_names[e].substr(2);
    double species_sum        = 0.0;
    double element_charge     = 0.0;
    for (size_t j = 0; j < tracer_names.size(); ++j) {
      const std::string &species = tracer_names[j];
      if (species.substr(0, 2) == "X_"
          || atomic_data.get_element_name(species) != element)
        continue;
      species_sum += input_vector[j + 2];
      element_charge += atomic_data.get_charge(species) * input_vector[j + 2];
    }
    element_charge += atomic_data.get_atomic_number(element)
                      * (input_vector[e + 2] - species_sum);
    total_charge += element_charge / atomic_data.get_atomic_mass(element);
  }
  return input_vector[0] * total_charge;
}

// Case 1 : HOp CX reaction
test_configuration setup_HOp_charge_exchange_test()
{
  /*
   * HOp CX relaxation test.
   * Reactions:
   *   O^+ + H^0 --> O^0 + H^+  (forward, coefficient k_forward)
   *   O^0 + H^+ --> O^+ + H^0  (reverse, coefficient k_reverse).
   * An electron is transferred between hydrogen and oxygen, conserving
   * the total number of nuclei of each element and the total charge.
   *
   * Initially, hydrogen is half ionised and oxygen is entirely O^+.
   * The total number densities are n_H = 1 cm^-3 and n_O = 1e-4 cm^-3,
   * with an initial temperature of 1000 K. Evolve one independent cell
   * with only charge exchange enabled.
   *
   * At equilibrium, the forward and reverse reaction rates balance:
   *   k_forward n(O^+) n(H^0) = k_reverse n(O^0) n(H^+).
   * Therefore, compare the final number-density ratio
   *   n(O^+) n(H^0) / [n(O^0) n(H^+)]
   * with k_reverse(T) / k_forward(T), evaluated at the final temperature.
   * Agreement within tolerance indicates CX equilibrium;
   * both reactions continue, but produce no net change in abundances.
   */
  spdlog::info(
      "Initialising Case(1) NEMO CX equilibrium test: O+ + H <-> O + H+");
  const double temperature      = 1000.0;                   // K
  const double finish_time      = 1.0e13;                   // seconds
  const double hydrogen_density = 1.0;                      // n_H = 1 cm^-3,
  const double oxygen_density = 1.0e-4 * hydrogen_density;  // n_O= n_H * 1e-4.
  const double hydrogen_mass  = 1.6738e-24;                 // g; NEMO m_H()
  const double oxygen_mass    = 2.6567628e-23;              // g; NEMO m_O()
  const double density =
      hydrogen_density * hydrogen_mass + oxygen_density * oxygen_mass;
  const double X_H = hydrogen_density * hydrogen_mass / density;
  const double X_O = oxygen_density * oxygen_mass / density;

  // Layout: density, pressure, elemental mass fractions, ion mass fractions.
  // H+ and O8+ are implicit, determined by elemental conservation.
  std::vector<std::string> tracer_names = {
      "X_H", "X_O", "H", "O", "O1+", "O2+", "O3+", "O4+", "O5+", "O6+", "O7+"};
  std::vector<double> initial_vector(2 + tracer_names.size(), 0.0);
  initial_vector[0] = density;    // Density
  initial_vector[1] = 1.0e-12;    // Pressure
  initial_vector[2] = X_H;        // H elemental mass fraction
  initial_vector[3] = X_O;        // O elemental mass fraction
  initial_vector[4] = 0.5 * X_H;  // Half neutral H, half H+.
  initial_vector[5] = 0.0;        // Neutral oxygen.
  initial_vector[6] = X_O;        // All oxygen initially in O1+ state.

  // Ratio of reverse to forward rate coefficients at T = 1000 K.
  // Rate-fit source: Kingdon & Ferland (1996), ApJS 106, 205--211.
  // https://doi.org/10.1086/192335
  // NEMO entries: ar-o1+h (forward), ar-o+h1 (reverse).
  // Numerical reference uses NEMO's interpolated tables, not a published ratio.
  const double reference_value = 5.1811635240e-01;  // k_reverse / k_forward
  spdlog::info(
      "Case 1: reference k_reverse / k_forward at T={:.0f} K = {:.10e}",
      temperature, reference_value);

  return {initial_vector, tracer_names,    temperature,
          finish_time,    reference_value, {}};
}

// Case 2: HNp CX relaxation towards equilibrium.
test_configuration setup_HNp_charge_exchange_test()
{
  // N^+ + H^0 <-> N^0 + H^+; initially half-ionised H and all N^+.
  spdlog::info(
      "Initialising Case(2) NEMO CX equilibrium test: N+ + H <-> N + H+");
  const double temperature      = 10000.0;  // K
  const double finish_time      = 5.0e13;   // seconds
  const double hydrogen_density = 1.0;      // cm^-3
  const double nitrogen_density = 1.0e-4 * hydrogen_density;
  const double hydrogen_mass    = 1.6738e-24;
  const double nitrogen_mass    = 2.325892e-23;
  const double density =
      hydrogen_density * hydrogen_mass + nitrogen_density * nitrogen_mass;
  const double X_H = hydrogen_density * hydrogen_mass / density;
  const double X_N = nitrogen_density * nitrogen_mass / density;
  std::vector<std::string> tracer_names = {"X_H", "X_N", "H",   "N",   "N1+",
                                           "N2+", "N3+", "N4+", "N5+", "N6+"};
  std::vector<double> initial_vector(2 + tracer_names.size(), 0.0);
  initial_vector[0] = density;
  initial_vector[1] = 1.0e-12;  // Replaced by Set_Temp before evolution.
  initial_vector[2] = X_H;
  initial_vector[3] = X_N;
  initial_vector[4] = 0.5 * X_H;
  initial_vector[6] = X_N;

  // Ratio of reverse to forward rate coefficients at T = 10000 K.
  // Rate-fit source: Kingdon & Ferland (1996), ApJS 106, 205--211.
  // https://doi.org/10.1086/192335
  // NEMO entries: ar-n1+h (forward), ar-n+h1 (reverse).
  // Numerical reference uses NEMO's interpolated tables, not a published ratio.
  const double reference_value = 1.5199216941e+00;
  spdlog::info(
      "Case 2: reference k_reverse / k_forward at T={:.0f} K = {:.10e}",
      temperature, reference_value);

  return {initial_vector, tracer_names,    temperature,
          finish_time,    reference_value, {}};
}

// Case 3: HCp CX relaxation towards equilibrium.
test_configuration setup_HCp_charge_exchange_test()
{
  // C^+ + H^0 <-> C^0 + H^+; initially half-ionised H and all C^+.
  spdlog::info(
      "Initialising Case(3) NEMO CX equilibrium test: C+ + H <-> C + H+");
  const double temperature      = 8000.0;  // K
  const double finish_time      = 1.0e16;  // seconds
  const double hydrogen_density = 1.0;     // cm^-3
  const double carbon_density   = 1.0e-4 * hydrogen_density;
  const double hydrogen_mass    = 1.6738e-24;
  const double carbon_mass      = 1.994374e-23;
  const double density =
      hydrogen_density * hydrogen_mass + carbon_density * carbon_mass;
  const double X_H = hydrogen_density * hydrogen_mass / density;
  const double X_C = carbon_density * carbon_mass / density;
  std::vector<std::string> tracer_names = {"X_H", "X_C", "H",   "C",  "C1+",
                                           "C2+", "C3+", "C4+", "C5+"};
  std::vector<double> initial_vector(2 + tracer_names.size(), 0.0);
  initial_vector[0] = density;
  initial_vector[1] = 1.0e-12;  // Replaced by Set_Temp before evolution.
  initial_vector[2] = X_H;
  initial_vector[3] = X_C;
  initial_vector[4] = 0.5 * X_H;
  initial_vector[6] = X_C;

  // Ratio of reverse to forward rate coefficients at T = 8000 K.
  // Rate-fit source: Kingdon & Ferland (1996), ApJS 106, 205--211.
  // https://doi.org/10.1086/192335
  // NEMO entries: ar-c1+h (forward), ar-c+h1 (reverse).
  // Numerical reference uses NEMO's interpolated tables, not a published ratio.
  const double reference_value = 1.3798769809e+03;
  spdlog::info(
      "Case 3: reference k_reverse / k_forward at T={:.0f} K = {:.10e}",
      temperature, reference_value);

  return {initial_vector, tracer_names,    temperature,
          finish_time,    reference_value, {}};
}

// Case 4: HeC2p CX relaxation towards equilibrium.
test_configuration setup_HeC2p_charge_exchange_test()
{
  /*
   * C^{2+} + He^0 --> C^+ + He^+ (forward)
   * C^+ + He^+ --> C^{2+} + He^0 (reverse).
   * Initially helium is half neutral and half He+, and all carbon is C2+.
   * At equilibrium, n(C2+) n(He) / [n(C+) n(He+)] = k_reverse/k_forward.
   * 8000 K lies within both rate fits' validity ranges.
   */
  spdlog::info(
      "Initialising Case(4) NEMO CX equilibrium test: C2+ + He <-> C+ + He+");
  const double temperature    = 8000.0;  // K
  const double finish_time    = 1.0e13;  // seconds
  const double helium_density = 1.0;     // cm^-3
  const double carbon_density = 1.0e-4 * helium_density;
  const double helium_mass    = 6.6464768e-24;
  const double carbon_mass    = 1.994374e-23;
  const double density =
      helium_density * helium_mass + carbon_density * carbon_mass;
  const double X_He = helium_density * helium_mass / density;
  const double X_C  = carbon_density * carbon_mass / density;
  std::vector<std::string> tracer_names = {"X_He", "X_C", "He",  "He1+", "C",
                                           "C1+",  "C2+", "C3+", "C4+",  "C5+"};
  std::vector<double> initial_vector(2 + tracer_names.size(), 0.0);
  initial_vector[0] = density;
  initial_vector[1] = 1.0e-12;  // Replaced by Set_Temp before evolution.
  initial_vector[2] = X_He;
  initial_vector[3] = X_C;
  initial_vector[4] = 0.5 * X_He;
  initial_vector[5] = 0.5 * X_He;
  initial_vector[8] = X_C;

  // Rate-fit source: Arnaud & Rothenflug (1985), A&AS 60, 425.
  // NEMO entries: ar-c2+he (forward), ar-c1+he1 (reverse).
  // Reference at 8000 K, calculated from the fits with NEMO's 300-point
  // temperature grid and linear interpolation; not a published ratio.
  const double reference_value = 3.4590722950e-06;
  spdlog::info(
      "Case 4: reference k_reverse / k_forward at T={:.0f} K = {:.10e}",
      temperature, reference_value);
  return {initial_vector, tracer_names,    temperature,
          finish_time,    reference_value, {}};
}


// Case 5: full-network CX evolution to a non-equilibrium state.
test_configuration setup_full_network_cx_test()
{
  spdlog::info(
      "Initialising Case(5) NEMO CX full-network non-equilibrium test");
  const double temperature = 8000.0;  // Initial temperature, K.
  const double finish_time = 1.0e7;   // One update to t1, seconds.
  // Include all elements represented in NEMO's active CX database.
  const std::vector<std::string> elements      = {"H",  "He", "C", "N", "O",
                                                  "Ne", "Si", "S", "Fe"};
  std::vector<double> elemental_mass_fractions = {
      0.7381,    0.2495736, 0.002368,  0.000693, 0.005739,
      0.0012579, 0.0006656, 0.0003096, 0.0012933};

  atomic_physics_database atomic_data;
  std::vector<std::string> tracer_names;
  std::vector<double> initial_vector(2, 0.0);
  // Set total hydrogen number density to 1 cm^-3.
  initial_vector[0] =
      atomic_data.get_atomic_mass("H") / elemental_mass_fractions[0];
  initial_vector[1] = 1.0e-12;  // Pressure
  for (size_t e = 0; e < elements.size(); ++e) {
    tracer_names.push_back("X_" + elements[e]);
    initial_vector.push_back(elemental_mass_fractions[e]);
  }
  // Populate every stage, including the implicit fully stripped stage,
  // equally by number within each element to provide available reactants.
  for (size_t e = 0; e < elements.size(); ++e) {
    const int Z = atomic_data.get_atomic_number(elements[e]);
    for (int q = 0; q < Z; ++q) {
      tracer_names.push_back(
          q == 0 ? elements[e] : elements[e] + std::to_string(q) + "+");
      initial_vector.push_back(elemental_mass_fractions[e] / (Z + 1));
    }
  }

  // Reference vector
  const std::vector<double> reference_vector = {
      // Density [g cm^-3], pressure [erg cm^-3]
      2.2677144018425685e-24,
      1.8505834520027823e-12,
      // Elemental mass fractions: H, He, C, N, O, Ne, Si, S, Fe
      0.7381,
      0.2495736,
      0.002368,
      0.000693,
      0.005739,
      0.0012579,
      0.0006656,
      0.0003096,
      0.0012933,
      // H: H0 (H+ is implicit)
      0.36904757566459245,
      // He: He0, He+ (He2+ is implicit)
      0.08319108673323329,
      0.0831913174263129,
      // C: C0--C5+ (C6+ is implicit)
      0.0003382856973836309,
      0.00033829735391446874,
      0.0003435802384958418,
      0.0003365662664086913,
      0.00033469901522593915,
      0.00033828571428571436,
      // N: N0--N6+ (N7+ is implicit)
      8.662492707613568e-05,
      8.689538142883452e-05,
      8.75345730610821e-05,
      8.618363942170177e-05,
      8.588647901224602e-05,
      8.662500000000002e-05,
      8.662500000000002e-05,
      // O: O0--O7+ (O8+ is implicit)
      0.0006379770685700943,
      0.0006405875813368865,
      0.0006472939894422271,
      0.000625676388403758,
      0.000636798305580369,
      0.0006376666666666669,
      0.0006376666666666669,
      0.0006376666666666669,
      // Ne: Ne0--Ne9+ (Ne10+ is implicit)
      0.00011435454545454549,
      0.0001143545515735244,
      0.00011723419195439352,
      0.00011493355070296612,
      0.0001108958875872979,
      0.00011435454545454549,
      0.00011435454545454549,
      0.00011435454545454549,
      0.00011435454545454549,
      0.00011435454545454549,
      // Si: Si0--Si13+ (Si14+ is implicit)
      4.437333333333335e-05,
      4.543587610631405e-05,
      4.341225579175121e-05,
      4.586815955382124e-05,
      4.277704188144692e-05,
      4.437333333333335e-05,
      4.437333333333335e-05,
      4.437333333333335e-05,
      4.437333333333335e-05,
      4.437333333333335e-05,
      4.437333333333335e-05,
      4.437333333333335e-05,
      4.437333333333335e-05,
      4.437333333333335e-05,
      // S: S0--S15+ (S16+ is implicit)
      1.8211763903472826e-05,
      1.8211766424321976e-05,
      1.8428197232135277e-05,
      1.8571875731652e-05,
      1.7635220237829726e-05,
      1.821176470588236e-05,
      1.821176470588236e-05,
      1.821176470588236e-05,
      1.821176470588236e-05,
      1.821176470588236e-05,
      1.821176470588236e-05,
      1.821176470588236e-05,
      1.821176470588236e-05,
      1.821176470588236e-05,
      1.821176470588236e-05,
      1.821176470588236e-05,
      // Fe: Fe0--Fe25+ (Fe26+ is implicit)
      4.790000000000002e-05,
      4.819757005284226e-05,
      4.8326594161617237e-05,
      4.717583578554055e-05,
      4.790000000000002e-05,
      4.790000000000002e-05,
      4.790000000000002e-05,
      4.790000000000002e-05,
      4.790000000000002e-05,
      4.790000000000002e-05,
      4.790000000000002e-05,
      4.790000000000002e-05,
      4.790000000000002e-05,
      4.790000000000002e-05,
      4.790000000000002e-05,
      4.790000000000002e-05,
      4.790000000000002e-05,
      4.790000000000002e-05,
      4.790000000000002e-05,
      4.790000000000002e-05,
      4.790000000000002e-05,
      4.790000000000002e-05,
      4.790000000000002e-05,
      4.790000000000002e-05,
      4.790000000000002e-05,
      4.790000000000002e-05,
  };

  return {
      initial_vector,
      tracer_names,
      temperature,
      finish_time,
      std::numeric_limits<double>::quiet_NaN(),
      reference_vector};
}

int main(void)
{
#ifdef NDEBUG
  spdlog::set_level(spdlog::level::info);
  spdlog::flush_on(spdlog::level::err);
#else
  spdlog::set_level(spdlog::level::trace);
  spdlog::flush_on(spdlog::level::trace);
#endif

  spdlog::info("NEMO CX test: git {} running, nth=1", git_version);

  struct nemo_options flags;
  flags.Min_Temperature         = 1.0e1;
  flags.Max_Temperature         = 1.0e9;
  flags.nemo_reltol             = 1.0e-6;
  flags.nemo_abstol             = 1.0e-12;
  flags.nemo_minerg             = 1.0e-20;
  flags.timestep_limiter        = 1;
  flags.update_energy           = 0;
  flags.cooling                 = 0;
  flags.collisional_ionization  = 0;
  flags.photoionization         = 0;
  flags.charge_exchange         = 1;
  flags.cosmicray_ionization    = 0;
  flags.radiative_recombination = 0;
  flags.onthespot_radiation     = 0;

  std::vector<struct nemo_star_properties> rs;
  std::vector<struct nemo_cell_rt_data> cell_rt;

  // Adiabatic index
  double eos_gamma = 5. / 3.;

  // Number of cell
  const int N_cell = 5;  // each test in each cell
  std::vector<std::vector<double> > inputs(N_cell);
  std::vector<std::vector<std::string> > tracer_names(N_cell);
  std::vector<double> temperatures(N_cell);
  std::vector<double> final_temperatures(N_cell);
  std::vector<double> finish_times(N_cell);
  std::vector<std::vector<double> > outputs(N_cell);

  // Case 1: O+ + H <-> O + H+.
  auto HOp_cx_reaction = setup_HOp_charge_exchange_test();
  inputs[0]            = HOp_cx_reaction.initial_vector;
  tracer_names[0]      = HOp_cx_reaction.tracer_names;
  temperatures[0]      = HOp_cx_reaction.temperature;
  finish_times[0]      = HOp_cx_reaction.finish_time;

  // Case 2: N+ + H <-> N + H+.
  auto HNp_cx_reaction = setup_HNp_charge_exchange_test();
  inputs[1]            = HNp_cx_reaction.initial_vector;
  tracer_names[1]      = HNp_cx_reaction.tracer_names;
  temperatures[1]      = HNp_cx_reaction.temperature;
  finish_times[1]      = HNp_cx_reaction.finish_time;

  // Case 3: C+ + H <-> C + H+.
  auto HCp_cx_reaction = setup_HCp_charge_exchange_test();
  inputs[2]            = HCp_cx_reaction.initial_vector;
  tracer_names[2]      = HCp_cx_reaction.tracer_names;
  temperatures[2]      = HCp_cx_reaction.temperature;
  finish_times[2]      = HCp_cx_reaction.finish_time;

  // Case 4: C2+ + He <-> C+ + He+.
  auto HeC2p_cx_reaction = setup_HeC2p_charge_exchange_test();
  inputs[3]              = HeC2p_cx_reaction.initial_vector;
  tracer_names[3]        = HeC2p_cx_reaction.tracer_names;
  temperatures[3]        = HeC2p_cx_reaction.temperature;
  finish_times[3]        = HeC2p_cx_reaction.finish_time;

  // Case 5: all-element CX evolution to t1, with one solver call.
  auto full_network_cx = setup_full_network_cx_test();
  inputs[4]            = full_network_cx.initial_vector;
  tracer_names[4]      = full_network_cx.tracer_names;
  temperatures[4]      = full_network_cx.temperature;
  finish_times[4]      = full_network_cx.finish_time;


  bool test_status = true;
  // Evolve each test cell independently.
  for (int cell = 0; cell < N_cell; ++cell) {

    const int N_state_vars = inputs[cell].size();
    const int n_tracers    = tracer_names[cell].size();
    outputs[cell].resize(N_state_vars);

    flags.nv_prim       = N_state_vars;
    flags.n_scalar      = n_tracers;
    flags.prim_dens_idx = 0;
    flags.prim_pres_idx = 1;

    // A separate solver allows each cell to use a different composition.
    class nemo MP(flags, tracer_names[cell], rs, eos_gamma);

    MP.Set_Temp(inputs[cell], temperatures[cell], eos_gamma);

    double x = 0.0;

    const int err = MP.update_single_cell(
        inputs[cell], cell_rt, outputs[cell], finish_times[cell], eos_gamma,
        &x);

    if (err) {
      spdlog::error("Error updating cell {}: code {}", cell + 1, err);
      return 1;
    }
    final_temperatures[cell] = MP.Temperature(outputs[cell], eos_gamma);
  }


  //////////////////////////////////////////////////////////////////////////
  //*************************************************************************
  // Case 1: HOp reaction check
  const bool HOp_temperature_passed =
      compare_scalar(final_temperatures[0], HOp_cx_reaction.temperature);
  spdlog::info(
      "Case 1: temperature {}: initial={:.0f} K, final={:.0f} K",
      HOp_temperature_passed ? "unchanged" : "changed",
      HOp_cx_reaction.temperature, final_temperatures[0]);
  if (!HOp_temperature_passed) {
    spdlog::error(
        "Case 1: final temperature {:.6e} K differs from expected {:.6e} K",
        final_temperatures[0], HOp_cx_reaction.temperature);
    test_status = false;
  }
  const std::vector<double> &HOp_output_vector = outputs[0];
  const double X_H0                            = HOp_output_vector[4];
  const double X_Hp                            = HOp_output_vector[2] - X_H0;
  const double X_O0                            = HOp_output_vector[5];
  const double X_Op                            = HOp_output_vector[6];
  // Density and atomic masses cancel in this number-density ratio.
  const double HOp_denominator = X_O0 * X_Hp;
  if (HOp_denominator > 0.0) {
    const double HOp_ratio = X_Op * X_H0 / HOp_denominator;
    spdlog::info("Case 1: n(O+) n(H) / [n(O) n(H+)] = {:.10e}", HOp_ratio);
    const bool HOp_passed =
        HOp_temperature_passed
        && compare_scalar(HOp_ratio, HOp_cx_reaction.reference_value);
    test_status = test_status && HOp_passed;
    spdlog::info(
        "NEMO CX equilibrium test: O+ + H <-> O + H+: {}",
        HOp_passed ? "PASSED" : "FAILED");
  }
  else {
    spdlog::error(
        "Case 1 exchange: cannot calculate ratio: denominator is not positive");
    test_status = false;
  }

  //*************************************************************************
  // Case 2: HNp ratio with its equilibrium reference.
  const bool HNp_temperature_passed =
      compare_scalar(final_temperatures[1], HNp_cx_reaction.temperature);
  spdlog::info(
      "Case 2: temperature {}: initial={:.0f} K, final={:.0f} K",
      HNp_temperature_passed ? "unchanged" : "changed",
      HNp_cx_reaction.temperature, final_temperatures[1]);
  if (!HNp_temperature_passed) {
    spdlog::error(
        "Case 2: final temperature differs from expected temperature");
    test_status = false;
  }
  const std::vector<double> &HNp_output_vector = outputs[1];
  const double HNp_X_H0                        = HNp_output_vector[4];
  const double HNp_X_Hp        = HNp_output_vector[2] - HNp_X_H0;
  const double X_N0            = HNp_output_vector[5];
  const double X_Np            = HNp_output_vector[6];
  const double HNp_denominator = X_N0 * HNp_X_Hp;
  if (HNp_denominator > 0.0) {
    const double HNp_ratio = X_Np * HNp_X_H0 / HNp_denominator;
    spdlog::info("Case 2: n(N+) n(H) / [n(N) n(H+)] = {:.10e}", HNp_ratio);
    const bool HNp_passed =
        HNp_temperature_passed
        && compare_scalar(HNp_ratio, HNp_cx_reaction.reference_value);
    test_status = test_status && HNp_passed;
    spdlog::info(
        "NEMO CX equilibrium test: N+ + H <-> N + H+: {}",
        HNp_passed ? "PASSED" : "FAILED");
  }
  else {
    spdlog::error(
        "Case 2: cannot calculate ratio: denominator is not positive");
    test_status = false;
  }

  //*************************************************************************
  // Case 3: HCp ratio with its equilibrium reference.
  const bool HCp_temperature_passed =
      compare_scalar(final_temperatures[2], HCp_cx_reaction.temperature);
  spdlog::info(
      "Case 3: temperature {}: initial={:.0f} K, final={:.0f} K",
      HCp_temperature_passed ? "unchanged" : "changed",
      HCp_cx_reaction.temperature, final_temperatures[2]);
  if (!HCp_temperature_passed) {
    spdlog::error(
        "Case 3: final temperature differs from expected temperature");
    test_status = false;
  }
  const std::vector<double> &HCp_output_vector = outputs[2];
  const double HCp_X_H0                        = HCp_output_vector[4];
  const double HCp_X_Hp        = HCp_output_vector[2] - HCp_X_H0;
  const double X_C0            = HCp_output_vector[5];
  const double X_Cp            = HCp_output_vector[6];
  const double HCp_denominator = X_C0 * HCp_X_Hp;
  if (HCp_denominator > 0.0) {
    const double HCp_ratio = X_Cp * HCp_X_H0 / HCp_denominator;
    spdlog::info("Case 3: n(C+) n(H) / [n(C) n(H+)] = {:.10e}", HCp_ratio);
    if (!std::isfinite(HCp_cx_reaction.reference_value)) {
      spdlog::warn("Case 3: PENDING reference value; test cannot pass yet");
      test_status = false;
    }
    else {
      const bool HCp_passed =
          HCp_temperature_passed
          && compare_scalar(HCp_ratio, HCp_cx_reaction.reference_value);
      test_status = test_status && HCp_passed;
      spdlog::info(
          "NEMO CX equilibrium test: C+ + H <-> C + H+: {}",
          HCp_passed ? "PASSED" : "FAILED");
    }
  }
  else {
    spdlog::error(
        "Case 3: cannot calculate ratio: denominator is not positive");
    test_status = false;
  }

  //*************************************************************************
  // Case 4: HeC2p equilibrium ratio, temperature and charge conservation.
  const std::vector<double> &HeC2p_output_vector = outputs[3];
  const double X_He0                             = HeC2p_output_vector[4];
  const double X_Hep                             = HeC2p_output_vector[5];
  const double HeC2p_X_Cp                        = HeC2p_output_vector[7];
  const double X_C2p                             = HeC2p_output_vector[8];
  const double HeC2p_denominator                 = HeC2p_X_Cp * X_Hep;
  const bool HeC2p_temperature_passed =
      compare_scalar(final_temperatures[3], HeC2p_cx_reaction.temperature);
  spdlog::info(
      "Case 4: temperature {}: initial={:.0f} K, final={:.0f} K",
      HeC2p_temperature_passed ? "unchanged" : "changed",
      HeC2p_cx_reaction.temperature, final_temperatures[3]);

  const bool HeC2p_charge_passed = compare_scalar(
      charge_density(HeC2p_output_vector, tracer_names[3]),
      charge_density(inputs[3], tracer_names[3]));
  spdlog::info(
      "Case 4: charge {}", HeC2p_charge_passed ? "conserved" : "not conserved");


  if (!HeC2p_temperature_passed) {
    spdlog::error(
        "Case 4: final temperature differs from expected temperature");
  }
  if (!HeC2p_charge_passed) {
    spdlog::error("Case 4: total charge is not conserved");
  }

  if (HeC2p_denominator > 0.0) {
    const double HeC2p_ratio = X_C2p * X_He0 / HeC2p_denominator;
    spdlog::info(
        "Case 4: n(C2+) n(He) / [n(C+) n(He+)] = {:.10e}", HeC2p_ratio);
    const bool HeC2p_passed =
        HeC2p_temperature_passed && HeC2p_charge_passed
        && compare_scalar(HeC2p_ratio, HeC2p_cx_reaction.reference_value);
    test_status = test_status && HeC2p_passed;
    spdlog::info(
        "NEMO CX equilibrium test: C2+ + He <-> C+ + He+: {}",
        HeC2p_passed ? "PASSED" : "FAILED");
  }
  else {
    spdlog::error(
        "Case 4: cannot calculate ratio: denominator is not positive");
    test_status = false;
  }

  //*************************************************************************
  // Case 5: NEQ reference and charge conservation.
  const std::vector<double> &full_network_output = outputs[4];
  const bool full_network_charge_passed          = compare_scalar(
      charge_density(full_network_output, tracer_names[4]),
      charge_density(inputs[4], tracer_names[4]));
  spdlog::info(
      "Case 5: charge {}",
      full_network_charge_passed ? "conserved" : "not conserved");
  const bool reference_passed =
      compare_vector(full_network_output, full_network_cx.reference_vector);
  const bool full_network_passed =
      reference_passed && full_network_charge_passed;
  test_status = test_status && full_network_passed;
  spdlog::info(
      "NEMO CX full-network non-equilibrium test: {}",
      full_network_passed ? "PASSED" : "FAILED");

  spdlog::info("NEMO CX tests: {}", test_status ? "ALL PASSED" : "FAILED");
  return test_status ? 0 : 1;
}
// *******************************************************************
