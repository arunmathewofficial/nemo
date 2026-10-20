/// test-nemo.cpp
///
/// For testing nemo module
///

#ifdef SPDLOG_FWD
#include <spdlog/fwd.h>
#endif
#include <spdlog/sinks/rotating_file_sink.h>
#include <spdlog/spdlog.h>

#include <fstream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

#ifdef NEMO_OMP
#include <omp.h>
#endif

#include "nemo-datastructs.h"
#include "nemo/nemo.h"
#include "nemo_git_version.h"

std::pair<std::vector<double>, std::vector<std::string> > setup_initial_state()
{
  // Initialize reference state vector with all ions zero (fully ionized state,
  // enforced by the conservation equations for each element)
  std::vector<double> reference(101, 0.0);
  reference[0] = 1.0e-24;  // density
  reference[1] = 1.0e-9;   // pressure (~1e5 K)

  // Elemental abundances (from Asplund 2010/Solar System)
  reference[2]  = 0.7381;     // H (72%)
  reference[3]  = 0.2495;     // He (26%)
  reference[4]  = 0.002368;   // C (0.25%)
  reference[5]  = 0.000693;   // N (0.15%)
  reference[6]  = 0.005739;   // O (0.15%)
  reference[7]  = 0.0012579;  // Ne (0.05%)
  reference[8]  = 0.0006656;  // Si (0.025%)
  reference[9]  = 0.0003096;  // S (0.015%)
  reference[10] = 0.0012933;  // Fe (0.01%)

  // Ion concentrations (all set to zero for fully ionized plasma)
  // H: 1 ion (H0) - index 11
  reference[11] = 0.0;  // H0

  // He: 2 ions (He0, He+) - indices 12-13
  reference[12] = 0.0;  // He0
  reference[13] = 0.0;  // He1+

  // C: 6 ions (C0-C6+) - indices 14-19
  reference[14] = 0.0;  // C0
  reference[15] = 0.0;  // C1+
  reference[16] = 0.0;  // C2+
  reference[17] = 0.0;  // C3+
  reference[18] = 0.0;  // C4+
  reference[19] = 0.0;  // C5+

  // N: 7 ions (N0-N7+) - indices 20-26
  reference[20] = 0.0;  // N0
  reference[21] = 0.0;  // N1+
  reference[22] = 0.0;  // N2+
  reference[23] = 0.0;  // N3+
  reference[24] = 0.0;  // N4+
  reference[25] = 0.0;  // N5+
  reference[26] = 0.0;  // N6+

  // O: 8 ions (O0-O8+) - indices 27-34
  reference[27] = 0.0;  // O0
  reference[28] = 0.0;  // O1+
  reference[29] = 0.0;  // O2+
  reference[30] = 0.0;  // O3+
  reference[31] = 0.0;  // O4+
  reference[32] = 0.0;  // O5+
  reference[33] = 0.0;  // O6+
  reference[34] = 0.0;  // O7+

  // Ne: 10 ions (Ne0-Ne9+) - indices 35-44
  reference[35] = 0.0;  // Ne0
  reference[36] = 0.0;  // Ne1+
  reference[37] = 0.0;  // Ne2+
  reference[38] = 0.0;  // Ne3+
  reference[39] = 0.0;  // Ne4+
  reference[40] = 0.0;  // Ne5+
  reference[41] = 0.0;  // Ne6+
  reference[42] = 0.0;  // Ne7+
  reference[43] = 0.0;  // Ne8+
  reference[44] = 0.0;  // Ne9+

  // Si: 14 ions (Si0-Si13+) - indices 45-58
  reference[45] = 0.0;  // Si0
  reference[46] = 0.0;  // Si1+
  reference[47] = 0.0;  // Si2+
  reference[48] = 0.0;  // Si3+
  reference[49] = 0.0;  // Si4+
  reference[50] = 0.0;  // Si5+
  reference[51] = 0.0;  // Si6+
  reference[52] = 0.0;  // Si7+
  reference[53] = 0.0;  // Si8+
  reference[54] = 0.0;  // Si9+
  reference[55] = 0.0;  // Si10+
  reference[56] = 0.0;  // Si11+
  reference[57] = 0.0;  // Si12+
  reference[58] = 0.0;  // Si13+

  // S: 16 ions (S0-S15+) - indices 59-74
  reference[59] = 0.0;  // S0
  reference[60] = 0.0;  // S1+
  reference[61] = 0.0;  // S2+
  reference[62] = 0.0;  // S3+
  reference[63] = 0.0;  // S4+
  reference[64] = 0.0;  // S5+
  reference[65] = 0.0;  // S6+
  reference[66] = 0.0;  // S7+
  reference[67] = 0.0;  // S8+
  reference[68] = 0.0;  // S9+
  reference[69] = 0.0;  // S10+
  reference[70] = 0.0;  // S11+
  reference[71] = 0.0;  // S12+
  reference[72] = 0.0;  // S13+
  reference[73] = 0.0;  // S14+
  reference[74] = 0.0;  // S15+

  // Fe: 26 ions (Fe0-Fe25+) - indices 75-100
  reference[75]  = 0.0;  // Fe0
  reference[76]  = 0.0;  // Fe1+
  reference[77]  = 0.0;  // Fe2+
  reference[78]  = 0.0;  // Fe3+
  reference[79]  = 0.0;  // Fe4+
  reference[80]  = 0.0;  // Fe5+
  reference[81]  = 0.0;  // Fe6+
  reference[82]  = 0.0;  // Fe7+
  reference[83]  = 0.0;  // Fe8+
  reference[84]  = 0.0;  // Fe9+
  reference[85]  = 0.0;  // Fe10+
  reference[86]  = 0.0;  // Fe11+
  reference[87]  = 0.0;  // Fe12+
  reference[88]  = 0.0;  // Fe13+
  reference[89]  = 0.0;  // Fe14+
  reference[90]  = 0.0;  // Fe15+
  reference[91]  = 0.0;  // Fe16+
  reference[92]  = 0.0;  // Fe17+
  reference[93]  = 0.0;  // Fe18+
  reference[94]  = 0.0;  // Fe19+
  reference[95]  = 0.0;  // Fe20+
  reference[96]  = 0.0;  // Fe21+
  reference[97]  = 0.0;  // Fe22+
  reference[98]  = 0.0;  // Fe23+
  reference[99]  = 0.0;  // Fe24+
  reference[100] = 0.0;  // Fe25+

  // Create tracer names vector (N_el - 2 elements)
  std::vector<std::string> tracer_names(99);  // 101 - 2 = 99 elements

  // Element names (X_ prefix)
  tracer_names[0] = "X_H";
  tracer_names[1] = "X_He";
  tracer_names[2] = "X_C";
  tracer_names[3] = "X_N";
  tracer_names[4] = "X_O";
  tracer_names[5] = "X_Ne";
  tracer_names[6] = "X_Si";
  tracer_names[7] = "X_S";
  tracer_names[8] = "X_Fe";

  // Ion names (element name + ion charge)
  // H: 1 ion (H0)
  tracer_names[9] = "H";

  // He: 2 ions (He0, He+)
  tracer_names[10] = "He";
  tracer_names[11] = "He1+";

  // C: 6 ions (C0-C6+)
  tracer_names[12] = "C";
  tracer_names[13] = "C1+";
  tracer_names[14] = "C2+";
  tracer_names[15] = "C3+";
  tracer_names[16] = "C4+";
  tracer_names[17] = "C5+";

  // N: 7 ions (N0-N7+)
  tracer_names[18] = "N";
  tracer_names[19] = "N1+";
  tracer_names[20] = "N2+";
  tracer_names[21] = "N3+";
  tracer_names[22] = "N4+";
  tracer_names[23] = "N5+";
  tracer_names[24] = "N6+";

  // O: 8 ions (O0-O8+)
  tracer_names[25] = "O";
  tracer_names[26] = "O1+";
  tracer_names[27] = "O2+";
  tracer_names[28] = "O3+";
  tracer_names[29] = "O4+";
  tracer_names[30] = "O5+";
  tracer_names[31] = "O6+";
  tracer_names[32] = "O7+";

  // Ne: 10 ions (Ne0-Ne9+)
  tracer_names[33] = "Ne";
  tracer_names[34] = "Ne1+";
  tracer_names[35] = "Ne2+";
  tracer_names[36] = "Ne3+";
  tracer_names[37] = "Ne4+";
  tracer_names[38] = "Ne5+";
  tracer_names[39] = "Ne6+";
  tracer_names[40] = "Ne7+";
  tracer_names[41] = "Ne8+";
  tracer_names[42] = "Ne9+";

  // Si: 14 ions (Si0-Si13+)
  tracer_names[43] = "Si";
  tracer_names[44] = "Si1+";
  tracer_names[45] = "Si2+";
  tracer_names[46] = "Si3+";
  tracer_names[47] = "Si4+";
  tracer_names[48] = "Si5+";
  tracer_names[49] = "Si6+";
  tracer_names[50] = "Si7+";
  tracer_names[51] = "Si8+";
  tracer_names[52] = "Si9+";
  tracer_names[53] = "Si10+";
  tracer_names[54] = "Si11+";
  tracer_names[55] = "Si12+";
  tracer_names[56] = "Si13+";

  // S: 16 ions (S0-S15+)
  tracer_names[57] = "S";
  tracer_names[58] = "S1+";
  tracer_names[59] = "S2+";
  tracer_names[60] = "S3+";
  tracer_names[61] = "S4+";
  tracer_names[62] = "S5+";
  tracer_names[63] = "S6+";
  tracer_names[64] = "S7+";
  tracer_names[65] = "S8+";
  tracer_names[66] = "S9+";
  tracer_names[67] = "S10+";
  tracer_names[68] = "S11+";
  tracer_names[69] = "S12+";
  tracer_names[70] = "S13+";
  tracer_names[71] = "S14+";
  tracer_names[72] = "S15+";

  // Fe: 26 ions (Fe0-Fe25+)
  tracer_names[73] = "Fe";
  tracer_names[74] = "Fe1+";
  tracer_names[75] = "Fe2+";
  tracer_names[76] = "Fe3+";
  tracer_names[77] = "Fe4+";
  tracer_names[78] = "Fe5+";
  tracer_names[79] = "Fe6+";
  tracer_names[80] = "Fe7+";
  tracer_names[81] = "Fe8+";
  tracer_names[82] = "Fe9+";
  tracer_names[83] = "Fe10+";
  tracer_names[84] = "Fe11+";
  tracer_names[85] = "Fe12+";
  tracer_names[86] = "Fe13+";
  tracer_names[87] = "Fe14+";
  tracer_names[88] = "Fe15+";
  tracer_names[89] = "Fe16+";
  tracer_names[90] = "Fe17+";
  tracer_names[91] = "Fe18+";
  tracer_names[92] = "Fe19+";
  tracer_names[93] = "Fe20+";
  tracer_names[94] = "Fe21+";
  tracer_names[95] = "Fe22+";
  tracer_names[96] = "Fe23+";
  tracer_names[97] = "Fe24+";
  tracer_names[98] = "Fe25+";

  return std::make_pair(reference, tracer_names);
}


// *******************************************************************
// check that result and reference are compatible.
bool compare_vectors(
    const std::vector<double> &a,
    const std::vector<double> &b,
    double relTol = 1e-8,
    double absTol = 1e-14)
{
  if (a.size() != b.size()) {
    throw std::invalid_argument("Vectors must have the same size");
  }

  for (size_t i = 0; i < a.size(); ++i) {
    double diff    = std::abs(a[i] - b[i]);
    double absSum  = std::abs(a[i]) + std::abs(b[i]);
    double relDiff = (absSum == 0.0) ? 0.0 : diff / absSum;
    if (diff > absTol && relDiff > relTol) {
      return false;
    }
  }
  return true;
}
// *******************************************************************


// *******************************************************************
int main(int argc, char **argv)
{
#ifdef NDEBUG
  spdlog::set_level(spdlog::level::info);
  spdlog::flush_on(spdlog::level::err);
#else
  spdlog::set_level(spdlog::level::trace);
  spdlog::flush_on(spdlog::level::trace);
#endif

#ifdef NEMO_OMP
  spdlog::info(
      "NEMO CIE test: git {} running, nth={}", git_version,
      omp_get_num_procs());
#else
  spdlog::info("NEMO CX-CIE test: git {} running, nth=1", git_version);
#endif


  // Microphysics settings
  struct which_physics wp;
  wp.coll_ionisation   = 1;
  wp.rad_recombination = 1;
  wp.charge_exchange   = 0;
  wp.cooling           = 0;
  wp.phot_ionisation   = 0;
  wp.update_erg        = 0;
  wp.MP_timestep_limit = 1;
  wp.MinTemperature    = 1.0e1;
  wp.MaxTemperature    = 1.0e9;

  struct rad_sources rs;
  rs.Nsources = 0;
  rs.sources.resize(0);
  rs.src_data.resize(0);

  double eos_gamma = 5. / 3.;
  // ############################
  // Get the vector length first
  std::pair<std::vector<double>, std::vector<std::string> > init =
      setup_initial_state();

  std::vector<double> initial_state = init.first;
  std::vector<string> tracer_names  = init.second;
  int N_el                          = initial_state.size();

  // nemo(Nel(state vec), Ntracer, dens-idx, pres-idx, tracer-names,
  // extra-physics, rad-sources, gamma)
  static class nemo *MP;
#ifdef NEMO_OMP
  #pragma omp threadprivate(MP)

  #pragma omp parallel
  {
    MP = new class nemo(
        N_el, N_el - 2, 0, 1, tracer_names.data(), &wp, &rs, eos_gamma);
  }
#else
  MP = new class nemo(
      N_el, N_el - 2, 0, 1, tracer_names.data(), &wp, &rs, eos_gamma);
#endif



  // Set number of state vectors as a variable
  int num_vectors = 100;  // This can be changed to any number

  // Create input vectors by copying the initial state
  std::vector<std::vector<double> > inputs(
      num_vectors, std::vector<double>(N_el, 0.0));
  for (int i = 0; i < num_vectors; i++) {
    inputs[i] = initial_state;
  }

  // Set up logarithmically spaced temperatures from 5000K to 5e8K
  double T_min = 10000.0;
  double T_max = 5e8;
  std::vector<double> temperatures(num_vectors);
  for (int i = 0; i < num_vectors; i++) {
    // Logarithmic spacing: T(i) = T_min * (T_max/T_min)^(i/(n-1))
    double log_ratio = std::log(T_max / T_min);
    temperatures[i]  = T_min * std::exp(i * log_ratio / (num_vectors - 1));
  }

  std::vector<std::vector<double> > outputs(
      num_vectors, std::vector<double>(N_el, 0.0));
  double x = 0.0;

  // Process each of the state vectors with temperature adjustment
#ifdef NEMO_OMP
  #pragma omp parallel
  {
    #pragma omp for
#endif
    for (int i = 0; i < num_vectors; i++) {
      int err = 0;
      // Apply temperature adjustment
      MP->Set_Temp(inputs[i].data(), temperatures[i], eos_gamma);

      // Run the update
      err = MP->update_single_cell(
          inputs[i].data(), 0, rs.src_data, 0, rs.src_data, outputs[i].data(),
          1e15, 5. / 3., 0, &x);
      if (err) {
        spdlog::error(
            "error updating cell {}: T = {:12.e}", i, temperatures[i]);
        spdlog::error("outputs: {}", outputs[i]);
        exit_nemo(err);
      }

      // Get temperature from output state
      double output_temperature = MP->Temperature(outputs[i].data(), eos_gamma);
      spdlog::info(
          "integrated vector {} with T_in= {:8.3e} K, T_out= {:8.3e} K", i,
          temperatures[i], output_temperature);
      outputs[i][1] = output_temperature;  // store T in pressure var.
    }
#ifdef NEMO_OMP
  }
#endif

  // Open file for writing
  std::ofstream outfile("output_states.txt");

  if (!outfile.is_open()) {
    spdlog::error("Failed to open output file");
    return 1;
  }
  // Write header line
  std::string header = "temperature";
  for (int i = 2; i < N_el; i++)
    header += " " + tracer_names[i - 2];
  outfile << header << "\n";
  outfile.setf(std::ios::scientific);
  outfile.precision(6);  // 6 decimal places in scientific notation
  for (int i = 0; i < num_vectors; i++) {
    for (int j = 1; j < N_el; j++) {
      outfile << outputs[i][j] << "  ";
    }
    outfile << "\n";
  }

  outfile.close();
  spdlog::info("Output saved to output_states.txt");

#ifdef NEMO_OMP
  #pragma omp parallel
  {
#endif
    delete MP;
#ifdef NEMO_OMP
  }
#endif

  return 0;
}
// *******************************************************************
