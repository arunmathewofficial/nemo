/// test-nemo.cpp
///
/// For testing nemo module
///

#ifdef SPDLOG_FWD
#include <spdlog/fwd.h>
#endif
#include <spdlog/sinks/rotating_file_sink.h>
#include <spdlog/spdlog.h>

#include <vector>
#include <stdexcept>

#include "nemo-datastructs.h"
#include "nemo/nemo.h"
#include "nemo_git_version.h"

// *******************************************************************
// check that result and reference are compatible.
bool compare_vectors(
      const std::vector<double>& a,
      const std::vector<double>& b,
      double relTol = 1e-8,
      double absTol = 1e-14)
{
    if (a.size() != b.size()) {
        throw std::invalid_argument("Vectors must have the same size");
    }

  for (size_t i = 0; i < a.size(); ++i) {
    double diff = std::abs(a[i] - b[i]);
    double absSum = std::abs(a[i]) + std::abs(b[i]);
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

  spdlog::info("NEMO tests: git {} running", git_version);

  std::string tr[5] = {"X_H", "X_He", "H", "He", "He1+"};
  struct which_physics wp;
  wp.dynamics          = 1;
  wp.raytracing        = 0;
  wp.cooling           = 0;
  wp.chemistry         = 1;
  wp.coll_ionisation   = 1;
  wp.rad_recombination = 1;
  wp.phot_ionisation   = 0;
  wp.charge_exchange   = 1;
  wp.sat_thermal_cond  = 1;
  wp.tc_strength       = 1;
  wp.update_erg        = 0;
  wp.MP_timestep_limit = 1;
  wp.MinTemperature    = 10.0;
  wp.MaxTemperature    = 1.0e9;
  wp.H_MassFrac        = 0.0;
  wp.Helium_MassFrac   = 0.0;
  wp.Metal_MassFrac    = 0.0;
  wp.compton_cool      = 0;
  wp.wind_acceleration = 0;
  wp.shocked_tracer    = 0;

  struct rad_sources rs;
  rs.Nsources = 0;
  rs.sources.resize(0);
  rs.src_data.resize(0);

  // nemo(Nel(state vec), Ntracer, dens-idx, pres-idx, tracer-names, extra-physics,
  // rad-sources, gamma)
  class nemo MP(7, 5, 0, 1, tr, &wp, &rs, 5.0 / 3.0);

  std::vector<double> input(7, 0.0);
  input[0]  = 1.0e-22;  // density
  input[1]  = 1.0e-9;   // pressure (~1e5 K)
  input[2] = 0.72;
  input[3] = 0.28;
  input[4] = 0.36;
  input[5] = 0.05;
  input[6] = 0.1;
  spdlog::info("input  [{:8.3e}", fmt::join(input, ", "));

  std::vector<double> output(7, 0.0);
  double x = 0.0;

  MP.update_single_cell(
      input.data(), 0, rs.src_data, 0, rs.src_data, output.data(), 1e14,
      5. / 3., 0, &x);

  spdlog::info("output [{:20.16e}]", fmt::join(output, ", "));

  std::vector<double> reference = {1.0e-22, 1.0e-09, 7.2e-01, 2.8e-01, 1.4590868581185208e-05, 2.1600536773812401e-04, 1.3033378001396104e-01};

  bool result = compare_vectors(output,reference);
  spdlog::info("comparison: {} {}", result, "(if passed, rel-err<1e-8 and abs-err<1e-14 for every element)");
  if (result) return 0;
  else return 1;
}
// *******************************************************************

