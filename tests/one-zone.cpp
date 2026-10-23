/// test-nemo.cpp
///
/// For testing nemo module
///

#ifdef SPDLOG_FWD
#include <spdlog/fwd.h>
#endif
#include <spdlog/sinks/rotating_file_sink.h>
#include <spdlog/spdlog.h>

#include <stdexcept>
#include <vector>

#include "nemo-datastructs.h"
#include "nemo/nemo.h"
#include "nemo_git_version.h"

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
int main(void)
{
#ifdef NDEBUG
  spdlog::set_level(spdlog::level::info);
  spdlog::flush_on(spdlog::level::err);
#else
  spdlog::set_level(spdlog::level::trace);
  spdlog::flush_on(spdlog::level::trace);
#endif

  spdlog::info("NEMO tests: git {} running", git_version);

  std::vector<std::string> tr({"X_H", "X_He", "H", "He", "He1+"});
  struct rad_sources rs;
  rs.Nsources = 0;
  rs.sources.resize(0);
  rs.src_data.resize(0);

  struct nemo_options flags;
  flags.Min_Temperature         = 1.0e1;
  flags.Max_Temperature         = 1.0e9;
  flags.nemo_reltol             = 1.0e-6;
  flags.nemo_abstol             = 1.0e-12;
  flags.nemo_minerg             = 1.0e-20;
  flags.timestep_limiter        = 1;
  flags.update_energy           = 0;
  flags.cooling                 = 0;
  flags.collisional_ionization  = 1;
  flags.photoionization         = 0;
  flags.charge_exchange         = 1;
  flags.cosmicray_ionization    = 1;
  flags.radiative_recombination = 1;
  flags.onthespot_radiation     = 0;

  // nemo(Nel(state vec), Ntracer, dens-idx, pres-idx, tracer-names,
  // nemo-flags, rad-sources, gamma)
  class nemo MP(7, 5, 0, 1, tr, flags, &rs, 5.0 / 3.0);

  std::vector<double> input(7, 0.0);
  input[0] = 1.0e-22;  // density
  input[1] = 1.0e-9;   // pressure (7.532e+04 K)
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
  spdlog::info("Temperature: {:8.3e}", MP.Temperature(output.data(), 5. / 3.));

  std::vector<double> reference = {
      1.0e-22,
      1.0e-09,
      7.2e-01,
      2.8e-01,
      1.4592800481802200e-05,
      4.8425499616648750e-05,
      1.3043365641742180e-01};

  bool result = compare_vectors(output, reference);
  spdlog::info("comparison: {:8.3e}", fmt::join(output, ", "));
  spdlog::info("reference:  {:8.3e}", fmt::join(reference, ", "));
  spdlog::info("result: {}", result);
  if (result)
    return 0;
  else
    return 1;
}
// *******************************************************************
