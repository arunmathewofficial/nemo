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

#include "nemo-datastructs.h"
#include "nemo/nemo.h"
#include "nemo_git_version.h"

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
  flags.timestep_limiter        = 1;
  flags.update_energy           = 1;
  flags.cooling                 = 1;
  flags.collisional_ionization  = 1;
  flags.photoionization         = 1;
  flags.charge_exchange         = 1;
  flags.cosmicray_ionization    = 1;
  flags.radiative_recombination = 1;


  // nemo(ndim, coord-sys, Nel(state vec), Ntracer, tracers, extra-physics,
  // rad-sources, gamma)
  class nemo MP(15, 5, 0, 1, tr, flags, &rs, 5.0 / 3.0);

  std::vector<double> input(15, 0.0);
  input[0]  = 1.0e-22;  // density
  input[1]  = 1.0e-8;   // pressure (~1e7 K)
  input[10] = 0.72;
  input[11] = 0.28;
  input[12] = 0.36;
  input[13] = 0.05;
  input[14] = 0.1;
  spdlog::info("input  [{:8.3e}", fmt::join(input, ", "));

  std::vector<double> output(15, 0.0);
  double x = 0.0;

  MP.update_single_cell(
      input.data(), 0, rs.src_data, 0, rs.src_data, output.data(), 1e10,
      5. / 3., 0, &x);

  spdlog::info("output [{:8.3e}]", fmt::join(output, ", "));

  return 0;
}
