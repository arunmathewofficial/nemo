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

  std::vector<struct nemo_star_properties> rs;
  std::vector<struct nemo_cell_rt_data> cell_rt;

  struct nemo_options flags;
  flags.Min_Temperature         = 1.0e1;
  flags.Max_Temperature         = 1.0e9;
  flags.nemo_reltol             = 1.0e-6;
  flags.nemo_abstol             = 1.0e-12;
  flags.nemo_minerg             = 1.0e-20;
  flags.nv_prim                 = 15;
  flags.n_scalar                = 5;
  flags.prim_dens_idx           = 0;
  flags.prim_pres_idx           = 1;
  flags.timestep_limiter        = 1;
  flags.update_energy           = 1;
  flags.cooling                 = 1;
  flags.collisional_ionization  = 1;
  flags.photoionization         = 0;
  flags.charge_exchange         = 1;
  flags.cosmicray_ionization    = 1;
  flags.radiative_recombination = 0;


  // nemo(ndim, coord-sys, Nel(state vec), Ntracer, tracers, extra-physics,
  // rad-sources, gamma)
  class nemo MP(flags, tr, rs, 5.0 / 3.0);

  std::vector<double> input(15, 0.0);
  input[0]  = 1.0e-22;  // density
  input[1]  = 1.0e-9;   // pressure (~1e7 K)
  input[10] = 0.72;
  input[11] = 0.28;
  input[12] = 0.36;
  input[13] = 0.05;
  input[14] = 0.1;
  spdlog::info("input  [{:8.3e}", fmt::join(input, ", "));

  std::vector<double> output(15, 0.0);
  double x = 0.0;

  MP.update_single_cell(input, cell_rt, output, 1e10, 5. / 3., &x);

  spdlog::info("output [{:8.3e}]", fmt::join(output, ", "));

  return 0;
}
