// Shows how to send NEMO's CVODE (SUNDIALS) error messages to their own log
// file, and how to print one line per step with the number of CVODE mxstep
// errors.  A hot, dense gas cell cools over 5 long steps.

#include <iostream>
#include <string>
#include <vector>

#include "nemo/nemo.h"  // everything in NEMO is in namespace NEMO

#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/spdlog.h>

int main()
{
  // NEMO sends CVODE error messages to the spdlog logger with this name.
  // Without it, they go to spdlog's default logger (the console).
  auto sundials_log =
      spdlog::basic_logger_mt(NEMO::sundials_logger_name, "sundials.log", true);
  // Alternative: following line prints sundials messages to the default logger
  // auto sundials_log = spdlog::default_logger();

  // State vector: density, pressure, then the tracers named below.
  std::vector<std::string> tracers = {
      "unrelated_mass_scalar", "X_H", "X_He", "H", "He", "He1+"};
  std::vector<NEMO::nemo_star_properties> sources;  // no radiation sources
  std::vector<NEMO::nemo_cell_rt_data> cell_rt;
  const double gamma = 5.0 / 3.0;

  NEMO::nemo_options flags;
  flags.Min_Temperature         = 1.0e1;
  flags.Max_Temperature         = 1.0e9;
  flags.nemo_reltol             = 1.0e-6;
  flags.nemo_abstol             = 1.0e-12;
  flags.nemo_minerg             = 1.0e-20;
  flags.nv_prim                 = 8;
  flags.n_scalar                = 6;
  flags.prim_dens_idx           = 0;
  flags.prim_pres_idx           = 1;
  flags.timestep_limiter        = 1;
  flags.update_energy           = true;
  flags.cooling                 = true;
  flags.collisional_ionization  = true;
  flags.photoionization         = false;
  flags.charge_exchange         = true;
  flags.cosmicray_ionization    = true;
  flags.radiative_recombination = true;
  flags.onthespot_radiation     = false;

  NEMO::nemo solver(flags, tracers, sources, gamma);

  // Hot, fully ionised gas: rho = 1e-18 g/cm^3, p = 1e-5 dyne/cm^2.
  std::vector<double> state = {1.0e-18, 1.0e-5, 1.0, 0.72, 0.28, 0.0, 0.0, 0.0};
  std::vector<double> output(state.size(), 0.0);
  double dt_min = 0.0;

  std::cout << "Initial temperature: " << solver.Temperature(state, gamma)
            << " K\n";
  for (int step = 1; step <= 5; ++step) {
    // Each step is 1e14 s, much longer than the cooling time.
    solver.update_single_cell(state, cell_rt, output, 1.0e14, gamma, &dt_min);
    state = output;
    std::cout << "step " << step << ": T = " << solver.Temperature(state, gamma)
              << " K\n";

    // NEMO recovers from mxstep errors by splitting the step, so it only
    // counts them (over all threads).  Log the count once per step.
    long int n = NEMO::get_and_reset_cvode_mxstep_count();
    if (n > 0) sundials_log->warn("step {}: {} CVODE mxstep errors", step, n);
  }

  // Write out the messages still in the buffer.  (Don't call
  // spdlog::shutdown() while the NEMO object exists: it removes the default
  // logger, which NEMO's destructor uses.)
  sundials_log->flush();
  std::cout << "CVODE messages are in sundials.log\n";
  return 0;
}
