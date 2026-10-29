#include <iostream>
#include <string>
#include <vector>

#include "nemo/nemo.h"

int main()
{
  // State vector: density, pressure, then the tracers named below.
  std::vector<std::string> tracers = {
      "unrelated_mass_scalar", "X_H", "X_He", "H", "He", "He1+"};
  std::vector<nemo_star_properties> sources;  // no radiation sources
  std::vector<nemo_cell_rt_data> cell_rt;
  const double gamma = 5.0 / 3.0;

  nemo_options flags;
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
  flags.update_energy           = false;
  flags.cooling                 = false;
  flags.collisional_ionization  = true;
  flags.photoionization         = false;
  flags.charge_exchange         = true;
  flags.cosmicray_ionization    = true;
  flags.radiative_recombination = true;
  flags.onthespot_radiation     = false;

  nemo solver(flags, tracers, sources, gamma);

  std::vector<double> input = {1.0e-22, 1.0e-9, 1.0,  0.72,
                               0.28,    0.36,   0.05, 0.1};
  std::vector<double> output(input.size(), 0.0);
  double dt_min = 0.0;

  // Evolve the cell for 1e14 seconds.
  solver.update_single_cell(input, cell_rt, output, 1.0e14, gamma, &dt_min);

  std::cout << "Temperature: " << solver.Temperature(output, gamma) << " K\n";
  for (size_t i = 0; i < output.size(); ++i)
    std::cout << "  output[" << i << "] = " << output[i] << "\n";
  return 0;
}
