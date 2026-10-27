
#include "common_functions.h"

// *******************************************************************
// Calculate optical depth to the inner edge of every shell.
void spherical_optical_depth(
    nemo &MP,
    const nemo_star_properties &star,
    const std::vector<std::vector<double> > &state,
    std::vector<std::vector<struct nemo_cell_rt_data> >
        &cell_rt  // vec of cell data
)
{
  int NTau = star.Nbins;
  // HARDCODED FOR JUST ONE RADIATION SOURCE
  // Each grid cell has a cell_rt vector of length 1 (for 1 source)
  for (auto el : cell_rt) {
    std::fill(el[0].tau.begin(), el[0].tau.end(), 0.0);
    std::fill(el[0].dtau.begin(), el[0].dtau.end(), 0.0);
  }
  std::vector<double> total_tau(NTau, 0.0);
  std::vector<double> shell_tau(NTau, 0.0);

  // Shells must be ordered from the star outward.
  for (std::size_t i = 0; i < state.size(); ++i) {
    cell_rt[i][0].tau = total_tau;  // Optical depth BEFORE this shell

    MP.get_dtau(state[i], cell_rt[i]);

    for (int b = 0; b < NTau; ++b) {
      total_tau[b] += cell_rt[i][0].dtau[b];
    }
  }
}
// *******************************************************************
