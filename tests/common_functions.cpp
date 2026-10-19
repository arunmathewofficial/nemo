
#include "common_functions.h"

// *******************************************************************
// Calculate optical depth to the inner edge of every shell.
void spherical_optical_depth(
    nemo &MP,
    const rad_src_info &star,
    const std::vector<std::vector<double> > &state,
    const std::vector<double> &thickness,
    std::vector<std::vector<double> > &optical_depth)
{
  optical_depth.assign(state.size(), std::vector<double>(star.NTau, 0.0));

  std::vector<double> total_tau(star.NTau, 0.0);
  std::vector<double> shell_tau(star.NTau, 0.0);

  // Shells must be ordered from the star outward.
  for (std::size_t i = 0; i < state.size(); ++i) {
    optical_depth[i] = total_tau;  // Optical depth BEFORE this shell

    MP.get_dtau(&star, thickness[i], state[i].data(), shell_tau.data());

    for (int b = 0; b < star.NTau; ++b) {
      total_tau[b] += shell_tau[b];
    }
  }
}
// *******************************************************************
