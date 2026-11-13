/// common_functions.cpp
///
/// purpose: define some commone functions that are used by multiple tests

#include "common_functions.h"

using namespace NEMO;

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
  // Every grid cell has a cell_rt vector with one entry per radiation
  // source.  All sources are at the centre, and use the same energy bins.
  const int NTau = star.Nbins;
  if (cell_rt.empty()) return;
  const std::size_t n_src = cell_rt[0].size();
  std::vector<std::vector<double> > total_tau(
      n_src, std::vector<double>(NTau, 0.0));

  // Shells must be ordered from the star outward.
  for (std::size_t i = 0; i < state.size(); ++i) {
    // Optical depth BEFORE this shell
    for (std::size_t src = 0; src < n_src; ++src)
      cell_rt[i][src].tau = total_tau[src];

    MP.get_dtau(state[i], cell_rt[i]);

    for (std::size_t src = 0; src < n_src; ++src) {
      for (int b = 0; b < NTau; ++b)
        total_tau[src][b] += cell_rt[i][src].dtau[b];
    }
  }
}
// *******************************************************************
