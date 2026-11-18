/**
 *  \file     photoionisation.cpp
 *  \author   Arun Mathew
 *  \date     03-07-2023
 **/
#include "photoionisation.h"

namespace NEMO {

photoionisation::photoionisation() {}
photoionisation::~photoionisation() {}

// SETUP PHOTO-IONISATION ################################################
void photoionisation::setup_photoionisation(
    const std::vector<std::vector<double> > &ebins,  // energy bins
    const int N_photo_species,                       // number of photo species
    const std::vector<std::vector<double> > &photo_xsection_table,  //
    const std::vector<std::vector<double> >
        &bin_weight_table,  // bin weight table
    const std::vector<std::vector<double> >
        &bin_luminosity,  // stellar spectrum of each source [src][bin]
    const std::vector<double> &photo_tracer_Eth  // photoionisation threshold
)
{

  spdlog::debug("Setting up photoionisation module");

  // replicating energy bin
  E_bins = ebins;
  // get the bin size
  bin_size = ebins.size();
  // make a replica of the following
  xsection_table   = photo_xsection_table;
  bin_frac_table   = bin_weight_table;
  bin_spectrum     = bin_luminosity;
  ionisation_thres = photo_tracer_Eth;

  // 1. Making the species threshold bin index table
  // This table contains the index of the bin where the threshold energy
  // lies in the energy bin vector for each photo-species.
  species_threshold_bin_index.resize(N_photo_species, -1);

  int b;
  for (int s = 0; s < N_photo_species; s++) {
    for (b = 0; b < bin_size; b++) {
      if (bin_frac_table[s][b] > 0.0) break;
    }
    species_threshold_bin_index[s] = b;
  }
  spdlog::debug(
      "photoionisation: species_threshold_bin_pos {}",
      species_threshold_bin_index);

  // Verifying the correctness of the bin indexing.
  auto it = std::find(
      species_threshold_bin_index.begin(), species_threshold_bin_index.end(),
      -1);
  if (it != species_threshold_bin_index.end()) {
    // If any entry in the species bin index table is zero, the
    // indexing has failed. It is necessary to reset the energy
    // bin properly.
    spdlog::error("photoionisation: species bin indexing failed.");
    throw nemo_error("photoionisation: species bin indexing failed");
  }
  else {
    spdlog::debug("photoionisation: species bin indexing successful");
  }

  // 2. Making the mean energy table
  mean_energy_table.resize(N_photo_species);
  double mean_bin_energy;
  for (int s = 0; s < N_photo_species; s++) {
    for (int b = 0; b < bin_size; b++) {

      // if the bin weight is zero then skip the iteration with zero mean value
      if (bin_frac_table[s][b] == 0.0) {
        mean_energy_table[s].push_back(0.0);
        continue;
      }
      // If the bin weight is <1, the ionisation threshold of this species
      // lies inside the bin, so the mean energy is taken between the
      // threshold and the upper edge of the bin; otherwise over the full
      // bin.  (E_bins is shared by all species and must not be changed.)
      const double e_min =
          (bin_frac_table[s][b] < 1.0) ? ionisation_thres[s] : E_bins[b][0];
      mean_bin_energy = 0.5 * (e_min + E_bins[b][1]);
      // append the mean bin energy to the table for the species
      mean_energy_table[s].push_back(mean_bin_energy * eV());
    }  // bin loop
  }  // species loop

  // 3. Tables used by the rate calculation: optical depth per unit column
  // density, 1/(mean photon energy), and the fraction of the photon energy
  // that goes into heating (the excess over the ionisation threshold).
  // They are zero in bins below the threshold of the species.
  tau_coeff.assign(N_photo_species, std::vector<double>(bin_size, 0.0));
  inv_energy.assign(N_photo_species, std::vector<double>(bin_size, 0.0));
  heat_frac.assign(N_photo_species, std::vector<double>(bin_size, 0.0));
  for (int s = 0; s < N_photo_species; s++) {
    for (int b = species_threshold_bin_index[s]; b < bin_size; b++) {
      tau_coeff[s][b]  = bin_frac_table[s][b] * xsection_table[s][b];
      inv_energy[s][b] = 1.0 / mean_energy_table[s][b];
      heat_frac[s][b] =
          1.0 - ionisation_thres[s] * eV() / mean_energy_table[s][b];
    }
  }

  // 4. Radiation field at the current cell, for each source (no flux
  // until set_cell_radiation() is called).
  cell_flux.assign(bin_spectrum.size(), std::vector<double>(bin_size, 0.0));
  cell_ds.assign(bin_spectrum.size(), 0.0);
}

// SET RADIATION FIELD AT THE CURRENT CELL ###############################
void photoionisation::set_cell_radiation(
    const int src,                       // source index
    const double vshell,                 // vshell
    const double ds,                     // ds
    const std::vector<double> &tau_bins  // Optical depth in each energy bin
)
{
  if (vshell <= 0.0 || ds <= 0.0) {
    cell_ds[src] = 0.0;
    return;
  }
  cell_ds[src] = ds;
  for (int b = 0; b < bin_size; b++)
    cell_flux[src][b] = bin_spectrum[src][b] * exp(-tau_bins[b]) / vshell;
}

// CALCULATING PHOTO-IONISATION AND PHOTO-HEATING RATES ##################
void photoionisation::calculate_photo_rates(
    const int species,                 // Species_identifier
    const double species_num_density,  // Species number density,
    const double element_num_density,  // Element number density,
    double &ion_rate,                  // photoionisation rate
    double &heat_rate                  // photoheating rate
) const
{
  // For each source and each energy bin b at or above the threshold of
  // the species, the energy absorbed by the species per unit volume and
  // time is
  //   A_b = L_b * exp(-tau_b) * (1 - exp(-dtau_b)) / Vshell
  // where tau_b is the optical depth from the source to the cell and
  // dtau_b = n_s * sigma_b * ds is the optical depth of this species
  // through the cell (photon-conserving form).  Then
  //   photoionisation rate = sum_b A_b / E_b / n_element
  //   photoheating rate    = sum_b A_b * (1 - E_threshold / E_b)
  // where E_b is the mean photon energy in the bin.
  ion_rate                      = 0.0;
  heat_rate                     = 0.0;
  const std::vector<double> &tc = tau_coeff[species];
  const std::vector<double> &ie = inv_energy[species];
  const std::vector<double> &hf = heat_frac[species];
  const int b0                  = species_threshold_bin_index[species];
  for (size_t src = 0; src < cell_ds.size(); src++) {
    if (cell_ds[src] <= 0.0) continue;
    const double column             = species_num_density * cell_ds[src];
    const std::vector<double> &flux = cell_flux[src];
    for (int b = b0; b < bin_size; b++) {
      const double dtau = column * tc[b];
      // if dtau << 1 then 1 - exp(-dtau) = dtau
      const double absorbed =
          flux[b] * ((dtau < 1.0E-02) ? dtau : (1.0 - exp(-dtau)));
      ion_rate += absorbed * ie[b];
      heat_rate += absorbed * hf[b];
    }
  }
  ion_rate /= element_num_density;
}

}  // namespace NEMO
