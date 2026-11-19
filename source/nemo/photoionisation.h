/**
 *  \file     photoionisation.h
 *  \author   Arun Mathew
 *  \date     03-07-2023
 **/

#ifndef NEMO_PHOTOIONISATION_H
#define NEMO_PHOTOIONISATION_H
#ifdef SPDLOG_FWD
#include <spdlog/fwd.h>
#endif
#include <spdlog/spdlog.h>
/* prevent clang-format reordering */
#include "nemo-fmt.h"

#include "nemo-constants.h"
#include "solver/interpolate.h"

namespace NEMO {

class photoionisation :
    virtual public interpolate_arrays,
    virtual public nemo_constants {
public:
  photoionisation();
  virtual ~photoionisation();

  /** \brief This function will adjust the number of photo-ionization
   * rate species tables and the photo-ionization heating rate species
   * tables based on the photo tracer list.
   **/
  void setup_photoionisation(
      const std::vector<std::vector<double> > &,  ///< energy bins
      const int,                                  ///< Number of photo species
      const std::vector<std::vector<double> > &,  ///< mean xsection table
      const std::vector<std::vector<double> > &,  ///< bin fraction table
      const std::vector<std::vector<double> >
          &,  ///< luminosity in each energy bin, for each source [src][bin]
      const std::vector<double> &  ///< photoionisation threshold
  );

  /// Set the radiation field of one source at the current cell.  This
  /// stores the attenuated energy flux per unit volume in each energy
  /// bin, L_b * exp(-tau_b) / Vshell, so that exp(-tau_b) is only
  /// calculated once per cell and not for every species and every call
  /// of the rate function.  A source with Vshell <= 0 or ds <= 0 does
  /// not contribute to the rates.
  void set_cell_radiation(
      const int,                   ///< source index
      const double,                ///< Vshell
      const double,                ///< ds, path length through the cell
      const std::vector<double> &  ///< optical depth to the cell in each bin
  );

  /// Photoionisation rate (per ion of the element, per second) and
  /// photoheating rate (erg/cm^3/s) of one species, summed over all
  /// sources and energy bins, for the radiation field set with
  /// set_cell_radiation().
  void calculate_photo_rates(
      const int,     ///< Species identifier
      const double,  ///< Species number density
      const double,  ///< Element number density
      double &,      ///< photoionisation rate
      double &       ///< photoheating rate
  ) const;

private:
  int bin_size;
  std::vector<std::vector<double> > E_bins;
  std::vector<std::vector<double> > mean_energy_table;

  /// This table record the index of the bin where the threshold energy
  /// lies for each photo-species.
  std::vector<int> species_threshold_bin_index;
  std::vector<std::vector<double> > xsection_table;
  std::vector<std::vector<double> > bin_frac_table;
  /// luminosity in each energy bin, for each source [src][bin]
  std::vector<std::vector<double> > bin_spectrum;
  std::vector<double> ionisation_thres;

  /// bin_frac_table * xsection_table: optical depth per unit column
  /// density of the species [species][bin]
  std::vector<std::vector<double> > tau_coeff;
  /// 1 / (mean photon energy) [species][bin]
  std::vector<std::vector<double> > inv_energy;
  /// fraction of the photon energy that heats the gas [species][bin]
  std::vector<std::vector<double> > heat_frac;

  /// attenuated energy flux per unit volume at the current cell,
  /// L_b * exp(-tau_b) / Vshell, [src][bin]
  std::vector<std::vector<double> > cell_flux;
  /// path length through the current cell for each source (0 = no flux)
  std::vector<double> cell_ds;
};

}  // namespace NEMO

#endif  // NEMO_PHOTOIONISATION_H
