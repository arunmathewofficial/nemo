/**
 *  \file     stellar_atmosphere_models.h
 *  \author   Arun Mathew
 *  \date     05-07-2023
 **/

#ifndef PION_STELLAR_ATMOSPHERE_MODELS_H
#define PION_STELLAR_ATMOSPHERE_MODELS_H

// logging
#ifdef SPDLOG_FWD
#include <spdlog/fwd.h>
#endif
#include <spdlog/spdlog.h>
/* prevent clang-format reordering */
#include <fmt/ranges.h>

// std lib
#include <fstream>
#include <string>
#include <sys/stat.h>
#include <variant>

// nemo
#include "nemo-constants.h"
#include "nemo-datastructs.h"

struct stellar_model_datatype {
  std::string info;
  std::vector<std::vector<double> > data;
};

class stellar_atmosphere_models : virtual public nemo_constants {

public:
  stellar_atmosphere_models();
  virtual ~stellar_atmosphere_models();

  /// define energy bins for stellar radiation
  void setup_stellar_atmosphere_models(
      struct nemo_star_properties &  ///< struct to write energy binning to
  );

  /** \brief This function calculate the luminosity number in different
   * energy bin for the given radiation source.
   **/
  int generate_bin_luminosity(
      struct nemo_star_properties
          &  ///< struct with all info on the radiation source
  );

  /// @brief  get number of energy bins used in the atmosphere models
  int get_N_Ebins() { return N_energy_bins; };

private:
  std::vector<std::vector<double> > energy_bins;
  int N_energy_bins;

  void stellar_atmosphere_model_database();

  // get stellar atmosphere model
  stellar_model_datatype get_stellar_atmosphere_model(int);

  // stellar atmosphere models;
  stellar_model_datatype atlas_ckp00_g45;  ///< O3V STAR
  stellar_model_datatype powr_mw_wc;       ///< WC
  stellar_model_datatype powr_mw_wne;      ///< WNE
  stellar_model_datatype blackbody;        ///< BLACKBODY
};

#endif  // PION_STELLAR_ATMOSPHERE_MODELS_H
