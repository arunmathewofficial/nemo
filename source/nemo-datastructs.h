#ifndef DATASTRUCTS_H
#define DATASTRUCTS_H

#include <array>
#include <string>
#include <vector>

/// List of options for NEMO setting out how it functions
struct nemo_options {
  double nemo_reltol;      ///< relative tolerance for ODE integration
  double nemo_abstol;      ///< absolute tolerance for ODE integration
  double nemo_minerg;      ///< absolute tolerance on Eint for ODE integration
  double Min_Temperature;  ///< enforce that T remains above this temperature
  double Max_Temperature;  ///< enforce that T remains below this temperature
  int timestep_limiter;    ///< integer flag that determines how to set
                           ///< timescales()
  bool update_energy;      ///< do/don't update the energy equation
  bool cooling;  ///< do/don't include radiative cooling in energy equation
  bool collisional_ionization;  ///< do/don't include coll.ion. reactions
  bool photoionization;         ///< do/don't include photo-ionization reactions
  bool charge_exchange;         ///< do/don't include charge-exchange reactions
  bool cosmicray_ionization;    ///< do/don't include CR ionization reactions
  bool radiative_recombination;  ///< do/don't include radiative recomb.
                                 ///< reactions
  bool onthespot_radiation;      ///< do/don't include recomb.rad. in
                                 ///< new ionizations
};

/// information for setting up radiation source properties
struct nemo_star_properties {
  double Luminosity;  ///< source luminosity in erg/s
  double Teff;        ///< stellar effective temperature (K).
  double Radius;      ///< stellar radius in Rsun
  double photo_emax;  ///< max. photon energy in bins.
  std::vector<double>
      bin_luminosity;  ///< luminosity in each energy bin (erg/s/bin)
  std::vector<std::vector<double> > Ebins;  ///< list of energy bins
  int id;                                   ///< id of source
  int SED_type;          ///< 7 = ATLAS OB, 8 = POWR WC, 9 = POWR WNE, 10 = BB
  int Nbins;             ///< number of energy bins for radiation
  std::string SED_info;  ///< SED info
};


/// Struct containing cell-data that we can use to calculate radiation
/// field in each energy bin.
struct nemo_cell_rt_data {
  double Vshell;  ///< Shell volume for discrete photo rates 4/3pi(R+^3-R-^3)
  double dS;      ///< Path length through cell.
  std::vector<double> tau;   ///< optical depth from star to near edge of cell.
  std::vector<double> dtau;  ///< optical depth through cell.
  int star_id;               ///< id of radiation source in list of sources
};

struct nemo_radiation_source {
  struct nemo_star_properties src;     ///< info on the radiation source
  struct nemo_cell_rt_data cell_data;  ///< info on RT to/through cell
};

/// Struct to track species participating in recombination; here each species
/// is the one which recombines, e.g. H+, He+, He2+,..., but not H0 or He0
struct nemo_recomb_species {
#ifdef NEMO_DEBUG
  std::string s_name;  ///< name of species (only for debugging)
#endif
  int element;   ///< position of parent element in element list
  int s_into;    ///< index in p_local species array that recombinations go into
  int s_from;    ///< index in p_local species array that recombinations go into
  int charge;    ///< charge of species
  bool top_ion;  ///< true if species is highest ionization stage tracked
};


#endif  // DATASTRUCTS_H
