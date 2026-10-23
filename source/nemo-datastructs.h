#ifndef DATASTRUCTS_H
#define DATASTRUCTS_H

#include <array>
#include <string>
#include <vector>

/*
#define NEMO_MAX_DIM 3


///
/// struct with lots of flags for what extra physics we are using.
///
/// All flags set to 0 (not used) by default, and the initial
/// conditions file can set whichever ones it wants, (or command line
/// overrides).  Some of these must be set if others are set...
/// i.e. they are not all independent.
///
struct which_physics {
  int dynamics;    ///< hydro, not relevant for nemo
  int raytracing;  ///< rt, not relevant for nemo
  int cooling;     ///< include radiative cooling?
  int chemistry;   ///< include chemistry?  YES
  int coll_ionisation;
  int rad_recombination;
  int phot_ionisation;
  int charge_exchange;
  int sat_thermal_cond;    ///< not relevant for nemo
  double tc_strength;      ///< not relevant for nemo
  int update_erg;          ///< update energy equation or not?
  int MP_timestep_limit;   ///< include rxn timescales in dt calc? not relevant
                           ///< for nemo
  double MinTemperature;   ///< set this to low value, in Kelvin
  double MaxTemperature;   ///< set this to high value, in Kelvin
  double H_MassFrac;       ///< not relevant for nemo
  double Helium_MassFrac;  ///< not relevant for nemo
  double Metal_MassFrac;   ///< not relevant for nemo
  int compton_cool;        ///< not relevant for nemo
  int wind_acceleration;   ///< not relevant for nemo
  int shocked_tracer;      ///< not relevant for nemo
  int dust_tracer;         ///< not relevant for nemo
};
*/

///
/// Radiation source struct.
///
struct rad_src_info {
  std::array<double, 3> pos;  ///< src position (cm) (not used by nemo)
  double strength;            ///< src luminosity (ergs/sec)
  double Rstar;               ///< stellar radius in solar radii.
  double Tstar;               ///< stellar effective temperature (K).
  int id;                     ///< src identifying number
  int type;                   ///< src type: == 1
  int update;                 ///< RT_UPDATE_EXPLICIT: == 2
  int at_infinity;            ///< set to 1 if source is at infinity.
  int ongrid;  ///< set to 1 if source is within the domain. (not used by nemo)
  int effect;  ///< 7 = ATLAS OBV, 8 = POWR WC, 9 = POWR WNE
  int NTau;    ///< number of energy bins
  int opacity_src;      ///< for NEMO, hardcoded to 10
  int opacity_var;      ///< (not used by nemo)
  std::string EvoFile;  ///< text file with info from stella evo code (optional)
};

/// List of all radiation sources.
struct rad_sources {
  int Nsources;
  std::vector<struct rad_src_info> sources;
  std::vector<struct rt_source_data> src_data;
};

/// Radiation Source data struct, used for passing info to microphysics classes.
#define MAX_TAU 100
struct rt_source_data {
  double
      Vshell;  ///< Shell volume for discrete photo-ionisation/-heating rates.
  double dS;   ///< Path length through cell.
  double strength[MAX_TAU];  ///< Luminosity (or flux if source at infinity).
  double Column[MAX_TAU];    ///< integral of quantities along LOS to near edge
                             ///< of cell.
  double DelCol[MAX_TAU];    ///< integral of quantities along LOS through cell.
  int id;                    ///< source id.
  int type;                  ///< diffuse-radiation or a real source.
  short unsigned int NTau;  ///< Number of LOS quantities traced for the source.
};

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
struct source {
  int SED_type;       ///< 7 = ATLAS OB, 8 = POWR WC, 9 = POWR WNE, 10 = BB
  double luminosity;  ///< source luminosity in erg/s
  double Rstar;       ///< stellar radius in cm.
  double Tstar;       ///< stellar effective temperature (K).
};

/// Struct containing cell-data that we can use to calculate radiation
/// field in each energy bin.
struct rad_flux_info {
  double Vshell;  ///< Shell volume for discrete photo rates 4/3pi(R+^3-R-^3)
  double dS;      ///< Path length through cell.
  std::vector<double> tau;   ///< optical depth from star to near edge of cell.
  std::vector<double> dtau;  ///< optical depth through cell.
};

struct nemo_rad_source {
  struct source s;         ///< info on the radiation source
  struct rad_flux_info f;  ///< info on RT to/through cell
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
