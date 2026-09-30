#ifndef DATASTRUCTS_H
#define DATASTRUCTS_H

#include <array>
#include <string>

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
  int dynamics; ///< hydro, not relevant for nemo
  int raytracing; ///< rt, not relevant for nemo
  int cooling;    ///< include radiative cooling?
  int chemistry; ///< include chemistry?  YES
  int coll_ionisation;
  int rad_recombination;
  int phot_ionisation;
  int charge_exchange;
  int sat_thermal_cond;   ///< not relevant for nemo
  double tc_strength;     ///< not relevant for nemo
  int update_erg;         ///< update energy equation or not?
  int MP_timestep_limit;  ///< include rxn timescales in dt calc? not relevant for nemo
  double MinTemperature;  ///< set this to low value, in Kelvin
  double MaxTemperature;  ///< set this to high value, in Kelvin
  double H_MassFrac;      ///< not relevant for nemo
  double Helium_MassFrac; ///< not relevant for nemo
  double Metal_MassFrac;  ///< not relevant for nemo
  int compton_cool;       ///< not relevant for nemo
  int wind_acceleration;  ///< not relevant for nemo
  int shocked_tracer;     ///< not relevant for nemo
  int dust_tracer;        ///< not relevant for nemo
};

///
/// Radiation source struct.
///
struct rad_src_info {
  std::array<double, NEMO_MAX_DIM> pos;  ///< src position (cm) (not used by nemo)
  double strength;  ///< src luminosity (ergs/sec)
  double Rstar;     ///< stellar radius in solar radii.
  double Tstar;     ///< stellar effective temperature (K).
  int id;           ///< src identifying number
  int type;         ///< src type: == 1
  int update;       ///< RT_UPDATE_EXPLICIT: == 2
  int at_infinity;  ///< set to 1 if source is at infinity.
  int ongrid;       ///< set to 1 if source is within the domain. (not used by nemo)
  int effect;       ///< 7 = ATLAS OBV, 8 = POWR WC, 9 = POWR WNE
  int NTau;         ///< number of energy bins
  int opacity_src;  ///< for NEMO, hardcoded to 10
  int opacity_var;  ///< (not used by nemo)
  std::string EvoFile; ///< text file with info from stella evo code (optional)
};


/// List of all radiation sources.
struct rad_sources {
  int Nsources;
  std::vector<struct rad_src_info> sources;
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



#endif // DATASTRUCTS_H
