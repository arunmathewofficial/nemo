///
/// \file NEMO.h
/// \authors Arun Mathew, Jonathan Mackey, Maggie Celeste Goulden
/// \date 2018.10
///
/// modifications:
/// - 2018.10 .

#ifndef NEMO_H
#define NEMO_H

/// Description:
/// A multi-ion chemistry solver with non-equilibrium ionization and
/// recombination of different species, and with elemental abundances
/// that are input for each solve (i.e. can vary within simulation).
///
/// The integration method uses the CVODES solver from the SUNDIALS
/// package by (Cohen, S. D., & Hindmarsh, A. C. 1996, Computers in
/// Physics, 10, 138) available from
///  https://computation.llnl.gov/casc/sundials/main.html
/// The method is backwards differencing (i.e. implicit) with Newton
/// iteration.
/// * This code is not yet working and should not be used.


/// Choose the cooling function
// #define MELLEMA
#define CHIANTI

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#ifdef SPDLOG_FWD
#include <spdlog/fwd.h>
#endif
#include <spdlog/spdlog.h>
/* prevent clang-format reordering */
#include "atomic_physics_database.h"
#include "charge_exchange.h"
#include "coll_ionise_recomb.h"
#include "constants.h"
#include "microphysics/cvode_integrator.h"
#include "microphysics/microphysics_base.h"
#include "sim_params.h"
#ifdef MELLEMA
#include "mellema_cooling.h"
#elif defined CHIANTI
#include "chianti_cooling.h"
#endif
#include "multi_ion_photo_xsections.h"
#include "photoionisation.h"
#include "stellar_atmosphere_models.h"
#include <algorithm>
#include <fstream>
#include <set>
#include <sys/stat.h>
#include <utility>

#define NEMO_RELTOL 1.0e-6   ///< relative-error tolerance
#define NEMO_ABSTOL 1.0e-12  ///< minimum neutral fraction I care.
#define NEMO_MINERG 1.0e-20  ///< Minimum internal energy density I care.

// ############################################################################
//  NEMO CLASS VARIABLES AND FUNCTIONS
///
/// Integrator for the microphysics equations for the non-equilibrium ion
/// fraction of various ions, and for the internal energy density.
///
class nemo :
    public microphysics_base,
    public cvode_solver,
    public coll_ionise_recomb,
    public charge_exchange,
#ifdef MELLEMA
    public mellema_cooling,
#elif defined CHIANTI
    public chianti_cooling,
#endif
    public stellar_atmosphere_models,
    public multi_ion_photo_xsections,
    public photoionisation {
public:
  // ****** CONSTRUCTOR ******
  nemo(
      const int,               ///< grid dimensions
      const int,               ///< Coordinate System flag
      const int,               ///< Total number of variables in state vector
      const int,               ///< Number of tracer variables in state vector.
      const std::string *,     ///< List of what the tracer variables mean.
      struct which_physics *,  ///< pointer to extra physics flags.
      struct rad_sources *,    ///< radiation sources.
      const double             ///< EOS Gamma
  );

  // ****** DESTRUCTOR ******
  ~nemo();

  // ****** CATALOGER ******


  std::vector<string> tracers_replica;  ///< microphysics tracers replica.
  int N_prim;              ///< Size of the primitive vector (prim_vec).
  int first_tracer_index;  /// Idx of the first tracer variable in the prim_vec.
  std::vector<vector<string> > all_species_list;  ///< Imported species list.
  std::vector<int>
      massfrac_idx;  ///< idx of elemental mass fraction in prim_vec.
  std::vector<std::string> element_tracers;  ///< Element tracers in prim_vec.
  std::unordered_map<std::string, int> element_position_index_map;
  std::vector<int> element_atomic_number;   ///< atomic number of elements.
  std::vector<double> element_atomic_mass;  ///< atomic mass of elements.
  int N_elements;  ///< Number of elements in the primitive-vector/tracer.
  std::vector<double> element_num_density;  ///< element number density vector.

  std::vector<string> local_tracers;  ///< Extracted local tracers.
  int N_species;                      ///< Number of species.
  std::vector<int> species_prim_idx;  ///< local species index in prim_vec.
  std::vector<int>
      species_local_idx;  ///< species idx in local vector (local_vec).
  std::vector<int>
      element_species_count;  ///< number of species belong to same elements.
  std::vector<int> species_element;  ///< species' element positions.
  std::unordered_map<std::string, int> species_prim_index_map;


  /// Cosmic-ray ionisation ===================================================
  std::vector<int>
      neutral_species_local_idx;  ///< neutral species idx in local_vec.
  double cosmicray_ionise_rate;
  double cosmic_ray_energy;

  /// Collisional ionisation ===================================================

  std::vector<string> collion_species;  ///< collisional ionisation species.
  std::vector<int> collion_species_charge;
  std::vector<int>
      upper_ion_local_idx;  /// local index of upper ionisation state if exist.
  std::vector<double>
      collion_species_IE;  /// coll-ionisation species ionisation energy (egs)


  /// Recombination ============================================================

  std::vector<string> recomb_species;  ///< recombination species.
  std::vector<int>
      recomb_species_element;  ///< recomb species element position.
  std::vector<int> recomb_species_charge;
  std::vector<int> recomb_top_ion_bool;

  /// Charge exchange reaction =================================================
  bool charge_exchange_enable;
  bool ionReact = true;
  bool recReact = false;
  std::vector<string> charge_exchange_ionise_species;
  std::vector<string> charge_exchange_recomb_species;
  std::vector<std::vector<std::pair<std::string, std::string> > >
      ionise_cx_reactions;
  std::vector<std::vector<std::pair<std::string, std::string> > >
      recomb_cx_reactions;
  std::vector<bool> has_ionise_cx_exist;
  std::vector<bool> has_recomb_cx_exist;
  std::vector<std::vector<int> > ionise_cx_second_reactant_element;
  std::vector<std::vector<int> > recomb_cx_second_reactant_element;
  std::vector<std::vector<int> > ionise_cx_second_reactant_charge;
  std::vector<std::vector<int> > recomb_cx_second_reactant_charge;

  std::vector<std::vector<double> >
      ionise_cx_energy;  ///< Charge exchange ionising reaction energy defect
  std::vector<std::vector<double> >
      recomb_cx_energy;  ///< Charge exchange recombining reaction energy defect


  /// Photo-ionisation =========================================================

  std::vector<string> photo_species;  ///< photoionisation species.
  int N_photo_species;  ///< Number of photo species in the tracer list.
  std::vector<int>
      photo_species_local_idx;  /// photo-species local vector index.
  std::vector<int>
      photo_upper_ion_local_idx;           /// photo-species upper local index.
  std::vector<int> photo_species_element;  /// Corresponding element's position.
  std::vector<int> photo_species_charge;   /// photo species charge.
  std::vector<int> photo_species_prim_idx;  /// photo-species prim-vec index.
  std::vector<double>
      photo_species_IE;  ///< Record the ionisation threshold of photo species

  std::vector<std::vector<double> > ebins;
  int N_ebins;
  double photo_Emax;

  std::vector<double>
      bin_luminosity;  ///< luminosity in each energy bin (ergs/s/bins)
  std::string SED_info;
  struct star_properties_datatype {
    int Type;              // identifier
    double Teff;           // Effective temperature of the star
    double Radius;         // Stellar radius
    double Luminosity;     // Total luminosity
    std::string SED_info;  // SED info
  };

  star_properties_datatype stellar_properties;


  /// Cooling ==================================================================
  std::vector<vector<string> > cooling_species;  ///< cooling species.
  std::ofstream cooling_outfile;                 /// CIE cooling function


  /// variable recording species number density for NEMO tracer list
  std::vector<std::vector<double> > species_number_density;



  int N_local;
  int E_index;



  /// Uses this steps to generate lookup tables
  double delta_logT;
  double delta_logne;
  std::vector<double> T_table;   ///< Temperature table
  std::vector<double> ne_table;  ///< electron number density table

  const double T_min;
  const double T_max;
  const int Num_temps;  // NB this needs to be const so can initialise arrays
                        // with it. If this is >=1e4, get stack overflow errors.

  const int Num_ne;
  const double ne_min;
  const double ne_max;



  // ****** GET RATES ******
  /// Function to return collisional ionisation rate of the given species for
  /// a given temperature.
  double collisional_ionisation_rate(
      const int,    // Species identifier
      const double  // Temperature
  );
  /// Function to return recombination rate.
  double recombination_rate(
      const int,     ///< Species identifier
      const double,  ///< Temperature
      double &       ///< Case A rate for H+,He+
  );

  /// Function to return ionising reaction rate.
  double ionise_charge_exchange_rate(
      const int,    // Elemental position of the first reactant
      const int,    // Second reactant identifier
      const double  // Temperature
  );

  /// Function to return recombining reaction rate.
  double recomb_charge_exchange_rate(
      const int,    // Elemental position of the first reactant
      const int,    // Second reactant identifier
      const double  // Temperature
  );

#ifdef MELLEMA
  /// Fucntion to return mellema colling rate.
  virtual double mellema_cooling_rate(
      const int,     //  species identifier
      const double,  // Temperature
      const double   // ne
  );
#endif

#ifdef CHIANTI
  /// Fucntion to return chianti colling rate.
  virtual double chianti_cooling_rate(
      const int,     //  species identifier
      const double,  // Temperature
      const double   // ne
  );
#endif


protected:
  /// Function will generate all the indexing necessary for NEMO module.
  void tracer_species_cataloger(
      const int,  ///< Number of variables in the primitive vector.
      const int,  ///< Number of tracer variables in the primitive vector.
      const std::string *  ///< List of tracer variables.
  );

  // END OF NEMO CLASS VARIABLES AND FUNCTIONS
  // ***************************************************************************



  // ############################################################################
  //  FUNCTIONS DERIVED FROM BASE CLASS FOLLOW
public:
  N_Vector y_in,    ///< current y-vector
      y_out;        ///< output y-vector
  int N_equations;  ///< number of equations i.e. number of species in Y.
  int N_extradata;  ///< number of elements in user-data array.

  ///
  /// calculate dy/dt for the vector of y-values.
  ///
  virtual int ydot(
      double,  ///< current time (probably not needed for rate equations)
      const N_Vector,  ///< current Y-value
      N_Vector,        ///< vector for Y-dot values
      const double *  ///< extra user-data vector, P, for evaluating ydot(y,t,p)
  );

  ///
  /// Calculate the Jacobian matrix d(dy_i/dt)/dy_j for a vector of y-values.
  ///
  virtual int Jacobian(
#if defined CVODE2
      int,  ///< N (not sure what this is for!)
#endif
      double,          ///< time, t
      const N_Vector,  ///< current Y-value
      const N_Vector,  ///< vector for Y-dot values
      const double
          *,    ///< extra user-data vector, P, for evaluating ydot(y,t,p)
      CVMatrix  ///< Jacobian matrix
  )
  {
    spdlog::info("Jacobian not implemented in NEMO!");
    return 1;
  }

  ///
  /// Get the number of extra parameters and the number of equations.
  ///
  virtual void get_problem_size(
      int *,  ///< number of equations
      int *   ///< number of parameters in user_data vector.
  );


protected:
  ///
  /// set the relative and absolute error tolerances
  ///
  virtual void get_error_tolerances(
      double *,  ///< relative error tolerance (single value)
      double[]   ///< absolute error tolerance (array)
  );
  // END OF FUNCTIONS DERIVED FROM BASE CLASS
  // ***************************************************************************



  // ############################################################################
  // ############################################################################


  ///
  /// Function for updating vectors according to species found in tracer list.
  ///
  // species_tracer_initialise

  ///
  /// A primitive vector is input, and lots of optional extra data in
  /// the last argument, and the ion fraction(s) and internal energy
  /// are updated by the requested timestep.  Results are saved to
  /// the destination vector, which can be the same pointer as the
  /// input vector.
  ///
  int TimeUpdateMP(
      const pion_flt *,  ///< Primitive Vector to be updated.
      pion_flt *,        ///< Destination Vector for updated values.
      const double,      ///< Time Step to advance by.
      const double,      ///< EOS gamma.
      const int,         ///< Switch for what type of integration to use.
      double *           ///< Vector of extra data (column densities, etc.).
  );

  ///
  ///
  ///
  int TimeUpdate_RTsinglesrc(
      const pion_flt *,  ///< Primitive Vector to be updated.
      pion_flt *,        ///< Destination Vector for updated values.
      const double,      ///< Time Step to advance by.
      const double,      ///< EOS gamma.
      const int,         ///< Switch for what type of integration to use.
      const double,      ///< flux in per unit length along ray (F/ds or L/dV)
      const double,      ///< path length ds through cell.
      const double,      ///< Optical depth to entry point of ray into cell.
      double *  ///< return optical depth through cell in this variable.
  )
  {
    spdlog::info("NEMO::TimeUpdate_RTsinglesrc() not implemented");
    return 1;
  }

  ///
  /// Update microphysics including photoionization and photoheating
  ///
  virtual int TimeUpdateMP_RTnew(
      const pion_flt *,  ///< Primitive Vector to be updated.
      const int,         ///< unused
      const std::vector<struct rt_source_data> &,  ///< unused
      const int,                                   ///< unused
      std::vector<struct rt_source_data> &,
      ///< list of ionising src column densities and source properties.
      pion_flt *,    ///< Destination Vector for updated values
      const double,  ///< Time Step to advance by.
      const double,  ///< EOS gamma.
      const int,     ///< unused
      double *       ///< any returned data (final temperature?).
  );

  /*
  ///
  /// Set the properties of a multifrequency ionising radiation source.
  ///
  int set_multifreq_source_properties(
      // int tempin,
      // double *
      const struct rad_src_info *,  ///< source data
      double *  ///< O/P source strength in different energy bins.
  );
   */

  ///
  /// Returns the gas temperature.  This is only needed for data output, so
  /// there is no need to make it highly optimized.
  /// - Is threadsafe.
  ///
  double Temperature(
      const pion_flt *,  ///< primitive vector
      const double       ///< eos gamma
  );

  ///
  /// Set the gas temperature to a specified value.
  /// Needed if you want this feature to set initial conditions.
  ///
  int Set_Temp(
      pion_flt *,    ///< primitive vector.
      const double,  ///< temperature
      const double   ///< eos gamma.
  );

  ///
  /// This returns the minimum timescale of the times flagged in the
  /// arguments.
  ///
  virtual double timescales(
      const pion_flt *,  ///< Current cell.
      const double,      ///< EOS gamma.
      const bool,        ///< set to 'true' if including cooling time.
      const bool,        ///< set to 'true' if including recombination time.
      const bool         ///< set to 'true' if including photo-ionsation time.
  );

  ///
  /// This returns the minimum timescale of all microphysical
  /// processes, including reaction times for each species and the
  /// total heating/cooling time for the gas.
  /// It requires the radiation field as an input, but this is not
  /// used for now so the vectors can be null.
  ///
  virtual double timescales_RT(
      const pion_flt *,  ///< Current cell.
      const int,         ///< Number of UV heating sources.
      const std::vector<struct rt_source_data> &,
      ///< list of UV-heating column densities and source properties.
      const int,  ///< number of ionising radiation sources.
      const std::vector<struct rt_source_data> &,
      ///< list of ionising src column densities and source properties.
      const double  ///< EOS gamma.
  );

  ///
  /// Initialise microphysics ionisation fractions to an equilibrium
  /// value. This is optionally used in the initial condition
  /// generator.  Not implemented here.
  ///
  int Init_ionfractions(
      pion_flt *,    ///< Primitive vector to be updated.
      const double,  ///< eos gamma.
      const double   ///< optional gas temperature to end up at.
  )
  {
    spdlog::info("NEMO::Init_ionfractions() not implemented");
    return 1;
  }

  ///
  /// Return index of tracer for a given string.
  ///
  // int Tr(const string  ///< name of tracer we are looking for.
  // );

  ///
  /// Get the total cooling rate.  This is for postprocessing the
  /// simulation data only -- IT IS NOT OPTIMISED FOR SPEED.
  ///
  virtual double total_cooling_rate(
      const pion_flt *,  ///< Current cell values.
      const int,         ///< Number of UV heating sources.
      const std::vector<struct rt_source_data> &,
      ///< list of UV-heating column densities and source properties.
      const int,  ///< number of ionising radiation sources.
      const std::vector<struct rt_source_data> &,
      ///< list of ionising src column densities and source properties.
      const double  ///< EOS gamma.
  )
  {
    spdlog::info("NEMO::total_cooling_rate() not implemented");
    return 1;
  }

  ///
  /// Get the total recombination rate for an ion, given the input
  /// state vector.
  ///
  virtual double get_recombination_rate(
      const int,         ///< ion index in tracer array (optional).
      const pion_flt *,  ///< input state vector (primitive).
      const double       ///< EOS gamma (optional)
  )
  {
    spdlog::info("NEMO::get_recombination_rate() not implemented");
    return 1;
  }

  ///
  /// Return the H mass fraction
  ///
  virtual double get_X_H();

  ///
  /// Get electron number density (cm^{-3})
  ///
  double get_n_elec(const pion_flt *  ///< primitive state vector.
  );

  ///
  /// Get electron number density (cm^{-3})
  ///
  double get_n_ion(
      std::string,      ///< ion name
      const pion_flt *  ///< primitive state vector.
  );

  ///
  /// Calculates the change in optical depth over the cell
  /// for use in the Raytracing module
  ///
  virtual void get_dtau(
      const rad_source *,  ///< pointer to radiation source struct
      const pion_flt,      ///< ds, thickness of the cell
      const pion_flt
          *,      ///< input primitive vector from grid cell (length nv_prim)
      pion_flt *  ///< output dtau vector
  );


  // ###########################################################################
  // ###########################################################################

protected:
  ///
  /// convert state vector from grid cell into local microphysics vector.
  ///
  virtual int convert_prim2local(
      const pion_flt *,       ///< primitive vector (nv_prim)
      std::vector<double> &,  ///< local vector [x(H0),E](n+1).
      int                     /// flags the function
  );

  ///
  /// Convert local microphysics vector into state vector for grid cell.
  /// This is the inverse of convert_prim2local.
  ///
  virtual int convert_local2prim(
      const std::vector<double> &,  ///< local (updated) vector [x(H0),E](n+1).
      const pion_flt *,             ///< input primitive vector (nv_prim)
      pion_flt *,                   ///< updated primitive vector (nv_prim)
      int                           /// flags the function
  );

  ///
  /// Set local variables (y_ion_frac, y_neutral_frac, ne) by calling
  /// get_y_ion_frac(), get_y_top_ion_frac, and get_ne().
  ///
  virtual void set_localvariables(
      const std::vector<double> &,  ///< local vector
      std::vector<double> &,        ///< y_ion_frac
      std::vector<double> &,        ///< y_neutral_frac
      double &                      ///< ne
  );

  /// Set y_ion_frac
  virtual void get_y_species_frac(
      const std::vector<double> &,  ///< local vector
      std::vector<double> &         ///< y_ion_frac
  );


  /// Set y_top_ion_frac
  virtual void get_y_top_ion_frac(
      const std::vector<double> &,  ///< y_species_frac
      std::vector<double> &         ///< y_top_ion_frac
  );


  /// Calculate electron number density
  virtual void get_ne(
      const std::vector<double> &,  ///< y_species_frac
      double &                      ///< electron number density
  );

  /// Correct local vector
  /// + check that ion fractions are >= ABSTOL.
  /// + check that energy is within temperature limits
  virtual int correct_localvector(
      std::vector<double> &,  /// local vector
      int                     /// flags the function
  );


  // ###########################################################################
  // ###########################################################################



  ///
  /// returns gas temperature according to E=nkT/(g-1) with n=1.1*nH*(1+x_in),
  /// appropriate for a gas with 10% Helium by number, and if He is singly
  /// ionised whenever H is.
  ///
  virtual double get_temperature(
      const std::vector<double> &,  ///< ion fractions
      const std::vector<double> &,  ///< element number densities
      const double                  ///< E_int (per unit volume)
  );

  ///
  /// Returns total number of particles.
  ///
  virtual double get_ntot(
      const std::vector<double> &,  ///< species fractions (input)
      const std::vector<double> &   ///< element number density (input)
  );

  ///
  /// Set the size of local vectors, and index them with integers.
  ///
  virtual void setup_local_vectors();


  ///
  /// This function is utilized to correct the input primitive vector of the
  /// grid. We utilize the sCMA (Simple Consistent Multi-fluid Advection, Plewa
  /// and Muller, Astron. Astrophys. 342, 179, 1999) scheme. The utility of this
  /// method requires the inputted mass fraction to satisfy the constraint that
  /// the sum of elemental mass fractions must equal unity.
  ///
  virtual void sCMA(
      pion_flt *,       ///< corrected vector.
      const pion_flt *  ///< input primitive vector.
  );

  const int ndim;          ///< Number of dimensions in grid.
  const double eos_gamma;  ///< EOS gamma for ideal gas.
  const int coord_sys;     ///< Coordinate System flag
  double gamma_minus_one;  ///< as named.
  double Min_NeutralFrac;  ///< minimum H0 fraction allowed (eps=1.0e-12)
  double Max_NeutralFrac;  ///< Maximum H0 fraction allowed (1-eps)


  ///***********************************************

  ///	Vector to store initial corrector values (used to modify flux according to
  /// sCMA)
  std::vector<pion_flt> prim_copy;


  int print_flag;


  // std::vector<double> elem_number_density;



  /// ===========================================================================
  ///  Photo-ionisation
  /// ===========================================================================


protected:
  /// Function will generate indexing necessary for photo-ionisation tracers
  /// void photo_species_cataloger();



  /// For each cell, ydot() needs to know the local radiation field.  This
  /// function takes all of the input radiation sources and calculates what
  /// ydot() needs to know for both heating and ionisation sources.
  ///
  void setup_radiation_source_parameters(
      // const pion_flt *,  ///< primitive input state vector.
      // vector<double> &,  ///< local input state vector (x_in,E_int)
      const int,  ///< number of ionising radiation sources.
      const std::vector<struct rt_source_data>
          &  ///< list of ionising src column densities and source properties.
  );

  double ds     = 0.0;           ///< Path length through the current cell
  double vshell = 0.0;           ///< Shell volume in the current cell
  std::vector<double> tau_bins;  ///< Optical depth in each energy bin



  /// Function to return photo-ionisation rate.
  double photoionisation_rate(
      const int,    // Species identifier
      const double  // Species number density
  );

  /// Function to return photo-heating rate.
  double photoheating_rate(
      const int,    // Species identifier
      const double  // Species number density
  );


  // int N_ion_srcs;    ///< No ionising sources --> 0, otherwise --> 1
  // int ion_src_type;  ///< Either RT_EFFECT_MFION or RT_EFFECT_PION_MONO.


  // Constant data in the cell, received from cell data.
  double mpv_nH;  ///< total hydrogen number density at current cell.
  int N_rsrc;     ///< number of radiation sources to include.
  std::vector<struct rt_source_data> rt_data;  ///< data on radiation sources.


  /// index matching photo_ion xsection to ion in local vector
  std::vector<std::vector<double> > y_ion_xsections;



  // #########################################################################
  // *********************** DEBUG  FUNCTIONS ********************************
  // #########################################################################

  /// This will print collisional ionisation or recombination rate for a
  /// different set of T points in the range T_min<T<T_max using already
  /// setup rate and slope lookup table.
  void print_CIR_LUT(
      const std::vector<string> &,  // tracer list
      const std::string info,
      const std::vector<double> &,  // T table
      // rate table
      const std::vector<std::vector<double> > &,
      // rate slope table
      const std::vector<std::vector<double> > &,
      // which rate (collisional or recombination )
      const std::string filename);

  /// ??
  /// ???
  /// ???
  void print_charge_exchange_LUT(
      const std::vector<std::vector<string> > &,  // tracer list
      const std::string info,
      const std::vector<double> &,  // T table
      // rate table
      const std::vector<std::vector<std::vector<double> > > &,
      // rate slope table
      const std::vector<std::vector<std::vector<double> > > &,
      // which rate (ionising or recombining )
      const std::vector<std::vector<double> > &,  // energy-defect table
      const std::string filename);

  /// ??
  /// ???
  /// ???
  void print_to_file_photoionisation(
      const std::string info,
      const std::vector<string> &,                // photo tracer list
      const std::vector<std::vector<double> > &,  // energy bins
      const std::vector<std::vector<double> > &,  // mean cross-section
      const std::string filename                  // filename
  );

  /// Debug tools: write bin luminosity to a file
  void print_stellar_properties(
      const star_properties_datatype &,           // stellar properties
      const std::vector<double> &,                // bin luminosity
      const std::vector<std::vector<double> > &,  // energy bin
      std::string                                 // file path
  );

  /// Print cooling rate look up tables to a file
  void print_cooling_table(
      std::vector<double> &,                // log10(ne)
      std::vector<double> &,                // log10(T)
      std::vector<std::vector<double> > &,  // Cooling rate
      const std::string,                    // Name of the ion
      std::string                           // file path
  );
};
#endif  // NEMO_H
