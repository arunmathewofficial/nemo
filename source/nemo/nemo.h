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

/// Choose the cooling function
// #define MELLEMA
#define CHIANTI

#ifdef SPDLOG_FWD
#include <spdlog/fwd.h>
#endif
#include <spdlog/spdlog.h>
/* prevent clang-format reordering */
#include "microphysics/cvode_integrator.h"
#include "nemo-constants.h"
#include "nemo-datastructs.h"
#include "nemo/atomic_physics_database.h"
#include "nemo/charge_exchange.h"
#include "nemo/coll_ionise_recomb.h"
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


// ############################################################################
//  NEMO CLASS VARIABLES AND FUNCTIONS
///
/// Integrator for the microphysics equations for the non-equilibrium ion
/// fraction of various ions, and for the internal energy density.
///
class nemo :
    virtual public cvode_solver,
    virtual public coll_ionise_recomb,
    virtual public charge_exchange,
#ifdef MELLEMA
    virtual public mellema_cooling,
#elif defined CHIANTI
    virtual public chianti_cooling,
#endif
    virtual public stellar_atmosphere_models,
    virtual public multi_ion_photo_xsections,
    virtual public photoionisation {

  // ############################################################################
  // ############################################################################
  // The public interface to NEMO
  // ############################################################################
public:
  /// CONSTRUCTOR: sets up the machinery for fast execution
  nemo(
      const int,  ///< Total number of variables in state vector
      const int,  ///< Number of passive-tracer vars in state vector
      const int,  ///< index of density variable in state vector
      const int,  ///< index of pressure variable in state vector
      const std::vector<std::string>,  ///< List of what the tracer variables
                                       ///< mean.
      const struct nemo_options,  ///< physics options to use with this NEMO
                                  ///< instance
      std::vector<nemo_star_properties>,  ///< radiation sources.
      const double                        ///< EOS Gamma
  );

  /// DESTRUCTOR: removes all local data
  virtual ~nemo();

  ///
  /// Update microphysics including photoionization and photoheating
  ///
  int update_single_cell(
      const double *,  ///< Primitive Vector to be updated.
      const std::vector<struct nemo_cell_rt_data>
          &,         ///< cell data for radiation source
      double *,      ///< Destination Vector for updated values
      const double,  ///< Time Step to advance by.
      const double,  ///< EOS gamma.
      double *       ///< any returned data (final temperature?).
  );

  ///
  /// Returns the gas temperature.  This is only needed for data output, so
  /// there is no need to make it highly optimized.
  /// - Is threadsafe.
  ///
  double Temperature(
      const double *,  ///< primitive vector
      const double     ///< eos gamma
  );

  ///
  /// Set the gas temperature to a specified value.
  /// Needed if you want this feature to set initial conditions.
  ///
  int Set_Temp(
      double *,      ///< primitive vector.
      const double,  ///< temperature
      const double   ///< eos gamma.
  );

  ///
  /// This returns the minimum timescale of all microphysical
  /// processes, including reaction times for each species and the
  /// total heating/cooling time for the gas.
  /// It requires the radiation field as an input, but the vectors can be null.
  /// Used for setting a simulation timestep, together with CFL condition.
  ///
  double timescales(
      const double *,                                 ///< Current cell.
      const std::vector<struct nemo_cell_rt_data> &,  ///< radiation field.
      const double                                    ///< EOS gamma.
  );

  /// Return the H mass fraction
  double get_X_H();

  /// Get electron number density (cm^{-3})
  double get_n_elec(const double *  ///< primitive state vector.
  );

  /// Get electron number density (cm^{-3})
  double get_n_ion(
      std::string,    ///< ion name
      const double *  ///< primitive state vector.
  );

  /// Calculates the change in optical depth over the cell
  /// for use in the Raytracing module
  void get_dtau(
      const std::vector<double> &,  // input primitive vector
      std::vector<struct nemo_cell_rt_data>
          &  ///< cell data (including dtau vec)
  );

  ///
  /// calculate dy/dt for the vector of y-values.
  ///
  int ydot(
      double,  ///< current time (probably not needed for rate equations)
      const N_Vector,  ///< current Y-value
      N_Vector,        ///< vector for Y-dot values
      const double *  ///< extra user-data vector, P, for evaluating ydot(y,t,p)
  );

  ///
  /// Calculate the Jacobian matrix d(dy_i/dt)/dy_j for a vector of y-values.
  ///
  int Jacobian(
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
    spdlog::info("Jacobian not yet implemented in NEMO!");
    return 1;
  }

  /// Get the number of extra parameters and the number of equations.
  void get_problem_size(
      int *,  ///< number of equations
      int *   ///< number of parameters in user_data vector.
  );

  /// Check that input stellar properties match those in nemo,
  /// and add energy binning info to the struct(s)
  void get_radiation_src_info(
      std::vector<struct nemo_star_properties> &  ///< list of stars
  );

  // ############################################################################
  // ############################################################################

protected:
  // --------------------------------------------------------------------------
  // static variables for book-keeping

  /// Struct with flags for which extra physics we are (not) doing.
  const struct nemo_options opt;
  std::vector<string> scalars_names;  ///< microphysics tracers replica.
  int N_prim;             ///< Size of the primitive vector (prim_vec).
  int prim_first_tracer;  ///< Idx of the first tracer variable in the prim_vec.
  int prim_dens_idx;      ///< index of density var in primitive vector
  int prim_pres_idx;      ///< index of pressure var in primitive vector
  std::vector<int>
      massfrac_species_prim_idx;  ///< local species index in prim_vec.
  std::vector<int>
      massfrac_el_prim_idx;  ///< idx of elemental mass fraction in prim_vec.

  std::vector<vector<string> > all_species_list;  ///< Imported species list.
  std::vector<std::string> element_tracers;  ///< Element tracers in prim_vec.
  std::unordered_map<std::string, int> element_position_index_map;
  std::vector<int> element_atomic_number;   ///< atomic number of elements.
  std::vector<double> element_atomic_mass;  ///< atomic mass of elements.
  int N_elements;  ///< Number of elements in the primitive-vector/tracer.

  std::vector<string> local_tracers;  ///< Extracted local tracers.
  int N_species;                      ///< Number of species.
  std::vector<int>
      species_local_idx;  ///< species idx in local vector (local_vec).
  std::vector<int>
      element_species_count;  ///< number of species belong to same elements.
  std::vector<int> species_element;  ///< species' element positions.
  std::unordered_map<std::string, int> species_prim_index_map;
  // --------------------------------------------------------------------------

  // --------------------------------------------------------------------------
  // imported from base class
  /// Struct with list of radiation sources
  std::vector<struct nemo_star_properties> RS;
  /// Struct with cell data for the radiation field
  std::vector<struct nemo_cell_rt_data> cell_rt;
  /// map tracer names to state-vec index
  std::map<std::string, int> tr_map;
  /// map element names to state-vec index
  std::map<std::string, int> el_map;
  std::vector<int> el_index;  ///< indices of elements (USED???)
  std::vector<int> tr_index;  ///< indices of tracers  (USED???)
  int n_el;                   ///< number of element tracers
  int nv_prim;                ///< number of vars in state vec
  int ntracer;                ///< number of tracer variables
  // --------------------------------------------------------------------------

  /// Cosmic-ray ionisation ===================================================
  std::vector<int>
      neutral_species_local_idx;  ///< neutral species idx in local_vec.
  double cosmicray_ionise_rate;
  double cosmic_ray_energy;
  // --------------------------------------------------------------------------

  /// Collisional ionisation ===================================================

  std::vector<string> collion_species;  ///< collisional ionisation species.
  std::vector<int> collion_species_charge;
  std::vector<int>
      upper_ion_local_idx;  /// local index of upper ionisation state if exist.
  std::vector<double>
      collion_species_IE;  /// coll-ionisation species ionisation energy (egs)
  // --------------------------------------------------------------------------

  /// Recombination ============================================================

  std::vector<string> recomb_species_list;  ///< recombination species.
  /// Track all species that participate in recombination.
  std::vector<struct nemo_recomb_species> rec_s;
  /// Track species reacting with the OTS radiation field from recombining ions
  std::vector<struct nemo_recomb_species> ots_s;
  // --------------------------------------------------------------------------



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
  // --------------------------------------------------------------------------

  /// Photo-ionisation =========================================================

  std::vector<string> photo_species;  ///< photoionisation species.
  int N_photo_species;  ///< Number of photo species in the tracer list.
  std::vector<int>
      photo_species_local_idx;  /// photo-species local vector index.
  std::vector<int>
      photo_upper_ion_local_idx;           /// photo-species upper local index.
  std::vector<int> photo_species_element;  /// Corresponding element's position.
  std::vector<int> photo_species_charge;   /// photo species charge.
  std::vector<int>
      photo_massfrac_species_prim_idx;  /// photo-species prim-vec index.
  std::vector<double>
      photo_species_IE;  ///< Record the ionisation threshold of photo species

  // std::vector<std::vector<double> > ebins;
  // int N_ebins;
  // double photo_Emax;

  // std::vector<double>
  //     bin_luminosity;  ///< luminosity in each energy bin (ergs/s/bins)
  // std::string SED_info;

  // struct nemo_star_properties stellar_properties;
  //  --------------------------------------------------------------------------

  /// Cooling ==================================================================
  std::vector<vector<string> > cooling_species;  ///< cooling species.
  std::ofstream cooling_outfile;                 /// CIE cooling function
  // --------------------------------------------------------------------------

  // --------------------------------------------------------------------------
  // cell data variables, that are specific to each update call
  std::vector<double> element_num_density;  ///< element number density vector.
  /// variable recording species number density for NEMO tracer list
  std::vector<std::vector<double> > species_number_density;
  // double ds     = 0.0;           ///< Path length through the current cell
  // double vshell = 0.0;           ///< Shell volume in the current cell
  // std::vector<double> tau_bins;  ///< Optical depth in each energy bin
  ///	Vector to store initial corrector values (used to modify flux according
  /// to
  /// sCMA)
  std::vector<double> prim_copy;
  // Constant data in the cell, received from cell data.
  double mpv_nH;  ///< total hydrogen number density at current cell.
  // int N_rsrc;     ///< number of radiation sources to include.
  // std::vector<struct rt_source_data> rt_data;  ///< data on radiation
  // sources.

  // --------------------------------------------------------------------------

  // --------------------------------------------------------------------------
  int N_local;  ///< number of variables to solve for: N_species + 1 (Eint)
  int E_index;
  int N_extradata;  ///< number of elements in user-data array.
  N_Vector y_in;    ///< current y-vector
  N_Vector y_out;   ///< output y-vector

  /// Uses this steps to generate lookup tables
  double delta_logT;
  double delta_logne;
  std::vector<double> T_table;   ///< Temperature table
  std::vector<double> ne_table;  ///< electron number density table

  const double tab_T_min;
  const double tab_T_max;
  const int tab_Num_temps;  // NB this needs to be const so can initialise
                            // arrays with it.

  const int tab_Num_ne;
  const double tab_ne_min;
  const double tab_ne_max;
  // --------------------------------------------------------------------------

  /// Function to return collisional ionisation rate of the given species for
  /// a given temperature.
  double collisional_ionisation_rate(
      const int,    // Species identifier
      const double  // Temperature
  );

  /// Function to return recombination rate.
  void recombination_rate(
      const int,     ///< Species identifier
      const double,  ///< Temperature
      double &,      ///< Case A rate for all species
      double &       ///< Case B rate for H+,He+,He2+, otherwise case A
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
  double mellema_cooling_rate(
      const int,     //  species identifier
      const double,  // Temperature
      const double   // ne
  );
#endif

#ifdef CHIANTI
  /// Fucntion to return chianti colling rate.
  double chianti_cooling_rate(
      const int,     //  species identifier
      const double,  // Temperature
      const double   // ne
  );
#endif

  /// Function will generate all the indexing necessary for NEMO module.
  void tracer_species_cataloger(
      const int,  ///< Number of variables in the primitive vector.
      const int,  ///< Number of tracer variables in the primitive vector.
      const std::vector<std::string> &  ///< List of tracer variables.
  );

  ///
  /// set the relative and absolute error tolerances
  ///
  void get_error_tolerances(
      double *,  ///< relative error tolerance (single value)
      double[]   ///< absolute error tolerance (array)
  );

  ///
  /// convert state vector from grid cell into local microphysics vector.
  ///
  int convert_prim2local(
      const double *,         ///< primitive vector (nv_prim)
      std::vector<double> &,  ///< local vector [x(H0),E](n+1).
      int                     /// flags the function
  );

  ///
  /// Convert local microphysics vector into state vector for grid cell.
  /// This is the inverse of convert_prim2local.
  ///
  int convert_local2prim(
      const std::vector<double> &,  ///< local (updated) vector [x(H0),E](n+1).
      const double *,               ///< input primitive vector (nv_prim)
      double *,                     ///< updated primitive vector (nv_prim)
      int                           /// flags the function
  );

  ///
  /// Set local variables (y_ion_frac, y_neutral_frac, ne) by calling
  /// get_y_ion_frac(), get_y_top_ion_frac, and get_ne().
  ///
  void set_localvariables(
      const std::vector<double> &,  ///< local vector
      std::vector<double> &,        ///< y_ion_frac
      std::vector<double> &,        ///< y_neutral_frac
      double &                      ///< ne
  );

  /// Set y_ion_frac
  void get_y_species_frac(
      const std::vector<double> &,  ///< local vector
      std::vector<double> &         ///< y_ion_frac
  );

  /// Set y_top_ion_frac
  void get_y_top_ion_frac(
      const std::vector<double> &,  ///< y_species_frac
      std::vector<double> &         ///< y_top_ion_frac
  );

  /// Calculate electron number density
  void get_ne(
      const std::vector<double> &,  ///< y_species_frac
      double &                      ///< electron number density
  );

  /// Correct local vector
  /// + check that ion fractions are >= ABSTOL.
  /// + check that energy is within temperature limits
  int correct_localvector(
      std::vector<double> &,  /// local vector
      int                     /// flags the function
  );

  ///
  /// returns gas temperature according to E=nkT/(g-1) with n=1.1*nH*(1+x_in),
  /// appropriate for a gas with 10% Helium by number, and if He is singly
  /// ionised whenever H is.
  ///
  double get_temperature(
      const std::vector<double> &,  ///< ion fractions
      const std::vector<double> &,  ///< element number densities
      const double                  ///< E_int (per unit volume)
  );

  ///
  /// Returns total number of particles.
  ///
  double get_ntot(
      const std::vector<double> &,  ///< species fractions (input)
      const std::vector<double> &   ///< element number density (input)
  );

  ///
  /// Set the size of local vectors, and index them with integers.
  ///
  void setup_local_vectors();

  ///
  /// This function is utilized to correct the input primitive vector of the
  /// grid. We utilize the sCMA (Simple Consistent Multi-fluid Advection, Plewa
  /// and Muller, Astron. Astrophys. 342, 179, 1999) scheme. The utility of this
  /// method requires the inputted mass fraction to satisfy the constraint that
  /// the sum of elemental mass fractions must equal unity.
  ///
  void sCMA(
      double *,       ///< corrected vector.
      const double *  ///< input primitive vector.
  );

  const double eos_gamma;  ///< EOS gamma for ideal gas.
  double gamma_minus_one;  ///< as named.
  double Min_NeutralFrac;  ///< minimum H0 fraction allowed (eps=1.0e-12)
  double Max_NeutralFrac;  ///< Maximum H0 fraction allowed (1-eps)

  int print_flag;

  /// For each cell, ydot() needs to know the local radiation field.  This
  /// function takes all of the input radiation sources and calculates what
  /// ydot() needs to know for both heating and ionisation sources.
  ///
  void setup_radiation_source_parameters(
      const std::vector<struct nemo_cell_rt_data>
          &  ///< cell data for radiation source
  );

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

  /// index matching photo_ion xsection to ion in local vector
  std::vector<std::vector<double> > y_ion_xsections;

  // #########################################################################
  // *********************** DEBUG  FUNCTIONS ********************************
  // #########################################################################

  /// This will print collisional ionisation or recombination rate for a
  /// different set of T points in the range tab_T_min<T<tab_T_max using already
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
      const nemo_star_properties &,  // stellar properties
      std::string                    // file path
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
