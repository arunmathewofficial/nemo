#pragma clang diagnostic push
#pragma ide diagnostic ignored "UnreachableCode"
/// \file nemo.cpp
/// \authors Arun Mathew, Jonathan Mackey, Maggie Celeste Goulden
/// \date 2018.10
///
/// Description:
/// - multi-species non-equilibrium chemistry solver.

///
/// The integration method uses the CVODES solver from the SUNDIALS
/// package by (Cohen, S. D., & Hindmarsh, A. C. 1996, Computers in
/// Physics, 10, 138) available from
///  https://computation.llnl.gov/casc/sundials/main.html
/// The method use backwards differencing (i.e. implicit) with Newton
/// iteration.
///
/// Modifications:
/// - 2018.10.09 JM: edited header
///

#include "nemo.h"
#ifdef SPDLOG_FWD
#include <spdlog/fwd.h>
#endif
#include <spdlog/spdlog.h>
/* prevent clang-format reordering */
#include <fmt/ranges.h>

#ifndef NDEBUG
#endif

using namespace std;

// #define NEMO_CIE_TEST
// #define NEMO_CIE_COOLING
// #define NEMO_DEBUG

// Timestep-limiting is important for making chemistry consistent
// with hydrodynamics.
// 0.25 is a good value for MPv3 (see Mackey,2012,A&A,539,A147)
// #define DTFRAC 0.25
#define DTFRAC 0.1

// Maximum charge of the reactant in charge exchange reactions
#define CX_MaxQ 4  // Note: CX_MaxQ can only take values 1, 2, 3, or 4

// #define OTS_HACK

// GET ERROR TOLERANCES #################################################
void nemo::get_error_tolerances(
    double *reltol,  ///< relative error tolerance.
    double atol[]    ///< absolute error tolerances
)
{
  *reltol = NEMO_RELTOL;
  for (int i = 0; i < N_equations - 1; i++)
    atol[i] = NEMO_ABSTOL;  ///< minimum neutral fraction I care about.
  atol[N_equations - 1] = NEMO_MINERG;  ///< E_int: for n=1.0, T=1.0e4, ==>
                                        ///< E=2.07e-12, so say 1e-17?
  return;
}



// GET PROBLEM SIZE  ######################################################
void nemo::get_problem_size(
    int *n_eqn,  ///< number of equations
    int *n_para  ///< number of parameters in user_data vector.
)
{
  *n_eqn  = N_equations;
  *n_para = N_extradata;
  return;
}

// CONSTRUCTOR ############################################################
nemo::nemo(
    const int nd,    ///< grid dimensions
    const int csys,  ///< Coordinate System flag
    const int nv,    ///< Total number of variables in state vector
    const int ntr,   ///< Number of tracer variables in state vector.
    const std::string *tracers,   ///< List of what the tracer variables mean.
    struct which_physics *ephys,  ///< pointer to extra physics flags.
    struct rad_sources *rsrcs,    ///< radiation sources.
    const double g                ///< EOS Gamma
    ) :
    microphysics_base(nv, ntr, tracers, ephys, rsrcs), coll_ionise_recomb(),
    charge_exchange(CX_MaxQ),
#ifdef MELLEMA
    mellema_cooling(),
#elif defined CHIANTI
    chianti_cooling(),
#else
#error "must define MELLEMA or CHIANTI"
#endif
    ndim(nd), eos_gamma(g), coord_sys(csys), T_min(1e1), T_max(1e9),
    Num_temps(300), Num_ne(100), ne_min(1.0), ne_max(1.0e6),
    stellar_atmosphere_models(), multi_ion_photo_xsections(), photoionisation()
{

  spdlog::debug(
      "nemo: Setting up tracer variables, assuming tracers"
      " are last {} variables in primitive vector",
      ntracer);

  // Setup stellar atmosphere models for stellar luminosity spectrum
  if (EP->phot_ionisation) {
    setup_stellar_atmosphere_models(ebins, N_ebins, photo_Emax);
    // Note: This will return energy bins, number of energy bins and maximum
    // energy in the energy bin vector.
  }

  // Enable the charge-exchange reaction only if both charge_exchange and
  // coll_ionisation are enabled in the parameter file.
  if (EP->charge_exchange && EP->coll_ionisation) {
    charge_exchange_enable = true;
    spdlog::info("Charge-exchange reaction(s) enabled");
  }
  else {
    charge_exchange_enable = false;
    spdlog::info("Charge-exchange reaction(s) disabled");
  }

  // Section: Invoking NEMO chemical Species cataloger ***********************
  tracer_species_cataloger(nv, ntr, tracers);

  // Section: Setting up LOOK-UP-TABLE variables ******************************
  // Setting up step size for log temperature (T) and log electron number
  // density (ne). This step size is used to generate look up tables. The same
  // is used to calculate temperature index and electron number density index to
  // get the rate for a specific value of T or ne. Setting up uniform step size
  // for log T.
  delta_logT = (log10(T_max) - log10(T_min)) / (Num_temps - 1);
  // Make temperature table for all LUTs.
  T_table.resize(Num_temps);
  for (int i = 0; i < Num_temps; i++) {
    // Temperature table for coll_ionise and recomb LUTs:
    T_table[i] = pow(10.0, (log10(T_min) + i * delta_logT));
  }
  // Setting up uniform step size for log ne.
  delta_logne = (log10(ne_max) - log10(ne_min)) / (Num_ne - 1);
  // Make electron number density table for all LUTs.
  ne_table.resize(Num_ne);
  for (int i = 0; i < Num_ne; i++) {
    ne_table[i] = pow(10, log10(ne_min) + i * delta_logne);
  }

  // Section: Cosmic-ray ionisation  *******************************************
  cosmicray_ionise_rate = 1.000E-17;  // cosmic ray ionisation rate in s^{-1}.
  // Following Goldsmith and Langer (1978), it is assumed that each ionisation
  // deposits 20 eV of energy into the gas.
  cosmic_ray_energy = 20.0 * pconst.eV();

  // Section: Collisional ionisation  ******************************************
  if (EP->coll_ionisation) {
    // setting up collisional ionisation
    setup_collisional_ionisation(collion_species);
    // Generating collisional ionisation rate lookup tables
    generate_collisional_ionisation_LUT(T_table);
#ifdef NEMO_DEBUG
    mkdir("NEMO-Tables", 0777);
    mkdir("NEMO-Tables/coll-recom-tables", 0777);
    spdlog::debug("nemo: writting collisional ionisation rate table into file");
    std::string collion_info = "-- File: Collisional ionisation rate table --\n"
                               "-- Generated by NEMO module PION v3.0 --\n";
    print_CIR_LUT(
        collion_species, collion_info, T_table, collisional_rate_table,
        collisional_slope_table,
        "NEMO-Tables/coll-recom-tables/coll-ionise-rate");
#endif
  }

  // Section: Recombination  ***************************************************
  if (EP->rad_recombination) {
    // setting up recombination
    setup_recombination(recomb_species);
    // Generating recombination rate lookup tables
    generate_recombination_LUT(T_table);
#ifdef NEMO_DEBUG
    mkdir("NEMO-Tables", 0777);
    mkdir("NEMO-Tables/coll-recom-tables", 0777);
    spdlog::debug("nemo: writting recombination rate table into txt file");
    std::string recomb_info = "-- File: Recombination rate table --\n"
                              "-- Generated by NEMO module PION v3.0 --\n";
    print_CIR_LUT(
        recomb_species, recomb_info, T_table, recombination_rate_table,
        recombination_slope_table, "NEMO-Tables/coll-recom-tables/recomb-rate");
#endif
  }

  // Section: Charge exchange reactions  ***************************************
  if (charge_exchange_enable) {
    // setting up charge-exchange reactions module
    setup_charge_exchange(ionise_cx_reactions, recomb_cx_reactions);
    // Generating charge-exchange rate lookup tables
    generate_charge_exchange_LUT(T_table);
    // Generate ionising charge-exchange reaction energy-defect table
    std::vector<std::vector<double> > ionise_cx_energy_eV;
    ionise_cx_energy_eV = import_energy_defect(ionise_cx_reactions);
    spdlog::debug("Ionising energy-defect table {}", ionise_cx_energy_eV);
    // converting energy defect in eV to erg
    ionise_cx_energy = ionise_cx_energy_eV;
    for (size_t i = 0; i < ionise_cx_energy.size(); i++)
      for (size_t j = 0; j < ionise_cx_energy[i].size(); j++)
        ionise_cx_energy[i][j] *= pconst.eV();
    // Generate recombining charge-exchange reaction energy-defect table
    std::vector<std::vector<double> > recomb_cx_energy_eV;
    recomb_cx_energy_eV = import_energy_defect(recomb_cx_reactions);
    spdlog::debug("Recombining energy-defect table {}", recomb_cx_energy_eV);
    // converting energy defect in eV to erg
    recomb_cx_energy = recomb_cx_energy_eV;
    for (size_t i = 0; i < recomb_cx_energy.size(); i++)
      for (size_t j = 0; j < recomb_cx_energy[i].size(); j++)
        recomb_cx_energy[i][j] *= pconst.eV();


#ifdef NEMO_DEBUG
    mkdir("NEMO-Tables", 0777);
    mkdir("NEMO-Tables/charge-exchange-tables", 0777);
    spdlog::debug(
        "nemo: writting charge-exchange reaction rate table into file");
    std::string ionise_cx_info =
        "-- File: Charge-exchange rate table for ionisation reactions. --\n"
        "-- Generated by NEMO module PION v3.0 --\n";
    print_charge_exchange_LUT(
        import_reactions_id(ionise_cx_reactions), ionise_cx_info, T_table,
        ionise_cx_rate_table, ionise_cx_slope_table, ionise_cx_energy_eV,
        "NEMO-Tables/charge-exchange-tables/cx-ionise-reactions");
    std::string recomb_cx_info =
        "-- File: Charge-exchange rate table for recombination reactions. --\n"
        "-- Generated by NEMO module PION v3.0 --\n";
    print_charge_exchange_LUT(
        import_reactions_id(recomb_cx_reactions), recomb_cx_info, T_table,
        recomb_cx_rate_table, recomb_cx_slope_table, recomb_cx_energy_eV,
        "NEMO-Tables/charge-exchange-tables/cx-recomb-reactions");
#endif
  }

  // Section: Photo-ionisation and photo-heating *******************************
  if (EP->phot_ionisation) {
    spdlog::debug(" nemo: Energy bins {}", ebins);
    spdlog::debug(" nemo: Number of energy bins {}", N_ebins);
    spdlog::debug(" nemo: Maximum photo-ionisation energy {} eV", photo_Emax);

    // Setup photo-ionisation x-section module: multi_ion_photo_xsection class.
    // This will determine the appropriate size for the x-section and bin
    // fraction tables. Additionally, it will return the photo tracer species'
    // ionisation threshold vector.
    photo_species_IE = setup_photo_xsections(photo_species);
    spdlog::debug(" nemo: ionisation threshold {}", photo_species_IE);

    // Generating mean cross-section and bin fraction table for the energy bins
    generate_ebin_mean_xsection(ebins);

#ifdef NEMO_DEBUG
    mkdir("NEMO-Tables/photo-tables", 0777);
    // print to file photo-ionisation cross-section tables
    spdlog::debug("nemo: writting photo cross-section table into file");
    const std::string xsection_info = "Mean photo-ionization cross-section"
                                      " for the photo-tracer species.";
    print_to_file_photoionisation(
        xsection_info, photo_species, ebins, photo_xsection_table,
        "NEMO-Tables/photo-tables/mean-photo-xsection");

    // print to file bin weight tables
    spdlog::debug(" nemo: writting bin weight table into file");
    const std::string weight_info = "Bin weight for the photo-tracer species.";
    print_to_file_photoionisation(
        weight_info, photo_species, ebins, bin_weight_table,
        "NEMO-Tables/photo-tables/bin-weights");
#endif

    // Identify the radiation sources
    N_rsrc = 0;
    for (int src = 0; src < RS->Nsources; src++) {
      // proceed if the source type is SINGLE
      if (RS->sources[src].type == RT_SRC_SINGLE) {

        // check number of Ntau equals to SED energy bins
        // Size of energy bin in parameter file should match the
        // energy bin hardcoded in SED.
        if (RS->sources[src].NTau != N_ebins) {
          spdlog::error(
              " nemo: Incompatible number of energy bins "
              "between SED and parameter file-{}.",
              RS->sources[src].NTau);
          exit(1);
        }
        N_rsrc++;
        // setting the size of the vector rt_data
        rt_data.resize(N_rsrc);
        // Generate the energy flux in different energy bin for the given
        // radiation source.
        int err = generate_bin_luminosity(
            RS->sources[src].effect, RS->sources[src].Tstar,
            RS->sources[src].Rstar, RS->sources[src].strength, bin_luminosity,
            SED_info);
        if (err) {
          spdlog::error(
              "Attempt to generate stellar spectrum for source-id: {}"
              " failed.",
              RS->sources[src].id);
          exit(1);
        }
      }
    }

    spdlog::debug("nemo: {} radiation source(s) detected.", N_rsrc);

    // Setup photo-ionisation module
    setup_photoionisation(
        ebins, N_photo_species, photo_xsection_table, bin_weight_table,
        bin_luminosity, photo_species_IE);

#ifdef NEMO_DEBUG
    mkdir("NEMO-Tables/stellar_model", 0777);
    stellar_properties = {
        RS->sources[0].type, RS->sources[0].Tstar, RS->sources[0].Rstar,
        RS->sources[0].strength, SED_info};
    // print to file photo-ionisation cross-section tables
    spdlog::debug(" nemo: writting bin luminosity into file");
    print_stellar_properties(
        stellar_properties, bin_luminosity, ebins,
        "NEMO-Tables/stellar_model/bin_luminosity");
#endif
  }


  // Section: Cooling Function *************************************************
#ifdef MELLEMA
  // MELLEMA COOLING FUNCTION ############
  setup_mellema_cooling(cooling_species);
  spdlog::debug(
      "nemo: Cooling species identifier = {}", cooling_species_identifier);
  generate_mellema_table(T_table, ne_table);
#ifdef NEMO_DEBUG
  mkdir("NEMO-Tables/cooling-tables", 0777);
  mkdir("NEMO-Tables/cooling-tables/mellema", 0777);
  spdlog::debug("nemo: writting Mellema cooling tables into txt files");
  int sct = 0;
  for (int i = 0; i < cooling_species.size(); i++) {
    for (int j = 0; j < cooling_species[i].size(); j++) {
      print_cooling_table(
          ne_table, T_table, mellema_table[sct].rate, cooling_species[i][j],
          "NEMO-Tables/cooling-tables/mellema/mellema_rates");
      sct++;
    }
  }
#endif
#endif

#ifdef CHIANTI
  // CHIANTI COOLING FUNCTION ############
  setup_chianti_cooling(cooling_species);
  spdlog::debug(
      "nemo: Cooling species identifier = {}", cooling_species_identifier);
  generate_chianti_table(T_table, ne_table);
#ifdef NEMO_DEBUG
  mkdir("NEMO-Tables/cooling-tables", 0777);
  mkdir("NEMO-Tables/cooling-tables/chianti", 0777);
  spdlog::debug("nemo: writting Chianti cooling tables into txt files");
  int sct = 0;
  for (size_t i = 0; i < cooling_species.size(); i++) {
    for (size_t j = 0; j < cooling_species[i].size(); j++) {
      print_cooling_table(
          ne_table, T_table, chianti_table[sct].rate, cooling_species[i][j],
          "NEMO-Tables/cooling-tables/chianti/chianti_rates");
      sct++;
    }
  }
#endif
#endif

#ifdef NEMO_CIE_COOLING
  // create a text file to write the cooling as the function of temperature.
  std::ofstream cooling_outfile("cooling_function.txt");
  cooling_outfile
      << "# Collisional ionisation equilibrium (CIE) cooling function for \n";
  cooling_outfile << "# all species in the cooling species list.\n";
  cooling_outfile
      << "# Column-1: Temperature (K), Other columns: Cooling function (erg cm^3/s).\n";
#endif



  // Initializing other important variables and vectors ************************
  // Set up local variables
  setup_local_vectors();
  // Set gamma - 1 from EOS gamma
  gamma_minus_one = eos_gamma - 1.0;
  // Set tolerance
  Min_NeutralFrac = NEMO_ABSTOL;
  Max_NeutralFrac = 1.0 - NEMO_ABSTOL;

  // Initializing CVODE Solver *************************************************
  // Initialise the CVODES solver memory etc.
  setup_cvode_solver_without_Jacobian();


  // CONSTRUCTOR END MESSAGE **********************************************
  spdlog::debug("nemo: Constructor finished and returning");
}



// #######################################################################
void nemo::setup_local_vectors()
{
  //
  // This is in a function so it can be replaced in an inherited class.
  //
  N_local     = N_species + 1;  // number of species plus energy
  N_extradata = 0;
  N_equations = N_local;
#if defined CVODE6 || defined CVODE7
  y_in  = N_VNew_Serial(N_equations, sunctx);
  y_out = N_VNew_Serial(N_equations, sunctx);
#else
  y_in  = N_VNew_Serial(N_equations);
  y_out = N_VNew_Serial(N_equations);
#endif
  E_index = N_species;
  return;
}



// NEMO DESTRUCTOR #########################################################
nemo::~nemo()
{
  // Free vector memory
  N_VDestroy_Serial(y_in);
  N_VDestroy_Serial(y_out);
}



// TRACER SPECIES CATALOGER ################################################
void nemo::tracer_species_cataloger(
    const int N_primitive, const int N_tracers, const std::string *tracers)
{

  // Beginning of cataloging tracers
  spdlog::debug("nemo: Cataloging chemical tracers");

  spdlog::debug(
      "nemo: Chemical tracers in primitive vector {} ",
      std::vector<string>(tracers, tracers + N_tracers));

  // Make a copy of tracer list
  std::vector<int> els(N_tracers, 0);
  tracers_replica.resize(N_tracers);
  for (int s = 0; s < N_tracers; s++)
    tracers_replica[s] = tracers[s];

  // Number of elements in the primitive vector.
  N_prim = N_primitive;
  prim_copy.resize(N_prim, 0.0);

  // Index of the first tracer variable in the primitive vector.
  first_tracer_index = N_prim - N_tracers;

  // Identifying the chemical elements and import all its species from *********
  // atomic-physics database.
  for (int s = 0; s < N_tracers; s++) {
    if (tracers_replica[s].substr(0, 2) == "X_") {
      tracers_replica[s].erase(0, 2);
      element_tracers.push_back(tracers_replica[s]);
      if (!atomic_physics_data.find_species(tracers_replica[s])) {
        spdlog::warn(
            "no species match for tracer: {}  {}", s, tracers_replica[s]);
        exit_nemo(1);
      }
      else
        all_species_list.push_back(
            atomic_physics_data.import_species(tracers_replica[s]));
      els[s] = 1;
    }
  }

  // get total number of chemical elements.
  N_elements = element_tracers.size();


  // Check if the condition for charge exchange reactions is met, ensuring
  // that the reaction involves exactly two reactants from two different
  // elemental species. If the number of elements is not equal to two or more,
  // pass an error message.
  if (charge_exchange_enable && N_elements == 1) {
    spdlog::error("Insufficient elements for charge exchange."
                  " Check config file.");
    exit(1);
  }

  // Sort element tracers according to position in period table.
  bool isElementOrdered = std::is_sorted(
      element_tracers.begin(), element_tracers.end(),
      [](const std::string &element_1, const std::string &element_2) {
        return atomic_physics_data.get_atomic_number(element_1)
               < atomic_physics_data.get_atomic_number(element_2);
      });

  if (!isElementOrdered) {
    spdlog::error(" NEMO: Elemental mass fraction are not arranged correctly"
                  " within the primitive vector.");
    exit_nemo(1);
  }

  int tr = 0;
  for (int e = 0; e < N_elements; e++) {
    while (els[tr] == 0)
      tr++;
    massfrac_idx.push_back(first_tracer_index + tr);
    tr++;
    element_atomic_number.push_back(
        atomic_physics_data.get_atomic_number(element_tracers[e]));
    element_atomic_mass.push_back(
        atomic_physics_data.get_atomic_mass(element_tracers[e]));
    element_position_index_map[element_tracers[e]] = e;
  }

  // Resizing element number density vector
  element_num_density.resize(N_elements);

  spdlog::debug(" nemo: Element tracers in prim-vec {}", element_tracers);
  spdlog::debug(" nemo: Elemental mass fraction prim idx {} ", massfrac_idx);
  spdlog::debug(
      " nemo: Element position index map {} ", element_position_index_map);
  spdlog::debug(
      " nemo: Element atomic number vector {} ", element_atomic_number);
  spdlog::debug(" nemo: Element atomic mass vector {} gm", element_atomic_mass);

  // Sort the imported species list according to their position in period table.
  std::sort(
      all_species_list.begin(), all_species_list.end(),
      [](const std::vector<string> &a, const std::vector<string> &b) {
        return atomic_physics_data.get_atomic_number(a[0])
               < atomic_physics_data.get_atomic_number(b[0]);
      });

  spdlog::debug(
      " nemo: Atomic-physics DB imported species list {}", all_species_list);
  // Removing mass fraction tracers from tracers_replica.
  vector<string>::iterator a = tracers_replica.begin();
  for (size_t v = 0; v < els.size(); v++) {
    a++;
    if (els[v] == 1) {
      tracers_replica.erase(a - 1);
      a--;  // we deleted an element so the iterator needs to go back one
    }
  }

  // Remove non-chemistry tracers from tracers_replica:
  int notchem = 0;
  a           = tracers_replica.end();
  vector<vector<string>::iterator> a2;
  for (int v = tracers_replica.size() - 1; v >= 0; v--) {
    int position = -1;
    for (size_t i = 0; i < all_species_list.size(); ++i) {
      auto it = find(
          all_species_list[i].begin(), all_species_list[i].end(),
          tracers_replica[v]);
      if (it != all_species_list[i].end()) {
        position = 1;
      }
    }
    if (position == -1) {
      notchem += 1;
      a = tracers_replica.begin() + v;
      a2.push_back(a);
    }
    else {
    }
  }
  for (int v = 0; v < a2.size(); v++) {
    spdlog::info(" NEMO: Removing non-chemistry tracers: {} {}", v, *(a2[v]));
    tracers_replica.erase(a2[v]);
  }

  // Now we have removed all non-species tracers, so copy this to local_tracers
  local_tracers = tracers_replica;
  if (local_tracers.size()
      != N_prim - N_elements - first_tracer_index - notchem) {
    spdlog::error(
        " NEMO: local tracers size mismatch with local variables"
        " in prim_vec {} ",
        local_tracers.size());
    exit_nemo(1);
  }

  // Check whether the local tracer is arranged in ascending order based on both
  // atomic number and charge.
  bool isSorted = std::is_sorted(
      local_tracers.begin(), local_tracers.end(),
      [](const std::string &species_1, std::string &species_2) {
        if (atomic_physics_data.get_atomic_number(species_1)
            == atomic_physics_data.get_atomic_number(species_2)) {
          return atomic_physics_data.get_charge(species_1)
                 < atomic_physics_data.get_charge(species_2);
        }
        return atomic_physics_data.get_atomic_number(species_1)
               < atomic_physics_data.get_atomic_number(species_2);
      });

  if (!isSorted) {
    spdlog::error(
        " NEMO: Incorrect arrangement of local tracers within the primitive vector.");
    exit_nemo(1);
  }

  spdlog::debug(
      " NEMO: Local tracers {} {} ", local_tracers.size(), local_tracers);
  N_species = local_tracers.size();

  // resizing all species number density 2D vector
  species_number_density.resize(all_species_list.size());
  for (int e = 0; e < all_species_list.size(); e++) {
    species_number_density[e].resize(all_species_list[e].size(), 0.0);
  }

  // title?? ############################################
  for (int s = 0; s < N_species; s++) {
    // get position and charge of the species by locating the species in
    // the atomic-physics imported species list, all_species_list.
    int position = -1, Z = -1, charge = -1;
    std::string element;
    for (size_t i = 0; i < all_species_list.size(); ++i) {
      auto it = find(
          all_species_list[i].begin(), all_species_list[i].end(),
          local_tracers[s]);
      if (it != all_species_list[i].end()) {
        position = i;
        charge   = distance(all_species_list[i].begin(), it);
        element  = atomic_physics_data.get_element_name(local_tracers[s]);
      }
    }
    if (position == -1) {
      spdlog::warn(
          "{} {} is not a known chemical species", s, local_tracers[s]);
      continue;
    }

    // get atomic number of the species.
    Z = atomic_physics_data.get_atomic_number(
        all_species_list[position][charge]);
    // Setting up species local and primitive index vector.
    // JM: loop through tracers until we find species
    int sp = 0;
    while (all_species_list[position][charge] != tracers[sp])
      sp++;

    species_prim_idx.push_back(first_tracer_index + sp);
    species_local_idx.push_back(s);

    // get element_species_count.
    element_species_count.resize(N_elements);
    element_species_count[position]++;

    // map tracer string to its primitive index
    species_prim_index_map[local_tracers[s]] = first_tracer_index + sp;

    // collisional ionisation species ******************************************
    if (EP->coll_ionisation) {
      if (charge < Z) {
        // appending collisional ionisation species.
        collion_species.push_back(local_tracers[s]);
        collion_species_charge.push_back(charge);

        // calculating upper ion local index.
        if (s + 1 < N_species
            && atomic_physics_data.get_element_name(local_tracers[s + 1])
                   == element)
          upper_ion_local_idx.push_back(s + 1);
        else
          upper_ion_local_idx.push_back(-1);

        // recording the position of corresponding element of the local species.
        species_element.push_back(position);
        // recording the ionisation energy (erg) of collisional ionisation
        // species.
        collion_species_IE.push_back(
            atomic_physics_data.get_ionisation_potential(local_tracers[s]));

        // record the position of neutral species in local vector.
        if (charge == 0) neutral_species_local_idx.push_back(s);
      }
      else {
        spdlog::error(
            " NEMO: Species cannot be collisionally ionized {} ",
            local_tracers[s]);
        exit_nemo(1);
      }
    }
    // end of collisional ionisation species ***********************************

    // recombination species ***************************************************
    if (EP->rad_recombination) {
      // Making recombination species list.
      // if hydrogen exist in the tracer list, then add H1+ ion.
      if (charge == 0 && charge == Z - 1) {
        recomb_species.push_back(all_species_list[position][Z]);
        recomb_species_element.push_back(position);
        recomb_species_charge.push_back(1);
        recomb_top_ion_bool.push_back(1);
      }
      // append all other species
      if (charge != 0) {  // avoiding neutral species
        recomb_species.push_back(local_tracers[s]);
        // Incorporating the highest ionization state of the respective element
        // into the recombination species list, as it undergoes recombination
        // with the electron. This chemical species is not included in the
        // primitive vector; thus, it is imported from the atomic physics
        // database. Note: In this step hydrogen ion is not included.
        recomb_species_element.push_back(position);
        recomb_species_charge.push_back(charge);
        recomb_top_ion_bool.push_back(0);
        if (charge == Z - 1) {
          recomb_species.push_back(all_species_list[position][Z]);
          recomb_species_element.push_back(position);
          recomb_species_charge.push_back(charge + 1);
          recomb_top_ion_bool.push_back(1);
        }
      }
    }
    // end of recombination species ********************************************

    // photo-ionisation species ************************************************
    if (EP->phot_ionisation) {
      // if the ionisation energy is less than the maximum photon energy,
      // then the tracer is qualified photo species
      if (charge < Z
          & atomic_physics_data.get_ionisation_energy(local_tracers[s])
                < photo_Emax) {
        photo_species.push_back(local_tracers[s]);
        // record the photo-species local index.
        photo_species_local_idx.push_back(s);
        // record the photo-species primitive index.
        photo_species_prim_idx.push_back(first_tracer_index + sp);

        // calculating upper ion local index for the current photo-species.
        if (s + 1 < N_species
            && atomic_physics_data.get_element_name(local_tracers[s + 1])
                   == element)
          photo_upper_ion_local_idx.push_back(s + 1);
        else
          photo_upper_ion_local_idx.push_back(-1);

        // recording the position of element of the photo species.
        photo_species_element.push_back(position);
        // recording the charge of each photo species;
        photo_species_charge.push_back(charge);
      }
    }
    // end of photo-ionisation species *****************************************
  }  // end of local tracers loop

  // Display species primitive and local indices
  spdlog::debug(
      " nemo: Found {} elements and {} species in primitive"
      " vector",
      N_elements, N_species);
  spdlog::debug(" NEMO: Local species primitive idx {} ", species_prim_idx);
  spdlog::debug(" NEMO: Species local idx {} ", species_local_idx);
  spdlog::debug(" NEMO: Species corresponding element {} ", species_element);
  spdlog::debug(" NEMO: Upper ion idx in local {} ", upper_ion_local_idx);
  spdlog::debug(
      " NEMO: Species count of same element {} ", element_species_count);
  spdlog::debug(
      " NEMO: Species-primitive index map {} ", species_prim_index_map);
  spdlog::debug(
      " NEMO: Neutral Species local idx {} ", neutral_species_local_idx);

  // Collisional species and their related indices.
  if (EP->coll_ionisation) {
    spdlog::debug(
        " NEMO: Collisional ionisation species {} {} ", collion_species.size(),
        collion_species);
    spdlog::debug(
        " NEMO: Collisional ionisation species element {} ", species_element);
    spdlog::debug(
        " NEMO: Collisional ionisation species charge {} ",
        collion_species_charge);
  }

  // Recombination species and their related indices.
  if (EP->rad_recombination) {
    spdlog::debug(
        " NEMO: Recombination species {} {} ", recomb_species.size(),
        recomb_species);
    spdlog::debug(
        " NEMO: Recombination species element {} ", recomb_species_element);
    spdlog::debug(
        " NEMO: Recombination species charge {} ", recomb_species_charge);
    spdlog::debug(
        " NEMO: Recombination top ion boolean {} ", recomb_top_ion_bool);
  }


  // Charge exchange reaction **************************************************
  if (charge_exchange_enable) {

    // Creating a reference species vector for the charge exchange reaction.
    // This vector will contain all possible species participating in the
    // ionisation-recombination charge exchange reaction.
    std::set<std::string> combined_set(
        collion_species.begin(), collion_species.end());
    combined_set.insert(recomb_species.begin(), recomb_species.end());
    std::vector<std::string> cx_ref_species_list(
        combined_set.begin(), combined_set.end());
    spdlog::debug(" NEMO: cx reference species set {}", cx_ref_species_list);

    // ionisation charge exchange reactions
    for (size_t s = 0; s < collion_species.size(); s++) {
      // Get the set of charge exchange reactant species participating
      // in ionisation reactions and the corresponding .
      charge_exchange_ionise_species.push_back(collion_species[s]);
      ionise_cx_reactions.push_back(
          get_cx_reactions(collion_species[s], ionReact, cx_ref_species_list));
    }

    // possible species participating in ionisation charge exchange reactions
    spdlog::debug(
        " NEMO: cx ionisation species {}", charge_exchange_ionise_species);
    // ionisation charge exchange reactions
    spdlog::debug(
        " NEMO: cx ionisation reactions list {}", ionise_cx_reactions);

    // recombining charge exchange reactions
    for (size_t s = 0; s < recomb_species.size(); s++) {
      // Get the set of charge exchange reactant species participating
      // in recombination reactions, limited by the maximum reactant charge.
      charge_exchange_recomb_species.push_back(recomb_species[s]);
      recomb_cx_reactions.push_back(
          get_cx_reactions(recomb_species[s], recReact, cx_ref_species_list));
    }

    spdlog::debug(
        " NEMO: cx recombination species {}", charge_exchange_recomb_species);
    spdlog::debug(
        " NEMO: cx recombination reactions list {}", recomb_cx_reactions);

    // Get position of the element and charge for the second reactants
    // for every ionisation charge-exchange reactions.
    ionise_cx_second_reactant_element.resize(ionise_cx_reactions.size());
    ionise_cx_second_reactant_charge.resize(ionise_cx_reactions.size());
    has_ionise_cx_exist.resize(ionise_cx_reactions.size(), false);
    for (size_t s = 0; s < ionise_cx_reactions.size(); s++) {
      for (size_t r = 0; r < ionise_cx_reactions[s].size(); r++) {
        // charge exchange reaction exist for this ionisable species
        has_ionise_cx_exist[s] = true;
        // get second reactant's elemental position in the tracer list
        ionise_cx_second_reactant_element[s].push_back(
            element_position_index_map[atomic_physics_data.get_element_name(
                ionise_cx_reactions[s][r].second)]);
        // get second reactant's charge
        ionise_cx_second_reactant_charge[s].push_back(
            atomic_physics_data.get_charge(ionise_cx_reactions[s][r].second));
      }
    }

    // spdlog::debug(
    //    " NEMO: cx exist for the ionisation species {}",
    //    has_ionise_cx_exist);
    spdlog::debug(
        " NEMO: ionise cx second reactant's elemental"
        " position {}",
        ionise_cx_second_reactant_element);
    spdlog::debug(
        " NEMO: ionise cx second reactant's charge {}",
        ionise_cx_second_reactant_charge);

    // get position of the element and charge for the second reactants
    // for every recombination charge-exchange reactions.
    recomb_cx_second_reactant_element.resize(recomb_cx_reactions.size());
    recomb_cx_second_reactant_charge.resize(recomb_cx_reactions.size());
    has_recomb_cx_exist.resize(ionise_cx_reactions.size(), false);
    for (size_t s = 0; s < recomb_cx_reactions.size(); s++) {
      for (size_t r = 0; r < recomb_cx_reactions[s].size(); r++) {
        // charge exchange reaction exist for this recombinable species
        has_recomb_cx_exist[s] = true;
        // get second reactant's elemental position in the tracer list
        recomb_cx_second_reactant_element[s].push_back(
            element_position_index_map[atomic_physics_data.get_element_name(
                recomb_cx_reactions[s][r].second)]);
        // get second reactant's charge
        recomb_cx_second_reactant_charge[s].push_back(
            atomic_physics_data.get_charge(recomb_cx_reactions[s][r].second));
      }
    }

    // spdlog::debug(
    //    " NEMO: cx exist for the recombination species {}",
    //    has_recomb_cx_exist);
    spdlog::debug(
        " NEMO: recomb cx second reactant's elemental"
        " position {}",
        recomb_cx_second_reactant_element);
    spdlog::debug(
        " NEMO: recomb cx second reactant's charge {}",
        recomb_cx_second_reactant_charge);
  }


  // Photo-ionisation species and their related indices. *********************
  if (EP->phot_ionisation) {
    // Calculate number of photo-species
    N_photo_species = photo_species.size();
    spdlog::debug(" NEMO: Photo-ionisation species {} ", photo_species);
    spdlog::debug(" NEMO: No of Photo-species {} ", N_photo_species);
    spdlog::debug(
        " NEMO: Photo-species local idx {} ", photo_species_local_idx);
    spdlog::debug(" NEMO: Photo-species prim idx {} ", photo_species_prim_idx);
    spdlog::debug(
        " NEMO: Photo upper ion local idx {} ", photo_upper_ion_local_idx);
    spdlog::debug(" NEMO: Photo-species element {} ", photo_species_element);
    spdlog::debug(" NEMO: Photo-species charge {} ", photo_species_charge);
  }
  // Cooling species and their related indices.
  cooling_species = all_species_list;
  if (EP->cooling) {
    spdlog::debug(" NEMO: Cooling species {} ", cooling_species);
  }
}


// GET ELECTRON NUMBER DENSITY ################################################
// Calculates electron number density from primitive vector
double nemo::get_n_elec(const pion_flt *prim  ///< primitive state vector array.
)
{
  int sct               = 0;
  double ne_per_density = 0.0;
  for (int e = 0; e < N_elements; e++) {  // loop over every element
    double top_ion_massfrac = prim[massfrac_idx[e]];
    double Ne_massfrac_product =
        0.0;  // product of Ne and species mass fraction.
    for (int q = 0; q < element_species_count[e]; q++) {
      // sum of product of number of electron and corresponding mass fraction.
      Ne_massfrac_product += q * prim[species_prim_idx[sct]];
      // calculating highest ion mass fraction
      top_ion_massfrac -= prim[species_prim_idx[sct]];
      sct++;
    }
    top_ion_massfrac = max(NEMO_ABSTOL, top_ion_massfrac);
    top_ion_massfrac = min(1.0 - NEMO_ABSTOL, top_ion_massfrac);

    ne_per_density +=
        (Ne_massfrac_product + element_atomic_number[e] * top_ion_massfrac)
        / element_atomic_mass[e];
  }
  return prim[RO] * ne_per_density;  // electron number density
}


// GET ION NUMBER DENSITY ##############################################
// Calculate given ion's number density from the primitive vector.
double nemo::get_n_ion(
    string s,             ///< ion name
    const pion_flt *prim  ///< primitive state vector.
)
{
  double species_num_density;
  int species_primitive_idx;
  auto it = species_prim_index_map.find(s);
  if (it != species_prim_index_map.end())
    species_primitive_idx = it->second;
  else {
    spdlog::error(" NEMO: Species {} do not found in tracer list", s);
    exit_nemo(1);
  }

  std::string element_name = atomic_physics_data.get_element_name(s);
  int element_idx          = element_position_index_map[element_name];

  return species_num_density = prim[RO] * prim[species_primitive_idx]
                               / element_atomic_mass[element_idx];
}


// ##################################################################
double nemo::get_X_H()
{
  spdlog::error("don't call get_X_H() in NEMO");
  exit_nemo(1);
  return 0.0;
}


// GET TEMPERATURE ##################################################
// Calculate temperature from ion fraction, element number density
// and internal energy using T = E*(gamma-1)/(k_B*n_total)
double nemo::get_temperature(
    const std::vector<double> &y_species_frac,
    const std::vector<pion_flt>
        &element_number_density,  // element number density
    const double E)
{
  return gamma_minus_one * E
         / (pconst.kB() * get_ntot(y_species_frac, element_number_density));
}



// GET TOTAL NUMBER DENSITY ##################################################
// Calculate total number density from y_ion fraction and element
// number density.
double nemo::get_ntot(
    const std::vector<double> &y_species_frac,  ///< y_ion_fraction
    const std::vector<pion_flt>
        &element_number_density  ///< element number density
)
{
  double ne = 0.0;
  get_ne(y_species_frac, ne);
  double n_tot = ne;

  for (int e = 0; e < N_elements; e++) {
    n_tot += element_number_density[e];
  }
  return n_tot;
}

// #############################################################################
// calculate local species fraction from local vector
void nemo::get_y_species_frac(
    const std::vector<double> &p_local, std::vector<double> &y_species_frac)
{
  for (int s = 0; s < N_species; s++) {
    y_species_frac[s] = p_local[s];
    // Set y_ion_frac within physically acceptable range.
    y_species_frac[s] = max(NEMO_ABSTOL, y_species_frac[s]);
    y_species_frac[s] = min(1.0 - NEMO_ABSTOL, y_species_frac[s]);
  }
  return;
}


// #############################################################################
// calculate top ion fraction from local species fraction
void nemo::get_y_top_ion_frac(
    const std::vector<double> &y_species_frac,
    std::vector<double> &y_top_ion_frac)
{
  // Calculate top ion fraction from the current local vector.

  int sct = 0;
  // loop over every element
  for (int e = 0; e < N_elements; e++) {
    y_top_ion_frac[e] = 1.0;
    // loop over every species in this element
    for (int q = 0; q < element_species_count[e]; q++) {
      y_top_ion_frac[e] -= y_species_frac[sct];
      sct++;
    }

    if (y_top_ion_frac[e] < -NEMO_RELTOL) {
      spdlog::error(
          "nemo: Negative species frac !!! E {} y-top-ion {}", e,
          y_top_ion_frac[e]);
    }

    // Ensure the top-ion frac is within the acceptable range
    y_top_ion_frac[e] = max(NEMO_ABSTOL, y_top_ion_frac[e]);
    y_top_ion_frac[e] = min(1.0 - NEMO_ABSTOL, y_top_ion_frac[e]);
  }
  return;
}


// #############################################################################
// SET ELECTRON NUMBER DENSITY
void nemo::get_ne(const std::vector<double> &y_species_frac, double &ne)
{
  // Calculate electron number density
  int sct = 0;
  // loop over every element
  for (int e = 0; e < N_elements; e++) {
    double electrons_per_ions   = 0.0;
    double element_top_ion_frac = 1.0;
    // loop over every species in THIS element
    for (int q = 0; q < element_species_count[e]; q++) {
      electrons_per_ions += q * y_species_frac[sct];
      element_top_ion_frac -= y_species_frac[sct];
      sct++;
    }
    element_top_ion_frac = max(NEMO_ABSTOL, element_top_ion_frac);
    element_top_ion_frac = min(1.0 - NEMO_ABSTOL, element_top_ion_frac);

    ne += (electrons_per_ions + element_atomic_number[e] * element_top_ion_frac)
          * element_num_density[e];
  }

  // Ensure electron number density is not negative.
  if (ne < 0.0) {
    spdlog::warn(
        "nemo: Negative ne = {:12.3e}, resetting to {:12.3e}", ne, 1.0e-10);
    ne = max(1.0e-10, ne);
    spdlog::warn("y_species_frac: {}", y_species_frac);
    exit_nemo(1);
  }

  return;
}


// ##################################################################
// GET OPTICAL DEPTH
// Calculate the optical depth of the current cell for individual energy bins
void nemo::get_dtau(
    const rad_source *,       ///< pointer to radiation source struct
    const pion_flt depth,     ///< ds, thickness of the cell
    const pion_flt *prim_in,  ///< input primitive vector
    pion_flt *dtau_vec        ///< output dtau vector
)
{

  // Calculate differential optical depth for each energy bin.
  for (int bin = 0; bin < N_ebins; bin++) {
    double bin_optical_depth = 0.0;
    for (int s = 0; s < N_photo_species; s++) {
      bin_optical_depth += prim_in[photo_species_prim_idx[s]]
                           * photo_xsection_table[s][bin]
                           * bin_weight_table[s][bin]
                           / element_atomic_mass[photo_species_element[s]];
      // if (prim_in[photo_species_prim_idx[0]]>0.71 && (s==0 || s==10)) {
      //  spdlog::info("s {}, b {}, dtau
      //  {:12.3e}",s,bin,bin_optical_depth*prim_in[RO] * depth);
      //}
    }
    bin_optical_depth *= prim_in[RO] * depth;
    dtau_vec[bin] = bin_optical_depth;
    // if (prim_in[photo_species_prim_idx[0]]>0.71) {
    //  spdlog::info("b {}, dtau {:12.3e}",bin,bin_optical_depth);
    //}
  }
}


// #############################################################################
// CONVERT PRIMITIVE VECTOR TO LOCAL VECTOR
int nemo::convert_prim2local(
    const pion_flt
        *prim_in,  ///< primitive vector from grid cell (length N_prim)
    std::vector<double> &p_local,  ///< converted local vector.
    int function_flag  /// < flag to say which function called this function
)
{
  // function flag 1 -> Set_Temp()
  // function flag 2 -> TimeUpdateMP_RTnew()
  // function flag 3 -> timescales_RT()
  // function flag 4 -> Temperature()

  // Make sure mass fraction in primitive vector is a self-consistent.
  sCMA(prim_copy.data(), prim_in);
  // spdlog::info("copy {} orig {}",prim_copy[PG],prim_in[PG]);

  // Set elemental number density from the current primitive vector.
  for (int e = 0; e < N_elements; e++) {  // loop over every element
    element_num_density[e] =
        prim_copy[RO] * prim_copy[massfrac_idx[e]] / element_atomic_mass[e];
  }

  // Set local vector species variables
  // Species fractions are determined by iterating over each species,
  // where the local_vec[i] is calculated as prim_vec[i] divided by X_e.
  // In this context, 'i' represents the mass fraction of species associated
  // with the element 'e' in the primitive vectors.
  int sct = 0;
  // loop over every element
  for (int e = 0; e < N_elements; e++) {
    // loop over every species in this element
    for (int q = 0; q < element_species_count[e]; q++) {
      if (prim_copy[massfrac_idx[e]] <= 0.0)
        p_local[sct] = 0.0;
      else
        p_local[sct] =
            prim_copy[species_prim_idx[sct]] / prim_copy[massfrac_idx[e]];
      p_local[sct] = min(1.0, p_local[sct]);
      p_local[sct] = max(0.0, p_local[sct]);
      sct++;
    }
  }

  // Set internal energy density in local vector.
  p_local[E_index] = prim_copy[PG] / (gamma_minus_one);
  if (!isfinite(p_local[E_index])) {
    spdlog::error(
        "nemo: Input pressure is not finite: {} {}", prim_copy[PG],
        prim_in[PG]);
    exit_nemo(4);
  }

  return 0;
}



// #############################################################################
// CONVERT LOCAL VECTOR TO PRIMITIVE VECTOR
int nemo::convert_local2prim(
    const std::vector<double> &p_local,
    const pion_flt *prim_in,  ///< input primitive vector
    pion_flt *prim_out,       ///< updated primitive vector
    int function_flag         /// < flag to say which function called me
)
{
  // copying input primitive vector to output primitive vector.
  for (int v = 0; v < N_prim; v++)
    prim_out[v] = prim_in[v];

  // Set pressure in the primitive vector.
  prim_out[PG] = p_local[E_index] * (gamma_minus_one);
  // resizing y_species vector.
  std::vector<double> y_species_frac(N_species, 0.0);
  // set primitive species variables and define y_species_frac
  int sct = 0;                            // species counter
  for (int e = 0; e < N_elements; e++) {  // loop over every element
    for (int q = 0; q < element_species_count[e]; q++) {
      prim_out[species_prim_idx[sct]] = p_local[sct] * prim_in[massfrac_idx[e]];
      y_species_frac[sct]             = p_local[sct];
      sct++;
    }
  }


  double T =
      get_temperature(y_species_frac, element_num_density, p_local[E_index]);
  // cout <<"current temp = " << T << "\n";
  if (T > 1.01 * EP->MaxTemperature) {
    spdlog::warn("nemo::convert_local2prim() High temperature encountered!!!");
    Set_Temp(prim_out, EP->MaxTemperature, 0);
    spdlog::warn(
        "T = {} K, eint = {} erg/cm^3, ...  limiting to T = {} K", T,
        p_local[E_index], EP->MaxTemperature);
    spdlog::warn(
        "prim_in {}", std::vector<pion_flt>(prim_in, prim_in + N_prim));
  }
  if (T < 0.99 * EP->MinTemperature) {
    spdlog::warn("nemo::convert_local2prim() LOW temperature encountered!!!");
    Set_Temp(prim_out, EP->MinTemperature, 0);
    spdlog::warn(
        "T = {} K, eint = {} erg/cm^3, ...  limiting to T = {} K", T,
        p_local[E_index], EP->MinTemperature);
  }

  return 0;
}



// #############################################################################
// SET LOCAL VARIABLES
// Set local variables from current local vector.
void nemo::set_localvariables(
    const std::vector<double> &p_local,
    std::vector<double> &y_species_frac,
    std::vector<double> &y_top_ion_frac,
    double &ne)
{
  // Set y species fraction from the local vector.
  get_y_species_frac(p_local, y_species_frac);

  // Set local neutral fraction from local ion fraction. TODO: remove this
  // set_y_neutral_frac(y_ion_frac, y_neutral_frac);

  // get local top ion fraction from the local y_species_frac;
  get_y_top_ion_frac(y_species_frac, y_top_ion_frac);

  // Calculate electron number density from local ion fraction.
  get_ne(y_species_frac, ne);

  return;
}



// ##################################################################
double nemo::Temperature(
    const pion_flt *prim,  ///< primitive vector
    const double           ///< eos gamma
)
{
  //
  // Check for negative pressure/density!  If either is found, return -1.0e99.
  if (prim[RO] <= 0.0 || prim[PG] <= 0.0) {
    // cout <<"nemo::Temperature() negative rho="<<pv[RO]<<" or
    // p="<<pv[PG]<<"\n";
    return -1.0e99;
  }
  if (!isfinite(prim[PG])) {
    spdlog::error(
        "NEMO Temperature() pressure is not finite: {}",
        std::vector<pion_flt>(prim, prim + nv_prim));
    exit_nemo(1);
  }

  // generate vector of (y_s, nE, Eint) to calculate Temperature  **************
  std::vector<double> p_local(N_local, 0.0);
  // p_local.resize(N_local);
  std::vector<double> y_species_frac(N_species, 0.0);
  // y_species_frac.resize(N_species);

  convert_prim2local(prim, p_local, 4);

  int sct = 0;
  for (int e = 0; e < N_elements; e++) {  // loop over every element
    // loop over every species in THIS element
    for (int q = 0; q < element_species_count[e]; q++) {
      if (prim[massfrac_idx[e]] <= 0.0) {
        y_species_frac[sct] = 0.0;
      }
      else {
        y_species_frac[sct] =
            prim[species_prim_idx[sct]] / prim[massfrac_idx[e]];
      }
      sct++;
    }
  }
  double T =
      get_temperature(y_species_frac, element_num_density, p_local[E_index]);
  // if (!isfinite(T)) {
  //   spdlog::error("Temperature {}",p_local);
  //   spdlog::error("{}",y_species_frac);
  //   spdlog::error("{}",element_num_density);
  //   exit_nemo(10);
  // }
  return get_temperature(y_species_frac, element_num_density, p_local[E_index]);
}



// ##################################################################
int nemo::Set_Temp(
    pion_flt *prim,  ///< primitive vector.
    const double T,  ///< temperature
    const double     ///< eos gamma.
)
{

  // Check for negative pressure.  If density<0 then we should bug
  // out because there is no way to set a temperature, but if p<0 we
  // can just overwrite it.

  if (prim[PG] <= 0.0 || !isfinite(prim[PG])) {
    spdlog::debug("nemo::Set_Temp() correcting negative pressure");
    prim[PG] = 1.0e-12;  // Need p>0 for prim-to-local conversion.
  }

  std::vector<double> p_local(N_local, 0.0);
  // p_local.resize(N_local);
  std::vector<double> y_species_frac(N_species, 0.0);
  // y_species_frac.resize(N_species);

  int err = convert_prim2local(prim, p_local, 1);

  // Determine y_ion_frac from the primitive vector
  int sct = 0;
  for (int e = 0; e < N_elements; e++) {  // loop over every element
    // loop over every species in THIS element
    for (int q = 0; q < element_species_count[e]; q++) {
      if (prim[massfrac_idx[e]] <= 0.0) {
        y_species_frac[sct] = 0.0;
      }
      else {
        y_species_frac[sct] =
            prim[species_prim_idx[sct]] / prim[massfrac_idx[e]];
      }
      sct++;
    }
  }

  // Determine internal energy using get_ntot
  p_local[E_index] = get_ntot(y_species_frac, element_num_density) * pconst.kB()
                     * T / gamma_minus_one;

  // Call convert_local2prim with the new local vector; this will
  // generate a new temperature value;
  err += convert_local2prim(p_local, prim, prim, 1);
  return err;
}



// ##################################################################
// TIME UPDATE MP
int nemo::TimeUpdateMP(
    const pion_flt *p_in,  ///< Primitive Vector to be updated.
    pion_flt *p_out,       ///< Destination Vector for updated values.
    const double dt,       ///< Time Step to advance by.
    const double,          ///< EOS gamma.
    const int,             ///< Switch for what type of integration to use.
    double *random_stuff   ///< Vector of extra data (column densities, etc.).
)
{
  //
  // Call the new update function, but with zero radiation sources.
  //
  std::vector<struct rt_source_data> temp;
  int err =
      TimeUpdateMP_RTnew(p_in, 0, temp, 0, temp, p_out, dt, 0, 0, random_stuff);
  return err;
}
// END OF TIME UPDATE MP ############################################



// ##################################################################
// TIME UPDATE MP RT NEW
int nemo::TimeUpdateMP_RTnew(
    const pion_flt *p_in,  ///< Primitive Vector to be updated.
    const int,             ///< unused.
    const std::vector<struct rt_source_data> &,  ///< unused.
    const int N_ion,  ///< number of ionising radiation sources.
    std::vector<struct rt_source_data> &ion_src,
    ///< list of ionising src column densities and source properties.
    pion_flt *p_out,      ///< Destination Vector for updated values
                          ///< (can be same as first Vector.
    const double dt,      ///< Time Step to advance by.
    const double,         ///< EOS gamma.
    const int,            ///< Switch (unused)
    double *random_stuff  ///< final temperature (debugging).
)
{
  int err = 0;
  std::vector<double> p_local(N_local, 0.0);
  err = convert_prim2local(p_in, p_local, 2);
  if (err) {
    spdlog::error("{}: {}", "Bad input state to nemo::TimeUpdateMP()", err);
    exit_nemo(1);
  }

  err = correct_localvector(p_local, 1);
  if (err) {
    spdlog::error("Bad input vector to TimeUpdateMP_RTnew");
    exit(1);
  }

  // issue is here, local vector P and NV_Ith_S(y_in, v) are different

  // setup_radiation_source_parameters(p_in, P, ion_src);
  setup_radiation_source_parameters(N_ion, ion_src);

  // Populates CVODE vector with initial conditions (input)
  for (int v = 0; v < N_local; v++) {
    NV_Ith_S(y_in, v) = p_local[v];
    // spdlog::info("v {} y_in {}", v, NV_Ith_S(y_in, v));
  }

  //============
  // Calculate y-dot[] to see if anything is changing significantly over dt
  // Here y_out is y_dot calculated from the ydot function.
  // err = ydot(0, y_in, y_out, 0);
  // if (err) {
  //  spdlog::error(
  //      "{}: {}", "dYdt() returned an error in nemo::TimeUpdateMP_RTnew()",
  //      err);
  //  exit(1);
  //}
  //============

  err = integrate_cvode_step(y_in, 0, 0.0, dt, y_out);
  if (err) {
    spdlog::error(
        "integration failed: nemo::TimeUpdateMP_RTnew() error {}", err);
    for (int v = 0; v < N_local; v++)
      p_local[v] = NV_Ith_S(y_in, v);
    spdlog::error("y_in: {}", p_local);
    for (int v = 0; v < N_local; v++)
      p_local[v] = NV_Ith_S(y_out, v);
    spdlog::error("y_out: {}", p_local);
    exit(1);
  }

  //
  // Now put the result into p_out[] and return.
  //
  for (int v = 0; v < N_local; v++)
    p_local[v] = NV_Ith_S(y_out, v);

  err = correct_localvector(p_local, 2);
  if (err) {
    spdlog::error("bad output vector from TimeUpdateMP_RTnew");
    exit(2);
  }
  err = convert_local2prim(p_local, p_in, p_out, 2);



#ifdef TEST_INF
  for (int v = 0; v < N_prim; v++) {
    if (!isfinite(p_in[v]) || !isfinite(p_out[v])) {
      spdlog::debug("NAN in MPv3 update: {}", v);
      spdlog::debug("Pin  : {}", std::vector<double>(p_in, p_in + N_prim));
      spdlog::debug("Pout : {}", std::vector<double>(p_out, p_out + N_prim));
      spdlog::debug("Ploc  : {}", p_local);
      // spdlog::error("{}: {}", "NAN in MPv3",P[2]);
      exit_nemo(v);
    }
  }
#endif

  return err;
}



// CORRECT LOCAL VECTOR (CLV) #####################################
int nemo::correct_localvector(
    std::vector<double> &p_local,  ///< check the local vector
    int function_flag              /// flags the function
)
{
  // function_flag 1 -> towards begining of TimeUpdateMP_RTnew
  // function_flag 2 -> towards the end of TimeUpdateMP_RTnew
  // function_flag 3 -> timescales_RT
  // function_flag 4 -> ydot


  int flag = 0;

  // *********************************************************************
  // Check if y_ion_fracs are within the acceptable range
  std::vector<double> y_species_frac(N_species);

  for (int s = 0; s < N_species; s++) {  // loop over ions in p_local
    y_species_frac[s] = p_local[s];
    y_species_frac[s] = max(NEMO_ABSTOL, y_species_frac[s]);
    y_species_frac[s] = min(1.0 - NEMO_ABSTOL, y_species_frac[s]);
  }
  // *********************************************************************


  // *********************************************************************
  // Scaling the y_species_frac if the total elemental fraction exceed 1.

  // 1. Calculating total ion fraction for each element.
  std::vector<double> y_frac_sum(N_elements, 0.0);
  int sct = 0;
  for (int e = 0; e < N_elements; e++) {
    for (int q = 0; q < element_species_count[e]; q++) {
      y_frac_sum[e] += y_species_frac[species_local_idx[sct]];
      sct++;
    }
  }


  // 2. Rescale the y_ion_frac if total ion fraction is > 1
  sct = 0;
  for (int e = 0; e < N_elements; e++) {
    // spdlog::info("correct localvec {}: checking el {}",function_flag, elem);
    if (1.0 - y_frac_sum[e] < NEMO_ABSTOL) {
      // if total ion fraction is nearly one or more than one,
      // the neutral fraction must be minumum (say, NEMO_ABSTOL)
      if (y_frac_sum[e] - 1.0 > 1.0e4 * NEMO_RELTOL && function_flag != 4) {
        spdlog::warn(
            "nemo::CorrLocalVec flag {}: too many species! element {}, y-1 = {:12.3e}",
            function_flag, e, y_frac_sum[e] - 1);
        if (e == 7)
          spdlog::warn(
              "y {}",
              std::vector<double>(
                  y_species_frac.begin() + 48, y_species_frac.begin() + 62));
        else
          spdlog::warn("y {}", y_species_frac);
      }
      // double total_fraction = y_total_ion_frac[elem] + 1.0e-30;
      for (int q = 0; q < element_species_count[e]; q++) {
        y_species_frac[sct] /= y_frac_sum[e] + NEMO_ABSTOL;  // total_fraction;
        sct++;
      }
    }
    else {
      sct += element_species_count[e];
    }
  }
  // **************************************************************

  // **************************************************************
  // Feeding the corrected y_ion_frac back in to p_local
  for (int s = 0; s < N_species; s++)
    p_local[s] = y_species_frac[s];
  // **************************************************************

  // **************************************************************
  // Check for negative pressure
  // Note: This shouldn't happen, so we output a warning) and set to
  // 10K if we find it.
  if (p_local[E_index] <= 0.0) {
    if (function_flag != 4)
      spdlog::warn(
          "nemo::CorrLocalVec flag: {} - Negative pressure input: Eint = {}, "
          "setting to {} K",
          function_flag, p_local[E_index], EP->MinTemperature);

    // reset the internal energy (requires using y_ion_frac in get_ntot)
    p_local[E_index] = get_ntot(y_species_frac, element_num_density)
                       * pconst.kB() * EP->MinTemperature / (gamma_minus_one);
  }
  // **************************************************************


  // **************************************************************
  // if the temperature is below Minimum Temperature, set it back to
  // Minimum Temperature by set the internal energy accordingly
  double T =
      get_temperature(y_species_frac, element_num_density, p_local[E_index]);

  if (T < 0.95 * EP->MinTemperature) {
    if (function_flag != 4)
      spdlog::debug(
          "nemo::CorrLocalVec flag: {} - Temperature = {} K, below T_MIN = {} K",
          function_flag, T, EP->MinTemperature);

    // Resetting the internal energy with Minimum Temperature
    p_local[E_index] = get_ntot(y_species_frac, element_num_density)
                       * pconst.kB() * EP->MinTemperature / (gamma_minus_one);
  }
  // **************************************************************

  return flag;
}
// END OF CORRECT LOCAL VECTOR ######################################



// ##################################################################
double nemo::timescales(
    const pion_flt *prim_in,  ///< Current cell state vector.
    const double,             ///< EOS gamma.
    const bool,               ///< set to 'true' if including cooling time.
    const bool,  ///< set to 'true' if including recombination time.
    const bool   ///< set to 'true' if including photo-ionsation time.
)
{
#ifdef NEMO_DEBUG
  if (RS->Nsources != 0) {
    spdlog::debug("WARNING: nemo::timescales() using non-RT version!");
  }
#endif  // NEMO_DEBUG
  std::vector<struct rt_source_data> temp;
  double tmin = timescales_RT(prim_in, 0, temp, 0, temp, 0.0);
  temp.clear();
  return tmin;
}



// TIME SCALE RT ###########################################################
/// This returns the minimum timescale of all microphysical processes,
/// including reaction times for each species and the total heating/cooling
/// time for the gas. It requires the radiation field as an input, so it has
/// substantially greater capability than the other timescales function.
/// Default setting is DT02, which limits by 0.25/ydot (and not by E/Edot)
///
double nemo::timescales_RT(
    const pion_flt *prim_in,  ///< Current cell state vector.
    const int N_heat,         ///< Number of UV heating sources.
    const std::vector<struct rt_source_data> &heat_src,
    ///< list of UV-heating column densities and source properties.
    const int N_ion,  ///< number of ionising radiation sources.
    const std::vector<struct rt_source_data> &ion_src,
    ///< list of ionising src column densities and source properties.
    const double  ///< EOS gamma.
)
{
  // spdlog::info("timestep limit: {}", EP->MP_timestep_limit);
  if (!EP->MP_timestep_limit) return HUGEVALUE;

  // spdlog::info("got past blockage 0");
  // Note: Minimum value of micro-physics trace variables in the primitive
  // vector is set to 1e-12.
  int err = 0;
  std::vector<double> p_local;
  p_local.resize(N_local, 0.0);

  // First convert to local variables.
  err = convert_prim2local(prim_in, p_local, 3);
  if (err) {
    spdlog::error("{}: {}", "Bad input state to nemo::timescales_RT()", err);
    exit_nemo(1);
  }

  // correct local vector
  err = correct_localvector(p_local, 3);

  // spdlog::debug("after correcting - p_local = {}", p_local);
  if (err) {
    spdlog::error("bad input vector to timescales_RT");
    exit_nemo(1);
  }

  // v runs over all species and internal energy
  for (int v = 0; v < N_local; v++) {
    NV_Ith_S(y_in, v) = p_local[v];
  }

  // Next set the radiation properties of the current cell.
  setup_radiation_source_parameters(N_ion, ion_src);

  // Now calculate y-dot[]...
  // ydot function returns RHS of ydot equation. Here y_out is RHS of ydot.
  err = ydot(0, y_in, y_out, 0);
  if (err) {
    spdlog::error(
        "{}: {}", "dYdt() returned an error in nemo::timescales_RT()", err);
    exit_nemo(1);
  }

  // spdlog::info("got past blockage 2");
  // double tt = Temperature(prim_in, 0);
  // spdlog::info("{:12.6e}  {:12.6e}", tt, NV_Ith_S(y_out,E_index));
  // spdlog::info("got past blockage 3");

  // And finally get the smallest timescale over which things are varying.
  // get the ionisation timescale, limited to dt = 0.25/|ydot|.
  double t = HUGEVALUE;
  std::vector<double> Y_dot;
  Y_dot.resize(N_local, 0.0);
  for (int v = 0; v < N_local; v++)
    Y_dot[v] = NV_Ith_S(y_out, v);

  for (int v = 0; v < N_equations; v++) {
    // t = min(t, DTFRAC / (fabs(NV_Ith_S(y_out, v)) + TINYVALUE));
    t = min(t, DTFRAC / (fabs(Y_dot[v]) + TINYVALUE));
  }

  return t;
}



// ##################################################################
void nemo::sCMA(
    pion_flt *corr,          ///< corrected vector.
    const pion_flt *prim_in  ///< input primitive vector from grid cell.
)
{
  // "Consistent multi-species advection" correction step from
  // Plewa and Muller, Astron. Astrophys. 342, 179 (1999).

  // This method ensures that the inputted mass fraction satisfies
  // the constraint that the sum of elemental mass fractions must equal unity
  // Similarly the sum of all species fractions of an element must not exceed
  // the total elemental abundance.


  // bool warn = false;
  std::vector<double> prim(N_prim);
  for (int i = 0; i < N_prim; i++)
    corr[i] = 1.0;
  for (int i = 0; i < N_prim; i++)
    prim[i] = prim_in[i];

  // first ensure elemental abundances are within allowed range
  for (int e = 0; e < N_elements; e++) {
    // if (prim[massfrac_idx[e]] < 0.0 || prim[massfrac_idx[e]] > 1.0) {
    //  spdlog::warn("element {} out of range: {}", e, prim[massfrac_idx[e]]);
    //  spdlog::info("prim: {}",std::vector<pion_flt>(prim,prim+N_prim));
    //  warn=true;
    //}
    // spdlog::info("e {}, el {}, val {:12.3e}, N_prim
    // {}",e,massfrac_idx[e],prim[massfrac_idx[e]], N_prim);
    prim[massfrac_idx[e]] = min(1.0, max(1.0e-20, prim[massfrac_idx[e]]));
  }
  // Then apply sCMA criterion to elements
  double element_massfrac = 0.0;
  for (int e = 0; e < N_elements; e++) {
    element_massfrac += prim[massfrac_idx[e]];
  }
  double element_correction = 1.0 / element_massfrac;

  element_massfrac = 0.0;  // test
  for (int e = 0; e < N_elements; e++) {
    prim[massfrac_idx[e]] *= element_correction;
    element_massfrac += prim[massfrac_idx[e]];  // test
  }
  if (!pconst.equalD(element_massfrac, 1.0))  // test
    spdlog::warn(
        "renormalisation of element mass fracs failed: {}",
        element_massfrac);  // test


  // now for each element, ensure species fractions of each ion are within
  // allowed range
  int sct = 0, s = 0;
  for (int e = 0; e < N_elements; e++) {
    for (int q = 0; q < element_species_count[e]; q++) {
      // if (prim[species_prim_idx[sct]] < 0.0 ||
      //    prim[species_prim_idx[sct]] > prim[massfrac_idx[e]]*(1.0+1e-12)) {
      //  spdlog::warn("species {},{} out of range: {}, {}", e,sct,
      //               prim[massfrac_idx[e]],prim[species_prim_idx[sct]]);
      //  spdlog::info("prim: {}",std::vector<pion_flt>(prim,prim+N_prim));
      //  warn=true;
      //}

      prim[species_prim_idx[sct]] =
          min(prim[massfrac_idx[e]], max(1.0e-20, prim[species_prim_idx[sct]]));
      sct++;
    }
  }

  // loop over every element
  sct                = 0;
  double species_sum = 0.0, species_corr = 1.0;
  for (int e = 0; e < N_elements; e++) {
    species_sum = 0.0;
    s           = sct;
    // loop over every species in this element to get sum
    for (int q = 0; q < element_species_count[e]; q++) {
      species_sum += prim[species_prim_idx[s]];
      s++;
    }
    // corrector is X_ei / sum(Y_ei,j), but we didn't count the top ion, so we
    // only correct if the sum is too large.
    species_corr = min(1.0, prim[massfrac_idx[e]] / species_sum);
    // apply corrector to each species.
    species_sum = 0.0;  // test
    for (int q = 0; q < element_species_count[e]; q++) {
      prim[species_prim_idx[sct]] *= species_corr;
      species_sum += prim[species_prim_idx[sct]];  // test
      sct++;
    }
    if (species_sum > prim[massfrac_idx[e]] * (1.0 + 1e-12))  // test
      spdlog::warn(
          "renormalisation of species mass fracs failed: {} {} {}", e,
          prim[massfrac_idx[e]], species_sum);  // test
  }

  for (int i = 0; i < N_prim; i++)
    corr[i] = prim[i];
  // if (warn) {
  //  spdlog::warn("before correction
  //  {}",std::vector<pion_flt>(prim_in,prim_in+N_prim)); spdlog::warn("after
  //  correction {}",std::vector<pion_flt>(corr,corr+N_prim));
  //}

  return;
}



// #######################################################################
// YDOT FUNCTION
// This function calculate the RHS of ydot equations.
int nemo::ydot(
    double,                ///< current time (UNUSED)
    const N_Vector y_now,  ///< current Y-value
    N_Vector y_dot,        ///< vector for Y-dot values
    const double *         ///< extra user-data vector (UNUSED)
)
{

  // SECTION: Determining the following from the current local vector **********
  // y_species_frac (ion fractions in the local vector)
  std::vector<double> y_species_frac;
  y_species_frac.resize(N_species, 0.0);
  // top ion fraction of elements in the tracer list.
  std::vector<double> y_top_ion_frac;
  y_top_ion_frac.resize(N_elements, 0.0);
  // ne (electron number density)
  double ne = 0.0;
  // local vector
  std::vector<double> p_local;
  p_local.resize(N_local, 0.0);
  // make a copy of y_now into p_local
  for (int v = 0; v < N_local; v++) {
    p_local[v] = NV_Ith_S(y_now, v);
  }

  correct_localvector(p_local, 4);

  // Set local quantities: y_ion_frac, y_neutral_frac and ne
  set_localvariables(p_local, y_species_frac, y_top_ion_frac, ne);
  // END OF SECTION ************************************************************


#ifdef NEMO_CIE_TEST
  // Restricting ne value not less than 1e-3 to perform CIE test.
  ne = max(1e-3, ne);
#endif

  // Initialise all elements of ydot vector to zero.
  for (int v = 0; v < N_equations; v++)
    NV_Ith_S(y_dot, v) = 0.0;

  // Initializing Edot to zero
  double Edot = 0.0;

  // Get internal energy from the current local vector (y_now)
  double E_in = p_local[E_index];

  // Calculate temperature from E_in, y_ion_frac and number density.
  double T = get_temperature(y_species_frac, element_num_density, E_in);
#ifdef CCCT1  // Iliev et al. (2006) test 1 imposes constant temperature.
  T = 1.0e4;
#endif

  // Restricting temperature within T_min and T_max
  if (T > T_max)
    T = T_max;
  else if (T < T_min)
    T = T_min;


  // SECTION: Species number density *******************************************
  int sct = 0;
  for (int e = 0; e < all_species_list.size(); e++) {
    // calculate the excess fraction corresponding to each element.
    double excess_elemental_frac = 0.0;
    std::vector<double> species_frac;
    for (int q = 0; q < all_species_list[e].size(); q++) {
      // Calculating the current number density of neutral atoms.
      if (q == element_atomic_number[e]) {
        species_number_density[e][q] =
            y_top_ion_frac[e] * element_num_density[e];
        excess_elemental_frac += y_top_ion_frac[e];
        species_frac.push_back(y_top_ion_frac[e]);
      }
      // Calculating the current number density of ions.
      else {
        species_number_density[e][q] =
            y_species_frac[sct] * element_num_density[e];
        excess_elemental_frac += y_species_frac[sct];
        species_frac.push_back(y_species_frac[sct]);
        sct++;
      }
    }

    if (fabs(excess_elemental_frac - 1.0) > 1e-7) {
      spdlog::error(
          "species counting is wrong {}: species {}, \n\t\t values {}",
          excess_elemental_frac, cooling_species[e], species_frac);
      if (fabs(excess_elemental_frac - 1.0) > 1e-5) {
        exit(1);
      }
    }
  }
  // END OF SECTION: Species number density ************************************


  // SECTION: Cosmic-ray Ionisation ********************************************
  if (EP->coll_ionisation) {
    double cosmicray_ydot = 0.0;
    for (int e = 0; e < N_elements; e++) {
      cosmicray_ydot =
          cosmicray_ionise_rate * p_local[neutral_species_local_idx[e]];

      // if the upper ionisation state exist then add cosmicray_ydot
      if (upper_ion_local_idx[neutral_species_local_idx[e]] != -1) {
        NV_Ith_S(y_dot, neutral_species_local_idx[e] + 1) += cosmicray_ydot;
      }
      // otherwise
      NV_Ith_S(y_dot, neutral_species_local_idx[e]) -= cosmicray_ydot;

      // Heating due to cosmic-ray ionisation of neutral
      // species neutral_species_local_idx[e]]
      Edot += cosmic_ray_energy * cosmicray_ydot * element_num_density[e];
    }
  }
  // END OF SECTION: Cosmic-ray Ionisation *************************************


  // SECTION: Collisional Ionisation *******************************************
  if (EP->coll_ionisation) {
    // Calculating RHS terms corresponding to collisional ionisation for
    // species in collion_species. These terms are calculated only once and
    // attached at relevant ydot equations.
    //
    // For example, for species s, the term ci_rate(s)y(s)ne is calculated
    // and substracted for ydot() equation while added to the ydot(s+1)
    // equation if the specie s+1 is the next higher ionisation level of
    // the species s.
    //
    // Since collisional ionization species share the same index as local
    // species, we leverage this similarity and utilize the local species index.

    // collisional ionisation rate and its corresponding ydot terms
    double coll_ionise_rate      = 0.0;
    double coll_ionise_ydot      = 0.0;
    double cx_ionise_energy_rate = 0.0;
    for (size_t s = 0; s < collion_species.size(); s++)  // loop over ci species
    {
      // colllisional ionisaation rate
      coll_ionise_rate = collisional_ionisation_rate(s, T);
      coll_ionise_ydot = coll_ionise_rate * p_local[s] * ne;

      // collisional ionisation term in ydot
      // if the upper ionisation state exist then add coll_ionise_ydot
      if (upper_ion_local_idx[s] != -1)
        NV_Ith_S(y_dot, s + 1) += coll_ionise_ydot;
      // otherwise
      NV_Ith_S(y_dot, s) -= coll_ionise_ydot;

      // Cooling due to collisional ionisation of species s.
      Edot -= collion_species_IE[s] * coll_ionise_ydot
              * element_num_density[species_element[s]];
    }
  }
  // END OF SECTION: Collisional Ionisation ************************************

  // SECTION: Recombination (Radiative + Dielectronic) *************************
  double rr_cool = 0.0;
  double cx_cool = 0.0;
  if (EP->rad_recombination) {
    // RHS terms corresponding to Recombination (R) are calculated for
    // all species. In this case, for specie s, the term rec_rate(s)y(s)ne is
    // subtracted for ydot(s-1) and the same terms is added for the ydot(s),
    // provided species s is a ion of same element. Whereas
    // if the species s is neutral, then only addition performed.

    // Radiative + Di-electronic recombination rate %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    double radiative_recomb_rate = 0.0;
    double radiative_recomb_ydot = 0.0;
    double caseA                 = 0.0;
    // std::vector<double> temp_rr(recomb_tracer_list.size());  //DEBUG

    for (size_t s = 0; s < recomb_species.size(); s++) {
      radiative_recomb_rate = recombination_rate(s, T, caseA);
      // radiative recombination terms
      // If s is a top level ion of an element then
      if (recomb_top_ion_bool[s]) {
        radiative_recomb_ydot =
            caseA * ne * y_top_ion_frac[recomb_species_element[s]];
        NV_Ith_S(y_dot, s) += radiative_recomb_ydot;
      }
      // otherwise
      else {
        radiative_recomb_ydot = caseA * ne * p_local[s + 1];
        NV_Ith_S(y_dot, s) += radiative_recomb_ydot;
        NV_Ith_S(y_dot, s + 1) -= radiative_recomb_ydot;
      }

      // Cooling due to radiative-recombination of this species
      // The factor 0.9 is from Hummer (1994,MNRAS,268,109), who found that
      // it should be nearly 1 at low temp, and around 0.9 at 1e4 K.
      // It should be lower at high temperature, but here free-free is
      // dominating anyway so it is not important.
      rr_cool = 0.9 * T * pconst.kB() * radiative_recomb_ydot
                * element_num_density[recomb_species_element[s]];
      Edot -= rr_cool;

      // on-the-spot approximation hack
#ifdef OTS_HACK
      // HACK!  Hard-code on-the-spot approx for Lexington benchmark.
      // on-the-spot approximation for He+: results in ionization of H
      double xsh = 0.0, xshe = 0.0, xshe2 = 0.0, fh = 0.0, fhe = 0.0,
             fhe2 = 0.0, ndotgamma = 0.0;
      if (s == 1) {
        // Recombination of He+ to ground state result in photons that can
        // ionize H0 and He0, with photon energy >24.6eV. x-sections from module
        // output, these variables are n_i * sigma_i for neutral H and neutral
        // He
        xsh = 9.551e-19 * element_num_density[recomb_species_element[0]]
              * p_local[0];
        xshe = 6.3372e-18 * element_num_density[recomb_species_element[1]]
               * p_local[1];
        fh  = xsh / (xsh + xshe);
        fhe = xshe / (xsh + xshe);
        ndotgamma =
            (caseA - radiative_recomb_rate) * p_local[s + 1] * ne
            * element_num_density[recomb_species_element[1]];  // photons per
                                                               // cm3 per s

        // H ionizations and heating
        NV_Ith_S(y_dot, 0) -=
            ndotgamma * fh / element_num_density[recomb_species_element[0]];
        Edot += ndotgamma * fh * (24.6 - 13.6) * pconst.eV();

        // He ionizations and heating
        NV_Ith_S(y_dot, 1) -=
            ndotgamma * fhe / element_num_density[recomb_species_element[1]];
        NV_Ith_S(y_dot, 2) +=
            ndotgamma * fhe / element_num_density[recomb_species_element[1]];
        Edot += 0.0;  // no net heating

        // 2.B - case B He+ recombinations can ionize H
        fh = 1.0;
        ndotgamma =
            radiative_recomb_rate * p_local[s + 1] * ne
            * element_num_density[recomb_species_element[1]];  // photons per
                                                               // cm3 per s
        NV_Ith_S(y_dot, 0) -=
            ndotgamma * fh / element_num_density[recomb_species_element[0]];
        Edot += ndotgamma * fh * (18.15 - 13.6) * pconst.eV();
      }
      else if (s == 2) {
        // Recombination of He2+ can ionize He+ and H.
        // ** Case 1 **
        // recombinations to ground state from 54.4eV
        xsh = 9.7296e-20 * element_num_density[recomb_species_element[0]]
              * p_local[0];
        xshe = 1.4245e-18 * element_num_density[recomb_species_element[1]]
               * p_local[1];
        xshe2 = 1.2841e-18 * element_num_density[recomb_species_element[2]]
                * p_local[2];
        fh   = xsh / (xsh + xshe + xshe2);
        fhe  = xshe / (xsh + xshe + xshe2);
        fhe2 = xshe2 / (xsh + xshe + xshe2);
        ndotgamma =
            (caseA - radiative_recomb_rate) * p_local[s + 1] * ne
            * element_num_density[recomb_species_element[2]];  // photons per
                                                               // cm3 per s
        // H ionizations and heating
        NV_Ith_S(y_dot, 0) -=
            ndotgamma * fh / element_num_density[recomb_species_element[0]];
        Edot += ndotgamma * fh * (54.4 - 13.6) * pconst.eV();
        // He ionizations and heating
        NV_Ith_S(y_dot, 1) -=
            ndotgamma * fhe / element_num_density[recomb_species_element[1]];
        NV_Ith_S(y_dot, 2) +=
            ndotgamma * fhe / element_num_density[recomb_species_element[1]];
        Edot += ndotgamma * fhe * (54.4 - 24.6) * pconst.eV();
        // He+ ionizations and heating
        NV_Ith_S(y_dot, 2) -=
            ndotgamma * fhe2 / element_num_density[recomb_species_element[2]];
        Edot += 0.0;  // re-ionization of He+ is energy neutral

        // ** Case B **
        // x-sections at 40.8 eV
        xsh = 2.4664e-19 * element_num_density[recomb_species_element[0]]
              * p_local[0];
        xshe = 2.738e-18 * element_num_density[recomb_species_element[1]]
               * p_local[1];
        fh  = xsh / (xsh + xshe);
        fhe = xshe / (xsh + xshe);
        ndotgamma =
            radiative_recomb_rate * p_local[s + 1] * ne
            * element_num_density[recomb_species_element[2]];  // photons per
                                                               // cm3 per s

        // H ionizations and heating
        NV_Ith_S(y_dot, 0) -=
            ndotgamma * fh / element_num_density[recomb_species_element[0]];
        Edot += ndotgamma * fh * (40.8 - 13.6) * pconst.eV();

        // He ionizations and heating
        NV_Ith_S(y_dot, 1) -=
            ndotgamma * fhe / element_num_density[recomb_species_element[1]];
        NV_Ith_S(y_dot, 2) +=
            ndotgamma * fhe / element_num_density[recomb_species_element[1]];
        Edot += ndotgamma * fhe * (40.8 - 24.6) * pconst.eV();
      }
#endif
    }
  }
  // END OF SECTION: Recombination *********************************************

  // SECTION: Charge Exchange Reactions ****************************************
  if (charge_exchange_enable) {

    // SUBSECTION: Charge Exchange Ionisation Reactions ************************
    for (size_t s = 0; s < collion_species.size(); s++)  // loop over ci species
    {
      if (has_ionise_cx_exist[s]) {
        double cx_ionise_rate        = 0.0;
        double cx_ionise_ydot        = 0.0;
        double cx_ionise_energy_rate = 0.0;
        // spdlog::debug(" s={} reaction {}", s, recomb_cx_reactions[s]);
        // Calculating the rate and energy rate coefficient for the
        // reaction of species s with other reactant r.
        for (size_t r = 0; r < ionise_cx_reactions[s].size(); r++) {
          double ionise_factor =
              species_number_density[ionise_cx_second_reactant_element[s][r]]
                                    [ionise_cx_second_reactant_charge[s][r]]
              * ionise_charge_exchange_rate(s, r, T);
          cx_ionise_rate += ionise_factor;
          cx_ionise_energy_rate += ionise_cx_energy[s][r] * ionise_factor;
        }

        // Calculating rhs of ydot for ionisation charge exchange reactions
        // of the species s.
        cx_ionise_ydot = cx_ionise_rate * p_local[s];
        // Inserting rhs
        // if the upper ionisation state exist for the species element
        // in the local vector then add cx_ionise_ydot
        if (upper_ion_local_idx[s] != -1) {
          NV_Ith_S(y_dot, s + 1) += cx_ionise_ydot;
        }
        // always subtract cx_ionise_ydot
        NV_Ith_S(y_dot, s) -= cx_ionise_ydot;

        // Cooling due to charge exchange ionisation reactions for the species
        // s.
        Edot += 0.5 * cx_ionise_energy_rate
                * species_number_density[species_element[s]]
                                        [collion_species_charge[s]];
      }
    }
    // END OF SUBSECTION: Charge Exchange Ionisation Reactions *****************

    // SUBSECTION: Charge Exchange Recombination Reactions *********************
    for (size_t s = 0; s < recomb_species.size(); s++) {

      if (has_recomb_cx_exist[s]) {
        double cx_recomb_rate        = 0.0;
        double cx_recomb_ydot        = 0.0;
        double cx_recomb_energy_rate = 0.0;
        // Calculating the rate and energy rate coefficient for the
        // reactions of species s with other reactant r.
        // spdlog::debug(" s={} reaction {}", s, recomb_cx_reactions[s]);
        for (size_t r = 0; r < recomb_cx_reactions[s].size(); r++) {
          // spdlog::debug("r={}, rate({}) = {}", r, T,
          // recombining_charge_exchange_rate(s, r, T));
          double recomb_factor =
              species_number_density[recomb_cx_second_reactant_element[s][r]]
                                    [recomb_cx_second_reactant_charge[s][r]]
              * recomb_charge_exchange_rate(s, r, T);
          cx_recomb_rate += recomb_factor;
          cx_recomb_energy_rate += recomb_cx_energy[s][r] * recomb_factor;
        }

        // Inserting the rhs of ydot due to recombination charge exchange
        // reactions of the species s.
        // If s is a top level ion of an element then
        if (recomb_top_ion_bool[s]) {
          cx_recomb_ydot =
              cx_recomb_rate * y_top_ion_frac[recomb_species_element[s]];
          NV_Ith_S(y_dot, s) += cx_recomb_ydot;
        }
        // otherwise
        else {
          cx_recomb_ydot = cx_recomb_rate * p_local[s + 1];
          NV_Ith_S(y_dot, s) += cx_recomb_ydot;
          NV_Ith_S(y_dot, s + 1) -= cx_recomb_ydot;
        }

        // Cooling due to charge exchange recombination reactions
        Edot += 0.5 * cx_recomb_energy_rate
                * species_number_density[recomb_species_element[s]]
                                        [recomb_species_charge[s]];
      }
    }
    // END OF SUBSECTION: Charge Exchange Recombination Reactions **************
  }
  // END OF SECTION: Charge Exchange Reactions *********************************

  // SECTION: Photoionisation **************************************************
  double photoheat_rate = 0.0;
  if (EP->phot_ionisation) {

    // Calculating RHS terms corresponding to photoionisation for species in
    // photo_tracer_list. These terms are calculated only once and attached
    // at relevant ydot equations.
    //
    // For example, for species s, if ion, the term pi_rate(s)y(s) is
    // subtracted for ydot(s) while the same terms is added for the
    // ydot(s+1) equation.
    // Whereas if the species s is neutral, then pi_rate(s)y(s) is added to
    // ydot(s+1) equation.

    double photoionise_rate = 0.0;
    double photoionise_ydot = 0.0;

    for (int s = 0; s < N_photo_species; s++)  // loop over ci species
    {

      // get photo species number density
      double species_num_density =
          species_number_density[photo_species_element[s]]
                                [photo_species_charge[s]];

      // get photo-ionisation rate for the species s
      photoionise_rate = photoionisation_rate(s, species_num_density);
      photoionise_ydot = photoionise_rate;

      // if the upper ionisation state exist then add photoionise_ydot
      if (photo_upper_ion_local_idx[s] != -1)
        NV_Ith_S(y_dot, photo_species_local_idx[s] + 1) += photoionise_ydot;
      // otherwise
      NV_Ith_S(y_dot, photo_species_local_idx[s]) -= photoionise_ydot;

      // Heating due to photo-ionisation of this species.
      photoheat_rate += photoheating_rate(s, species_num_density);
      // if (photoheating_rate(s, species_num_density) < 0.0) {
      //  spdlog::info("photoheating rate: {} {:12.6e} {:12.6e}",s,
      //  species_num_density, photoheating_rate(s, species_num_density));
      //}
      // if (p_local[0] > 0.99) {
      //  spdlog::info("s {} {:12.3e} {:12.3e}, pir {:12.6e}", s,
      //  p_local[photo_species_local_idx[s]], species_num_density,
      //  photoionise_rate);
      //}
    }
    Edot += photoheat_rate;
    // spdlog::info("photoheat_rate: {:12.6e}",photoheat_rate);
  }
  // END OF SECTION: Photoionisation *******************************************


  // SECTION: Cooling Function *************************************************
  // Edot source term associated with Cooling Function is calculated.
  double Lambda         = 0.0;
  double L              = 0.0;
  double number_density = 0.0;

#ifdef NEMO_CIE_COOLING
  cooling_outfile.open("cooling_function.txt", ios::app);
  cooling_outfile << T << "  ";
#endif

  int cooling_sct = 0;
  for (size_t e = 0; e < cooling_species.size(); e++) {
    for (size_t q = 0; q < cooling_species[e].size(); q++) {
      number_density = species_number_density[e][q];

#ifdef MELLEMA
      // mellema cooling rate
      L = mellema_cooling_rate(cooling_sct, T, ne);
#elif defined CHIANTI
      // chianti cooling rate
      L = chianti_cooling_rate(cooling_sct, T, ne);
#endif

      Lambda += ne * number_density * L;
      // if (ne * number_density * L > 4e-23) {
      //  spdlog::info("{} {:.2f} {:.2f} {:12.3e} {:12.3e}
      //  {:12.3e}",cooling_sct, T, ne, number_density, L, ne * number_density *
      //  L); spdlog::info("State vec {}, {}",y_species_frac, y_top_ion_frac);
      //}

#ifdef NEMO_CIE_COOLING
      // cooling function of individual species ne*number_density*L
      cooling_outfile << ne * number_density * L << "  ";
#endif

      cooling_sct++;
    }
  }

#ifdef NEMO_CIE_COOLING
  cooling_outfile << ne << "  " << element_num_density[0] << "\n";
  cooling_outfile.close();
#endif

  // Updating Edot due to cooling function.
  Edot -= Lambda;
  // spdlog::info("{:12.6e}  {:12.6e}",photoheat_rate, Lambda);

  // END OF SECTION: Cooling Function ******************************************


  //
  // We want to limit cooling as we approach the minimum temperature, so we
  // scale the rate to linearly approach zero as we reach T_min.
  //
  if (Edot < 0.0 && T < 2.0 * EP->MinTemperature) {
    Edot = min(0.0, (Edot) * (T - EP->MinTemperature) / EP->MinTemperature);
  }

  // Switch to toggle the time derivative of energy density on and off.
  if (!EP->update_erg) Edot = 0.0;

  // Updating the ydot vector element corresponding to the index E_index.
  NV_Ith_S(y_dot, E_index) = Edot;

  /*
  for (int v = 0; v < N_species; v++) {
    if (y_species_frac[v] < NEMO_ABSTOL && NV_Ith_S(y_dot, v) < 0.0)
      NV_Ith_S(y_dot, v) = 0.0;
    if (y_species_frac[v] > 1.0 - NEMO_ABSTOL && NV_Ith_S(y_dot, v) > 0.0)
      NV_Ith_S(y_dot, v) = 0.0;
  }
   */

  // for (int v = 0; v < N_local; v++) p_local[v] = NV_Ith_S(y_dot, v);
  // spdlog::info("ydot: {}",p_local);

  return 0;
}
// END OF YDOT FUNCTION ####################################################



// GET COLLISIONAL IONISATION RATE FROM LUT ################################
double nemo::collisional_ionisation_rate(
    const int species_index,  // species index (first index)
    const double T            // temperature
)
{

  int T_index = int((log10(T) - log10(T_min)) / delta_logT);
  double dT   = T - T_table[T_index];

  double rate = collisional_rate_table[species_index][T_index]
                + dT * collisional_slope_table[species_index][T_index];

  return rate;
}



// GET RECOMBINATION RATE ###################################################
double nemo::recombination_rate(
    const int species_index,  // species identifier
    const double T,           // temperature
    double &caseA             // case A rate (For H+ and He+)
)
{

  int T_index = int((log10(T) - log10(T_min)) / delta_logT);
  double dT   = T - T_table[T_index];

  caseA = caseA_recombination_rate_table[species_index][T_index]
          + dT * caseA_recombination_slope_table[species_index][T_index];

  return recombination_rate_table[species_index][T_index]
         + dT * recombination_slope_table[species_index][T_index];
}

// GET IONISING CHARGE EXCHANGE RATE ###########################################
double nemo::ionise_charge_exchange_rate(
    const int element,   // elemental position of first reactant
    const int reaction,  // second reactant identifier
    const double T       // temperature
)
{
  // calculating temperature index and dT
  int T_index = int((log10(T) - log10(T_min)) / delta_logT);
  double dT   = T - T_table[T_index];
  double rate = ionise_cx_rate_table[element][reaction][T_index]
                + dT * ionise_cx_slope_table[element][reaction][T_index];
  return rate;
}


// GET RECOMBINING CHARGE EXCHANGE RATE ########################################
double nemo::recomb_charge_exchange_rate(
    const int element,   // elemental position of first reactant
    const int reaction,  // second reactant identifier
    const double T       // temperature
)
{
  // calculating temperature index and dT
  int T_index = int((log10(T) - log10(T_min)) / delta_logT);
  double dT   = T - T_table[T_index];
  // calculating rate
  double rate = recomb_cx_rate_table[element][reaction][T_index]
                + dT * recomb_cx_slope_table[element][reaction][T_index];
  return rate;
}

// GET PHOTO-IONISATION RATE  ###############################################
double nemo::photoionisation_rate(
    const int species_idx,            // photo species index (first index)
    const double species_num_density  // species number density
)
{
  // spdlog::info("index = {}, n_H = {}, n_s = {}", species_index,
  // element_num_density, species_num_density);
  double elem_num_density =
      element_num_density[photo_species_element[species_idx]];
  return calculate_photoionisation_rate(
      species_idx, species_num_density, elem_num_density, vshell, ds, tau_bins);
}

// GET PHOTO-HEATING RATE  ###############################################
double nemo::photoheating_rate(
    const int species_index,          // photo species index (first index)
    const double species_num_density  // species number density
)
{
  return calculate_photoheating_rate(
      species_index, species_num_density, vshell, ds, tau_bins);
}



#ifdef MELLEMA
// GET MELLEMA COOLING RATE #################################################
double nemo::mellema_cooling_rate(
    const int location,  // species identifier (database location)
    double T,            // Temperature
    double ne            // electron number density
)
{

  // 1. Set up log10(T) and log10(ne) indices to access LUT.
  // A conditional statement is added to access the Mellema Cooling table only
  // if ne and T range lie within the table in the database of
  // Mellema_/// rate.cpp

  if (T > T_max) T = T_max;
  if (T < T_min) T = T_min;
  if (ne < ne_min) ne = ne_min;

  // Note: Since the Mellema cooling rates are given in log scales,
  // i.e, L(Log T, Log ne), we calculate the differential change dLog_T
  // (step) and dLog_ne along with their indices. Delta_Log_T and
  // Delta_Log_ne is a constant step size calculated in the
  // Mellema_cooling_rate.cpp

  // Get temperature vector index
  int T_index = int((log10(T) - log10(T_min)) / delta_logT);
  // Get the difference in temperature.
  double dT = T - T_table[T_index];

  // Get electron number density vector index.
  int ne_index = int((log10(ne) - log10(ne_min)) / delta_logne);
  // Get the difference in electron number density.
  double dne = ne - ne_table[ne_index];

  // cout<<"T= "<<T<< ", "<<"T_index= "<<T_Index<<", "<<"T= "<<Vec_T[T_Index]
  // << endl; cout<<"ne= "<<ne<< ", "<<"ne_index= "<<ne_Index<<", "<<"ne=
  // "<<Vec_ne[ne_Index] << endl;

  // Calculating Mellema Cooling Rate.
  double L = 0.0;
  // Calculating cooling from the LookUp Table (LUT) with the location of
  // each Ion and atom in the database.
  L = mellema_table[location].rate[T_index][ne_index]
      + dT * mellema_table[location].T_rateslope[T_index][ne_index]
      + dne * mellema_table[location].ne_rateslope[T_index][ne_index];
  return L;
}
#endif


#ifdef CHIANTI
// GET CHIANTI COOLING RATE ###################################################
double nemo::chianti_cooling_rate(
    const int location,  // species identifier (database location)
    double T,            // Temperature
    double ne            // electron number density
)
{

  // 1. Set up log10(T) and log10(ne) indices to access LUT.
  // A conditional statement is added to access the chinati cooling table only
  // if ne and T range lie within the table in the database of
  // Mellema_/// rate.cpp

  if (T > T_max) T = T_max;
  if (T < T_min) T = T_min;
  if (ne < ne_min) ne = ne_min;

  // Note: Since the chianti cooling rates are given in log scales,
  // i.e, L(Log T, Log ne), we calculate the differential change dLog_T
  // (step) and dLog_ne along with their indices. Delta_Log_T and
  // Delta_Log_ne is a constant step size calculated in the
  // chianti_cooling_rate.cpp

  // Get temperature vector index
  int T_index = int((log10(T) - log10(T_min)) / delta_logT);
  // Get the difference in temperature.
  double dT = T - T_table[T_index];

  // Get electron number density vector index.
  int ne_index = int((log10(ne) - log10(ne_min)) / delta_logne);
  // Get the difference in electron number density.
  double dne = ne - ne_table[ne_index];

  // cout<<"T= "<<T<< ", "<<"T_index= "<<T_Index<<", "<<"T= "<<Vec_T[T_Index]
  // << endl; cout<<"ne= "<<ne<< ", "<<"ne_index= "<<ne_Index<<", "<<"ne=
  // "<<Vec_ne[ne_Index] << endl;

  // Calculating chianti cooling rate.
  double L = 0.0;
  // Calculating cooling from the LookUp Table (LUT) with the location of
  // each Ion and atom in the database.
  L = chianti_table[location].rate[T_index][ne_index]
      + dT * chianti_table[location].T_rateslope[T_index][ne_index]
      + dne * chianti_table[location].ne_rateslope[T_index][ne_index];

  return L;
}
#endif


// SET RADIATION SOURCE PARAMETERS ############################################
void nemo::setup_radiation_source_parameters(
    const int N_ion, const std::vector<struct rt_source_data> &ion_src
    ///< list of ionising src column densities and source properties.
)
{

  if (ion_src.size() != static_cast<unsigned int>(N_ion)) {
    spdlog::error(
        "{}: {}",
        "Timescales: N_ionising_srcs doesn't match vector size in NEMO",
        ion_src.size());
    exit(1);
  }

  // proceed if there are any ionizing sources present.
  if (N_ion > 0) {
    // IMPO: this is hard coded for a single source.
    // path length through the current cell
    ds = ion_src[0].dS;
    // shell volume in the current cell
    vshell = ion_src[0].Vshell;
    // optical depth in each energy bin at the edge of the cell
    tau_bins.resize(N_ebins, 0.0);
    for (int b = 0; b < N_ebins; b++)
      tau_bins[b] = ion_src[0].Column[b];
  }
  // spdlog::debug("ds {:12.6e}, Vsh {:12.6e}, col {}", ds, vshell, tau_bins);
}



// ############################################################################
// ***************************** END OF NEMO *********************************
// ############################################################################



// ############################################################################
//  DEBUG  FUNCTIONS

// ############################################################################
// Print collisional ionisation and recombination LUT to file
void nemo::print_CIR_LUT(
    const std::vector<string> &tracer_list,
    const std::string info,
    const std::vector<double> &X,
    const std::vector<std::vector<double> > &Fx,
    const std::vector<std::vector<double> > &Fx_slope,
    const std::string filename)
{


  std::ofstream outfile(filename + ".txt");
  outfile << info;
  outfile << "------ Attributes ------\n";
  outfile << "# First column: Temperature (K) \n";
  outfile << "# Subsequent columns: Rates (cm^3/s)\n";
  outfile << "# No of species: " << tracer_list.size() << endl;
  outfile << "# No of temperature points: 200\n";
  outfile << "------- Species --------\n";
  outfile << "# Species NEMO Names: ";
  for (int j = 0; j < tracer_list.size(); j++) {
    outfile << "\t" << tracer_list[j];
  }
  outfile << "\n";

  outfile << "-------- Data ----------\n";

  if (!outfile.is_open()) {
    spdlog::error("NEMO {} : {}", "couldn't open outfile", 1);
    exit_nemo(1);
  }

  outfile.setf(ios_base::scientific);
  outfile.precision(4);

  double x_min      = X.front();
  double x_max      = X.back();
  double delta_logx = log10(X[1]) - log10(X[0]);

  // Make a new set of X array, call it X_new
  int N_points = 200;  // No of new x-points, set to 200
  double delta = (log10(x_max) - log10(x_min)) / (N_points - 1);
  std::vector<double> X_new;
  for (int i = 0; i < N_points; i++) {
    double logx_new = log10(x_min) + i * delta;
    X_new.push_back(pow(10.0, logx_new));
  }

  for (int k = 0; k < X_new.size(); k++) {
    int x_index = int((log10(X_new[k]) - log10(x_min)) / delta_logx);
    double dx   = X_new[k] - X[x_index];
    outfile << X_new[k];
    for (int i = 0; i < Fx.size(); i++) {
      outfile << "  " << Fx[i][x_index] + dx * Fx_slope[i][x_index];
    }
    outfile << "\n";
  }
}

// ############################################################################
// Print charge exchange reaction rate LUT to file
void nemo::print_charge_exchange_LUT(
    const std::vector<std::vector<string> > &reactions,
    const std::string info,
    const std::vector<double> &T,
    const std::vector<std::vector<std::vector<double> > > &rate,
    const std::vector<std::vector<std::vector<double> > > &rate_slope,
    const std::vector<std::vector<double> > &energy_defect,
    const std::string filename)
{
  // Caculating number of reactions
  int N_reaction = 0;
  std::string reaction_list;  // Corrected variable name
  for (size_t e = 0; e < reactions.size(); e++) {
    for (size_t r = 0; r < reactions[e].size(); r++) {
      reaction_list += reactions[e][r] + "\t";
      N_reaction++;
    }
  }
  std::ofstream outfile(filename + ".txt");
  outfile << info;

  outfile << "------ Reactions ------\n";
  outfile << "# CX Reactions: ";
  outfile << reaction_list;
  outfile << "\n";

  outfile << "------ Energy ------\n";
  outfile << "# Energy-defect (eV): ";
  for (size_t i = 0; i < energy_defect.size(); i++)
    for (size_t j = 0; j < energy_defect[i].size(); j++)
      outfile << energy_defect[i][j] << "  ";
  outfile << "\n";

  outfile << "------ Attributes ------\n";
  outfile << "# First column: Temperature (K) \n";
  outfile << "# Subsequent columns: Reaction rates (cm^3/s)\n";
  outfile << "# No of reactions: " << N_reaction << endl;
  outfile << "# No of temperature points: 200\n";

  outfile << "------ Data ------\n";

  if (!outfile.is_open()) {
    spdlog::error("NEMO {} : {}", "couldn't open outfile", 1);
    exit_nemo(1);
  }

  outfile.setf(ios_base::scientific);
  outfile.precision(4);

  double temp_min       = T.front();
  double temp_max       = T.back();
  double delta_log_temp = log10(T[1]) - log10(T[0]);

  // Make a new set of temperature array, call it temp_new.
  int N_points = 400;  // No of new temperature-points, set to 200.
  double delta = (log10(temp_max) - log10(temp_min)) / (N_points - 1);
  std::vector<double> temp_new;
  for (int i = 0; i < N_points; i++) {
    double log_temp_new = log10(T_min) + i * delta;
    temp_new.push_back(pow(10.0, log_temp_new));
  }

  for (size_t k = 0; k < temp_new.size(); k++) {
    int temp_index =
        int((log10(temp_new[k]) - log10(temp_min)) / delta_log_temp);
    double dx = temp_new[k] - T[temp_index];
    outfile << temp_new[k];
    for (size_t s = 0; s < rate.size(); s++) {
      for (size_t r = 0; r < rate[s].size(); r++) {
        outfile << "  "
                << rate[s][r][temp_index] + dx * rate_slope[s][r][temp_index];
      }
    }
    outfile << "\n";
  }
}

// ############################################################################
// Use the following function to print photo-ionisation mean cross-section
// and bin_fraction to a txt file
void nemo::print_to_file_photoionisation(
    const std::string info,
    const std::vector<string> &species_list,
    const std::vector<std::vector<double> > &X,
    const std::vector<std::vector<double> > &Fx,
    const std::string filename)
{
  std::ofstream outfile(filename + ".txt");
  if (!outfile.is_open()) {
    spdlog::error("NEMO {} : {}", "couldn't open outfile", 1);
    exit(1);
  }

  outfile << "#FILE: Generated by PION v3.0 - NEMO Module \n";
  outfile << "#DATA DESCRIPTION: " + info + "\n";
  // Column names
  outfile << "#ATTRIBUTE LABEL: Bin_Min\tBin_Max";
  for (int j = 0; j < species_list.size(); j++) {
    outfile << "\t" << species_list[j];
  }
  outfile << "\n";
  // column units
  outfile << "#UNIT: eV\teV";
  size_t found = info.find("cross-section");
  if (found != std::string::npos) {
    for (int j = 0; j < species_list.size(); j++)
      outfile << "\t"
              << "cm^2";
  }
  else {
    for (int j = 0; j < species_list.size(); j++)
      outfile << "\t"
              << "1";
  }
  outfile << "\n";


  outfile.setf(ios_base::scientific);
  outfile.precision(4);
  outfile << "#DATA: \n";
  for (int k = 0; k < X.size(); k++) {
    outfile << X[k][0] << "\t" << X[k][1];
    for (int j = 0; j < Fx.size(); j++) {
      outfile << "\t" << Fx[j][k];
    }
    outfile << "\n";
  }
}

// ###########################################################################
// Print cooling tables to file
void nemo::print_cooling_table(
    std::vector<double> &logne,  // log10(ne)
    std::vector<double> &logT,   // log10(T)
    std::vector<std::vector<double> > &cooling_rate,
    const std::string ion_name,
    std::string filename)
{

  std::ofstream outfile(filename + "_" + ion_name + ".txt");
  if (!outfile.is_open()) {
    spdlog::error("NEMO {} : {}", "couldn't open outfile", 1);
    exit(1);
  }

  outfile << "#FILE: Generated by PION 3.0.0 - NEMO Module \n";
  outfile << "#DATA DESCRIPTION: " + ion_name + " cooling rate table\n";
  outfile << "#ATTRIBUTE LABEL: Temperature";  // Column names

  for (int j = 0; j < logne.size(); j++) {
    outfile << "\t"
            << "Rate";
  }
  outfile << endl;
  outfile << "#UNIT: K";
  for (int j = 0; j < logne.size(); j++) {
    outfile << "\t"
            << "ergs cm^3/s";
  }
  outfile << endl;

  outfile << "#DATA: \n";
  outfile << "#Log(ne): ";
  for (int j = 0; j < logne.size(); j++) {
    outfile << "\t" << log10(logne[j]);
  }
  outfile << endl;

  for (size_t j = 0; j < cooling_rate.size(); j++) {
    outfile << log10(logT[j]) << ",";
    for (size_t k = 0; k < cooling_rate[j].size(); k++) {
      outfile << log10(cooling_rate[j][k]);
      if (k == cooling_rate[j].size() - 1) {
        outfile << endl;
      }
      else
        outfile << ",";
    }
  }

  outfile.close();
}


// #########################################################################
// Print stellar attributes and binned luminosity to file
void nemo::print_stellar_properties(
    const star_properties_datatype &star_attributes,
    const std::vector<double> &binned_luminosity,
    const std::vector<std::vector<double> > &energy_bins,
    std::string filename)
{
  filename =
      filename + "_srctype" + std::to_string(star_attributes.Type) + ".txt";
  ofstream file;
  file.open(filename);
  file << "#FILE: Generated by PION v3.0 - Stellar Atmosphere Model Module\n\n";
  file << "#COMMENT: Source id " << star_attributes.Type << endl;
  file << "#COMMENT: Teff " << star_attributes.Teff << " K" << endl;
  file << "#COMMENT: Radius " << star_attributes.Radius << " solar radii"
       << endl;
  file << "#COMMENT: Luminosity " << star_attributes.Luminosity << " ergs/s"
       << endl;
  file << star_attributes.SED_info << endl;
  file << "#UNIT: eV\teV\tegrs/s/bins" << endl;
  file << "#NAME: Bin_Min\tBin_Max\tBin_Luminosity" << endl;
  file << "#DATA:" << endl;
  for (size_t i = 0; i < binned_luminosity.size(); i++) {
    file << energy_bins[i][0] << "\t" << energy_bins[i][1] << "\t"
         << binned_luminosity[i];
    if (i == binned_luminosity.size() - 1)
      break;
    else
      file << endl;
  }
  file << '\n';
  file << '\n';
  file << "#COMMENT: Total fraction of Luminosity in the bins ";
  double tot_frac_L = 0.0;
  for (size_t i = 0; i < bin_luminosity.size(); i++) {
    tot_frac_L += bin_luminosity[i] / star_attributes.Luminosity;
  }
  file << tot_frac_L << endl;
  file.close();
}

#pragma clang diagnostic pop
