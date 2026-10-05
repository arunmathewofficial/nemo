/**
 *  \file      charge_exchange.h
 *  \brief     Charge exchange reaction class definition.
 *  \author    Arun Mathew
 *  \date      10-03-2024
 **/

#include "charge_exchange.h"
#include <algorithm>  // for std::find
#include <string>     // for std::string
#include <vector>     // for std::vector

// #############################################################################
// Constructor
charge_exchange::charge_exchange(const int max_reactant_charge)
{
  // assigning the maximum reactant charge.
  q_max = max_reactant_charge;
  // initializing charge exchange database.
  charge_exchange_database();
}

// #############################################################################
// Destructor
charge_exchange::~charge_exchange() {}

// #############################################################################
void charge_exchange::setup_charge_exchange(
    const std::vector<std::vector<std::pair<std::string, std::string> > >
        &ionise_reaction_list,
    const std::vector<std::vector<std::pair<std::string, std::string> > >
        &recomb_reaction_list)
{
  // Number of ionisation and recombination species
  // Note: The outer size of the 2D array gives the number of species
  N_ionise_species = ionise_reaction_list.size();
  N_recomb_species = recomb_reaction_list.size();

  // Resizing ionisation charge exchange rate and slope tables
  ionise_cx_rate_table.resize(N_ionise_species);
  ionise_cx_slope_table.resize(N_ionise_species);
  for (int s = 0; s < N_ionise_species; s++) {
    ionise_cx_rate_table[s].resize(ionise_reaction_list[s].size());
    ionise_cx_slope_table[s].resize(ionise_reaction_list[s].size());
  }

  // Resizing recombination charge exchange rate and slope tables
  recomb_cx_rate_table.resize(N_recomb_species);
  recomb_cx_slope_table.resize(N_recomb_species);
  for (int s = 0; s < N_recomb_species; s++) {
    recomb_cx_rate_table[s].resize(recomb_reaction_list[s].size());
    recomb_cx_slope_table[s].resize(recomb_reaction_list[s].size());
  }

  // make copies of ionisation reactions and recombination reactions
  ionise_reactions = ionise_reaction_list;
  recomb_reactions = recomb_reaction_list;
}

// #############################################################################
std::vector<std::vector<double> > charge_exchange::import_energy_defect(
    const std::vector<std::vector<std::pair<std::string, std::string> > >
        &reaction_list)
{
  // Export the energy defect of each reaction in the reaction list to the
  // NEMO module.
  int N_species = reaction_list.size();
  std::vector<std::vector<double> > charge_exchange_energies(N_species);
  for (size_t s = 0; s < N_species; s++) {
    for (const auto &species_reaction : reaction_list[s]) {
      // Retrieve the reaction energy defect, however get actual key to find
      // the reaction
      AR_parameters &reaction = AR_reaction[get_actual_key(species_reaction)];
      charge_exchange_energies[s].push_back(reaction.dE);
    }
  }
  // Note: Energies are returned in electron volts (eV).
  return charge_exchange_energies;
}

// #############################################################################
bool charge_exchange::reactant_exists_in_reference(
    const std::string &species, const std::vector<std::string> &reference)
{
  return std::find(reference.begin(), reference.end(), species)
         != reference.end();
}

// #############################################################################
std::vector<std::pair<std::string, std::string> >
charge_exchange::get_cx_reactions(
    const std::string &species,
    const bool type,
    const std::vector<std::string> &reference  // reference species list
)
{
  // All reactions corresponding to the given species
  std::vector<std::pair<std::string, std::string> > species_reactions;
  // Search for all reactions with the given species as the first reactant,
  // specific type, and the first reactant charge less than or equal to qmax.
  for (const auto &reaction : AR_reaction) {
    if (reaction.first.first == species && reaction.second.type == type
        && atomic_physics_data.get_charge(reaction.first.first) <= q_max
        && reactant_exists_in_reference(reaction.first.second, reference)) {
      species_reactions.push_back(reaction.first);
    }
    // Search for all reactions with the given species as the second reactant,
    // specific type, and the first reactant charge less than or equal to qmax.
    if (reaction.first.second == species && reaction.second.type == !type
        && atomic_physics_data.get_charge(reaction.first.first) <= q_max
        && reactant_exists_in_reference(reaction.first.first, reference)) {
      std::pair<std::string, std::string> swapped_key =
          std::make_pair(reaction.first.second, reaction.first.first);
      species_reactions.push_back(swapped_key);
    }
  }
  return species_reactions;
}

// #############################################################################
std::pair<std::string, std::string> charge_exchange::get_actual_key(
    const std::pair<std::string, std::string> &reaction_key)
{
  // This function returns the actual key of the reaction. The bookkeeping
  // process for the charge exchange reaction flips the key position, if
  // necessary, to calculate the ionization and recombination rates. However,
  // after swapping the key positions, the reaction may no longer exist for
  // that key in the database. Therefore, it is important to select the correct
  // reaction for the modified keys while generating lookup tables.
  std::pair<std::string, std::string> actual_key;

  // Retrieve the reaction corresponding to the key from the database.
  AR_parameters &reaction = AR_reaction[reaction_key];
  // If the reaction ID returns a non-empty string, it indicates that the
  // reaction exists and that the inputted key matches the actual key. In
  // this case, the function will return the given key.
  if (!reaction.reaction_id.empty()) {
    actual_key = std::make_pair(reaction_key.first, reaction_key.second);
  }
  // If the ID returns an empty string, swap the reactants' order
  // in the key, retrieve the reaction from the database, and push
  // back the corresponding reaction ID.
  else {
    std::pair<std::string, std::string> swapped_key =
        std::make_pair(reaction_key.second, reaction_key.first);
    AR_parameters &swapped_reaction = AR_reaction[swapped_key];
    // Check for the reaction in the database, if exist return the swapped key
    if (!swapped_reaction.reaction_id.empty()) {
      actual_key = std::make_pair(swapped_key.first, swapped_key.second);
    }
  }
  return actual_key;
}

// #############################################################################
std::vector<std::vector<std::string> > charge_exchange::import_reactions_id(
    std::vector<std::vector<std::pair<std::string, std::string> > >
        &reaction_list)
{
  // The outer size of the 2D array gives the number of species
  int N_species = reaction_list.size();
  std::vector<std::vector<std::string> > reactions_identifier(N_species);
  // loop over collisional/recombination species
  for (size_t s = 0; s < N_species; s++) {
    // loop over reactions of each species
    for (const auto &species_reaction : reaction_list[s]) {
      // Retrieve the reaction corresponding to the key, however get the
      // actual key before finding the reaction
      AR_parameters &reaction = AR_reaction[get_actual_key(species_reaction)];
      reactions_identifier[s].push_back(reaction.reaction_id);
    }
  }
  return reactions_identifier;
}

// #############################################################################
double charge_exchange::get_charge_exchange_rate(
    const std::pair<std::string, std::string> &reaction_key, const double T)
{
  double rate             = 0.0;
  const auto &reaction_it = AR_reaction.find(reaction_key);
  if (reaction_it != AR_reaction.end()) {
    const auto &reaction = reaction_it->second;
    // Note: We use the charge transfer reaction rates between neutral hydrogen
    // and elements with a parent ion charge q ranging from 1 to 4,
    // obtained from an analytic fit formula presented by J. B. Kingdon
    // and G. J Ferland (ApJS, 106, 206, 1996). The recombination rate
    // coefficient is represented via a four-parameter fit formula,
    // a*(T/10**4K)**b*[1+c*exp(d*T/10**4K)] cm**3/s, while the ionization
    // reaction rate coefficients are expressed as a product of a fitting
    // function and the Boltzmann factor exp(-dE_fit/kT). The value dE_fit is
    // non-negative, as referenced in J. B. Kingdon & G. J. Ferland, ApJS,
    // 106:205 (1996).

    // For charge transfer recombination reactions with He, we use the fit
    // formula from M. Arnaud & R. Rothenflug (1985, A&AS, 60, 425), which is
    // the same four-parameter formula used for recombination reactions with H.
    // However, the ionization charge exchange reaction rate coefficient is
    // given by a slightly different formula: rate = a*(T/10**4K)**b
    // *exp(-c*T/10**4K)*exp(-dE_fit/kT) cm**3/s.

    // SECTION 1: Charge-exchange recombination reaction with H and He
    if (!reaction.type) {
      rate = reaction.a * pow(T * 1.00E-04, reaction.b)
             * (1.00 + reaction.c * exp(reaction.d * T * 1.00E-04));
    }

    // SECTION 2:. Charge-exchange ionization reaction with H1+
    if (reaction_key.second == "H1+" && reaction.type) {
      rate = reaction.a * pow(T * 1.00E-04, reaction.b)
             * (1.00 + reaction.c * exp(reaction.d * T * 1.00E-04))
             * exp(-reaction.dE_fit * eV() / (T * kB()));
    }

    // SECTION 3: Charge-exchange ionisation reaction with He1+
    if (reaction_key.second == "He1+" && reaction.type) {
      rate = reaction.a * pow(T * 1.00E-04, reaction.b)
             * exp(-reaction.c * T * 1.00E-04)
             * exp(-reaction.dE_fit * eV() / (T * kB()));
    }

    // Note: To ensure rates have acceptable values, we modify the rate
    // equation for temperatures outside the range of T_min and T_max.
    if (T < reaction.T_min)
      rate *= pow(exp(-reaction.T_min / T) / 2.71828, 3.0);
    if (T > reaction.T_max)
      rate *= pow(exp(-T / reaction.T_max) / 2.71828, 3.0);
  }
  return rate;
}

// #############################################################################
void charge_exchange::generate_charge_exchange_LUT(
    const std::vector<double> &T_table)
{

  spdlog::info("Generating charge-exchange LUTs");

  // std::pair<std::string, std::string> key = std::make_pair("C1+", "H");
  // spdlog::debug("rate = {}", get_charge_exchange_rate(key, 5.0E+6));
  // exit(1);

  // SECTION: LUT for Ionisation charge-exchange reactions *********************
  // 1. Generate ionisation charge-exchange reaction rate tables.
  for (int s = 0; s < N_ionise_species; s++) {
    for (size_t r = 0; r < ionise_reactions[s].size(); r++) {
      for (int k = 0; k < T_table.size(); k++) {
        double ionise_rate = get_charge_exchange_rate(
            get_actual_key(ionise_reactions[s][r]), T_table[k]);
        ionise_cx_rate_table[s][r].push_back(ionise_rate);
      }
    }
  }

  //  2.Generate charge-exchange rate slope table for ionisation reactions.
  for (int s = 0; s < N_ionise_species; s++) {
    for (size_t r = 0; r < ionise_reactions[s].size(); r++) {
      for (int k = 0; k < T_table.size() - 1; k++) {
        ionise_cx_slope_table[s][r].push_back(
            (ionise_cx_rate_table[s][r][k + 1] - ionise_cx_rate_table[s][r][k])
            / (T_table[k + 1] - T_table[k]));
      }
      ionise_cx_slope_table[s][r].push_back(0.0);
    }
  }
  // END OF SECTION: LUT for Ionisation charge-exchange reactions **************

  // SECTION: LUT for recombination charge-exchange reactions ******************
  // 1. Generate charge-exchange rate tables for recombination reactions.
  for (int s = 0; s < N_recomb_species; s++) {
    for (size_t r = 0; r < recomb_reactions[s].size(); r++) {
      for (int k = 0; k < T_table.size(); k++) {
        double recomb_rate = get_charge_exchange_rate(
            get_actual_key(recomb_reactions[s][r]), T_table[k]);
        recomb_cx_rate_table[s][r].push_back(recomb_rate);
      }
    }
  }

  // 2. Generate charge-exchange rate slope table for recombination reactions.
  for (int s = 0; s < N_recomb_species; s++) {
    for (size_t r = 0; r < recomb_reactions[s].size(); r++) {
      for (int k = 0; k < T_table.size() - 1; k++) {
        recomb_cx_slope_table[s][r].push_back(
            (recomb_cx_rate_table[s][r][k + 1] - recomb_cx_rate_table[s][r][k])
            / (T_table[k + 1] - T_table[k]));
      }
      recomb_cx_slope_table[s][r].push_back(0.0);
    }
  }
  // END OF SECTION: LUT for recombination charge-exchange reactions ***********
}

// #############################################################################
void charge_exchange::charge_exchange_database()
{

  spdlog::info("Setting up charge exchange database");

  // This function serves as the database for charge exchange reactions.
  // The reactions included in the database follow the M. Arnaud &
  // R. Rothenflug (AR for short) fit formula.

  // The reference for the energy defect (dE) for each reaction is provided
  // above each entry.
  //
  // The dE values missing in the literature are crudely
  // calculated using the equation dE = I(X^{n-1}) - I(A) for the reaction
  // X^n + A --> X^{n-1} + A^x, where I represents the ionisation potential
  // of the species.

  // [RECOMBINATION REACTION WITH H *******************************************]
  // The fit parameters and energy defect for the recombination charge exchange
  // reaction with H1 are obtained from J. B. Kingdon and G. J. Ferland
  // (ApJS, 106, 206, 1996).

  // [dE NOT CLEAR] Ref: B. Zygelman & A. Dalgarno, ApJ, 365:239 (1990)
  // [USES dE FROM] Ref: B.Friedman et al., J.Phys.B:At.Mol.Opt,50:115202 (2017)
  AR_reaction[{"He1+", "H"}] = AR_parameters{
      "ar-he1+h", 0, 7.47E-15, 2.06, 9.93, -3.89, 6.00e3, 1.00e5, 0.0, -8.83};

  // [dE NOT CLEAR] Ref: B.W. West & N.F. Lane, Phy. Rev. A, 26:3164 (1982)
  // [USES dE FROM] Ref: B.Friedman et al., J.Phys.B:At.Mol.Opt,50:115202 (2017)
  AR_reaction[{"He2+", "H"}] = AR_parameters{
      "ar-he2+h", 0, 1.00E-14, 0.00, 0.00, 0.00, 1.00e3, 1.00e7, 0.0, 40.82};

  // [USES dE FROM] Ref: S.E. Butler & A. Dalgarno, A&A, 85:144 (1980)
  // [RATE CHECK WITH] Ref: S.E. Butler & A. Dalgarno, A&A, 85:144 (1980) [OK]
  AR_reaction[{"C1+", "H"}] = AR_parameters{
      "ar-c1+h", 0, 1.76E-18, 8.33, 4278.78, -6.41, 3.00e3, 1.00e4, 0.0, -2.33};

  // [RATE CHECK WITH] Ref: S.E. Butler, T. G Heil & A. Dalgarno, ApJ, 241, 442
  // (1980) [OK]
  // [USES dE FROM] Ref: B.Friedman et al., J.Phys.B:At.Mol.Opt,50:115202 (2017)
  AR_reaction[{"C2+", "H"}] = AR_parameters{
      "ar-c2+h", 0, 1.67E-13, 2.79, 304.72, -4.07, 5.00e3, 5.00e4, 0.0, 10.78};

  // [RATE CHECK WITH] Ref: S.E. Butler, T. G Heil & A. Dalgarno, ApJ, 241, 442
  // (1980) [OK]
  // [USES dE FROM] Ref: B.Friedman et al., J.Phys.B:At.Mol.Opt,50:115202 (2017)
  AR_reaction[{"C3+", "H"}] = AR_parameters{
      "ar-c3+h", 0, 3.25E-09, 0.21, 0.19, -3.29, 1.00e3, 1.00e5, 0.0, 21.60};

  // [RATE CHECK WITH] Ref: M. Gargaud et al., J Phys. B: At. Mol. Phys. 14:
  // 2259 (1981) [OK] [USES dE FROM] Ref: B.Friedman et al.,
  // J.Phys.B:At.Mol.Opt,50:115202 (2017)
  AR_reaction[{"C4+", "H"}] =
      AR_parameters{"ar-c4+h", 0,      332.46E-09, -0.11, -9.95e-1,
                    -1.58e-3,  1.00e1, 1.00e5,     0.0,   13.34};

  // [RATE CHECK WITH] Ref: S.E. Butler & A. Dalgarno, ApJ, 234:767 (1979) [OK]
  // [USES dE FROM] Ref: S.E. Butler & A. Dalgarno, ApJ, 234:767 (1979)
  AR_reaction[{"N1+", "H"}] = AR_parameters{
      "ar-n1+h", 0, 1.01E-12, -0.29, -0.92, -8.38, 1.00e2, 5.00e4, 0.0, 0.94};

  // [RATE CHECK WITH] Ref: S.E. Butler, T. G Heil & A. Dalgarno, ApJ, 241, 442
  // (1980), B. Herrero et al., J Phys. B: At. Mol. Opt. Phys. 28:711 (1995)[OK]
  // [USES dE FROM] Ref: B.Friedman et al., J.Phys.B:At.Mol.Opt,50:115202 (2017)
  AR_reaction[{"N2+", "H"}] = AR_parameters{
      "ar-n2+h", 0, 3.05E-10, 0.60, 2.65, -0.93, 1.00e3, 1.00e5, 0.0, 11.95};

  // [RATE CHECK WITH] Refs: S.E. Butler, T. G Heil & A. Dalgarno, ApJ, 241, 442
  // (1980), M. Gargaud et al., J Phys. B: At. Mol. Phys. 14: 2259 (1981) [OK]
  AR_reaction[{"N3+", "H"}] = AR_parameters{
      "ar-n3+h", 0, 4.54e-09, 0.57, -0.65, -0.89, 1.00e1, 1.00e5, 0.0, 0.00};

  // [CHECK PENDING]
  AR_reaction[{"N4+", "H"}] = AR_parameters{
      "ar-n4+h", 0, 3.28e-09, 0.52, -0.52, -0.19, 1.00e3, 3.00e5, 0.0, 0.00};

  // [RATE CHECK WITH] Ref: G.B. Field & G. Steigman, ApJ, 166:59 (1971)[OK]
  // [USES dE FROM] Ref: S.E. Butler & J.C. Raymond, ApJ, 240:680 (1980)
  AR_reaction[{"O1+", "H"}] = AR_parameters{
      "ar-o1+h", 0, 1.04e-09, 3.15e-2, -0.16, -9.73, 1.00e1, 1.00e4, 0.0, 0.02};

  // [RATE CHECK WITH] Refs: S.E. Butler, T. G Heil & A. Dalgarno, ApJ, 241, 442
  // (1980), P. Honvault et al., A&A, 302:931 (1995) [OK]
  // [USES dE FROM] Ref: B.Friedman et al., J.Phys.B:At.Mol.Opt,50:115202 (2017)
  AR_reaction[{"O2+", "H"}] = AR_parameters{
      "ar-o2+h", 0, 1.04e-09, 0.27, 2.02, -5.92, 1.00e2, 1.00e5, 0.0, 16.50};

  // [RATE CHECK WITH]  Ref: S.E. Butler, T. G Heil & A. Dalgarno, ApJ, 241, 442
  // (1980) [OK]
  AR_reaction[{"O3+", "H"}] = AR_parameters{
      "ar-o3+h", 0, 3.98e-09, 0.26, 0.56, -2.62, 1.00e3, 5.00e4, 0.00, 0.00};

  // [RATE CHECK WITH] Ref: S.E. Butler & A. Dalgarno, ApJ, 241, 838 (1980) [OK]
  AR_reaction[{"O4+", "H"}] = AR_parameters{
      "ar-o4+h", 0, 2.52e-10, 0.63, 2.08, -4.16, 1.00e3, 3.00e4, 0.0, 0.00};

  // [RATE CHECK WITH] Refs: S.E. Butler, T. G Heil & A. Dalgarno, ApJ, 241, 442
  // (1980) [OK]. [COMMENT] The reaction is radiative and will proceed more
  // rapidly with a rate coefficient between 10^{-15} and 10^{-12} cm^3/s.
  AR_reaction[{"Ne2+", "H"}] = AR_parameters{
      "ar-ne2+h", 0, 1.00e-14, 0.00, 0.00, 0.00, 5.00e3, 5.00e4, 0.0, 0.00};

  // [RATE CHECK WITH] Refs: S.E. Butler, T. G Heil & A. Dalgarno, ApJ, 241, 442
  // (1980) [OK]
  AR_reaction[{"Ne3+", "H"}] =
      AR_parameters{"ar-ne3+h", 0,      14.73e-9, 4.52e-2, -0.84,
                    -0.31,      5.00e3, 5.00e4,   0.0,     0.00};

  // [RATE CHECK WITH] Ref: S.E. Butler & A. Dalgarno, ApJ, 241, 838 (1980) [OK]
  AR_reaction[{"Ne4+", "H"}] = AR_parameters{
      "ar-ne4+h", 0, 6.47e-9, 0.54, 3.59, -5.22, 1.00e3, 3.00e4, 0.0, 8.60};

  // [RATE CHECK WITH] Ref: M. Gargaud, R. McCarroll & P. Valiron, A&A, 106:197
  // (1982) [OK]
  // [USES dE FROM] Ref: R. McCarrol & P. Vairon, A&A, 53:83 (1979)
  AR_reaction[{"Si2+", "H"}] = AR_parameters{
      "ar-si2+h", 0, 1.23e-9, 0.24, 3.17, 4.18e-3, 1.00e1, 1.00e6, 0.0, 2.74};

  // [RATE CHECK WITH] Ref: S.E. Butler & A. Dalgarno, ApJ, 241, 838 (1980) [OK]
  AR_reaction[{"Si3+", "H"}] =
      AR_parameters{"ar-si3+h", 0,      4.90e-10, -8.74e-2, -0.36,
                    -0.79,      1.00e3, 3.00e4,   0.0,      0.00};

  // [RATE CHECK PENDING]
  AR_reaction[{"Si4+", "H"}] = AR_parameters{
      "ar-si4+h", 0, 7.58E-09, 0.37, 1.06, -4.09, 1.00e3, 5.00e4, 0.0, 0.0};

  // [USES dE FROM] Ref: S.E. Butler & A. Dalgarno, A&A, 85:144 (1980)
  // [RATE CHECK WITH] Ref: S.E. Butler & A. Dalgarno, A&A, 85:144 (1980) [OK]
  AR_reaction[{"S1+", "H"}] = AR_parameters{
      "ar-s1+h", 0, 3.82e-16, 11.10, 2.57e4, -8.22, 3.00e3, 1.00e4, 0.0, -3.24};

  // [RATE CHECK WITH] Ref: S.E. Butler & A. Dalgarno, ApJ, 241, 838 (1980) [OK]
  AR_reaction[{"S2+", "H"}] = AR_parameters{
      "ar-s2+h", 0, 1.00e-14, 0.00, 0.00, 0.00, 1.00e3, 3.00e4, 0.0, 0.00};

  // [RATE CHECK WITH] Ref: S.E. Butler & A. Dalgarno, ApJ, 241, 838 (1980) [OK]
  AR_reaction[{"S3+", "H"}] = AR_parameters{
      "ar-s3+h", 0, 2.29e-9, 4.02e-2, 1.59, -6.06, 1.00e3, 3.00e4, 0.0, 0.00};

  // [RATE CHECK WITH] Ref: S.E. Butler & A. Dalgarno, ApJ, 241, 838 (1980) [OK]
  AR_reaction[{"S4+", "H"}] = AR_parameters{
      "ar-s4+h", 0, 6.44e-9, 0.13, 2.69, -5.69, 1.00e3, 3.00e4, 0.0, 0.00};

  // [RATE CHECK WITH] Ref: D.A Neufeld & A. Dalgarno, Phys.Rev.A, 35, 3142
  // (1987) [OK]
  // [USES dE FROM] Ref: D.A Neufeld & A. Dalgarno, Phys.Rev.A, 35, 3142 (1987)
  AR_reaction[{"Fe2+", "H"}] = AR_parameters{
      "ar-fe2+h", 0, 1.26e-9, 7.72e-2, -0.41, -7.31, 1.00e3, 1.00e5, 0.0, 2.56};

  // [RATE CHECK PENDING]
  AR_reaction[{"Fe3+", "H"}] = AR_parameters{
      "ar-fe3+h", 0, 3.42e-9, 0.51, -2.06, -8.99, 1.00e3, 1.00e5, 0.0, 0.00};
  //****************************************************************************

  // [IONISATION REACTION WITH H1+ ********************************************]
  // The fit parameters and energy defect for the ionization charge exchange
  // reaction with H1+ are obtained from J. B. Kingdon and G. J. Ferland
  // (ApJS, 106, 206, 1996). The reference provides Delta E in units of Kelvin,
  // which I have converted to eV for convenience.

  // [USES dE_fit FROM] Ref: J.B. Kingdon & G.J. Ferland, ApJS, 106:205 (1996)
  AR_reaction[{"C", "H1+"}] =
      AR_parameters{"ar-c+h1", 1,        1.00E-14, 0.00E+00, 0.00E+00,
                    0.00E+00,  1.00E+03, 1.00E+04, 0.0,      2.33};

  // [USES dE_fit FROM] Ref: J.B. Kingdon & G.J. Ferland, ApJS, 106:205 (1996)
  // [USES dE FROM] Using Symmetry
  AR_reaction[{"N", "H1+"}] =
      AR_parameters{"ar-n+h1", 1,      4.55E-12, -2.9E-01, -0.92,
                    -8.38,     1.00e2, 5.00e4,   0.936,    -0.94};

  // [USES dE_fit FROM] Ref: J.B. Kingdon & G.J. Ferland, ApJS, 106:205 (1996)
  // [USES dE FROM] Using Symmetry
  AR_reaction[{"O", "H1+"}] =
      AR_parameters{"ar-o+h1", 1,      7.40E-11, 4.7E-01, 24.37,
                    -0.74,     1.00e1, 1.00e4,   0.0198,  -0.02};

  // [USES dE_fit FROM] Ref: J.B. Kingdon & G.J. Ferland, ApJS, 106:205 (1996)
  // [USES dE FROM] Ref: S.E. Butler & J.C. Raymond, ApJ, 240:680 (1980)
  AR_reaction[{"Si1+", "H1+"}] =
      AR_parameters{"ar-si1+h1", 1,      4.10E-10, 0.24,  3.17,
                    4.18e-3,     5.00e3, 1.00e6,   2.739, -2.74};

  // [USES dE_fit FROM] Ref: J.B. Kingdon & G.J. Ferland, ApJS, 106:205 (1996)
  AR_reaction[{"S", "H1+"}] = AR_parameters{
      "ar-s+h1", 1, 1.00e-14, 0.00, 0.00, 0.00, 1.00e3, 1.00e4, 0.0, 0.0};

  // [USES dE_fit FROM] Ref: J.B. Kingdon & G.J. Ferland, ApJS, 106:205 (1996)
  AR_reaction[{"Fe1+", "H1+"}] =
      AR_parameters{"ar-fe1+h1", 1,      2.10e-9, 7.72e-2, -0.41,
                    -7.31,       1.00e4, 1.00e5,  2.590,   0.00};
  //****************************************************************************

  // [RECOMBINATION REACTIONS WITH He *****************************************]
  // The fit parameters for the recombination charge exchange reaction with He
  // are taken from M. Arnaud & R. Rothenflug Astron. Astrophys. Supple. Ser.
  // 60, 425 (1985). In this case, the fit formula is same as in the
  // recombination charge exchange reaction used above.

  // [USES dE FROM] Refs: S.E. Butler, T. G Heil & A. Dalgarno, ApJ, 241, 442
  // (1980),  G. Steigman, ApJ, 199, 336 (1975). Delta E = - 0.2.
  // Important Note: The upper and lower temperature limits of 10^4 K, as
  // reported in M. Arnaud & R. Rothenflug, Astron. Astrophys. Suppl. Ser. 60,
  // 425 (1985), are inappropriate for calculating the rate. We replace the
  // lower temperature limit with 300 K, which is suitable since this reaction
  // becomes important above this temperature.
  AR_reaction[{"C2+", "He"}] =
      AR_parameters{"ar-c2+he", 0,        1.00E-10, 0.00E+00, 0.00E+00,
                    0.00E+00,   3.00E+02, 1.00E+04, 0.0,      -0.2};

  // [RATE CHECK WITH] Ref: S.E. Butler & A. Dalgarno, ApJ, 241, 838 (1980) [OK]
  // [USES dE FROM] Ref: B.Friedman et al., J.Phys.B:At.Mol.Opt,50:115202 (2017)
  AR_reaction[{"C3+", "He"}] =
      AR_parameters{"ar-c3+he", 0,        4.60E-11, 2.00E+00, 0.00E+00,
                    0.00E+00,   1.00E+03, 3.00E+04, 0.0,      16.81};

  // [RATE CHECK WITH] Ref: S.E. Butler & A. Dalgarno, ApJ, 241, 838 (1980) [OK]
  // [USES dE FROM] Ref: B.Friedman et al., J.Phys.B:At.Mol.Opt,50:115202 (2017)
  AR_reaction[{"C4+", "He"}] =
      AR_parameters{"ar-c4+he", 0,        1.00E-14, 0.00E+00, 0.00E+00,
                    0.00E+00,   1.00E+03, 3.00E+04, 0.0,      31.89};

  // [RATE CHECK WITH] Ref: S.E. Butler & A. Dalgarno, ApJ, 241, 838 (1980) [OK]
  AR_reaction[{"N2+", "He"}] =
      AR_parameters{"ar-n2+he", 0,        3.30E-10, 2.90E-01, 1.30E+00,
                    -4.50E+00,  1.00E+03, 3.00E+04, 0.0,      0.00};

  // [RATE CHECK WITH] Ref: S.E. Butler & A. Dalgarno, ApJ, 241, 838 (1980) [OK]
  AR_reaction[{"N3+", "He"}] =
      AR_parameters{"ar-n3+he", 0,        1.50E-10, 0.00E+00, 0.00E+00,
                    0.00E+00,   1.00E+03, 3.00E+04, 0.0,      0.00};

  // [RATE CHECK WITH] Ref: S.E. Butler & A. Dalgarno, ApJ, 241, 838 (1980) [OK]
  AR_reaction[{"N4+", "He"}] =
      AR_parameters{"ar-n4+he", 0,        1.70E-09, 0.00E+00, 2.50E+00,
                    -3.70E+01,  1.00E+03, 3.00E+04, 0.0,      0.00};

  // [RATE CHECK PENDING]
  // [USES dE FROM] Ref: B.Friedman et al., J.Phys.B:At.Mol.Opt,50:115202 (2017)
  AR_reaction[{"O2+", "He"}] =
      AR_parameters{"ar-o2+he", 0,        2.00E-10, 9.50E-01, 0.00E+00,
                    0.00E+00,   5.00E+03, 5.00E+04, 0.0,      10.52};

  // [RATE CHECK WITH] Ref: S.E. Butler & A. Dalgarno, ApJ, 241, 838 (1980) [OK]
  // [USES dE FROM] Ref: B.Friedman et al., J.Phys.B:At.Mol.Opt,50:115202 (2017)
  AR_reaction[{"O3+", "He"}] =
      AR_parameters{"ar-o3+he", 0,        1.00E-09, 0.00E+00, 1.25E+00,
                    -5.80E+00,  1.00E+03, 3.00E+04, 0.0,      25.00};

  // [RATE CHECK WITH] Ref: S.E. Butler & A. Dalgarno, ApJ, 241, 838 (1980) [OK]
  AR_reaction[{"O4+", "He"}] =
      AR_parameters{"ar-o4+he", 0,        6.40E-10, 0.00E+00, 2.00E+00,
                    -5.50E+00,  1.00E+03, 3.00E+04, 0.0,      0.00};

  // [RATE CHECK WITH] Ref: S.E. Butler & A. Dalgarno, ApJ, 241, 838 (1980) [OK]
  // [USES dE FROM] Ref: B.Friedman et al., J.Phys.B:At.Mol.Opt,50:115202 (2017)
  AR_reaction[{"Ne2+", "He"}] =
      AR_parameters{"ar-ne2+he", 0,        1.00E-14, 0.00E+00, 0.00E+00,
                    0.00E+00,    1.00E+03, 3.00E+04, 0.0,      16.37};

  // [RATE CHECK WITH] Ref: S.E. Butler & A. Dalgarno, ApJ, 241, 838 (1980) [OK]
  AR_reaction[{"Ne3+", "He"}] =
      AR_parameters{"ar-ne3+he", 0,        1.00E-14, 5.10E-01, 0.00E+00,
                    0.00E+00,    1.00E+03, 3.00E+04, 0.0,      0.00};

  // [RATE CHECK WITH] Ref: S.E. Butler & A. Dalgarno, ApJ, 241, 838 (1980) [OK]
  AR_reaction[{"Ne4+", "He"}] =
      AR_parameters{"ar-ne4+he", 0,        1.70E-09, 5.20E-01, 3.30E+00,
                    -5.30E+00,   1.00E+03, 3.00E+04, 0.00,     0.00};

  // [RATE CHECK WITH] Ref: S.E. Butler & A. Dalgarno, ApJ, 241, 838 (1980) [OK]
  // [USES dE FROM] Using Symmetry
  AR_reaction[{"Si3+", "He"}] =
      AR_parameters{"ar-si3+he", 0,        9.50E-10, 7.50E-01, 0.00E+00,
                    0.00E+00,    1.00E+03, 3.00E+04, 0.0,      8.88};

  // [RATE CHECK WITH] Ref: S.E. Butler & A. Dalgarno, ApJ, 241, 838 (1980) [OK]
  // dE Ref: S. L. Balliunas and S. E. Butler, ApJ, 235:L45 (1980)
  AR_reaction[{"Si4+", "He"}] =
      AR_parameters{"ar-si4+he", 0,        1.20E-09, 0.00E+00, 0.00E+00,
                    0.00E+00,    1.00E+03, 3.00E+04, 0.0,      0.00};

  // [RATE CHECK WITH] Ref: S.E. Butler & A. Dalgarno, ApJ, 241, 838 (1980) [OK]
  AR_reaction[{"S3+", "He"}] =
      AR_parameters{"ar-s3+he", 0,        1.10E-09, 5.60E-01, 0.00E+00,
                    0.00E+00,   1.00E+03, 3.00E+04, 0.0,      0.00};

  // [RATE CHECK WITH] Ref: S.E. Butler & A. Dalgarno, ApJ, 241, 838 (1980) [OK]
  AR_reaction[{"S4+", "He"}] =
      AR_parameters{"ar-s4+he", 0,        7.60E-13, 3.20E-01, 3.40E+00,
                    -5.25E+00,  1.00E+03, 3.00E+04, 0.0,      0.00};
  //****************************************************************************

  // [IONISATION REACTIONS WITH He1+ ******************************************]
  // The fit parameters and energy defect for the ionisation charge exchange
  // reaction with He1+ are taken from https://www.pa.uky.edu/~verner/chtr.html.
  // The fit formula for these reactions is given in M. Arnaud & R. Rothenflug
  // Astron. Astrophys. Supple. Ser. 60, 425 (1985).

  // [USES dE FROM] Ref: S.E. Butler, T. G Heil & A. Dalgarno, ApJ, 241, 442
  // (1980), Delta E = -6.29 [ENDOTHERMIC REACTIONS]
  // [USES dE_fit FROM] Ref: M. Arnaud & R. Rothenflug, A&AS, 60:425 (1985)
  // [RATE CHECK WITH] Ref: S.E. Butler, T. G Heil & A. Dalgarno, ApJ, 241, 442
  // (1980) [OK]
  AR_reaction[{"C1+", "He1+"}] =
      AR_parameters{"ar-c1+he1", 1,        5.00E-12, 2.00E+00, 7.00E-02,
                    0.00,        5.00E+03, 5.00E+04, 6.29,     -6.29};

  // [USES dE_fit FROM] Ref: M. Arnaud & R. Rothenflug, A&AS, 60:425 (1985)
  // [RATE CHECK WITH] Ref: S.E. Butler & A. Dalgarno, ApJ, 241, 838 (1980) [OK]
  AR_reaction[{"N1+", "He1+"}] =
      AR_parameters{"ar-n1+he1", 1,        3.70E-12, 2.10E+00, 6.30E-02,
                    0.00,        1.00E+04, 3.00E+05, 1.44,     0.00};

  // [USES dE_fit FROM] Ref: S.E. Butler & A. Dalgarno, ApJ, 241, 838 (1980)
  // [RATE CHECK WITH] Ref: S.E. Butler & A. Dalgarno, ApJ, 241, 838 (1980) [OK]
  AR_reaction[{"Si1+", "He1+"}] =
      AR_parameters{"ar-si1+he1", 1,        1.50E-10, 2.40E-01, 0.00E+00,
                    0.00,         1.00E+04, 3.00E+05, 6.91,     0.00};

  // [USES dE_fit FROM] Ref: S.L. Balliunas & S.E. Butler, ApJ, 235:L45 (1980)
  // [RATE CHECK WITH] Ref: S.E. Butler & A. Dalgarno, ApJ, 241, 838 (1980) [OK]
  // [USES dE FROM] Ref: S.E. Butler & J.C. Raymond, ApJ, 240:680 (1980)
  AR_reaction[{"Si2+", "He1+"}] =
      AR_parameters{"ar-si2+he1", 1,        1.15E-09, 4.40E-01, 0.00E+00,
                    0.00,         1.00E+04, 3.00E+05, 8.88,     -8.88};

  // [USES dE_fit FROM] Ref: S.E. Butler & A. Dalgarno, ApJ, 241, 838 (1980)
  // [RATE CHECK WITH] Ref: S.E. Butler & A. Dalgarno, ApJ, 241, 838 (1980) [OK]
  AR_reaction[{"S1+", "He1+"}] =
      AR_parameters{"ar-s1+he1", 1,        2.80E-11, 1.20E+00, 3.60E-02,
                    0.00,        1.00E+04, 3.00E+05, 9.17,     0.00};

  // [USES dE_fit FROM] Ref: S.E. Butler & A. Dalgarno, ApJ, 241, 838 (1980)
  // [RATE CHECK WITH] Ref: S.E. Butler & A. Dalgarno, ApJ, 241, 838 (1980) [OK]
  AR_reaction[{"S2+", "He1+"}] =
      AR_parameters{"ar-s2+he1", 1,        1.40E-11, 1.60E+00, 4.60E-02,
                    0.00,        1.00E+04, 3.00E+05, 10.47,    0.00};
  // ***************************************************************************

  /** The following reactions are not included **/
  // HYDROGEN ################################################################
  UMIST_reaction_parameters[{"H", "He1+"}] = UMIST_datatype{
      "UMIST-H+He1+", 1.20e-15, 0.25, 0.00, 10, 41000, 10.99};  // #502 [1]
  // Reaction: H1+ + He -> H + He1+, Not in Database

  // Reaction: H + C1+ - > H1+ + C,  Not in Database
  // Reaction: H1+ + C - > H + C1+,  Not in Database

  // Reaction: H + N1+ - > H1+ + N,  Not in Database
  // Reaction: H1+ + N - > H + N1+,  Not in Database

  UMIST_reaction_parameters[{"H", "O1+"}] = UMIST_datatype{
      "UMIST-H+O1+", 5.66e-10, 0.36, -8.60, 10, 41000, 0.02};  // 503, [2]
  UMIST_reaction_parameters[{"H1+", "O"}] = UMIST_datatype{
      "UMIST-H1++O", 6.86e-10, 0.26, 224.30, 10, 41000, -0.02};  // #417, [3]

  // Reaction: H + Ne1+ -> H1+ + Ne, Not in Database
  // Reaction: H1+ + Ne -> H + Ne1+, Not in Database

  // Reaction: H + Si1+ -> H1+ + Si, Not in Database
  UMIST_reaction_parameters[{"H1+", "Si"}] = UMIST_datatype{
      "UMIST-H1++Si", 9.90e-10, 0.00, 0.00, 10, 41000, 0.12};  // #429, [4]

  // Reaction: H + S1+ - > H1+ + S,  Not in Database
  UMIST_reaction_parameters[{"H1+", "S"}] = UMIST_datatype{
      "UMIST-H1++S", 1.30e-9, 0.00, 0.00, 10, 41000};  // #426, [5]

  // Reaction: H + Fe1+ -> H1+ + Fe, Not in Database
  UMIST_reaction_parameters[{"H1+", "Fe"}] = UMIST_datatype{
      "UMIST-H1++Fe", 7.40e-9, 0.00, 0.00, 10, 41000};  // #385, [6]

  // HELIUM ##################################################################
  // Reaction: He + H1+ -> He1+ + H, Already in Database.

  // Reaction: He + C1+ -> He1+ + C, Not in Database
  UMIST_reaction_parameters[{"He1+", "C"}] = UMIST_datatype{
      "UMIST-He1++C", 6.30e-15, 0.75, 0.00, 10, 300};  // #524, [7]

  // Reaction: He + N1+ - > He1+ + N,  Not in Database
  // Reaction: He1+ + N - > He + N1+,  Not in Database

  // Reaction: He + O1+ - > H1+ + O,   Not in Database
  // Reaction: He1+ + O - > He + O1+,  Not in Database

  // Reaction: He + Ne1+ -> He1+ + Ne, Not in Database
  // Reaction: He1+ + Ne -> He + Ne1+, Not in Database

  // Reaction: He + Si1+ -> He1+ + Si, Not in Database
  UMIST_reaction_parameters[{"He1+", "Si"}] = UMIST_datatype{
      "UMIST-He1++Si", 3.30e-9, 0.00, 0.00, 10, 41000};  // #539 [8]

  // Reaction: He + S1+ - > He1+ + S,  Not in Database
  // Reaction: He1+ + S ->  He + S1+,  Not in Database

  // Reaction: He + Fe1+ -> He1+ + Fe, Not in Database
  // Reaction: He1+ + Fe -> He + Fe1+, Not in Database

  // CARBON ##################################################################
  // Reaction: C + H1+ -> C1+ + H, Not in Database
  // Reaction: C1+ + H -> C + H1+, Not in Database

  // Reaction: C + He1+ -> C+ + He,  Already in Database
  // Reaction: C1+ + He -> C + He1+, Not in Database

  // Reaction: C + N1+ - > C1+ + N,  Not in Database
  // Reaction: C1+ + N - > C + N1+,  Not in Database

  // Reaction: C + O1+ - > C1+ + O,  Not in Database
  // Reaction: C1+ + O - > C + O1+,  Not in Database

  // Reaction: C + Ne1+ -> C1+ + Ne, Not in Database
  // Reaction: C1+ + Ne -> C + Ne1+, Not in Database

  // Reaction: C + Si1+ -> C1+ + Si, Not in Database
  UMIST_reaction_parameters[{"C1+", "Si"}] = UMIST_datatype{
      "UMIST-C1++Si", 2.10e-9, 0.00, 0.00, 10, 41000};  // #203, [9]

  // Reaction: C + S1+ -> C1+ + S, Not in Database
  UMIST_reaction_parameters[{"C1+", "S"}] = UMIST_datatype{
      "UMIST-C1++S", 5.00e-11, 0.00, 0.00, 10, 4999};  // #717, [10]
  UMIST_reaction_parameters[{"C1+", "S"}] = UMIST_datatype{
      "UMIST-C1++S", 5.54e-11, 0.86, 680.70, 5000, 50000};  // #717, [10]

  // Reaction: C + Fe1+ -> C1+ + Fe, Not in Database
  UMIST_reaction_parameters[{"C1+", "Fe"}] = UMIST_datatype{
      "UMIST-C1++Fe", 2.60e-9, 0.00, 0.00, 10, 41000};  // #184 [11]

  // NITROGEN ################################################################
  // Reaction: N + H1+ - > N1+ + H,  Not in Database
  // Reaction: N1+ + H - > N + H1+,  Not in Database

  // Reaction: N + He1+ -> N1+ + He, Not in Database
  // Reaction: N1+ + He -> N + He1+, Not in Database

  // Reaction: N + C1+ - > N1+ + C,  Not in Database
  // Reaction: N1+ + C - > N + C1+,  Not in Database

  // Reaction: N + O1+ - > N1+ + O,  Not in Database
  // Reaction: N1+ + O - > N + O1+,  Not in Database

  // Reaction: N + Ne1+ -> N1+ + Ne, Not in Database
  // Reaction: N1+ + Ne -> N + Ne1+, Not in Database

  // Reaction: N + Si1+ -> N1+ + Si, Not in Database
  // Reaction: N1+ + Ne -> N + Ne1+, Not in Database

  // Reaction: N + S1+ - > N1+ + S,  Not in Database
  // Reaction: N1+ + Ne -> N + Ne1+, Not in Database

  // Reaction: N + Fe1+ -> N1+ + Fe, Not in Database
  UMIST_reaction_parameters[{"N1+", "Fe"}] = UMIST_datatype{
      "UMIST-N1++Fe", 1.50e-9, 0.00, 0.00, 10, 41000};  // #561, [12]

  // OXYGEN ##################################################################
  // Reaction: O + H1+ -> O+ + H, Already in Database
  // Reaction: O1+ + H -> O + H+, Already in Database

  // Reaction: O + He1+ -> O1+ + He, Not in Database
  // Reaction: O1+ + He -> O + He1+, Not in Database

  // Reaction: O + C1+ - > O1+ + C,  Not in Database
  // Reaction: O1+ + C - > O + C1+,  Not in Database

  // Reaction: O + N1+ - > O1+ + N,  Not in Database
  // Reaction: O1+ + N - > O + N1+,  Not in Database

  // Reaction: O + Ne1+ -> O1+ + Ne, Not in Database
  // Reaction: O1+ + Ne -> O + Ne1+, Not in Database

  // Reaction: O + Si1+ -> O1+ + Si, Not in Database
  // Reaction: O1+ + Si -> O + Si1+, Not in Database

  // Reaction: O + S1+ - > O1+ + S,  Not in Database
  // Reaction: O1+ + S - > O + S1+,  Not in Database

  // Reaction: O + Fe1+ -> O1+ + Fe, Not in Database
  UMIST_reaction_parameters[{"O1+", "Fe"}] = UMIST_datatype{
      "UMIST-O1++Fe", 2.90e-9, 0.00, 0.00, 10, 41000};  // #665, [13]

  // NEON ####################################################################
  // Reaction: Ne + H1+ -> Ne1+ + H, Not in Database
  // Reaction: Ne1+ + H -> Ne + H1+, Not in Database

  // Reaction: Ne + He1+ -> Ne1+ + He, Not in Database
  // Reaction: Ne1+ + He -> Ne + He1+, Not in Database

  // Reaction: Ne + C1+ - > Ne1+ + C,  Not in Database
  // Reaction: Ne1+ + C - > Ne + C1+,  Not in Database

  // Reaction: Ne + N1+ - > Ne1+ + N,  Not in Database
  // Reaction: Ne1+ + N - > Ne + N1+,  Not in Database

  // Reaction: Ne + O1+ - > Ne1+ + O,  Not in Database
  // Reaction: Ne1+ + O - > Ne + O1+,  Not in Database

  // Reaction: Ne + Si1+ -> Ne1+ + Si, Not in Database
  // Reaction: Ne1+ + Si -> Ne + Si1+, Not in Database

  // Reaction: Ne + S1+ - > Ne1+ + S,  Not in Database
  // Reaction: Ne1+ + S - > Ne + S1+,  Not in Database

  // Reaction: O + Fe1+ -> Ne1+ + Fe,  Not in Database
  // Reaction: Ne1+ + Fe -> Ne + Fe1+, Not in Database

  // SILICON #################################################################
  // Reaction: Si + H+ -> Si1+ + H, Already in Database.
  // Reaction: Si1+ + H -> Si + H1+, Not in Database

  // Reaction: Si + He+ -> Si1+ + He, Already in Database.
  // Reaction: Si1+ + He -> Si + He1+, Not in Database

  // Reaction: Si + C+ -> Si1+ + C, Already in Database.
  // Reaction: Si1+ + C - > Si + C1+,  Not in Database

  // Reaction: Si + N1+ - > Si1+ + N,  Not in Database
  // Reaction: Si1+ + N - > Si + N1+,  Not in Database

  // Reaction: Si + O1+ - > Si1+ + O,  Not in Database
  // Reaction: Si1+ + O - > Si + O1+,  Not in Database

  // Reaction: Si + Ne1+ -> Si1+ + Ne, Not in Database
  // Reaction: Si1+ + Ne -> Si + Ne1+, Not in Database

  UMIST_reaction_parameters[{"Si", "S1+"}] = UMIST_datatype{
      "UMIST-Si+S1+", 1.60e-9, 0.00, 0.00, 10, 41000};  // #733, [8]
  // Reaction: Si1+ + S - > Si + S1+,  Not in Database

  // Reaction: Si + Fe1+ -> Si1+ + Fe,  Not in Database
  UMIST_reaction_parameters[{"Si1+", "Fe"}] = UMIST_datatype{
      "UMIST-Si1++Fe", 1.90e-9, 0.00, 0.00, 10, 41000};  // #725, [9]

  // SULFUR ##################################################################
  // Reaction: S + H1+ -> S+ + H, Already in Database.
  // Reaction: S1+ + H - > S + H1+,  Not in Database

  // Reaction: S + He1+ -> S1+ + He, Not in Database
  // Reaction: S1+ + He -> S + He1+, Not in Database

  // Reaction: S + C1+ -> S1+ + C, Already in Database
  // Reaction: S1+ + C - > S + C1+,  Not in Database

  // Reaction: S + N1+ - > S1+ + N,  Not in Database
  // Reaction: S1+ + N - > S + N1+,  Not in Database

  // Reaction: S + O1+ - > S1+ + O,  Not in Database
  // Reaction: S1+ + O - > S + O1+,  Not in Database

  // Reaction: S + Ne1+ -> S1+ + Ne, Not in Database
  // Reaction: S1+ + Ne -> S + Ne1+, Not in Database

  // Reaction: S + Si1+ -> S1+ + Si, Not in Database
  // Reaction: S1+ + Si -> S + Si1+, Already in Database

  // Reaction: S + Fe1+ -> S1+ + Fe, Not in Database
  UMIST_reaction_parameters[{"S1+", "Fe"}] = UMIST_datatype{
      "UMIST-S1++Fe", 1.80e-10, 0.00, 0.00, 10, 41000};  // #714, [14]

  // IRON ####################################################################
  // Reaction: Fe + H1+ -> Fe1+ + H, Already in Database
  // Reaction: Fe1+ + H -> Fe + H1+,     Not in Database

  // Reaction: Fe + He1+ -> Fe1+ + He, Not in Database
  // Reaction: Fe1+ + He -> Fe + He1+, Not in Database

  // Reaction: Fe + C1+ -> Fe1+ + C, Already in Database
  // Reaction: Fe1+ + C - > Fe + C1+,  Not in Database

  // Reaction: Fe + N1+ -> Fe1+ + N, Already in Database
  // Reaction: Fe1+ + N - > Fe + N1+,  Not in Database

  // Reaction: Fe + O1+ -> Fe1+ + O, Already in Database
  // Reaction: Fe1+ + O - > Fe + O1+,  Not in Database

  // Reaction: Fe + Ne1+ -> Fe1+ + Ne, Not in Database
  // Reaction: Fe1+ + Ne -> Fe + Ne1+, Not in Database

  // Reaction: Fe + Si1+ -> Fe1+ + Si, Already in Database
  // Reaction: Fe1+ + Si -> Fe + Si1+, Not in Database

  // Reaction: Fe + S1+ -> Fe1+ + S, Already in Database
  // Reaction: Fe1+ + S - > Fe + S1+,  Not in Database
}
