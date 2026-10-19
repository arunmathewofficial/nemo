/**
 *  \file      charge_exchange.h
 *  \brief     Charge exchange reaction class.
 *  \author    Arun Mathew
 *  \date      10-03-2024
 **/

#ifndef PION_CHARGE_EXCHANGE_H
#define PION_CHARGE_EXCHANGE_H

#include <spdlog/spdlog.h>
/* prevent clang-format reordering */
#include "atomic_physics_database.h"
#include <fmt/ranges.h>
#include <map>
#include <variant>

// UMIST charge exchange datatype
struct UMIST_datatype {
  std::string reaction_id;
  double alpha;  ///< Charge exchange fit parameter alpha (cm^3/s)
  double beta;   ///< Charge exchange fit parameter beta  (dimensionless)
  double gamma;  ///< Charge exchange fit parameter gamma X (K)
  double T_min;
  double T_max;
  double Delta_E;  ///< Energy-defect of reactions
};

// Arnaud and Rothenflug charge exchange rate parameters
struct AR_parameters {
  std::string reaction_id;
  bool type;
  double a;
  double b;
  double c;
  double d;
  double T_min;
  double T_max;
  double dE_fit;  /// Engery defect for fit formula for H1+ (10^4 K) and for
                  /// He1+ (eV)
  double dE;      ///< Energy-defect of reactions (eV)
};

class charge_exchange : virtual public atomic_physics_database {
public:
  /** \brief Constructor: Setup charge exchange database. **/
  charge_exchange(const int);
  /** \brief Destructor */
  virtual ~charge_exchange();

  void setup_charge_exchange(
      const std::vector<std::vector<std::pair<std::string, std::string> > >
          &,  ///< relevant ionising reactions
      const std::vector<std::vector<std::pair<std::string, std::string> > >
          &  ///< relevant recombining reactions
  );

  /** \brief This function exports element list from NEMO module. **/
  void generate_charge_exchange_LUT(
      const std::vector<double> &  ///< Temperature table
  );

  /** \brief  **/
  bool reactant_exists_in_reference(
      const std::string &,              ///< species
      const std::vector<std::string> &  ///< cx reference species list
  );

  /** \brief  **/
  std::vector<std::pair<std::string, std::string> > get_cx_reactions(
      const std::string &,  ///< ionisation/recombination species list
      const bool,           ///< type of reaction
      const std::vector<std::string> &  ///< cx reference species list
  );

  /** \brief  **/
  std::pair<std::string, std::string> get_actual_key(
      const std::pair<std::string, std::string> &);

  /** \brief  **/
  std::vector<std::vector<std::string> > import_reactions_id(
      std::vector<std::vector<std::pair<std::string, std::string> > > &);

  /** \brief  **/
  std::vector<std::vector<double> > import_energy_defect(
      const std::vector<std::vector<std::pair<std::string, std::string> > > &);

  std::vector<std::vector<std::pair<std::string, std::string> > >
      ionise_reactions;
  std::vector<std::vector<std::pair<std::string, std::string> > >
      recomb_reactions;

  std::vector<std::vector<std::vector<double> > >
      ionise_cx_rate_table;  ///< Charge exchange rate LUT
  std::vector<std::vector<std::vector<double> > >
      ionise_cx_slope_table;  ///< Charge exchange rate slope LUT
  std::vector<std::vector<std::vector<double> > >
      recomb_cx_rate_table;  ///< Charge exchange rate LUT
  std::vector<std::vector<std::vector<double> > >
      recomb_cx_slope_table;  ///< Charge exchange rate slope LUT

  // PRIVATE METHODS AND VARIABLES *************************
private:
  int N_ionise_species;
  int N_recomb_species;
  int q_max;

  /// Function dedicated to managing reaction parameters.
  void charge_exchange_database();
  /// Map reaction parameters to the reactants (species symbols).
  std::map<std::pair<std::string, std::string>, UMIST_datatype>
      UMIST_reaction_parameters;

  std::map<std::pair<std::string, std::string>, AR_parameters> AR_reaction;

  double get_charge_exchange_rate(
      const std::pair<std::string, std::string> &,  ///< Reaction
      const double                                  ///< Temperature
  );
};
#endif  // PION_CHARGE_EXCHANGE_H
