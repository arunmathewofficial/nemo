/**
 *  \file     atomic_physics_database.h
 *  \author   Arun Mathew
 *  \date     12/12/2022
 **/

#include "atomic_physics_database.h"
#include <algorithm>

// Global instance of class
class atomic_physics_database atomic_physics_data;

// ##################################################################
// Constructor
atomic_physics_database::atomic_physics_database()
{
  database();
}

// ##################################################################
// Destructor
atomic_physics_database::~atomic_physics_database() {}

// ##################################################################
// Find the species in the database.
bool atomic_physics_database::find_species(const std::string species)
{

  auto it = species_attributes.find(species);
  if (it != species_attributes.end())
    return true;
  else {
    spdlog::error("Species not found in database {}", species);
    return false;
  }
}

// ##################################################################
// Export species of an element from the database.
std::vector<std::string> atomic_physics_database::import_species(
    const std::string element_name)
{
  // get all the species corresponding to an element.
  std::vector<std::string> element_species_set;
  std::string key = "";
  for (const auto &pair : species_attributes) {
    if (pair.second.element == element_name) {
      key = pair.first;
      element_species_set.push_back(key);
    }
  }

  // Sort the obtained string species vector according to their charge
  // Comparator function to compare species based on charges
  auto charge_comparator =
      [this](const std::string &species_1, const std::string &species_2) {
        return get_charge(species_1) < get_charge(species_2);
      };
  // Sorting the vector using the charge comparator
  std::sort(
      element_species_set.begin(), element_species_set.end(),
      charge_comparator);

  // if the element is not in the database return error.
  if (key == "") {
    spdlog::error("Element not found in database {}", element_name);
    exit_nemo(1);
  }

  return element_species_set;
}

// ##################################################################
// Get atomic element from the database
std::string atomic_physics_database::get_element_name(const std::string species)
{
  return species_attributes[species].element;
}

// ##################################################################
// Get atomic mass from the database
int atomic_physics_database::get_atomic_number(const std::string species)
{
  return species_attributes[species].atomic_number;
}

// ##################################################################
// Get atomic mass from the database
double atomic_physics_database::get_atomic_mass(const std::string species)
{
  return species_attributes[species].mass;
}

// ##################################################################
// Get the charge of the species
int atomic_physics_database::get_charge(const std::string species)
{
  return species_attributes[species].charge;
}

// ##################################################################
// Get number of electron released by the species
int atomic_physics_database::get_electron_num(const std::string species)
{
  return species_attributes[species].charge;
}

// ##################################################################
// Get ionisation potential in units of erg
double atomic_physics_database::get_ionisation_potential(
    const std::string species)
{
  return species_attributes[species].ionise_E * eV();
}

// ##################################################################
// Get ionisation potential in units of eV
double atomic_physics_database::get_ionisation_energy(const std::string species)
{
  return species_attributes[species].ionise_E;
}

// ##################################################################
// Atomic physics database
void atomic_physics_database::database()
{
  spdlog::debug("Setting up atomic physics database");

  // Ionisation energies are taken from NIST Atomic Spectra
  // Database (ver.5.10), Kramida, A., Ralchenko, Yu., Reader,
  // J., and NIST ASD Team (2022).
  // Available at DOI: https://doi.org/10.18434/T4W30F
  // Ionisation energies are given in the units of eV.

  // Entry: H
  species_attributes["H"] = species_datatype{"H", 1, 0, 1.6738e-24, 1.35984e+1};

  // Entry: H1+
  species_attributes["H1+"] = species_datatype{
      "H", 1, 1,
      0,       // mass of H1+ is set equal to zero.
      1.0e-99  // small dummy value.
  };

  // Entry: He
  species_attributes["He"] =
      species_datatype{"He", 2, 0, 6.6464768e-24, 2.45874e+1};

  // Entry: He1+
  species_attributes["He1+"] = species_datatype{"He", 2, 1, 0, 5.44178e+1};

  // Entry: He2+
  species_attributes["He2+"] = species_datatype{"He", 2, 2, 0, 1.0e-99};

  // Entry: C
  species_attributes["C"] =
      species_datatype{"C", 6, 0, 1.994374e-23, 1.12603e+1};

  // Entry: C1+
  species_attributes["C1+"] = species_datatype{"C", 6, 1, 0, 2.43831e+1};

  // Entry: C2+
  species_attributes["C2+"] = species_datatype{"C", 6, 2, 0, 4.78878e+1};

  // Entry: C3+
  species_attributes["C3+"] = species_datatype{"C", 6, 3, 0, 6.44935e+1};

  // Entry: C4+
  species_attributes["C4+"] = species_datatype{"C", 6, 4, 0, 3.92091e+2};

  // Entry: C5+
  species_attributes["C5+"] = species_datatype{"C", 6, 5, 0, 4.89993e+2};

  // Entry: C6+
  species_attributes["C6+"] = species_datatype{"C", 6, 6, 0, 1.0e-99};

  // Entry: N
  species_attributes["N"] =
      species_datatype{"N", 7, 0, 2.325892e-23, 1.45341e+1};

  // Entry: N1+
  species_attributes["N1+"] = species_datatype{"N", 7, 1, 0, 2.96013e+1};

  // Entry: N2+
  species_attributes["N2+"] = species_datatype{"N", 7, 2, 0, 4.74453e+1};

  // Entry: N3+
  species_attributes["N3+"] = species_datatype{"N", 7, 3, 0, 7.74735e+1};

  // Entry: N4+
  species_attributes["N4+"] = species_datatype{"N", 7, 4, 0, 9.78901e+1};

  // Entry: N5+
  species_attributes["N5+"] = species_datatype{"N", 7, 5, 0, 5.52067e+2};

  // Entry: N6+
  species_attributes["N6+"] = species_datatype{"N", 7, 6, 0, 6.67046e+2};

  // Entry: N7+
  species_attributes["N7+"] = species_datatype{"N", 7, 7, 0, 1.0e-99};

  // Entry: O
  species_attributes["O"] =
      species_datatype{"O", 8, 0, 2.6567628e-23, 1.36181e+1};

  // Entry: O1+
  species_attributes["O1+"] = species_datatype{"O", 8, 1, 0, 3.51211e+1};

  // Entry: O2+
  species_attributes["O2+"] = species_datatype{"O", 8, 2, 0, 5.49355e+1};

  // Entry: O3+
  species_attributes["O3+"] = species_datatype{"O", 8, 3, 0, 7.74135e+1};

  // Entry: O4+
  species_attributes["O4+"] = species_datatype{"O", 8, 4, 0, 1.13899e+2};

  // Entry: O5+
  species_attributes["O5+"] = species_datatype{"O", 8, 5, 0, 1.38119e+2};

  // Entry: O6+
  species_attributes["O6+"] = species_datatype{"O", 8, 6, 0, 7.39327e+2};

  // Entry: O7+
  species_attributes["O7+"] = species_datatype{"O", 8, 7, 0, 8.71410e+2};

  // Entry: O8+
  species_attributes["O8+"] = species_datatype{"O", 8, 8, 0, 1.e-99};

  // Entry: Ne
  species_attributes["Ne"] =
      species_datatype{"Ne", 10, 0, 3.3509177e-23, 2.15645e+1};

  // Entry: Ne1+
  species_attributes["Ne1+"] = species_datatype{"Ne", 10, 1, 0, 4.09630e+1};

  // Entry: Ne2+
  species_attributes["Ne2+"] = species_datatype{"Ne", 10, 2, 0, 6.34233e+1};

  // Entry: Ne3+
  species_attributes["Ne3+"] = species_datatype{"Ne", 10, 3, 0, 9.71900e+1};

  // Entry: Ne4+
  species_attributes["Ne4+"] = species_datatype{"Ne", 10, 4, 0, 1.26247e+2};

  // Entry: Ne5+
  species_attributes["Ne5+"] = species_datatype{"Ne", 10, 5, 0, 1.57934e+2};

  // Entry: Ne6+
  species_attributes["Ne6+"] = species_datatype{"Ne", 10, 6, 0, 2.07271e+2};

  // Entry: Ne7+
  species_attributes["Ne7+"] = species_datatype{"Ne", 10, 7, 0, 2.39097e+2};

  // Entry: Ne8+
  species_attributes["Ne8+"] = species_datatype{"Ne", 10, 8, 0, 1.19581e+3};

  // Entry: Ne9+
  species_attributes["Ne9+"] = species_datatype{"Ne", 10, 9, 0, 1.36220e+3};

  // Entry: Ne10+
  species_attributes["Ne10+"] = species_datatype{"Ne", 10, 10, 0, 1.0e-99};

  // Entry: Si
  species_attributes["Si"] =
      species_datatype{"Si", 14, 0, 4.6637066e-23, 8.15168e+0};

  // Entry: Si1+
  species_attributes["Si1+"] = species_datatype{"Si", 14, 1, 0, 1.63459e+1};

  // Entry: Si2+
  species_attributes["Si2+"] = species_datatype{"Si", 14, 2, 0, 3.34930e+1};

  // Entry: Si3+
  species_attributes["Si3+"] = species_datatype{"Si", 14, 3, 0, 4.51418e+1};

  // Entry: Si4+
  species_attributes["Si4+"] = species_datatype{"Si", 14, 4, 0, 1.66767e+2};

  // Entry: Si5+
  species_attributes["Si5+"] = species_datatype{"Si", 14, 5, 0, 2.05279e+2};

  // Entry: Si6+
  species_attributes["Si6+"] = species_datatype{"Si", 14, 6, 0, 2.46570e+2};

  // Entry: Si7+
  species_attributes["Si7+"] = species_datatype{"Si", 14, 7, 0, 3.03590e+2};

  // Entry: Si8+
  species_attributes["Si8+"] = species_datatype{"Si", 14, 8, 0, 3.51280e+2};

  // Entry: Si9+
  species_attributes["Si9+"] = species_datatype{"Si", 14, 9, 0, 4.01380e+2};

  // Entry: Si10+
  species_attributes["Si10+"] = species_datatype{"Si", 14, 10, 0, 4.76273e+2};

  // Entry: Si11+
  species_attributes["Si11+"] = species_datatype{"Si", 14, 11, 0, 5.23415e+2};

  // Entry: Si12+
  species_attributes["Si12+"] = species_datatype{"Si", 14, 12, 0, 2.43766e+3};

  // Entry: Si13+
  species_attributes["Si13+"] = species_datatype{"Si", 14, 13, 0, 2.67318e+3};

  // Entry: Si14+
  species_attributes["Si14+"] = species_datatype{"Si", 14, 14, 0, 1.0e-99};

  // Entry: S
  species_attributes["S"] =
      species_datatype{"S", 16, 0, 5.3245181e-23, 1.03600e+1};

  // Entry: S1+
  species_attributes["S1+"] = species_datatype{"S", 16, 1, 0, 2.33379e+1};

  // Entry: S2+
  species_attributes["S2+"] = species_datatype{"S", 16, 2, 0, 3.48600e+1};

  // Entry: S3+
  species_attributes["S3+"] = species_datatype{"S", 16, 3, 0, 4.72220e+1};

  // Entry: S4+
  species_attributes["S4+"] = species_datatype{"S", 16, 4, 0, 7.25945e+1};

  // Entry: S5+
  species_attributes["S5+"] = species_datatype{"S", 16, 5, 0, 8.80529e+1};

  // Entry: S6+
  species_attributes["S6+"] = species_datatype{"S", 16, 6, 0, 2.80954e+2};

  // Entry: S7+
  species_attributes["S7+"] = species_datatype{"S", 16, 7, 0, 3.28794e+2};

  // Entry: S8+
  species_attributes["S8+"] = species_datatype{"S", 16, 8, 0, 3.79840e+2};

  // Entry: S9+
  species_attributes["S9+"] = species_datatype{"S", 16, 9, 0, 4.47700e+2};

  // Entry: S10+
  species_attributes["S10+"] = species_datatype{"S", 16, 10, 0, 5.04550e+2};

  // Entry: S11+
  species_attributes["S11+"] = species_datatype{"S", 16, 11, 0, 5.64410e+2};

  // Entry: S12+
  species_attributes["S12+"] = species_datatype{"S", 16, 12, 0, 6.51960e+2};

  // Entry: S13+
  species_attributes["S13+"] = species_datatype{"S", 16, 13, 0, 7.06994e+2};

  // Entry: S14+
  species_attributes["S14+"] = species_datatype{"S", 16, 14, 0, 3.22378e+3};

  // Entry: S15+
  species_attributes["S15+"] = species_datatype{"S", 16, 15, 0, 3.49419e+3};

  // Entry: S16+
  species_attributes["S16+"] = species_datatype{"S", 16, 16, 0, 1.0e-99};

  // Entry: Fe
  species_attributes["Fe"] =
      species_datatype{"Fe", 26, 0, 9.2732796e-23, 7.90247e+0};

  // Entry: Fe1+
  species_attributes["Fe1+"] = species_datatype{"Fe", 26, 1, 0, 1.61992e+1};

  // Entry: Fe2+
  species_attributes["Fe2+"] = species_datatype{"Fe", 26, 2, 0, 3.06510e+1};

  // Entry: Fe3+
  species_attributes["Fe3+"] = species_datatype{"Fe", 26, 3, 0, 5.49100e+1};

  // Entry: Fe4+
  species_attributes["Fe4+"] = species_datatype{"Fe", 26, 4, 0, 7.50000e+1};

  // Entry: Fe5+
  species_attributes["Fe5+"] = species_datatype{"Fe", 26, 5, 0, 9.89850e+1};

  // Entry: Fe6+
  species_attributes["Fe6+"] = species_datatype{"Fe", 26, 6, 0, 1.24967e+2};

  // Entry: Fe7+
  species_attributes["Fe7+"] = species_datatype{"Fe", 26, 7, 0, 1.51060e+2};

  // Entry: Fe8+
  species_attributes["Fe8+"] = species_datatype{"Fe", 26, 8, 0, 2.33600e+2};

  // Entry: Fe9+
  species_attributes["Fe9+"] = species_datatype{"Fe", 26, 9, 0, 2.62100e+2};

  // Entry: Fe10+
  species_attributes["Fe10+"] = species_datatype{"Fe", 26, 10, 0, 2.90900e+2};

  // Entry: Fe11+
  species_attributes["Fe11+"] = species_datatype{"Fe", 26, 11, 0, 3.30800e+2};

  // Entry: Fe12+
  species_attributes["Fe12+"] = species_datatype{"Fe", 26, 12, 0, 3.61000e+2};

  // Entry: Fe13+
  species_attributes["Fe13+"] = species_datatype{"Fe", 26, 13, 0, 3.92200e+2};

  // Entry: Fe14+
  species_attributes["Fe14+"] = species_datatype{"Fe", 26, 14, 0, 4.56200e+2};

  // Entry: Fe15+
  species_attributes["Fe15+"] = species_datatype{"Fe", 26, 15, 0, 4.89312e+2};

  // Entry: Fe16+
  species_attributes["Fe16+"] = species_datatype{"Fe", 26, 16, 0, 1.26270e+3};

  // Entry: Fe17+
  species_attributes["Fe17+"] = species_datatype{"Fe", 26, 17, 0, 1.35780e+3};

  // Entry: Fe18+
  species_attributes["Fe18+"] = species_datatype{"Fe", 26, 18, 0, 1.46000e+3};

  // Entry: Fe19+
  species_attributes["Fe19+"] = species_datatype{"Fe", 26, 19, 0, 1.57560e+3};

  // Entry: Fe20+
  species_attributes["Fe20+"] = species_datatype{"Fe", 26, 20, 0, 1.68700e+3};

  // Entry: Fe21+
  species_attributes["Fe21+"] = species_datatype{"Fe", 26, 21, 0, 1.79840e+3};

  // Entry: Fe22+
  species_attributes["Fe22+"] = species_datatype{"Fe", 26, 22, 0, 1.95040e+3};

  // Entry: Fe23+
  species_attributes["Fe23+"] = species_datatype{"Fe", 26, 23, 0, 2.04576e+3};

  // Entry: Fe24+
  species_attributes["Fe24+"] = species_datatype{"Fe", 26, 24, 0, 8.82819e+3};

  // Entry: Fe25+
  species_attributes["Fe25+"] = species_datatype{"Fe", 26, 25, 0, 9.27769e+3};

  // Entry: Fe26+
  species_attributes["Fe26+"] = species_datatype{"Fe", 26, 26, 0, 1.0e-99};
}

// **************************************************************
