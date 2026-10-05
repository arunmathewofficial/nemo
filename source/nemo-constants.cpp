///
/// \file nemo-constants.cpp
/// \author Jonathan Mackey
/// \date 2026.10.01
///
/// Class with inline functions for many physics and astronomy
/// constants.  All constants are in cgs units unless stated
/// otherwise.
///
/// Modifications:
/// - modified version of PION constants.cpp from 2015

#ifdef SPDLOG_FWD
#include <spdlog/fwd.h>
#endif
#include <spdlog/spdlog.h>
/* prevent clang-format reordering */

#include "nemo-constants.h"

using namespace std;

// ##################################################################
// ##################################################################

nemo_constants::nemo_constants()
nemo_constants::nemo_constants():
  c_cgs(units::constants::c.value_as(units::precise::cgs::cm / units::precise::s)),
  kB_cgs(units::constants::k.value_as(units::precise::cgs::erg / units::precise::K)),
  h_cgs(units::constants::h.value_as(units::precise::cgs::erg * units::precise::s)),
  StefanBoltzmann_cgs(units::constants::s.value_as(units::precise::cgs::erg / (units::precise::cgs::cm * units::precise::cgs::cm * units::precise::cgs::cm * units::precise::K * units::precise::K * units::precise::K * units::precise::K))),
  eV_cgs(units::constants::e.value_as(units::precise::cgs::erg)),
  ThomsonXSection_cgs(6.6524587321e-25), // This constant doesn't have units in LLNL, keep as is
  uamu_cgs(1.660539040e-24), // Keep as is
  m_e_cgs(9.1093837015e-28), // Keep as is
  m_p_cgs(1.672621898e-24), // Keep as is
  m_H_cgs(1.6738e-24), // Keep as is
  m_He_cgs(6.6464768e-24), // Keep as is
  m_C_cgs(1.994374e-23), // Keep as is
  m_N_cgs(2.325892e-23), // Keep as is
  m_O_cgs(2.6567628e-23), // Keep as is
  G_cgs(6.67408e-8), // Keep as is
  Eth_H_cgs(2.178710264e-11), // Keep as is
  NuTh_H_cgs(3.28808819e+15), // Keep as is
  K_per_eV_cgs(1.16045221e4), // Keep as is
  year_cgs(3.1558150e7), // Keep as is
  asec_per_rad_cgs(206264.806), // Keep as is
  sqasec_per_sr_cgs(4.2545250225e10), // Keep as is
  Msun_cgs(1.9891e33), // Keep as is
  Lsun_cgs(3.839e33), // Keep as is
  Rsun_cgs(6.955e10), // Keep as is
  AU_cgs(1.49597870700e+13), // Keep as is
  parsec_cgs(3.0856776e18) // Keep as is
{}




// ##################################################################
// ##################################################################

nemo_constants::~nemo_constants() {}


// ##################################################################
// ##################################################################



bool nemo_constants::equalD(const double a, const double b)
{
  if (a == b) {
    return (true);
  }
  if (fabs(a) + fabs(b) < TINYVALUE) {
#ifdef DEBUG3
    spdlog::debug(
        "tiny numbers in equalD(a,b); a,b <1.e-100... a={}, b={}; returning true.\n",
        a, b);
#endif
    return (true);
  }
  if ((fabs(a - b) / (fabs(a) + fabs(b) + TINYVALUE)) < SMALLVALUE) {
    return (true);
  }
  else {
    return (false);  // false is zero.
  }
}



// ##################################################################
// ##################################################################



double nemo_constants::pow_fast(double a, double b)
{
  return exp(b * log(a));
}


// ##################################################################
// ##################################################################
