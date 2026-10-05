///
/// \file nemo-constants.h
/// \author Jonathan Mackey
/// \date 2026.09.30
///
/// Class with inline functions for many physics and astronomy
/// constants.  All constants are in cgs units unless stated
/// otherwise.  Adapted from PION constants.h file.
///
///
/// Modifications:
/// - 

#ifndef NEMO_CONSTANTS_H
#define NEMO_CONSTANTS_H

#ifdef SPDLOG_FWD
#include <spdlog/fwd.h>
#endif
#include <spdlog/spdlog.h>
/* prevent clang-format reordering */
#include <fmt/ranges.h>

#ifndef INTEL
#include <cmath>
#else
#include <mathimf.h>
#endif

#include <units/units.hpp>

///
/// Physical constants are from CODATA 2014, converted to CGS units.
/// http://dx.doi.org/10.5281/zenodo.22826
/// Masses of elements from http://www.ciaaw.org/atomic-weights.htm
///
class nemo_constants {
private:
  // Cached CGS values for efficiency
  const double c_cgs;           ///< Speed of light in cm/s
  const double kB_cgs;          ///< Boltzmann constant in erg/K
  const double h_cgs;           ///< Planck constant in erg*s
  const double StefanBoltzmann_cgs; ///< Stefan-Boltzmann constant in erg/(cm^3*K^4)
  const double eV_cgs;          ///< Electron volt in erg
  const double ThomsonXSection_cgs; ///< Thomson cross section in cm^2
  const double uamu_cgs;        ///< Unified atomic mass unit in g
  const double m_e_cgs;         ///< Electron mass in g
  const double m_p_cgs;         ///< Proton mass in g
  const double m_H_cgs;         ///< Hydrogen mass in g
  const double m_He_cgs;        ///< Helium mass in g
  const double m_C_cgs;         ///< Carbon mass in g
  const double m_N_cgs;         ///< Nitrogen mass in g
  const double m_O_cgs;         ///< Oxygen mass in g
  const double G_cgs;           ///< Gravitational constant in cm^3/(g*s^2)
  const double Eth_H_cgs;       ///< Ionization potential of H in erg
  const double NuTh_H_cgs;      ///< Ionization potential of H in Hz
  const double K_per_eV_cgs;    ///< Kelvin per eV
  const double year_cgs;        ///< Year in seconds
  const double asec_per_rad_cgs;///< Arcseconds per radian
  const double sqasec_per_sr_cgs;///< Square arcseconds per steradian
  const double Msun_cgs;        ///< Solar mass in g
  const double Lsun_cgs;        ///< Solar luminosity in erg/s
  const double Rsun_cgs;        ///< Solar radius in cm
  const double AU_cgs;          ///< Astronomical unit in cm
  const double parsec_cgs;      ///< Parsec in cm

public:
  nemo_constants(
  );
  ~nemo_constants();

  //
  // Mathematics
  //
  inline double ln10() { return 2.3025850929940459; }
  inline double pi() { return 3.14159265358979324; }
  inline double halfpi() { return 0.5 * pi(); }
  inline double root4pi() { return 3.5449077018110318; }
  inline double sqrt2() { return 1.4142135623730950; }
  //
  // Physics
  //
  inline double c() { return c_cgs; }
  inline double kB() { return kB_cgs; }
  inline double h() { return h_cgs; }
  inline double StefanBoltzmannConst() { return StefanBoltzmann_cgs; }
  inline double eV() { return eV_cgs; }
  inline double ThomsonXSection() { return ThomsonXSection_cgs; }
  ///
  /// Unified atomic mass unit.
  ///
  inline double uamu() { return uamu_cgs; }
  ///
  /// Electron mass (CODATA 2018) in grams
  ///
  inline double m_e() { return m_e_cgs; }
  ///
  /// proton mass (CODATA 2014) in grams
  ///
  inline double m_p() { return m_p_cgs; }
  ///
  /// Hydrogen mass (IUPAC value is 1.0080+/-0.0002 * uamu).  in grams
  ///
  inline double m_H() { return m_H_cgs; }
  ///
  /// Helium mass (CIAWW online table) in grams
  ///
  inline double m_He() { return m_He_cgs; }

  ///
  /// Carbon mass (CIAWW online table) in grams
  /// 0.5*(12.0096 + 12.0116) uamu
  ///
  inline double m_C() { return m_C_cgs; }

  ///
  /// Nitrogen mass 0.5*(14.00643 + 14.00728) uamu in grams
  /// (CIAWW online table)
  ///
  inline double m_N() { return m_N_cgs; }

  ///
  /// Oxygen mass (CIAWW online table) in grams
  /// 0.5*(15.99903 + 15.99977) uamu
  ///
  inline double m_O() { return m_O_cgs; }

  /// Gravitational constant (CODATA 2014)
  ///
  inline double G() { return G_cgs; }
  ///
  /// Ionization potential of H (erg)
  ///
  inline double Eth_H() { return Eth_H_cgs; }
  ///
  /// Ionization potential of H, in Hertz.
  ///
  inline double NuTh_H() { return NuTh_H_cgs; }
  //
  // Unit conversion
  //
  inline double K_per_eV() { return K_per_eV_cgs; }
  inline double year() { return year_cgs; }
  inline double asec_per_rad() { return asec_per_rad_cgs; }
  inline double sqasec_per_sr() { return sqasec_per_sr_cgs; }
  //
  // Astronomy
  //
  inline double Msun() { return Msun_cgs; }
  inline double Lsun() { return Lsun_cgs; }
  inline double Rsun() { return Rsun_cgs; }
  inline double AU() { return AU_cgs; }
  inline double parsec() { return parsec_cgs; }

  ///
  /// tests if two doubles are equal to an accuracy of 1e-8.
  ///
  bool equalD(const double, const double);
  /// Function to replace pow(a, b) - exp(b*log(a)) is twice as fast
  double pow_fast(
      double,  ///< a
      double   ///< b
  );
};


//-------------------------------------------------------
// Defines for different sized numbers.
//-------------------------------------------------------
#define SMALLVALUE 1.0e-12
#define MACHINEACCURACY 5.e-16
#define TINYVALUE 1.0e-100
#define VERY_TINY_VALUE 1.0e-200
#define VERY_LARGE_VALUE 1.0e100
#define HUGEVALUE 1.0e+200
#define ONE_PLUS_EPS (1.0 + SMALLVALUE)
#define ONE_MINUS_EPS (1.0 - SMALLVALUE)

// Effect of source:
#define RT_EFFECT_MP10_ATLAS_O3V 7   ///< ATLAS models
#define RT_EFFECT_MP10_POWR_WC 8     ///< POWR WC models
#define RT_EFFECT_MP10_POWR_WNE 9    ///< POWR WNE models
#define RT_EFFECT_MP10_BLACKBODY 10  ///< BLACKBODY

///
/// Global function to exit NEMO
///
inline int exit_nemo(const int i)
{
  spdlog::error("NEMO encountered an error and will now exit: code = {}", i);
  spdlog::shutdown();  // should flush the logfile buffer to file.
  exit(i);
}

//--------------------------------------------------------
/// \brief Primitive Variables Enum
/// The code requires that the velocities follow each other in order,
/// and the same for the Magnetic Field.  I'm not sure if it will work
/// without them being in that order, but I'm writing the code assuming
/// they are in order.
///
enum primitive {
  RO = 0,  ///< Density
  PG = 1,  ///< Gas pressure
  VX = 2,  ///< x-velocity
  VY = 3,  ///< y-velocity
  VZ = 4,  ///< z-velocity
  BX = 5,  ///< x-magnetic field
  BY = 6,  ///< y-magnetic field
  BZ = 7,  ///< z-magnetic field
  SI = 8   ///< The GLM scalar field, psi.
};
//--------------------------------------------------------



#endif  // NEMO_CONSTANTS_H
