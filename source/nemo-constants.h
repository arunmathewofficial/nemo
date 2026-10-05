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

///
/// Physical constants are from CODATA 2014, converted to CGS.
/// http://dx.doi.org/10.5281/zenodo.22826
/// Masses of elements from http://www.ciaaw.org/atomic-weights.htm
///
class nemo_constants {
public:
  nemo_constants();
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
  inline double c() { return 2.99792458e+10; }
  inline double kB() { return 1.38064852e-16; }
  inline double h() { return 6.62606957e-27; }
  inline double StefanBoltzmannConst() { return 5.670367e-5; }
  inline double eV() { return 1.6021766208e-12; }
  inline double ThomsonXSection() { return 6.6524587321e-25; }
  ///
  /// Unified atomic mass unit.
  ///
  inline double uamu() { return 1.660539040e-24; }
  ///
  /// Electron mass (CODATA 2018) in grams
  ///
  inline double m_e() { return 9.1093837015e-28; }
  ///
  /// proton mass (CODATA 2014) in grams
  ///
  inline double m_p() { return 1.672621898e-24; }
  ///
  /// Hydrogen mass (IUPAC value is 1.0080+/-0.0002 * uamu).  in grams
  ///
  inline double m_H() { return 1.6738e-24; }
  ///
  /// Helium mass (CIAWW online table) in grams
  ///
  inline double m_He() { return 6.6464768e-24; }

  ///
  /// Carbon mass (CIAWW online table) in grams
  /// 0.5*(12.0096 + 12.0116) uamu
  ///
  inline double m_C() { return 1.994374e-23; }

  ///
  /// Nitrogen mass 0.5*(14.00643 + 14.00728) uamu in grams
  /// (CIAWW online table)
  ///
  inline double m_N() { return 2.325892e-23; }

  ///
  /// Oxygen mass (CIAWW online table) in grams
  /// 0.5*(15.99903 + 15.99977) uamu
  ///
  inline double m_O() { return 2.6567628e-23; }

  /// Gravitational constant (CODATA 2014)
  ///
  inline double G() { return 6.67408e-8; }
  ///
  /// Ionization potential of H (erg)
  ///
  inline double Eth_H() { return 2.178710264e-11; }
  ///
  /// Ionization potential of H, in Hertz.
  ///
  inline double NuTh_H() { return 3.28808819e+15; }
  //
  // Unit conversion
  //
  inline double K_per_eV() { return 1.16045221e4; }
  inline double year() { return 3.1558150e7; }
  inline double asec_per_rad() { return 206264.806; }
  inline double sqasec_per_sr() { return 4.2545250225e10; }
  //
  // Astronomy
  //
  inline double Msun() { return 1.9891e33; }
  inline double Lsun() { return 3.839e33; }
  inline double Rsun() { return 6.955e10; }
  inline double AU() { return 1.49597870700e+13; }
  inline double parsec() { return 3.0856776e18; }

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
