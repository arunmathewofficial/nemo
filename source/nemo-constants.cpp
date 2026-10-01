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

nemo_constants::nemo_constants() {}


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
