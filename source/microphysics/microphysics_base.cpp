/// \file microphysics_base.cpp
/// \author Jonathan Mackey
/// \date 2015.03.26

#include "microphysics_base.h"

#ifdef SPDLOG_FWD
#include <spdlog/fwd.h>
#endif
#include <spdlog/spdlog.h>
/* prevent clang-format reordering */
#include <fmt/ranges.h>

using namespace std;



// ##################################################################
// ##################################################################

microphysics_base::microphysics_base(
    const int nv,           ///< Total number of variables in state vector
    const int ntr,          ///< Number of tracer variables in state vector
    const std::string *tr,  ///< List of tracer variable names.
    struct which_physics *ephys,  ///< which physics to calculate.
    struct rad_sources *rsrcs     ///< radiation sources.
    ) : nv_prim(nv), ntracer(ntr)
{
  EP = ephys;
  RS = rsrcs;

#ifndef NDEBUG
  spdlog::debug(
      "microphysics_base:: setting up EP and RS: {}\t{}", fmt::ptr(EP),
      fmt::ptr(RS));
#endif

  // set up a map of tracer names to indices in state vec.
  int offset = nv - ntr;
  n_el       = 0;
  for (int i = 0; i < ntr; i++) {
    // spdlog::warn("tri {} {}",i,tr[i]);
    tr_map[tr[i]] = i + offset;
    tr_index.push_back(i + offset);
    // second map just for elements.
    if (tr[i].substr(0, 2) == "X_") {
      el_map[tr[i]] = i + offset;
      el_index.push_back(i + offset);
      n_el++;
    }
  }
  // spdlog::warn("el_index {}",el_index);
  // spdlog::warn("tr_index {}",tr_index);

#ifdef DEBUG_MP
  // for (auto elem : el_map) { // c++11 only (not working with intel compiler
  // 2018r4)
  for (std::map<std::string, int>::iterator it = el_map.begin();
       it != el_map.end(); ++it) {
    spdlog::debug("Elements Map: {} {}", it->first, it->second);
  }
  // for (auto elem : tr_map) { // c++11 only (not working with intel compiler
  // 2018r4)
  for (std::map<std::string, int>::iterator it = tr_map.begin();
       it != tr_map.end(); ++it) {
    spdlog::debug("Tracers Map: {} {}", it->first, it->second);
  }
#endif
}

// ##################################################################
// ##################################################################

// #define DEBUG_SCMA

void microphysics_base::sCMA(
    pion_flt *corrected,  ///< output corrected vector
    const pion_flt *p_in  ///< input primitive vector (nv_prim)
)
{
  // "Consistent multi-species advection" correction step from
  // Plewa and Muller, Astron. Astrophys. 342, 179 (1999).
#ifdef DEBUG_SCMA
  int print_flag = 0;
  double val     = 0.0;
#endif  // DEBUG_SCMA

  for (int i = 0; i < nv_prim; i++)
    corrected[i] = p_in[i];
  // enforce all tracers to be within [0,1]
  for (int v = 0; v < ntracer; v++)
    corrected[tr_index[v]] = min(1.0, max(0.0, p_in[tr_index[v]]));

  // sum over all elements and add up mass fraction, correcting
  // state vector to be within allowed limits for mass fraction.
  pion_flt corr = 0.0;
  for (int v = 0; v < n_el; v++)
    corr += corrected[el_index[v]];
  corr = 1.0 / corr;  // correction factor for elements

#ifdef DEBUG_SCMA
  if (print_flag)
    spdlog::warn(
        "Corrected 1: {}",
        std::vector<pion_flt>(corrected, corrected + nv_prim));
#endif  // DEBUG_SCMA

  // renormalise element values to sum to 1.
  for (int v = 0; v < n_el; v++)
    corrected[el_index[v]] *= corr;

#ifdef DEBUG_SCMA
  if (print_flag)
    spdlog::warn(
        "Corrector 2: {}",
        std::vector<pion_flt>(corrected, corrected + nv_prim));
#endif  // DEBUG_SCMA

#ifdef DEBUG_SCMA
  if (print_flag) {
    spdlog::warn("------------  sCMA  ---------------");
    spdlog::warn("p_in     : {}", std::vector<pion_flt>(p_in, p_in + nv_prim));
    spdlog::warn(
        "Corrected: {}", std::vector<pion_flt>(corrected, corrected + nv_prim));
  }
#endif  // DEBUG_SCMA
}

// ##################################################################
// ##################################################################

int microphysics_base::Tr(const string s)
{
  if (tr_map.find(s) == tr_map.end())
    return -1;
  else
    return tr_map[s];
}

// ##################################################################
// ##################################################################
