///
/// \file microphysics_base.h
/// \author Jonathan Mackey
/// \date 2011.01.14
///
/// \description
/// Microphysics base class, from which all microphysics classes are
/// derived.  This defines the interface between the microphysics
/// class and the rest of the simulation code.
///
/// - 2026.09.30

#ifndef MICROPHYSICS_BASE_H
#define MICROPHYSICS_BASE_H

#include "nemo-datastructs.h"
#include <map>
#include <string>
#include <vector>
// #define MP_DEBUG

// sets constant temperature for Iliev et al. 2006 test 1 in NEMO, and sets
// monochromatic radiation in 13.6-14.5 eV band
// #define CCCT1


/// virtual base/interface class for in-cell microphysics update.
class microphysics_base {
public:
  ///
  /// Constructor sets up parameters used by derived classes.
  ///
  microphysics_base(
      const int,               ///< Total number of variables in state vector
      const int,               ///< Number of tracer variables in state vector
      const std::string *,     ///< List of tracer variable names.
      struct which_physics *,  ///< which physics to calculate.
      struct rad_sources *     ///< radiation sources.
  );
  virtual ~microphysics_base() {}  ///< non-virtual destructor.

  /// Non-RT microphysics update, so cooling and heating and chemistry.
  ///
  /// This uses various integration methods to update the elements and the
  /// internal energy.
  ///
  virtual int TimeUpdateMP(
      const double *,  ///< Primitive Vector to be updated.
      double *,        ///< Destination Vector for updated values.
      const double,      ///< Time Step to advance by.
      const double,      ///< EOS gamma.
      const int,  ///< Switch for what type of integration to use. (0=adaptive
                  ///< RK5, 1=adaptive Euler,2=onestep o4-RK)
      double *    ///< final temperature.
      ) = 0;

  ///
  /// If doing ray-tracing, the tracer can call this function to
  /// integrate the microphysics variables forward one timestep
  /// given an external radiation flux.
  ///
  virtual int TimeUpdate_RTsinglesrc(
      const double *,  ///< Primitive Vector to be updated.
      double *,        ///< Destination Vector for updated values.
      const double,      ///< Time Step to advance by.
      const double,      ///< EOS gamma.
      const int,  ///< Switch for what type of integration to use. (0=adaptive
                  ///< RK5, 1=adaptive Euler,2=onestep o4-RK)
      const double,  ///< flux in per unit length along ray (F/ds or L/dV)
      const double,  ///< path length ds through cell.
      const double,  ///< Optical depth to entry point of ray into cell.
      double *       ///< return optical depth through cell in this variable.
      ) = 0;

  /// Initialise microphysics ionisation fractions to an equilibrium value.
  virtual int Init_ionfractions(
      double *,    ///< Primitive vector to be updated.
      const double,  ///< eos gamma.
      const double   ///< optional gas temperature to end up at.
                     ///< (negative means use pressure)
      ) = 0;

  ///
  /// This takes a copy of the primitive vector and advances it in time over
  /// the step requested, and at the end copies the updated vector into the
  /// destination vector.  For fully local microphysics but WITH radiative
  /// transfer, where the column densities for diffuse and direct radiation
  /// are included as parameters.  The input list of column densities is
  /// ordered in a list of UV heating sources and a list of ionising sources.
  ///
  virtual int TimeUpdateMP_RTnew(
      const double *,  ///< Primitive Vector to be updated.
      const int,         ///< Number of UV heating sources.
      const std::vector<struct rt_source_data> &,
      ///< list of UV-heating column densities and source properties.
      const int,  ///< number of ionising radiation sources.
      std::vector<struct rt_source_data> &,
      ///< list of ionising src column densities and source properties.
      double *,  ///< Destination Vector for updated values
      ///< (can be same as first Vector.
      const double,  ///< Time Step to advance by.
      const double,  ///< EOS gamma.
      const int,     ///< Switch for what type of integration to use.
      ///< (0=adaptive RK5, 1=adaptive Euler,2=onestep o4-RK)
      double *  ///< final temperature (not strictly needed).
  )
  {
    return 1;
  }

  /// Given a string, return the index of the tracer that it refers to.
  ///
  /// \retval int>0 success: int=tracer value
  /// \retval int<0 failure: string doesn't correspond to any tracer.
  ///
  virtual int Tr(const std::string  ///< name of tracer we are looking for.
  );

  /// Sets the pressure in p-vec to be such that temperature is
  /// what you want it to be.
  virtual int Set_Temp(
      double *,    ///< primitive vector.
      const double,  ///< temperature.
      const double   ///< eos gamma.
      ) = 0;

  ///
  /// Returns the gas temperature (not very optimized though)
  /// Assumes primitive vector is in cgs units.
  ///
  virtual double Temperature(
      const double *,  ///< primitive vector
      const double       ///< eos gamma
      ) = 0;

  ///
  /// This returns the minimum timescale of the times flagged in the
  /// arguments.  Time is returned in seconds.
  ///
  virtual double timescales(
      const double *,  ///< Current cell.
      const double,      ///< EOS gamma.
      const bool,        ///< set to true if including cooling time.
      const bool,        ///< set to true if including recombination time.
      const bool         ///< set to true if including photo-ionsation time.
      ) = 0;

  ///
  /// This returns the minimum timescale of all microphysical processes,
  /// including reaction times for each species and the total heating/cooling
  /// time for the gas. It requires the radiation field as an input, so it has
  /// substantially greater capability than the other timescales function.
  ///
  virtual double timescales_RT(
      const double *,  ///< Current cell.
      const int,         ///< Number of UV heating sources.
      const std::vector<struct rt_source_data> &,
      ///< list of UV-heating column densities and source properties.
      const int,  ///< number of ionising radiation sources.
      const std::vector<struct rt_source_data> &,
      ///< list of ionising src column densities and source properties.
      const double  ///< EOS gamma.
      ) = 0;

  ///
  /// Set the properties of a multifrequency ionising radiation source.
  ///
  virtual int set_multifreq_source_properties(
      const struct rad_src_info *,  ///< source data
      double *  ///< source strength in different energy bins.
  )
  {
    return -999;
  }

  ///
  /// Get the total cooling rate.  This is for postprocessing the
  /// simulation data only -- IT IS NOT OPTIMISED FOR SPEED.  Note
  /// that this is dE/dt, in erg/cm^3/s, so it is the difference
  /// between heating and cooling.
  ///
  virtual double total_cooling_rate(
      const double *,  ///< Current cell values.
      const int,         ///< Number of UV heating sources.
      const std::vector<struct rt_source_data> &,
      ///< list of UV-heating column densities and source properties.
      const int,  ///< number of ionising radiation sources.
      const std::vector<struct rt_source_data> &,
      ///< list of ionising src column densities and source properties.
      const double  ///< EOS gamma.
  )
  {
    return -999.999;
  }

  ///
  /// Get the total recombination rate for an ion, given the input
  /// state vector.
  ///
  virtual double get_recombination_rate(
      const int,         ///< ion index in tracer array (optional).
      const double *,  ///< input state vector (primitive).
      const double       ///< EOS gamma (optional)
      ) = 0;

  ///
  /// Updates the corrector vector
  /// according to sCMA (simple Consistent Multi-fluid Advection,
  /// Plewa + Muller, 1999).
  /// Used for modifying tracer fluxes, normalising state vec.
  ///
  virtual void sCMA(
      double *,       ///< output corrector vector
      const double *  ///< input primitive vector from grid cell
  );

  ///
  /// Get optical depth for a range of frequencies based on the
  /// local abundances of elements and species in the input primitive
  /// vector.
  ///
  virtual void get_dtau(
      const rad_src_info *,  ///< pointer to radiation source struct
      const double,      ///< ds, thickness of the cell
      const double *,    ///< input primitive vector from cell
      double *           ///< output dtau vector
  ) 
  {
    return;
  }

  ///
  /// Get the ionisation cross section for an atom/ion at its
  /// threshold frequency.
  ///
  virtual double get_th_xsection(const int  ///< integer identifier for the ion.
  )
  {
    return -1.0e99;
  }

  ///
  /// Return number density of a given element.
  ///
  virtual double get_n_el(
      const double *,  ///< primitive state vector.
      const int          ///< integer identifier for the element.
  )
  {
    return -1.0e99;
  }

  ///
  /// Get electron number density (cm^{-3})
  ///
  virtual double get_n_elec(const double *  ///< primitive state vector.
  )
  {
    return -1.0e99;
  }

  ///
  /// Get ion number density (cm^{-3})
  ///
  virtual double get_n_ion(
      std::string,      ///< ion name
      const double *  ///< primitive state vector.
  )
  {
    return -1.0e99;
  }

  ///
  /// Get H+ number density (cm^{-3})
  ///
  virtual double get_n_Hplus(const double *  ///< primitive state vector.
  )
  {
    return -1.0e99;
  }

  ///
  /// Get neutral H number density (cm^{-3})
  ///
  virtual double get_n_Hneutral(const double *  ///< primitive state vector.
  )
  {
    return -1.0e99;
  }

  ///
  /// Return the H mass fraction
  ///
  virtual inline double get_X_H() { return 0.7154; }

protected:
  /// Struct with flags for which extra physics we are (not) doing.
  struct which_physics *EP;
  /// Struct with list of radiation sources
  struct rad_sources *RS;
  /// map tracer names to state-vec index
  std::map<std::string, int> tr_map;
  /// map element names to state-vec index
  std::map<std::string, int> el_map;
  std::vector<int> el_index;  ///< indices of elements
  std::vector<int> tr_index;  ///< indices of tracers
  int n_el;                   ///< number of element tracers

  int nv_prim;  ///< number of vars in state vec
  int ntracer;  ///< number of tracer variables
};

#endif  // MICROPHYSICS_BASE_H
