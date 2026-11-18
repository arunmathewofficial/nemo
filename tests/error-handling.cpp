/// error-handling.cpp
///
/// Checks that NEMO reports errors to the caller instead of exiting:
///  - the constructor throws NEMO::nemo_error for an invalid tracer list;
///  - update_single_cell() returns non-zero for a bad input state and for
///    radiation-source data that do not match the sources;
///  - other functions throw NEMO::nemo_error for bad input;
///  - a valid update still works after a failed one.
/// Returns 0 if all checks pass, 1 otherwise.
///

#include <cmath>
#include <cstdio>
#include <limits>
#include <string>
#include <vector>

#include "nemo/nemo.h"
using namespace NEMO;

static int nfail = 0;

static void check(bool ok, const char *what)
{
  printf("  %s: %s\n", ok ? "pass" : "FAIL", what);
  if (!ok) nfail++;
}

static nemo_options make_options()
{
  nemo_options o;
  o.Min_Temperature         = 1.0e1;
  o.Max_Temperature         = 1.0e9;
  o.nemo_reltol             = 1.0e-6;
  o.nemo_abstol             = 1.0e-12;
  o.nemo_minerg             = 1.0e-20;
  o.nv_prim                 = 7;
  o.n_scalar                = 5;
  o.prim_dens_idx           = 0;
  o.prim_pres_idx           = 1;
  o.timestep_limiter        = 1;
  o.update_energy           = true;
  o.cooling                 = true;
  o.collisional_ionization  = true;
  o.photoionization         = false;
  o.charge_exchange         = true;
  o.cosmicray_ionization    = false;
  o.radiative_recombination = true;
  o.onthespot_radiation     = false;
  return o;
}

int main()
{
  spdlog::set_level(spdlog::level::off);  // errors below are expected
  const double g = 5.0 / 3.0;
  std::vector<nemo_star_properties> no_src;

  // Constructor: elements out of order (X_He before X_H).
  {
    std::vector<std::string> tr = {"X_He", "X_H", "H", "He", "He1+"};
    bool thrown                 = false;
    try {
      nemo bad(make_options(), tr, no_src, g);
    }
    catch (const nemo_error &e) {
      thrown = (e.code() != 0);
    }
    check(thrown, "constructor throws nemo_error for unordered elements");
  }

  // Constructor: charge exchange with only one element.
  {
    nemo_options o              = make_options();
    o.nv_prim                   = 4;
    o.n_scalar                  = 2;
    std::vector<std::string> tr = {"X_H", "H"};
    bool thrown                 = false;
    try {
      nemo bad(o, tr, no_src, g);
    }
    catch (const nemo_error &) {
      thrown = true;
    }
    check(thrown, "constructor throws nemo_error for 1-element charge exch.");
  }

  std::vector<std::string> tr = {"X_H", "X_He", "H", "He", "He1+"};
  nemo mp(make_options(), tr, no_src, g);
  std::vector<nemo_cell_rt_data> rt;
  const std::vector<double> good = {1.0e-22, 1.0e-12, 0.72, 0.28,
                                    0.36,    0.05,    0.1};
  std::vector<double> out(good.size(), 0.0);
  double dt_min = 0.0;

  // Non-finite pressure: update returns an error code.
  {
    std::vector<double> bad = good;
    bad[1]                  = std::numeric_limits<double>::quiet_NaN();
    int err = mp.update_single_cell(bad, rt, out, 1.0e10, g, &dt_min);
    check(err != 0, "update_single_cell returns non-zero for NaN pressure");
  }

  // Functions without an error code throw.
  {
    bool thrown = false;
    try {
      mp.get_n_ion("C2+", good);
    }
    catch (const nemo_error &) {
      thrown = true;
    }
    check(thrown, "get_n_ion throws nemo_error for unknown ion");

    thrown = false;
    try {
      std::vector<double> bad = good;
      bad[1]                  = std::numeric_limits<double>::quiet_NaN();
      mp.timescales(bad, rt, g);
    }
    catch (const nemo_error &) {
      thrown = true;
    }
    check(thrown, "timescales throws nemo_error for NaN pressure");
  }

  // With photoionisation: cell data must have one entry per source, with
  // one optical depth per energy bin.
  {
    nemo_options o    = make_options();
    o.photoionization = true;
    std::vector<nemo_star_properties> rs(1);
    rs[0].SED_type   = RT_EFFECT_MP10_BLACKBODY;
    rs[0].Luminosity = 3.0e39;
    rs[0].Teff       = 4.0e4;
    rs[0].Radius     = 18.0;
    rs[0].SED_info   = "";
    rs[0].id         = 0;
    nemo pmp(o, tr, rs, g);
    pmp.get_radiation_src_info(rs);

    std::vector<nemo_cell_rt_data> rt2(2);  // two entries for one source
    for (auto &c : rt2) {
      c.Vshell = 1.0e54;
      c.dS     = 1.0e17;
      c.tau.assign(rs[0].Nbins, 0.0);
      c.dtau.assign(rs[0].Nbins, 0.0);
      c.star_id = 0;
    }
    int err = pmp.update_single_cell(good, rt2, out, 1.0e10, g, &dt_min);
    check(err != 0, "update_single_cell returns non-zero for extra source");

    std::vector<nemo_cell_rt_data> rt1(1, rt2[0]);
    rt1[0].tau.resize(1);  // too few optical depths
    err = pmp.update_single_cell(good, rt1, out, 1.0e10, g, &dt_min);
    check(err != 0, "update_single_cell returns non-zero for short tau");

    bool thrown = false;
    try {
      pmp.timescales(good, rt2, g);
    }
    catch (const nemo_error &) {
      thrown = true;
    }
    check(thrown, "timescales throws nemo_error for extra source");

    thrown = false;
    try {
      std::vector<nemo_star_properties> rs2(2, rs[0]);
      pmp.get_radiation_src_info(rs2);
    }
    catch (const nemo_error &) {
      thrown = true;
    }
    check(thrown, "get_radiation_src_info throws for wrong source count");

    rt1[0].tau.assign(rs[0].Nbins, 0.0);
    err         = pmp.update_single_cell(good, rt1, out, 1.0e10, g, &dt_min);
    bool finite = true;
    for (double v : out)
      finite = finite && std::isfinite(v);
    check(err == 0 && finite, "photoionisation update succeeds afterwards");
  }

  // The solver still works after the failures above.
  {
    int err     = mp.update_single_cell(good, rt, out, 1.0e10, g, &dt_min);
    bool finite = true;
    for (double v : out)
      finite = finite && std::isfinite(v);
    check(err == 0 && finite, "valid update succeeds after failed calls");
  }

  printf("%s: %d failure(s)\n", nfail ? "FAILED" : "PASSED", nfail);
  return nfail ? 1 : 0;
}
