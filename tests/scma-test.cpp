/// scma-test.cpp
///
/// Unit test of nemo::sCMA(), the consistent multi-species advection
/// correction (Plewa & Muller, 1999, A&A, 342, 179).  Checks that:
///  - element mass fractions are renormalised to sum to 1;
///  - the fraction of each element in each species (X_s/X_e) is unchanged
///    by the renormalisation (so no spurious ions are created);
///  - the species of an element are limited to sum to <= X_e;
///  - negative mass fractions are clipped.
/// Returns 0 if all checks pass, 1 otherwise.
///

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "nemo/nemo.h"
using namespace NEMO;

// state vector: rho, p, vx, vy, vz, X_H, X_He, H, He, He1+
static const int iXH = 5, iXHe = 6, iH = 7, iHe = 8, iHe1 = 9, NV = 10;

static int nfail = 0;

static void check(bool ok, const char *what, double a, double b)
{
  if (!ok) {
    printf("  FAIL: %s (%.15e vs %.15e)\n", what, a, b);
    nfail++;
  }
}

static bool close(double a, double b, double rtol = 1.0e-12)
{
  return std::fabs(a - b) <= rtol * std::fabs(b) + 1.0e-30;
}

// Input with element sum `sum`, H+ fraction yHp, He+ fraction yHep and
// He2+ fraction yHepp: check that output elements sum to 1 and the
// ionisation fractions of each element are preserved.
static void test_renorm(
    nemo &n, double sum, double yHp, double yHep, double yHepp)
{
  printf(
      "element sum %.3f, y(H+)=%.2f y(He+)=%.2f y(He2+)=%.2f\n", sum, yHp, yHep,
      yHepp);
  double XH = 0.7 * sum, XHe = 0.3 * sum;
  std::vector<double> p(NV, 0.0), c(NV, 0.0);
  p[0]    = 1.0e-22;
  p[1]    = 1.0e-12;
  p[iXH]  = XH;
  p[iXHe] = XHe;
  p[iH]   = XH * (1.0 - yHp);
  p[iHe]  = XHe * (1.0 - yHep - yHepp);
  p[iHe1] = XHe * yHep;
  n.sCMA(c, p);

  check(
      close(c[iXH] + c[iXHe], 1.0), "elements sum to 1", c[iXH] + c[iXHe], 1.0);
  check(
      close(c[iXH] / c[iXHe], 0.7 / 0.3), "element ratio preserved",
      c[iXH] / c[iXHe], 0.7 / 0.3);
  // fractions of each element in each species
  check(
      close(c[iH] / c[iXH], 1.0 - yHp), "y(H0) preserved", c[iH] / c[iXH],
      1.0 - yHp);
  if (1.0 - yHep - yHepp > 0.0)
    check(
        close(c[iHe] / c[iXHe], 1.0 - yHep - yHepp), "y(He0) preserved",
        c[iHe] / c[iXHe], 1.0 - yHep - yHepp);
  if (yHep > 0.0)
    check(
        close(c[iHe1] / c[iXHe], yHep), "y(He+) preserved", c[iHe1] / c[iXHe],
        yHep);
  // top ions (untracked) preserved
  check(
      close((c[iXH] - c[iH]) / c[iXH] + 1.0, yHp + 1.0), "y(H+) preserved",
      (c[iXH] - c[iH]) / c[iXH], yHp);
  check(
      close((c[iXHe] - c[iHe] - c[iHe1]) / c[iXHe] + 1.0, yHepp + 1.0),
      "y(He2+) preserved", (c[iXHe] - c[iHe] - c[iHe1]) / c[iXHe], yHepp);
}

int main(void)
{
  std::vector<std::string> tr({"X_H", "X_He", "H", "He", "He1+"});
  std::vector<struct nemo_star_properties> rs;

  struct nemo_options o;
  o.Min_Temperature         = 1.0e1;
  o.Max_Temperature         = 1.0e9;
  o.nemo_reltol             = 1.0e-6;
  o.nemo_abstol             = 1.0e-12;
  o.nemo_minerg             = 1.0e-20;
  o.nv_prim                 = NV;
  o.n_scalar                = 5;
  o.prim_dens_idx           = 0;
  o.prim_pres_idx           = 1;
  o.timestep_limiter        = 1;
  o.update_energy           = 0;
  o.cooling                 = 0;
  o.collisional_ionization  = 1;
  o.photoionization         = 0;
  o.charge_exchange         = 0;
  o.cosmicray_ionization    = 0;
  o.radiative_recombination = 1;
  o.onthespot_radiation     = 0;
  nemo n(o, tr, rs, 5.0 / 3.0);

  // renormalisation must preserve ionisation fractions
  for (double sum : {1.0, 0.95, 1.05, 0.5, 1.5}) {
    test_renorm(n, sum, 0.0, 0.0, 0.0);
    test_renorm(n, sum, 0.3, 0.2, 0.1);
    test_renorm(n, sum, 0.9, 0.5, 0.4);
  }

  // species summing to more than their element are scaled down to X_e
  {
    printf("species sum > element\n");
    std::vector<double> p = {1e-22, 1e-12, 0, 0, 0, 0.7, 0.3, 0.8, 0.2, 0.2};
    std::vector<double> c(NV, 0.0);
    n.sCMA(c, p);
    check(close(c[iH], c[iXH]), "H <= X_H", c[iH], c[iXH]);
    check(
        close(c[iHe] + c[iHe1], c[iXHe]), "He + He1+ == X_He", c[iHe] + c[iHe1],
        c[iXHe]);
    check(
        close(c[iHe] / c[iHe1], 1.0), "He/He1+ ratio preserved",
        c[iHe] / c[iHe1], 1.0);
  }

  // negative values are clipped
  {
    printf("negative mass fractions\n");
    std::vector<double> p = {1e-22, 1e-12, 0,   0,     0,
                             1.0,   -0.01, 0.5, -0.01, 0.0};
    std::vector<double> c(NV, 0.0);
    n.sCMA(c, p);
    for (int v = iXH; v < NV; v++)
      check(c[v] > 0.0 && c[v] <= 1.0, "mass fraction in (0,1]", c[v], 0.0);
    check(
        close(c[iXH] + c[iXHe], 1.0), "elements sum to 1", c[iXH] + c[iXHe],
        1.0);
  }

  if (nfail)
    printf("sCMA test: %d checks FAILED\n", nfail);
  else
    printf("sCMA test: all checks passed\n");
  return nfail ? 1 : 0;
}
