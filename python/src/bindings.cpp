///
/// \file bindings.cpp
///
/// Python bindings for NEMO (module nemo_plasma._core), using pybind11.
///
/// Each method of the Python class Nemo wraps one public function of
/// NEMO::nemo.  The methods take a state vector of shape (nv,) for one cell,
/// or an array of shape (..., nv) for many cells, which are processed in C++
/// without holding the GIL.  Radiation-field data (for photoionisation) are
/// NumPy arrays with one row per source.
///
/// A Nemo object is not shared between threads by NEMO, so each Python
/// object holds a mutex: calls from several Python threads on the same
/// object are serialised, and different objects run in parallel.
///

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <mutex>
#include <string>
#include <vector>

#include "nemo/nemo.h"

namespace py = pybind11;

namespace {

/// Contiguous double array, converted from any array-like input.
using darray = py::array_t<double, py::array::c_style | py::array::forcecast>;

/// Python exception type for NEMO::nemo_error (set at module import, and
/// never released, so that nothing is destroyed after the interpreter).
py::handle nemo_error_type;

[[noreturn]] void raise_value_error(const std::string &msg)
{
  throw py::value_error(msg);
}

/// Options with the same defaults as PION.  nv_prim and n_scalar default
/// to -1, meaning "density, pressure and then the tracers".
NEMO::nemo_options default_options()
{
  NEMO::nemo_options o;
  o.nemo_reltol             = 1.0e-6;
  o.nemo_abstol             = 1.0e-12;
  o.nemo_minerg             = 1.0e-20;
  o.Min_Temperature         = 10.0;
  o.Max_Temperature         = 1.0e9;
  o.nv_prim                 = -1;
  o.n_scalar                = -1;
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

/// Shape of a state-vector array (..., nv): number of cells and the shape
/// of the leading dimensions.
struct cells {
  py::ssize_t n;                  ///< number of cells
  std::vector<py::ssize_t> lead;  ///< leading dimensions (empty: 1 cell)
};

cells check_prim(const darray &prim, const int nv, const char *name)
{
  if (prim.ndim() < 1 || prim.shape(prim.ndim() - 1) != nv)
    raise_value_error(
        std::string(name)
        + ": last dimension must have length nv_prim = " + std::to_string(nv));
  cells c;
  c.n = 1;
  for (py::ssize_t d = 0; d < prim.ndim() - 1; d++) {
    c.lead.push_back(prim.shape(d));
    c.n *= prim.shape(d);
  }
  return c;
}

/// Output array with shape lead + extra.
darray make_out(const cells &c, std::vector<py::ssize_t> extra)
{
  std::vector<py::ssize_t> shape = c.lead;
  shape.insert(shape.end(), extra.begin(), extra.end());
  return darray(shape);
}

/// Broadcast a scalar or per-cell array (shape lead) to n values.
std::vector<double> per_cell(const darray &a, const cells &c, const char *name)
{
  std::vector<double> v(c.n);
  if (a.size() == 1) {
    std::fill(v.begin(), v.end(), *a.data());
    return v;
  }
  bool ok = (a.ndim() == static_cast<py::ssize_t>(c.lead.size()));
  for (size_t d = 0; ok && d < c.lead.size(); d++)
    ok = (a.shape(d) == c.lead[d]);
  if (!ok)
    raise_value_error(
        std::string(name)
        + ": must be a scalar or have the shape of the state array without "
          "its last dimension");
  std::copy(a.data(), a.data() + c.n, v.begin());
  return v;
}

/// Radiation-field input for many cells: tau and dtau (..., n_src, n_bins),
/// ds and vshell (..., n_src).
struct radiation {
  bool present = false;
  int n_src    = 0;
  int n_bins   = 0;
  std::vector<double> tau, dtau, ds, vshell;

  /// Fill rt (one entry per source) for cell i.
  void fill(py::ssize_t i, std::vector<NEMO::nemo_cell_rt_data> &rt) const
  {
    if (!present) return;
    for (int s = 0; s < n_src; s++) {
      auto &r        = rt[s];
      const size_t k = (static_cast<size_t>(i) * n_src + s);
      r.star_id      = s;
      r.dS           = ds[k];
      r.Vshell       = vshell[k];
      std::copy(
          tau.begin() + k * n_bins, tau.begin() + (k + 1) * n_bins,
          r.tau.begin());
      std::copy(
          dtau.begin() + k * n_bins, dtau.begin() + (k + 1) * n_bins,
          r.dtau.begin());
    }
  }
};

/// Copy array `a` with expected shape lead + extra into a vector.
std::vector<double> copy_shaped(
    const darray &a,
    const cells &c,
    const std::vector<py::ssize_t> &extra,
    const char *name)
{
  std::vector<py::ssize_t> want = c.lead;
  want.insert(want.end(), extra.begin(), extra.end());
  bool ok = (a.ndim() == static_cast<py::ssize_t>(want.size()));
  for (size_t d = 0; ok && d < want.size(); d++)
    ok = (a.shape(d) == want[d]);
  if (!ok) {
    std::string s = "(";
    for (size_t d = 0; d < want.size(); d++)
      s += std::to_string(want[d]) + (d + 1 < want.size() ? ", " : "");
    raise_value_error(std::string(name) + ": expected shape " + s + ")");
  }
  return std::vector<double>(a.data(), a.data() + a.size());
}

/// Find a tracer's position in the state vector.
int tracer_index(
    const std::vector<std::string> &tracers,
    const NEMO::nemo_options &o,
    const std::string &name)
{
  auto it = std::find(tracers.begin(), tracers.end(), name);
  if (it == tracers.end())
    raise_value_error("index: no tracer named '" + name + "'");
  return o.nv_prim - o.n_scalar + static_cast<int>(it - tracers.begin());
}

// ##########################################################################

class py_nemo {
public:
  py_nemo(
      NEMO::nemo_options o,
      std::vector<std::string> tracers,
      std::vector<NEMO::nemo_star_properties> sources,
      const double gamma) :
      opt(o), tracer_names(std::move(tracers)), srcs(std::move(sources)),
      eos_gamma(gamma)
  {
    const int ntr = static_cast<int>(tracer_names.size());
    if (opt.n_scalar < 0) opt.n_scalar = ntr;
    if (opt.nv_prim < 0) opt.nv_prim = 2 + opt.n_scalar;
    // NEMO does not check these, and would read out of bounds.
    if (opt.n_scalar != ntr)
      raise_value_error(
          "options.n_scalar (" + std::to_string(opt.n_scalar)
          + ") must equal the number of tracer names (" + std::to_string(ntr)
          + ")");
    const int nfirst = opt.nv_prim - opt.n_scalar;
    if (nfirst < 2)
      raise_value_error(
          "options.nv_prim must leave room for density and pressure before "
          "the tracers");
    if (opt.prim_dens_idx < 0 || opt.prim_dens_idx >= nfirst
        || opt.prim_pres_idx < 0 || opt.prim_pres_idx >= nfirst
        || opt.prim_dens_idx == opt.prim_pres_idx)
      raise_value_error(
          "options.prim_dens_idx and prim_pres_idx must be different and "
          "before the tracers");
    if (!(gamma > 1.0)) raise_value_error("gamma must be > 1");
    if (opt.photoionization && srcs.empty())
      raise_value_error("photoionization needs at least one source");
    for (size_t s = 0; s < srcs.size(); s++) {
      srcs[s].id         = static_cast<int>(s);
      srcs[s].photo_emax = 0.0;
      srcs[s].Nbins      = 0;
      srcs[s].SED_info   = "";
    }

    mp.reset(new NEMO::nemo(opt, tracer_names, srcs, eos_gamma));
    if (opt.photoionization) {
      mp->get_radiation_src_info(srcs);
      n_bins = srcs.empty() ? 0 : srcs[0].Nbins;
    }
    else {
      srcs.clear();  // NEMO ignores sources without photoionisation
    }
  }

  // --------------------------------------------------------------------
  // properties

  const NEMO::nemo_options &options() const { return opt; }
  const std::vector<std::string> &tracers() const { return tracer_names; }
  const std::vector<NEMO::nemo_star_properties> &sources() const
  {
    return srcs;
  }
  double gamma() const { return eos_gamma; }
  int nv_prim() const { return opt.nv_prim; }
  int nbins() const { return n_bins; }
  int index(const std::string &name) const
  {
    return tracer_index(tracer_names, opt, name);
  }

  // --------------------------------------------------------------------
  // update_single_cell

  darray update(
      const darray &prim,
      const darray &dt_in,
      const py::object &tau,
      const py::object &dtau,
      const py::object &ds,
      const py::object &vshell,
      const std::string &on_error)
  {
    if (on_error != "raise" && on_error != "nan")
      raise_value_error("on_error must be 'raise' or 'nan'");
    const cells c          = check_prim(prim, opt.nv_prim, "update");
    std::vector<double> dt = per_cell(dt_in, c, "dt");
    radiation rad          = get_radiation(c, tau, dtau, ds, vshell, true);
    darray out             = make_out(c, {opt.nv_prim});

    const double *pin  = prim.data();
    double *pout       = out.mutable_data();
    const int nv       = opt.nv_prim;
    py::ssize_t n_fail = 0, first_fail = -1;
    int first_code = 0;
    {
      py::gil_scoped_release release;
      std::lock_guard<std::mutex> lock(mtx);
      std::vector<double> vin(nv), vout(nv);
      std::vector<NEMO::nemo_cell_rt_data> rt = make_rt(rad);
      double unused                           = 0.0;
      for (py::ssize_t i = 0; i < c.n; i++) {
        std::copy(pin + i * nv, pin + (i + 1) * nv, vin.begin());
        rad.fill(i, rt);
        int err =
            mp->update_single_cell(vin, rt, vout, dt[i], eos_gamma, &unused);
        if (err) {
          if (n_fail == 0) {
            first_fail = i;
            first_code = err;
          }
          n_fail++;
          std::fill(
              vout.begin(), vout.end(),
              std::numeric_limits<double>::quiet_NaN());
        }
        std::copy(vout.begin(), vout.end(), pout + i * nv);
      }
    }
    if (n_fail > 0 && on_error == "raise") {
      std::string msg =
          "update failed for " + std::to_string(n_fail) + " cell(s)";
      if (!c.lead.empty())
        msg += "; first failure at flat index " + std::to_string(first_fail);
      msg += " (NEMO error code " + std::to_string(first_code)
             + ": 1 = bad input, 2 = ODE integration failed, 3 = bad "
               "output)";
      throw NEMO::nemo_error(msg, first_code);
    }
    return out;
  }

  // --------------------------------------------------------------------
  // Temperature, get_n_elec, get_n_ion, timescales: one value per cell

  darray temperature(const darray &prim)
  {
    return per_cell_value(prim, "temperature", [this](const auto &v) {
      return mp->Temperature(v, eos_gamma);
    });
  }

  darray electron_density(const darray &prim)
  {
    return per_cell_value(prim, "electron_density", [this](const auto &v) {
      return mp->get_n_elec(v);
    });
  }

  darray ion_density(const std::string &ion, const darray &prim)
  {
    return per_cell_value(prim, "ion_density", [this, &ion](const auto &v) {
      return mp->get_n_ion(ion, v);
    });
  }

  darray timescale(
      const darray &prim,
      const py::object &tau,
      const py::object &dtau,
      const py::object &ds,
      const py::object &vshell)
  {
    const cells c     = check_prim(prim, opt.nv_prim, "timescale");
    radiation rad     = get_radiation(c, tau, dtau, ds, vshell, true);
    darray out        = make_out(c, {});
    const double *pin = prim.data();
    double *pout      = out.mutable_data();
    const int nv      = opt.nv_prim;
    run_released([&]() {
      std::vector<double> v(nv);
      std::vector<NEMO::nemo_cell_rt_data> rt = make_rt(rad);
      for (py::ssize_t i = 0; i < c.n; i++) {
        std::copy(pin + i * nv, pin + (i + 1) * nv, v.begin());
        rad.fill(i, rt);
        pout[i] = mp->timescales(v, rt, eos_gamma);
      }
    });
    return out;
  }

  // --------------------------------------------------------------------
  // Set_Temp and sCMA: one state vector per cell

  darray set_temperature(const darray &prim, const darray &T_in)
  {
    const cells c         = check_prim(prim, opt.nv_prim, "set_temperature");
    std::vector<double> T = per_cell(T_in, c, "T");
    darray out            = make_out(c, {opt.nv_prim});
    const double *pin     = prim.data();
    double *pout          = out.mutable_data();
    const int nv          = opt.nv_prim;
    int err_code          = 0;
    run_released([&]() {
      std::vector<double> v(nv);
      for (py::ssize_t i = 0; i < c.n && !err_code; i++) {
        std::copy(pin + i * nv, pin + (i + 1) * nv, v.begin());
        err_code = mp->Set_Temp(v, T[i], eos_gamma);
        std::copy(v.begin(), v.end(), pout + i * nv);
      }
    });
    if (err_code) throw NEMO::nemo_error("set_temperature failed", err_code);
    return out;
  }

  darray scma(const darray &prim)
  {
    const cells c     = check_prim(prim, opt.nv_prim, "scma");
    darray out        = make_out(c, {opt.nv_prim});
    const double *pin = prim.data();
    double *pout      = out.mutable_data();
    const int nv      = opt.nv_prim;
    run_released([&]() {
      std::vector<double> v(nv), w(nv);
      for (py::ssize_t i = 0; i < c.n; i++) {
        std::copy(pin + i * nv, pin + (i + 1) * nv, v.begin());
        mp->sCMA(w, v);
        std::copy(w.begin(), w.end(), pout + i * nv);
      }
    });
    return out;
  }

  // --------------------------------------------------------------------
  // get_dtau: optical depth through each cell in each energy bin

  darray dtau(const darray &prim, const darray &ds_in)
  {
    if (!opt.photoionization)
      raise_value_error("dtau needs photoionization to be on");
    const cells c          = check_prim(prim, opt.nv_prim, "dtau");
    std::vector<double> ds = per_cell(ds_in, c, "ds");
    darray out             = make_out(c, {n_bins});
    const double *pin      = prim.data();
    double *pout           = out.mutable_data();
    const int nv = opt.nv_prim, nb = n_bins;
    run_released([&]() {
      std::vector<double> v(nv);
      std::vector<NEMO::nemo_cell_rt_data> rt(1);
      rt[0].star_id = 0;
      rt[0].Vshell  = 0.0;
      rt[0].tau.assign(nb, 0.0);
      rt[0].dtau.assign(nb, 0.0);
      for (py::ssize_t i = 0; i < c.n; i++) {
        std::copy(pin + i * nv, pin + (i + 1) * nv, v.begin());
        rt[0].dS = ds[i];
        mp->get_dtau(v, rt);
        std::copy(rt[0].dtau.begin(), rt[0].dtau.end(), pout + i * nb);
      }
    });
    return out;
  }

private:
  NEMO::nemo_options opt;
  std::vector<std::string> tracer_names;
  std::vector<NEMO::nemo_star_properties> srcs;
  double eos_gamma;
  int n_bins = 0;
  std::unique_ptr<NEMO::nemo> mp;
  std::mutex mtx;

  /// Run f without the GIL and with this object locked.  A nemo_error
  /// thrown by f is passed on after the GIL is taken back.
  template<typename F>
  void run_released(F f)
  {
    py::gil_scoped_release release;
    std::lock_guard<std::mutex> lock(mtx);
    f();
  }

  template<typename F>
  darray per_cell_value(const darray &prim, const char *name, F f)
  {
    const cells c     = check_prim(prim, opt.nv_prim, name);
    darray out        = make_out(c, {});
    const double *pin = prim.data();
    double *pout      = out.mutable_data();
    const int nv      = opt.nv_prim;
    run_released([&]() {
      std::vector<double> v(nv);
      for (py::ssize_t i = 0; i < c.n; i++) {
        std::copy(pin + i * nv, pin + (i + 1) * nv, v.begin());
        pout[i] = f(v);
      }
    });
    return out;
  }

  std::vector<NEMO::nemo_cell_rt_data> make_rt(const radiation &rad) const
  {
    std::vector<NEMO::nemo_cell_rt_data> rt(rad.n_src);
    for (auto &r : rt) {
      r.tau.assign(rad.n_bins, 0.0);
      r.dtau.assign(rad.n_bins, 0.0);
    }
    return rt;
  }

  /// Read and check the radiation-field arrays.  With photoionisation all
  /// four are required; without it they must not be given.
  radiation get_radiation(
      const cells &c,
      const py::object &tau,
      const py::object &dtau,
      const py::object &ds,
      const py::object &vshell,
      const bool)
  {
    radiation r;
    const bool any =
        !tau.is_none() || !dtau.is_none() || !ds.is_none() || !vshell.is_none();
    if (!opt.photoionization) {
      if (any)
        raise_value_error(
            "tau, dtau, ds and vshell are only used with photoionization");
      return r;
    }
    if (tau.is_none() || dtau.is_none() || ds.is_none() || vshell.is_none())
      raise_value_error(
          "photoionization is on: tau, dtau, ds and vshell must be given");
    r.present                         = true;
    r.n_src                           = static_cast<int>(srcs.size());
    r.n_bins                          = n_bins;
    const std::vector<py::ssize_t> sb = {r.n_src, r.n_bins}, s1 = {r.n_src};
    r.tau    = copy_shaped(darray::ensure(tau), c, sb, "tau");
    r.dtau   = copy_shaped(darray::ensure(dtau), c, sb, "dtau");
    r.ds     = copy_shaped(darray::ensure(ds), c, s1, "ds");
    r.vshell = copy_shaped(darray::ensure(vshell), c, s1, "vshell");
    return r;
  }
};

}  // namespace

// ##########################################################################

PYBIND11_MODULE(_core, m)
{
  m.doc() = "Python bindings for NEMO (non-equilibrium multi-ion solver).";
  m.attr("__version__") = NEMO_PY_VERSION;

  nemo_error_type =
      py::exception<NEMO::nemo_error>(m, "NemoError", PyExc_RuntimeError)
          .release();
  py::register_exception_translator([](std::exception_ptr p) {
    try {
      if (p) std::rethrow_exception(p);
    }
    catch (const NEMO::nemo_error &e) {
      py::object err   = nemo_error_type(e.what());  // NemoError(msg)
      err.attr("code") = e.code();
      PyErr_SetObject(nemo_error_type.ptr(), err.ptr());
    }
  });

  m.attr("SED_ATLAS_OB")  = RT_EFFECT_MP10_ATLAS_O3V;
  m.attr("SED_POWR_WC")   = RT_EFFECT_MP10_POWR_WC;
  m.attr("SED_POWR_WNE")  = RT_EFFECT_MP10_POWR_WNE;
  m.attr("SED_BLACKBODY") = RT_EFFECT_MP10_BLACKBODY;

  m.def(
      "set_log_level",
      [](const std::string &level) {
        auto l = spdlog::level::from_str(level);
        if (l == spdlog::level::off && level != "off")
          raise_value_error("unknown log level '" + level + "'");
        spdlog::set_level(l);
      },
      py::arg("level"),
      "Set the level of NEMO's log messages (to stderr): 'trace', 'debug', "
      "'info', 'warning', 'error', 'critical' or 'off'.");

  // Options ----------------------------------------------------------------
  py::class_<NEMO::nemo_options>(
      m, "Options",
      "Physics and solver options (NEMO::nemo_options).  The defaults are "
      "PION's.  nv_prim and n_scalar = -1 mean: density and pressure, "
      "followed by the tracers.")
      .def(
          py::init([](const py::kwargs &kw) {
            py::object o = py::cast(default_options());
            for (auto item : kw) {
              if (!py::hasattr(o, item.first))
                throw py::type_error(
                    "Options: unknown option '"
                    + py::str(item.first).cast<std::string>() + "'");
              py::setattr(o, item.first, item.second);
            }
            return o.cast<NEMO::nemo_options>();
          }),
          "Options(**kwargs): PION's defaults, with any field set by "
          "keyword.")
      .def(
          "__repr__",
          [](const py::object &self) {
            std::string s       = "Options(";
            const char *names[] = {
                "reltol",
                "abstol",
                "minerg",
                "min_temperature",
                "max_temperature",
                "nv_prim",
                "n_scalar",
                "prim_dens_idx",
                "prim_pres_idx",
                "timestep_limiter",
                "update_energy",
                "cooling",
                "collisional_ionization",
                "photoionization",
                "charge_exchange",
                "cosmicray_ionization",
                "radiative_recombination",
                "onthespot_radiation"};
            bool first = true;
            for (const char *n : names) {
              s += (first ? "" : ", ") + std::string(n) + "="
                   + py::repr(self.attr(n)).cast<std::string>();
              first = false;
            }
            return s + ")";
          })
      .def_readwrite("reltol", &NEMO::nemo_options::nemo_reltol)
      .def_readwrite("abstol", &NEMO::nemo_options::nemo_abstol)
      .def_readwrite("minerg", &NEMO::nemo_options::nemo_minerg)
      .def_readwrite("min_temperature", &NEMO::nemo_options::Min_Temperature)
      .def_readwrite("max_temperature", &NEMO::nemo_options::Max_Temperature)
      .def_readwrite("nv_prim", &NEMO::nemo_options::nv_prim)
      .def_readwrite("n_scalar", &NEMO::nemo_options::n_scalar)
      .def_readwrite("prim_dens_idx", &NEMO::nemo_options::prim_dens_idx)
      .def_readwrite("prim_pres_idx", &NEMO::nemo_options::prim_pres_idx)
      .def_readwrite("timestep_limiter", &NEMO::nemo_options::timestep_limiter)
      .def_readwrite("update_energy", &NEMO::nemo_options::update_energy)
      .def_readwrite("cooling", &NEMO::nemo_options::cooling)
      .def_readwrite(
          "collisional_ionization", &NEMO::nemo_options::collisional_ionization)
      .def_readwrite("photoionization", &NEMO::nemo_options::photoionization)
      .def_readwrite("charge_exchange", &NEMO::nemo_options::charge_exchange)
      .def_readwrite(
          "cosmicray_ionization", &NEMO::nemo_options::cosmicray_ionization)
      .def_readwrite(
          "radiative_recombination",
          &NEMO::nemo_options::radiative_recombination)
      .def_readwrite(
          "onthespot_radiation", &NEMO::nemo_options::onthespot_radiation);

  // Radiation sources -------------------------------------------------------
  py::class_<NEMO::nemo_star_properties>(
      m, "Source",
      "An ionising point source (NEMO::nemo_star_properties).  After the "
      "Nemo object is created, Nemo.sources holds copies with the energy "
      "bins filled in.")
      .def(
          py::init(
              [](double luminosity, double teff, double radius, int sed_type) {
                NEMO::nemo_star_properties s;
                s.Luminosity = luminosity;
                s.Teff       = teff;
                s.Radius     = radius;
                s.SED_type   = sed_type;
                s.photo_emax = 0.0;
                s.Nbins      = 0;
                s.id         = 0;
                return s;
              }),
          py::arg("luminosity"), py::arg("teff"), py::arg("radius"),
          py::arg("sed_type") = RT_EFFECT_MP10_BLACKBODY,
          "luminosity in erg/s, teff in K, radius in Rsun, sed_type one of "
          "the SED_* constants.")
      .def_readwrite("luminosity", &NEMO::nemo_star_properties::Luminosity)
      .def_readwrite("teff", &NEMO::nemo_star_properties::Teff)
      .def_readwrite("radius", &NEMO::nemo_star_properties::Radius)
      .def_readwrite("sed_type", &NEMO::nemo_star_properties::SED_type)
      .def_readonly("sed_info", &NEMO::nemo_star_properties::SED_info)
      .def_readonly("n_bins", &NEMO::nemo_star_properties::Nbins)
      .def_readonly("photo_emax", &NEMO::nemo_star_properties::photo_emax)
      .def_readonly("energy_bins", &NEMO::nemo_star_properties::Ebins)
      .def_readonly(
          "bin_luminosity", &NEMO::nemo_star_properties::bin_luminosity);

  // Solver -----------------------------------------------------------------
  const char *rt_doc =
      "With photoionization, also give for each cell and source: tau and "
      "dtau, shape (..., n_sources, n_bins), the optical depth from the "
      "source to the cell and through the cell; ds and vshell, shape "
      "(..., n_sources), the path length through the cell (cm) and the "
      "shell volume 4/3 pi (r+^3 - r-^3) (cm^3).";

  py::class_<py_nemo>(
      m, "Nemo",
      "NEMO solver for one set of tracers and options (NEMO::nemo).\n\n"
      "State vectors have nv_prim entries: density (g/cm^3), pressure "
      "(dyne/cm^2), any other variables, then the tracers (mass fractions) "
      "in the order of the tracer list.  Methods accept one state vector "
      "(shape (nv_prim,)) or many (shape (..., nv_prim)).")
      .def(
          py::init([](std::vector<std::string> tracers,
                      const py::object &options,
                      std::vector<NEMO::nemo_star_properties> sources,
                      double gamma) {
            NEMO::nemo_options o = options.is_none() ?
                                       default_options() :
                                       options.cast<NEMO::nemo_options>();
            return new py_nemo(
                o, std::move(tracers), std::move(sources), gamma);
          }),
          py::arg("tracers"), py::kw_only(), py::arg("options") = py::none(),
          py::arg("sources") = std::vector<NEMO::nemo_star_properties>(),
          py::arg("gamma")   = 5.0 / 3.0,
          "tracers: names of the tracer variables (e.g. ['X_H', 'X_He', "
          "'H', 'He', 'He1+']); options: an Options (default Options()); "
          "sources: list of Source, used if options.photoionization; gamma: "
          "adiabatic index.")
      .def(
          "__repr__",
          [](const py_nemo &n) {
            std::string t;
            for (const auto &s : n.tracers())
              t += (t.empty() ? "'" : ", '") + s + "'";
            return "Nemo(tracers=[" + t
                   + "], nv_prim=" + std::to_string(n.nv_prim())
                   + ", n_sources=" + std::to_string(n.sources().size()) + ")";
          })
      .def_property_readonly("options", &py_nemo::options)
      .def_property_readonly("tracers", &py_nemo::tracers)
      .def_property_readonly("sources", &py_nemo::sources)
      .def_property_readonly("gamma", &py_nemo::gamma)
      .def_property_readonly("nv_prim", &py_nemo::nv_prim)
      .def_property_readonly(
          "n_bins", &py_nemo::nbins, "number of photon-energy bins")
      .def(
          "index", &py_nemo::index, py::arg("tracer"),
          "Position of a tracer in the state vector.")
      .def(
          "update", &py_nemo::update, py::arg("prim"), py::arg("dt"),
          py::kw_only(), py::arg("tau") = py::none(),
          py::arg("dtau") = py::none(), py::arg("ds") = py::none(),
          py::arg("vshell") = py::none(), py::arg("on_error") = "raise",
          (std::string("Advance the ionisation (and energy) of each cell "
                       "by dt seconds (scalar or per cell), and return the "
                       "new state.  If the update of a cell fails, raise "
                       "NemoError (on_error='raise') or set that cell to "
                       "NaN (on_error='nan').  ")
           + rt_doc)
              .c_str())
      .def(
          "temperature", &py_nemo::temperature, py::arg("prim"),
          "Gas temperature (K).")
      .def(
          "set_temperature", &py_nemo::set_temperature, py::arg("prim"),
          py::arg("T"),
          "Return the state with the pressure set for temperature T (K; "
          "scalar or per cell).")
      .def(
          "timescale", &py_nemo::timescale, py::arg("prim"), py::kw_only(),
          py::arg("tau") = py::none(), py::arg("dtau") = py::none(),
          py::arg("ds") = py::none(), py::arg("vshell") = py::none(),
          (std::string("Shortest microphysics timescale (s), for limiting a "
                       "simulation timestep.  ")
           + rt_doc)
              .c_str())
      .def(
          "electron_density", &py_nemo::electron_density, py::arg("prim"),
          "Electron number density (cm^-3).")
      .def(
          "ion_density", &py_nemo::ion_density, py::arg("ion"), py::arg("prim"),
          "Number density (cm^-3) of a tracked ion, e.g. 'H', 'He1+'.")
      .def(
          "dtau", &py_nemo::dtau, py::arg("prim"), py::arg("ds"),
          "Optical depth in each energy bin through a path length ds (cm; "
          "scalar or per cell); shape (..., n_bins).  Needs "
          "photoionization.")
      .def(
          "scma", &py_nemo::scma, py::arg("prim"),
          "Return the state with the mass fractions made consistent (sCMA, "
          "Plewa & Muller 1999): elements sum to 1, and the ions of each "
          "element to at most the element.");
}
