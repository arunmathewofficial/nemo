"""Python versions of NEMO's C++ regression tests (tests/*), compared with
the same reference data and tolerances.

The reference data are read from tests/data in the NEMO repository; the
tests are skipped if that directory is not there (e.g. in an unpacked
sdist).  Tests marked "slow" take 10 s or more each (about 2 minutes in
all); skip them with `pytest -m "not slow"`.
"""

from pathlib import Path

import numpy as np
import pytest

import nemo_plasma as nemo

REPO_DATA = Path(__file__).resolve().parents[2] / "tests" / "data"
LOCAL_DATA = Path(__file__).resolve().parent / "data"

needs_data = pytest.mark.skipif(
    not REPO_DATA.is_dir(), reason="needs tests/data from the NEMO repository")

ELEMENTS = ["H", "He", "C", "N", "O", "Ne", "Si", "S", "Fe"]
ATOMIC_NUMBER = {
    "H": 1,
    "He": 2,
    "C": 6,
    "N": 7,
    "O": 8,
    "Ne": 10,
    "Si": 14,
    "S": 16,
    "Fe": 26
}
# Masses (g) from NEMO's atomic_physics_database.
ATOMIC_MASS = {
    "H": 1.6738e-24,
    "He": 6.6464768e-24,
    "C": 1.994374e-23,
    "N": 2.325892e-23,
    "O": 2.6567628e-23,
    "Ne": 3.3509177e-23,
    "Si": 4.6637066e-23,
    "S": 5.3245181e-23,
    "Fe": 9.2732796e-23
}
# Solar-like element mass fractions used by several tests.
SOLAR = [
    0.7381, 0.2495736, 0.002368, 0.000693, 0.005739, 0.0012579, 0.0006656,
    0.0003096, 0.0012933
]

PC = 3.085677581e18  # cm


def compare_data(new, ref, rel_tol, abs_tol):
    """The test of tests/compare_data.py: each value must agree within
    abs_tol, or within rel_tol relative to the mean of the two values."""
    new, ref = np.asarray(new), np.asarray(ref)
    assert new.shape == ref.shape
    diff = np.abs(new - ref)
    mean = 0.5 * (np.abs(new) + np.abs(ref))
    rel = np.divide(diff, mean, out=np.zeros_like(diff), where=mean != 0)
    bad = (diff > abs_tol) & (rel > rel_tol)
    if bad.any():
        idx = np.argwhere(bad)[:5]
        detail = ", ".join(
            f"{tuple(int(j) for j in i)}: {new[tuple(i)]:.6e} vs "
            f"{ref[tuple(i)]:.6e}" for i in idx)
        pytest.fail(f"{bad.sum()} values differ from the reference: {detail}")


def read_table(name):
    """Header names and data of a reference file with a one-line header."""
    path = REPO_DATA / name
    with open(path) as f:
        header = f.readline().split()
    return header, np.loadtxt(path, skiprows=1, ndmin=2)


def all_ion_tracers(elements):
    """Element tracers, then all ions (neutral to Z-1) of each element."""
    tracers = ["X_" + e for e in elements]
    for e in elements:
        tracers += [e] + [f"{e}{q}+" for q in range(1, ATOMIC_NUMBER[e])]
    return tracers


def element_of(ion):
    return ion.rstrip("+0123456789")


def charge_of(ion):
    q = ion[len(element_of(ion)):].rstrip("+")
    return int(q) if q else 0


def charge_density(prim, tracers):
    """Charge density / e (cm^-3), as in tests/cx-test."""
    total = 0.0
    for e_name in [t[2:] for t in tracers if t.startswith("X_")]:
        X_e = prim[2 + tracers.index("X_" + e_name)]
        ions = [
            i for i, t in enumerate(tracers)
            if not t.startswith("X_") and element_of(t) == e_name
        ]
        charge = sum(charge_of(tracers[i]) * prim[2 + i] for i in ions)
        charge += ATOMIC_NUMBER[e_name] * (X_e - sum(prim[2 + i]
                                                     for i in ions))
        total += charge / ATOMIC_MASS[e_name]
    return prim[0] * total


# ---------------------------------------------------------------------------
# Collisional ionisation equilibrium (tests/cie-test)


@pytest.mark.slow
@needs_data
def test_cie():
    header, ref = read_table("cie-test-output.txt")
    tracers = header[1:]
    prim0 = np.zeros(2 + len(tracers))
    prim0[0], prim0[1] = 1.0e-24, 1.0e-9
    prim0[2:11] = [
        0.7381, 0.2495, 0.002368, 0.000693, 0.005739, 0.0012579, 0.0006656,
        0.0003096, 0.0012933
    ]
    opt = nemo.Options(update_energy=False,
                       cooling=False,
                       charge_exchange=True,
                       cosmicray_ionization=True,
                       onthespot_radiation=True)
    s = nemo.Nemo(tracers, options=opt)
    T = np.geomspace(1.0e4, 5.0e8, 100)
    prim = s.set_temperature(np.tile(prim0, (100, 1)), T)
    out = s.update(prim, 1.0e15)
    new = np.column_stack([s.temperature(out), out[:, 2:]])
    compare_data(new, ref, 3.0e-5, 1.0e-10)


# ---------------------------------------------------------------------------
# Isochoric radiative cooling at fixed ionisation (tests/isochoric-cooling)


@needs_data
def test_isochoric_cooling():
    ref = np.loadtxt(REPO_DATA / "isochoric_cooling_reference.txt")
    t_ref, T_ref, dTdt = ref.T
    tracers = ["X_" + e for e in ELEMENTS]
    prim0 = [3.847988437033092e-22, 1.0e-9] + SOLAR
    for e, X in zip(ELEMENTS, SOLAR):
        Z = ATOMIC_NUMBER[e]
        tracers += [e] + [f"{e}{q}+" for q in range(1, Z)]
        prim0 += [X / (Z + 1)] * Z
    prim0 = np.array(prim0)
    opt = nemo.Options(min_temperature=1.0e3,
                       reltol=1.0e-8,
                       timestep_limiter=0,
                       collisional_ionization=False,
                       charge_exchange=False,
                       radiative_recombination=False)
    s = nemo.Nemo(tracers, options=opt)
    # Each interval starts from the reference temperature at its start.
    prim = np.tile(prim0, (len(t_ref) - 1, 1))
    prim[:, 1] *= T_ref[:-1] / s.temperature(prim0)
    interval = np.diff(t_ref)
    out = s.update(prim, interval)
    np.testing.assert_allclose(out[:, 0], prim0[0], rtol=1e-8, atol=0)
    np.testing.assert_allclose(out[:, 2:], prim[:, 2:], rtol=1e-8, atol=0)
    T = s.temperature(out)
    dT = np.abs(T - T_ref[1:])
    time_error = dT / np.abs(dTdt[1:]) / interval
    ok = (dT <= 1.0e-6 * T_ref[1:]) | (time_error <= 1.0e-5)
    assert ok.all(), f"intervals failing: {np.flatnonzero(~ok) + 1}"


# ---------------------------------------------------------------------------
# Non-equilibrium isochoric cooling from CIE (tests/isochoric-nei)


@needs_data
def test_isochoric_nei():
    ref = np.loadtxt(REPO_DATA / "isochoric_nei_reference.txt")
    tracers = all_ion_tracers(ELEMENTS)
    ions = tracers[len(ELEMENTS):]
    n_ions = len(ions)
    t, T_ref, dTdt = ref[:, 0], ref[:, 1], ref[:, 2]
    y_ref, dydt = ref[:, 3:3 + n_ions], ref[:, 3 + n_ions:]
    X_el = np.array([SOLAR[ELEMENTS.index(element_of(i))] for i in ions])

    opt = nemo.Options(min_temperature=1.0e3,
                       reltol=1.0e-8,
                       timestep_limiter=0,
                       cosmicray_ionization=False)
    s = nemo.Nemo(tracers, options=opt)
    n = len(t) - 1
    prim = np.zeros((n, 2 + len(tracers)))
    prim[:, 0] = 3.847988437033092e-22
    prim[:, 1] = 1.0e-9
    prim[:, 2:2 + len(ELEMENTS)] = SOLAR
    prim[:, 2 + len(ELEMENTS):] = y_ref[:-1] * X_el
    prim[:, 1] *= T_ref[:-1] / s.temperature(prim)
    interval = np.diff(t)
    out = s.update(prim, interval)

    def agrees(value, ref_value, dref_dt, rel_tol, abs_tol):
        err = np.abs(value - ref_value)
        ok = err <= rel_tol * np.abs(ref_value) + abs_tol
        with np.errstate(divide="ignore", invalid="ignore"):
            time_error = np.where(dref_dt != 0,
                                  err / np.abs(dref_dt) / interval[:, None],
                                  np.inf)
        return np.isfinite(value) & (ok | (time_error <= 1.0e-5))

    T = s.temperature(out)[:, None]
    assert agrees(T, T_ref[1:, None], dTdt[1:, None], 1e-6, 0.0).all()
    y = out[:, 2 + len(ELEMENTS):] / X_el
    checked = y_ref[1:] >= 1.0e-6
    ok = agrees(y, y_ref[1:], dydt[1:], 1e-6, 1e-10) | ~checked
    assert ok.all(), (
        f"ions failing: {sorted({ions[j] for j in np.argwhere(~ok)[:, 1]})}")


# ---------------------------------------------------------------------------
# Charge exchange (tests/cx-test)

CX_OPTIONS = dict(update_energy=False,
                  cooling=False,
                  collisional_ionization=False,
                  charge_exchange=True,
                  radiative_recombination=False)


def cx_two_element(el1, el2, T, t_end, start):
    """el1 (H or He) with n = 1 cm^-3, el2 with n = 1e-4 cm^-3.  start
    gives the initial mass fractions of the ions."""
    rho = ATOMIC_MASS[el1] + 1.0e-4 * ATOMIC_MASS[el2]
    X1, X2 = ATOMIC_MASS[el1] / rho, 1.0e-4 * ATOMIC_MASS[el2] / rho
    tracers = ["X_" + el1, "X_" + el2]
    for e in (el1, el2):
        tracers += [e] + [f"{e}{q}+" for q in range(1, ATOMIC_NUMBER[e])]
    prim = np.zeros(2 + len(tracers))
    prim[:4] = [rho, 1.0e-12, X1, X2]
    for name, frac in start.items():
        X = X1 if element_of(name) == el1 else X2
        prim[2 + tracers.index(name)] = frac * X
    s = nemo.Nemo(tracers, options=nemo.Options(**CX_OPTIONS))
    prim = s.set_temperature(prim, T)
    out = s.update(prim, t_end)
    assert s.temperature(out) == pytest.approx(T, rel=2e-6)
    return s, tracers, prim, out


@pytest.mark.parametrize("ion, T, t_end, ratio",
                         [("O", 1000.0, 1.0e13, 5.1811635240e-01),
                          ("N", 1.0e4, 5.0e13, 1.5199216941e+00),
                          ("C", 8000.0, 1.0e16, 1.3798769809e+03)])
def test_cx_equilibrium_with_hydrogen(ion, T, t_end, ratio):
    """X+ + H <-> X + H+ reaches n(X+) n(H) / [n(X) n(H+)] = k_rev/k_fwd."""
    s, tr, _, out = cx_two_element("H", ion, T, t_end, {
        "H": 0.5,
        ion + "1+": 1.0
    })
    X_H0 = out[2 + tr.index("H")]
    X_Hp = out[2 + tr.index("X_H")] - X_H0
    X_0 = out[2 + tr.index(ion)]
    X_p = out[2 + tr.index(ion + "1+")]
    assert X_p * X_H0 / (X_0 * X_Hp) == pytest.approx(ratio, rel=2e-6)


def test_cx_equilibrium_c2p_helium():
    """C2+ + He <-> C+ + He+, and charge is conserved."""
    s, tr, prim, out = cx_two_element("He", "C", 8000.0, 1.0e13, {
        "He": 0.5,
        "He1+": 0.5,
        "C2+": 1.0
    })
    v = {n: out[2 + tr.index(n)] for n in ("He", "He1+", "C1+", "C2+")}
    r = v["C2+"] * v["He"] / (v["C1+"] * v["He1+"])
    assert r == pytest.approx(3.4590722950e-06, rel=2e-6)
    assert charge_density(out, tr) == pytest.approx(charge_density(prim, tr),
                                                    rel=2e-6)


def test_cx_full_network():
    """All ions of 9 elements, one update with only charge exchange."""
    tracers = all_ion_tracers(ELEMENTS)
    prim = [ATOMIC_MASS["H"] / SOLAR[0], 1.0e-12] + SOLAR
    for e, X in zip(ELEMENTS, SOLAR):
        Z = ATOMIC_NUMBER[e]
        prim += [X / (Z + 1)] * Z
    prim = np.array(prim)
    s = nemo.Nemo(tracers, options=nemo.Options(**CX_OPTIONS))
    prim = s.set_temperature(prim, 8000.0)
    out = s.update(prim, 1.0e7)
    ref = np.loadtxt(LOCAL_DATA / "cx_full_network_reference.txt")
    # compare_scalar of cx-test: rel. to the sum, or abs 1e-14
    diff = np.abs(out - ref)
    assert ((diff <= 1e-14) | (diff <= 1e-6 *
                               (np.abs(out) + np.abs(ref)))).all()
    assert charge_density(out, tracers) == pytest.approx(charge_density(
        prim, tracers),
                                                         rel=2e-6)


# ---------------------------------------------------------------------------
# Photoionised spherical shells (tests/test4-photoion-eq, tests/lexington)


def evolve_shells(s, prim, n_cells, t_end=3.156e13, ds_factor=None):
    """Evolve n_cells spherical shells (0 to 5.76 pc) around the sources of
    s, as in tests/common_functions.cpp and the two tests.  ds_factor gives
    a factor (1 or 0) per source for the path length through each shell."""
    n_src, nb = len(s.sources), s.n_bins
    if ds_factor is None:
        ds_factor = np.ones(n_src)
    r = np.linspace(0.0, 5.76 * PC, n_cells + 1)
    dr = r[1] - r[0]
    ds = np.outer(np.ones(n_cells), dr * np.asarray(ds_factor))
    vshell = 4.0 * np.pi / 3.0 * (r[1:]**3 - r[:-1]**3)
    vshell = np.column_stack([vshell] * n_src)
    vshell[:, np.asarray(ds_factor) == 0] = 1.0e200
    state = np.tile(prim, (n_cells, 1))
    time = 0.0
    while time < t_end:
        # optical depth through each shell, and from the source to it
        dtau = np.stack([s.dtau(state, ds[:, k]) for k in range(n_src)],
                        axis=1)
        tau = np.cumsum(dtau, axis=0) - dtau
        rt = dict(tau=tau, dtau=dtau, ds=ds, vshell=vshell)
        remaining = t_end - time
        dt = min(remaining, s.timescale(state, **rt).min())
        dt = min(max(dt, 0.1 * time), remaining)
        state = s.update(state, dt, **rt)
        time += dt
    out = state.copy()
    out[:, 1] = s.temperature(state)
    out[:, 0] = 0.5 * (r[1:] + r[:-1])
    return out, tau[:, 0, :]


HII_ELEMENTS = ["H", "He", "C", "N", "O", "Ne", "S"]
HII_X = [
    7.0933e-01, 2.8373e-01, 1.8725e-03, 3.9726e-04, 3.7452e-03, 7.1647e-04,
    2.0492e-04
]


def hii_initial(tracers, rho, p):
    prim = np.zeros(2 + len(tracers))
    prim[:2] = rho, p
    prim[2:2 + len(HII_X)] = HII_X
    for name, value in {
            "H": 7.0933e-01,
            "He": 2.8373e-01,
            "C": 1.0e-3,
            "C1+": 0.8725e-3,
            "N": 3.9726e-04,
            "O": 3.7452e-03,
            "Ne": 7.1647e-04,
            "S": 2.0492e-04
    }.items():
        prim[2 + tracers.index(name)] = value
    return prim


def star(luminosity=3.07e39, radius=18.6517, teff=4.0e4):
    return nemo.Source(luminosity, teff, radius, nemo.SED_BLACKBODY)


@pytest.mark.slow
@needs_data
def test_photoionisation_equilibrium():
    """tests/test4-photoion-eq: 32 shells, photoionisation and
    recombination only, isothermal."""
    header, ref = read_table("pieq_test_refdata.txt")
    tracers = header[2:header.index("[plus")]
    assert tracers == all_ion_tracers(HII_ELEMENTS)
    opt = nemo.Options(update_energy=False,
                       cooling=False,
                       collisional_ionization=False,
                       photoionization=True,
                       charge_exchange=False,
                       onthespot_radiation=True)
    s = nemo.Nemo(tracers, options=opt, sources=[star()])
    out, tau = evolve_shells(s, hii_initial(tracers, 2.36e-22, 2.43e-10), 32)
    compare_data(np.column_stack([out, tau]), ref, 1.0e-4, 1.0e-9)


def hii40_solver(tracers, sources):
    opt = nemo.Options(min_temperature=1.0e3,
                       photoionization=True,
                       cosmicray_ionization=True,
                       onthespot_radiation=True)
    return nemo.Nemo(tracers, options=opt, sources=sources)


HII40_RHO = 2.3589999999999998e-22


@pytest.mark.slow
@needs_data
def test_hii40_lexington():
    """tests/lexington (run-hii40-test.sh): 20 shells, all processes."""
    header, ref = read_table("hii40_test_output.txt")
    tracers = header[2:]
    s = hii40_solver(tracers, [star()])
    out, _ = evolve_shells(s, hii_initial(tracers, HII40_RHO, 1.0e-12), 20)
    compare_data(out, ref, 3.0e-5, 1.0e-10)


@pytest.mark.slow
def test_hii40_two_sources():
    """tests/lexington (run-hii40-2src.sh): the star split into two at the
    same place gives the same result, and a second source with no flux
    changes nothing."""
    tracers = all_ion_tracers(HII_ELEMENTS)
    prim = hii_initial(tracers, HII40_RHO, 1.0e-12)
    single, _ = evolve_shells(hii40_solver(tracers, [star()]), prim, 20)
    half = star(3.07e39 * 0.5, 18.6517 / np.sqrt(2.0))
    split, _ = evolve_shells(hii40_solver(tracers, [half, half]), prim, 20)
    compare_data(split, single, 1.0e-4, 1.0e-11)
    teff, radius = 3.0e4, 10.0
    lum = 4.0 * np.pi * 5.670374e-5 * teff**4 * (radius * 6.957e10)**2
    other = star(lum, radius, teff)
    absent, _ = evolve_shells(hii40_solver(tracers, [star(), other]),
                              prim,
                              20,
                              ds_factor=[1, 0])
    compare_data(absent, single, 1.0e-12, 1.0e-20)
