"""Basic tests of the nemo_plasma bindings."""

import threading

import numpy as np
import pytest

import nemo_plasma as nemo

TRACERS = ["X_H", "X_He", "H", "He", "He1+"]
PRIM = np.array([1.0e-22, 1.0e-12, 0.72, 0.28, 0.36, 0.05, 0.1])


def test_minimal_example():
    """Same set-up and result as the C++ examples/minimal."""
    opt = nemo.Options(update_energy=False,
                       cooling=False,
                       cosmicray_ionization=True)
    tr = ["unrelated_mass_scalar"] + TRACERS
    s = nemo.Nemo(tr, options=opt)
    p = np.array([1.0e-22, 1.0e-9, 1.0, 0.72, 0.28, 0.36, 0.05, 0.1])
    out = s.update(p, 1.0e14)
    np.testing.assert_allclose(out[5:], [1.45928e-05, 4.84255e-05, 0.130434],
                               rtol=1e-5)
    np.testing.assert_array_equal(out[:5], p[:5])


def test_options_and_properties():
    o = nemo.Options(cooling=False, reltol=1e-8)
    assert o.cooling is False and o.reltol == 1e-8
    assert "cooling=False" in repr(o)
    with pytest.raises(TypeError):
        nemo.Options(not_an_option=1)
    s = nemo.Nemo(TRACERS)
    assert s.nv_prim == 7
    assert s.options.n_scalar == 5
    assert s.tracers == TRACERS
    assert s.index("X_H") == 2 and s.index("He1+") == 6
    with pytest.raises(ValueError):
        s.index("C")


def test_many_cells_match_one_cell():
    s = nemo.Nemo(TRACERS)
    one = s.update(PRIM, 1.0e10)
    grid = np.tile(PRIM, (3, 4, 1))
    many = s.update(grid, 1.0e10)
    assert many.shape == (3, 4, 7)
    np.testing.assert_array_equal(many, np.broadcast_to(one, many.shape))
    np.testing.assert_array_equal(s.temperature(grid),
                                  np.full((3, 4), s.temperature(PRIM)))
    # dt per cell
    dts = np.array([1.0e9, 1.0e10])
    two = s.update(np.tile(PRIM, (2, 1)), dts)
    np.testing.assert_array_equal(two[1], one)


def test_derived_quantities():
    s = nemo.Nemo(TRACERS)
    T = 1.0e4
    p = s.set_temperature(PRIM, T)
    assert s.temperature(p) == pytest.approx(T, rel=1e-10)
    assert p[0] == PRIM[0] and np.all(p[2:] == PRIM[2:])
    ne = s.electron_density(PRIM)
    assert ne > 0
    m_H = 1.6735575e-24
    assert s.ion_density("H", PRIM) == pytest.approx(PRIM[0] * PRIM[4] / m_H,
                                                     rel=1e-3)
    assert s.timescale(PRIM) > 0
    # sCMA renormalises the element fractions to sum to 1
    bad = PRIM.copy()
    bad[2:4] *= 1.1
    fixed = s.scma(bad)
    assert fixed[2] + fixed[3] == pytest.approx(1.0)


def test_errors():
    with pytest.raises(nemo.NemoError) as e:
        nemo.Nemo(["X_He", "X_H", "H", "He", "He1+"])  # wrong order
    assert e.value.code != 0
    with pytest.raises(ValueError):
        nemo.Nemo(TRACERS, options=nemo.Options(n_scalar=3))
    s = nemo.Nemo(TRACERS)
    with pytest.raises(ValueError):
        s.update(PRIM[:-1], 1.0e10)  # wrong length
    bad = PRIM.copy()
    bad[1] = np.nan
    nemo.set_log_level("off")
    try:
        with pytest.raises(nemo.NemoError) as e:
            s.update(bad, 1.0e10)
        assert e.value.code == 1
        out = s.update(np.array([PRIM, bad]), 1.0e10, on_error="nan")
        assert np.all(np.isfinite(out[0])) and np.all(np.isnan(out[1]))
        with pytest.raises(nemo.NemoError):
            s.timescale(bad)
        with pytest.raises(nemo.NemoError):
            s.ion_density("C2+", PRIM)
    finally:
        nemo.set_log_level("warning")
    # still works afterwards
    assert np.all(np.isfinite(s.update(PRIM, 1.0e10)))


def test_photoionisation():
    src = nemo.Source(3.0e39, 4.0e4, 18.0, nemo.SED_BLACKBODY)
    s = nemo.Nemo(TRACERS,
                  options=nemo.Options(photoionization=True),
                  sources=[src])
    nb = s.n_bins
    assert nb > 0 and s.sources[0].n_bins == nb
    assert len(s.sources[0].energy_bins) == nb
    p = np.array([1.0e-22, 1.0e-12, 0.72, 0.28, 0.72, 0.28, 0.0])
    r, dr = 3.0e18, 1.0e16
    rt = dict(
        tau=np.zeros((1, nb)),
        dtau=s.dtau(p, dr)[None, :],
        ds=np.array([dr]),
        vshell=np.array([4 / 3 * np.pi * ((r + dr / 2)**3 - (r - dr / 2)**3)]),
    )
    out = s.update(p, 3.0e11, **rt)
    assert out[4] / out[2] < 1e-4  # hydrogen ionised
    assert s.temperature(out) > 5000.0
    # radiation data are required with photoionisation, and wrong shapes fail
    with pytest.raises(ValueError):
        s.update(p, 1.0e10)
    with pytest.raises(ValueError):
        s.update(p, 1.0e10, **{**rt, "tau": np.zeros((2, nb))})
    # many cells
    n = 5
    grid = np.tile(p, (n, 1))
    rtg = {k: np.repeat(v[None], n, axis=0) for k, v in rt.items()}
    np.testing.assert_array_equal(s.update(grid, 3.0e11, **rtg),
                                  np.tile(out, (n, 1)))


def test_threads():
    s = nemo.Nemo(TRACERS)
    grid = np.tile(PRIM, (200, 1))
    ref = s.update(grid, 1.0e12)
    results = [None] * 4

    def work(i):
        results[i] = s.update(grid, 1.0e12)  # same object: serialised

    th = [threading.Thread(target=work, args=(i, )) for i in range(4)]
    for t in th:
        t.start()
    for t in th:
        t.join()
    for r in results:
        np.testing.assert_array_equal(r, ref)
