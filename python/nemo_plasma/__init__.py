"""Python interface to NEMO, a non-equilibrium multi-ion solver for the
ionisation and thermal state of astrophysical plasmas.

Mathew, Mackey, Celeste, Haworth & Mellema (2025), A&A 695, A73.

Example::

    import numpy as np
    import nemo_plasma as nemo

    tracers = ["X_H", "X_He", "H", "He", "He1+"]
    solver = nemo.Nemo(tracers, options=nemo.Options(cooling=False))
    prim = np.array([1e-22, 1e-12, 0.72, 0.28, 0.36, 0.05, 0.1])
    new = solver.update(prim, dt=1e10)
    print(solver.temperature(new))

State vectors are density (g/cm^3), pressure (dyne/cm^2), then the tracer
mass fractions in the order given.  All methods accept one cell (shape
(nv_prim,)) or many (shape (..., nv_prim)).
"""

from ._core import (
    SED_ATLAS_OB,
    SED_BLACKBODY,
    SED_POWR_WC,
    SED_POWR_WNE,
    Nemo,
    NemoError,
    Options,
    Source,
    __version__,
    set_log_level,
)

# NEMO logs at "info" level by default, which is too chatty for Python use.
set_log_level("warning")

__all__ = [
    "Nemo",
    "NemoError",
    "Options",
    "Source",
    "SED_ATLAS_OB",
    "SED_BLACKBODY",
    "SED_POWR_WC",
    "SED_POWR_WNE",
    "set_log_level",
    "__version__",
]
