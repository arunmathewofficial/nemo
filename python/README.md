# nemo-plasma

Python interface to [NEMO](https://gitlab.git.nrw/jmackey1/nemo), the
Non-Equilibrium Multi-ion mOdule for ionised astrophysical plasmas.  NEMO
integrates the ionisation state of many elements and ions, and the gas
temperature, with the CVODE solver from SUNDIALS.

Method paper: Mathew, A., Mackey, J., Celeste, M., Haworth, T. J., &
Mellema, G. (2025), A&A, 695, A73,
[doi:10.1051/0004-6361/202452373](https://doi.org/10.1051/0004-6361/202452373).

## Installation

```bash
pip install nemo-plasma
```

The wheels contain everything that NEMO needs (SUNDIALS is built in); no MPI
or other libraries are required.  To build from source you need a C++14
compiler; CMake and the other dependencies are downloaded during the build.

## Example

```python
import numpy as np
import nemo_plasma as nemo

tracers = ["X_H", "X_He", "H", "He", "He1+"]   # elements, then ions
solver = nemo.Nemo(tracers, options=nemo.Options(cooling=False))

# density (g/cm^3), pressure (dyne/cm^2), then the tracers (mass fractions)
prim = np.array([1e-22, 1e-12, 0.72, 0.28, 0.36, 0.05, 0.1])
new = solver.update(prim, dt=1e10)       # advance by 1e10 s
print(solver.temperature(new), solver.electron_density(new))

# Many cells at once: any array of shape (..., nv_prim)
grid = np.tile(prim, (100, 100, 1))
new_grid = solver.update(grid, dt=1e10)   # shape (100, 100, 7)
```

Tracer names: elements are `X_H`, `X_He`, `X_C`, ... in order of atomic
number, followed by the ions of each element from the neutral upwards
(`H`, `He`, `He1+`, `C`, `C1+`, ...); the highest ionisation stage of each
element is not listed.  Element mass fractions must sum to 1.  See the
[NEMO README](https://gitlab.git.nrw/jmackey1/nemo) for details.

Errors (bad input, failed integration) raise `nemo_plasma.NemoError`, which
has an error `code`.  `nemo_plasma.set_log_level("info")` shows NEMO's log
messages.

A `Nemo` object can be used from several Python threads (calls are
serialised); for parallel work use one object per thread.  The computation
runs without the GIL.

## Tests

From a clone of the NEMO repository:

```bash
pip install ".[test]"
pytest -m "not slow"   # about 5 s
pytest                 # all tests, about 2 minutes
```

The tests in `python/tests/test_reference.py` are the Python versions of
NEMO's C++ regression tests, and compare with the same reference data
(`tests/data`) and tolerances.
