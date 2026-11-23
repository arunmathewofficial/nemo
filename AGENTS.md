# NEMO: notes for coding agents

NEMO (Non-Equilibrium Multi-ion mOdule) integrates the ionisation and thermal
state of a plasma with many elements and ions, using SUNDIALS CVODE. It is
built as a static library (`libnemo.a`) and used by PION. Method paper:
Mathew et al. (2025), A&A 695, A73. See `README.md` for the API and the rules
for tracer names and mass fractions.

## Formatting (checked by CI)

- Format every C++ file you change with `tools/format.sh <file>`. It uses
  **clang-format-18**, on purpose: it matches the CI static-analysis stage.
  PION uses clang-format-19, so don't mix them up.
- Use the script, not clang-format directly: it hides `#pragma omp` lines from
  clang-format.
- Format only the files you changed (`tools/run_format.sh` formats the whole
  tree).

## Building, testing and installing

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=mpicxx \
  -DNEMO_OMP=ON -DCMAKE_INSTALL_PREFIX="$HOME/.local"
cmake --build build -j 8
make -C build/tests test      # the tests CI runs
cmake --install build         # installs libnemo.a, headers and CMake files
```

- CI builds with `NEMO_OMP=ON`, runs the tests in `tests/`, then installs NEMO
  and builds `examples/minimal` against the installed copy to check
  `find_package(NEMO)` and the `NEMO::nemo` target.
- Unit tests are small programs in `tests/` registered with `add_test`. Some
  tests compare against reference data in `tests/data`; if a physics change
  alters results on purpose, update the reference data in the same commit and
  say why.
- PION links the **installed** library. After a change, reinstall it and
  rebuild PION to test the change there.
- The Python package `nemo-astro` (pybind11, in `python/`) wraps the same
  library: `pip install ".[test]"`, then `pytest -m "not slow"` (about 5 s)
  or `pytest` (all tests). If you change the public C++ API, update the
  Python bindings and their tests too. Releases to PyPI follow
  `python/RELEASING.md`; never upload to PyPI or TestPyPI without the
  maintainers' go-ahead (uploads can't be undone).

## Code conventions

- Everything public is in the C++ namespace `NEMO`. Keep new classes, functions
  and global data inside it, so they can't clash with code that uses NEMO.
- `nemo` objects are not shared between threads: the calling code creates one
  per OpenMP thread. Keep member functions free of shared mutable state.
- Log through spdlog's default logger. The calling program owns it (PION
  can redirect it to a file), so don't create loggers or change logger
  settings inside the library.
- The library never exits the program. `update_single_cell()` and
  `Set_Temp()` return a non-zero error code (e.g. on a CVODE failure); other
  functions, including the constructor, throw `NEMO::nemo_error`
  (`nemo-constants.h`). Log the error with spdlog before reporting it.
  `exit_nemo()` is only for programs that use NEMO (tests, examples), never
  for the library itself.
- C++14.

## Physics details that are easy to get wrong

- Tracer order: elements (`X_H`, `X_He`, ...) by atomic number, then each
  element's ions from neutral upwards. The highest ionisation stage is not
  tracked; it is the element fraction minus the tracked ions.
- Element mass fractions must sum to 1; an element's tracked ions must not
  sum to more than the element. `sCMA()` renormalises after advection and
  must keep each element's ionisation fractions unchanged
  (`tests/scma-test.cpp`).
- `Temperature()` is capped at `Max_Temperature`, so a temperature check will
  not catch an unphysical input state; check density and pressure directly.

## Git

- Development happens in the private GitLab repository
  (https://gitlab.git.nrw/jmackey1/nemo). Branch from `dev` and merge back
  into `dev` through a merge request (fast-forward merges, so keep branches
  rebased on their target). `main` only moves at a release, and fixes that
  are needed before the next release.
- The public GitHub repository (https://github.com/arunmathewofficial/nemo)
  receives only `main` and release tags, at a release. Don't push anything
  else there, and don't merge pull requests on GitHub: bring their commits
  into a GitLab branch, keeping the original author (see `CONTRIBUTING.md`).
- If a change alters the public API (`nemo.h`, `nemo-datastructs.h`), PION's
  `source/microphysics/nemo_interface.cpp` has to change as well. Merge the
  NEMO change first.
