# nemo
Non-Equilibrium Multi-ion mOdule for ionised plasmas

NEMO: 
+ [gitlab](https://gitlab.git.nrw/jmackey1/nemo)
+ [github](https://github.com/jmackey-astro/nemo) 
+ [Description](https://arunmathewofficial.github.io/software/nemo/)
+ [Methods paper](https://doi.org/10.1051/0004-6361/202452373): Mathew, A., Mackey, J., Celeste, M., Haworth, T. J., & Mellema, G. (2025), *"A multi-ion non-equilibrium solver for ionised astrophysical plasmas with arbitrary elemental abundances," Astronomy & Astrophysics, 695, A73.*

# Requirements:

+ CMake
+ MPI wrappers for C++ compiler (e.g. openmpi-bin), and C++ compiler with OpenMP support (e.g. g++)
+ Boost developement library (libboost-dev, >1.73)
+ SUNDIALS development library (libsundials-dev)
+ spdlog (libspdlog-dev)

# Installing and using NEMO in your own code

NEMO is built as a static library, `libnemo.a`.  To install it, together with
its headers and CMake package files, into e.g. `~/.local`:
```bash
cmake -S . -B build \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_CXX_COMPILER=mpicxx \
    -DCMAKE_INSTALL_PREFIX="$HOME/.local"
cmake --build build -j 4
cmake --install build
```

In your own project's `CMakeLists.txt`:
```cmake
find_package(NEMO 1.0 REQUIRED)
target_link_libraries(my-program PRIVATE NEMO::nemo)
```
and configure it with `-DCMAKE_PREFIX_PATH="$HOME/.local"`.  Then
`#include "nemo/nemo.h"` in your code.  `NEMO::nemo` brings in all of the
include paths, compile definitions and libraries that NEMO needs (SUNDIALS,
Boost, MPI, spdlog).

See [examples/minimal](examples/minimal) for a complete example project.

# Using NEMO to integrate the ionisation and thermal state of a plasma

NEMO requires as inputs a vector containing gas density (g/cm3), pressure (dyne/cm2) and mass fractions of elements and ions (collectively called tracer variables).  The constructor takes as an argument a `struct nemo_constants` (declared in `nemo-datastructs.h`) where the relevant fields are:
+ nv_prim: length of primitive state vec
+ n_scalar: number of passive scalars in state vec
+ prim_dens_idx: index of density var (g/cm3) in state vec (usuallly 0)
+ prim_pres_idx: index of pressure var (dyne/cm2) in state vec (usually 1)

The constructor also takes as an argument a list of the tracer variable names (as a `std::vector<std::string>` vector of strings).  This should include all tracers in the vector, whether or not they are to be used by NEMO, because NEMO will take the last `n_scalar` elements in the state vector and scan through the names and only select the ones that correspond to valid elements or ions.

Important points:
+ The primitive state vector should be a STL vector `std::vector<double>`
+ The tracer variables are required to be at the end of the state vector, and so the first tracer is at index `nv_prim-n_scalar`.
+ Tracer names should be in the same order as they appear in the primitive state vector.
+ Elements must be ordered by increasing atomic number in the list of tracer variables.
+ Ions must be ordered by increasing charge and increasing atomic number in the list of tracer variables.
+ Tracer names should be formatted as follows:
    + elements begin with `X_`, e.g. `X_H`, `X_He`, `X_S`, `X_Si` etc.
    + ionic species are distinguished by their charge and NEMO explicitly tracks the neutral species and the highest ionisation stage (bare ion) is not solved explicitly but is solved via the conservation equation from the other species.
    + neutral species are named simply with the element, e.g., `H`, `He`, `S`, `Si` etc.
    + Charged species are named e.g. `H1+`, `Si1+`, `Si2+`, `Si3+`, etc.
    + Any tracer variable that does not match these patterns is assumed to be irrelevant for NEMO and is ignored.
+ Tracer names should be unique and not conflict with other variables in the state vector.
+ Each element should have a sequential list of ions from the neutral up to the 2nd last ionisation stage.
+ The sum of elemental mass fractions should be 1.
+ The sum of ion mass fractions for each element may not exceed the mass fraction of the corresponding element.

An example of a valid tracer list is:
```
std::vector<double> primitive_vec = {1.0e-22, 1.0e-9, 1.0, 0.72, 0.28, 0.36, 0.05, 0.1};
std::vector<std::string> tracers = {"unrelated-mass-scalar", "X_H", "X_He", "H", "He", "He1+"};
```
