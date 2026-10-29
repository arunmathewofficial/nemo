# Minimal NEMO example (by Claude AI)

A standalone project showing how to use an installed copy of NEMO from CMake.
It sets up NEMO for a single cell with H and He, evolves it for 1e14 seconds
with collisional ionisation, charge exchange, cosmic-ray ionisation and
radiative recombination, and prints the resulting temperature and state vector.

First install NEMO (see the top-level [README](../../README.md)), e.g. into
`~/.local`.  Then, from this directory:
```bash
cmake -S . -B build -DCMAKE_PREFIX_PATH="$HOME/.local"
cmake --build build
./build/nemo-example
```

The output should be close to:
```
Temperature: 74897.4 K
  output[0] = 1e-22
  output[1] = 1e-09
  output[2] = 1
  output[2] = 0.72
  output[3] = 0.28
  output[4] = 1.45928e-05
  output[5] = 4.84255e-05
  output[6] = 0.130434
```

The only NEMO-specific CMake lines are in [CMakeLists.txt](CMakeLists.txt):
```cmake
find_package(NEMO 1.0 REQUIRED)
target_link_libraries(nemo-example PRIVATE NEMO::nemo)
```
No MPI compiler wrapper is needed: `NEMO::nemo` passes on the MPI include
paths and libraries.
