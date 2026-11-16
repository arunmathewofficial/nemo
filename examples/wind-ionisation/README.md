# NEMO example: wind tracers in photoionisation equilibrium (by Claude AI)

PION sets the gas inside the boundary of a wind source (radius
`WIND_<k>_radius`) to fixed values, including the ion mass fractions
`WIND_<k>_TR<n>`.  If the star's own radiation keeps the gas in a different
ionisation state, the cells just outside the boundary change their
ionisation state abruptly (e.g. C2+ at the boundary, C3+ in the first
cells outside it).

`wind-ionisation` reads a PION parameter file and, for each wind source,
evolves one cell just outside the boundary with NEMO until the ion fractions
no longer change.  It then prints all `WIND_<k>_TR<n>` lines of that wind
source, which can replace the block in the parameter file:

- the element mass fractions (`X_El`).  If they do not add up to 1, X_He
  is changed by the difference (or, without He, the most abundant element),
  and its ions are scaled with it;
- the ions, in equilibrium;
- the other tracers (e.g. wind and shock tracers), unchanged.

The cell is set up like PION's:

- density `Mdot / (4 pi r^2 v)` at `r = WIND_<k>_radius + ds/2`, with
  `v = vinf` (or the beta = 1 law if `WIND_<k>_acceleration` is 1);
- temperature as PION sets it in the wind (adiabatic expansion from
  `WIND_<k>_temp` at `WIND_<k>_Rstr`, but at least `EP_Min_Temperature`),
  kept fixed;
- the radiation of the source at the position of the star, unattenuated,
  because PION treats the gas inside a wind boundary as transparent;
- cell size `ds` of the finest grid level that contains the star (it
  matters for very optically thick ions: NEMO's rates are photon
  conserving);
- the NEMO options from the `EP_` parameters.

Options:

| option | effect |
|---|---|
| `--wind <k>` | only wind source k |
| `--ds <cm>` | cell size |
| `--T <K>` | temperature |
| `--thermal` | evolve the temperature too (thermal equilibrium) |
| `--all-sources` | also include the other radiation sources, unattenuated, at their distance from the wind boundary (an upper limit) |

Build it against an installed NEMO (1.3 or later) and run it on the example
parameter file, a WC Wolf-Rayet star whose wind tracers put all carbon in
C2+:
```bash
cmake -S . -B build -DCMAKE_PREFIX_PATH="$HOME/.local"
cmake --build build
./build/wind-ionisation example_wc_wind.txt
```

The output should be close to:
```
# Wind source 0 at x = (0.0000e+00, 0.0000e+00, 0.0000e+00) cm: boundary radius 1.5000e+13 cm, cell size 2.7344e+12 cm
#   rho = 1.4403e-15 g/cm^3, n_e = 2.1116e+08 cm^-3, T = 5.0000e+03 K (fixed), radiation of source 0
#   converged after 24 steps (3.367e+04 s), largest change in the last step 4.1e-11
#   element mass fractions add up to 1
#   ion           input  equilibrium
#   C2+      1.0000e+00   3.5710e-03
#   C3+      2.5000e-12   9.9643e-01
#   He1+     1.0000e+00   9.9999e-01
#   O2+      1.0000e+00   9.9994e-01
# New tracer values for the parameter file (all tracers of wind source 0):
WIND_0_TR0   5.500000000000e-01   # X_He (was 5.500000e-01)
WIND_0_TR1   5.675437212396e-09   # He (was 1.000000e-12)
WIND_0_TR2   5.499954212658e-01   # He1+ (was 5.500000e-01)
WIND_0_TR3   4.000000000000e-01   # X_C (was 4.000000e-01)
...
WIND_0_TR6   1.428405124833e-03   # C2+ (was 4.000000e-01)
WIND_0_TR7   3.985711415467e-01   # C3+ (was 1.000000e-12)
...
```
The ion table lists ions with more than 1e-4 of their element, as fractions
of the element.  The new tracer values are mass fractions.  Ion mass
fractions below 1e-12 are set to 1e-12, and the most abundant ion is reduced
so that the ions of an element add up to at most (1 - 1e-11) times the
element, so that they never add up to more than the element even after
rounding of the printed values.

This is approximate: the cells just outside the boundary also see the
radiation of other stars, attenuated by the gas in between, and their
temperature is set by the hydrodynamics and heating and cooling.  Compare
with the ion fractions in those cells in a snapshot.
