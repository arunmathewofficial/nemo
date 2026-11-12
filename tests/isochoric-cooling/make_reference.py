#!/usr/bin/env python3
"""Write the reference data for NEMO's isochoric-cooling test.

The gas (composition and density as in isochoric-cooling.cpp) cools from
1e9 K at constant density, with fixed ion fractions and n_e = 100 cm^-3.
This script repeats NEMO's cooling-rate calculation, using the tables and
constants in NEMO's own source code:

1. log10(Lambda) of each ion is interpolated bilinearly in (log T, log n_e)
   from the CHIANTI tables in source/nemo/chianti_cooling.cpp onto NEMO's
   lookup grid (sizes from source/nemo/nemo.cpp).
2. Between lookup-grid nodes, NEMO interpolates linearly in T and n_e.

At fixed n_e the total loss rate L(T) = n_e * sum(n_i Lambda_i) is then
linear in T between grid nodes, so the cooling time
    t = n_tot k_B / (gamma - 1) * integral dT / L(T)
is integrated exactly, interval by interval.

The output has 50 temperatures, log-spaced from 1e9 to 1e4 K, the exact
times at which the gas reaches them, and dT/dt there (the test uses it to
convert a temperature error into a time error).

Run it again whenever the cooling tables or NEMO's interpolation change.
From the top-level NEMO directory:
    python3 tests/isochoric-cooling/make_reference.py . \
        tests/data/isochoric_cooling_reference.txt

usage: make_reference.py <NEMO source dir> <output file>
"""

import re
import sys
from pathlib import Path

import numpy as np
from scipy.interpolate import RegularGridInterpolator

nemo_src = Path(sys.argv[1])
out_file = Path(sys.argv[2])

# --- Test setup (tests/isochoric-cooling/isochoric-cooling.cpp) ------------
elements = ['H', 'He', 'C', 'N', 'O', 'Ne', 'Si', 'S', 'Fe']
mass_fractions = dict(H=0.7381,
                      He=0.2495736,
                      C=0.002368,
                      N=0.000693,
                      O=0.005739,
                      Ne=0.0012579,
                      Si=0.0006656,
                      S=0.0003096,
                      Fe=0.0012933)
rho = 3.847988437033092e-22  # g/cm^3, gives n_e = 100 cm^-3
T_start = 1.0e9
gamma = 5.0 / 3.0

# --- Constants from NEMO's source ------------------------------------------
KB = float(
    re.search(r'kB\(\)\s*\{\s*return\s*([0-9.eE+-]+)',
              (nemo_src / 'source/nemo-constants.h').read_text()).group(1))
apd = (nemo_src / 'source/nemo/atomic_physics_database.cpp').read_text()
Z, mass = {}, {}
for el in elements:
    m = re.search(
        r'species_attributes\["%s"\]\s*=\s*species_datatype\{\s*"%s",'
        r'\s*(\d+),\s*0,\s*([0-9.eE+-]+)' % (el, el), apd)
    Z[el], mass[el] = int(m.group(1)), float(m.group(2))

# --- Lookup-grid size from nemo.cpp -----------------------------------------
nemo_cpp = (nemo_src / 'source/nemo/nemo.cpp').read_text()
n_T = int(re.search(r'tab_Num_temps\((\d+)\)', nemo_cpp).group(1))
n_ne = int(re.search(r'tab_Num_ne\((\d+)\)', nemo_cpp).group(1))
T_min = float(re.search(r'tab_T_min\(([0-9.eE+-]+)\)', nemo_cpp).group(1))
T_max = float(re.search(r'tab_T_max\(([0-9.eE+-]+)\)', nemo_cpp).group(1))
m = re.search(r'tab_ne_max\(([^)]*\)?)\)', nemo_cpp).group(1)
ne_max = eval(m.replace('std::pow', 'pow'), {'pow': pow})
ne_min = float(re.search(r'tab_ne_min\(([0-9.eE+-]+)\)', nemo_cpp).group(1))

dlogT = (np.log10(T_max) - np.log10(T_min)) / (n_T - 1)
T_tab = 10.0**(np.log10(T_min) + np.arange(n_T) * dlogT)
dlogne = (np.log10(ne_max) - np.log10(ne_min)) / (n_ne - 1)
ne_tab = 10.0**(np.log10(ne_min) + np.arange(n_ne) * dlogne)

# --- CHIANTI tables (log10 Lambda on logT 1..9 x logne 0, 0.5, ...) -------
cc = (nemo_src / 'source/nemo/chianti_cooling.cpp').read_text()
names = {
    int(i): n
    for i, n in re.findall(r'species\[(\d+)\]\.Name\s*=\s*"([^"]+)"', cc)
}
tables = {}
for i, body in re.findall(r'species\[(\d+)\]\.rate\s*=\s*\{(.*?)\};', cc,
                          re.S):
    rows = re.findall(r'\{([^{}]*)\}', body)
    tables[names[int(i)]] = np.array(
        [[float(x) for x in r.split(',') if x.strip()] for r in rows])
n_cols = next(iter(tables.values())).shape[1]
logT_data = np.round(1.0 + 0.1 * np.arange(81), 1)
logne_data = np.round(0.5 * np.arange(n_cols), 1)

# --- Number densities --------------------------------------------------------
n_tot, ne = 0.0, 0.0
n_el = {}
for el in elements:
    n_el[el] = rho * mass_fractions[el] / mass[el]
    n_tot += n_el[el] * (1.0 + Z[el] / 2.0)  # equal fractions: Z/2 electrons
    ne += n_el[el] * Z[el] / 2.0

# --- Total loss rate at the T nodes, for this n_e ---------------------------
j = min(max(int((np.log10(ne) - np.log10(ne_min)) / dlogne), 0), n_ne - 2)
fn = (ne - ne_tab[j]) / (ne_tab[j + 1] - ne_tab[j])
logT_nodes = np.clip(np.log10(T_tab), logT_data[0], logT_data[-1])


def lut(name, jj):
    """NEMO's lookup table for one ion at n_e node jj, at all T nodes."""
    f = RegularGridInterpolator((logT_data, logne_data), tables[name])
    ln = np.clip(np.log10(ne_tab[jj]), logne_data[0], logne_data[-1])
    return 10.0**f(np.column_stack((logT_nodes, np.full(n_T, ln))))


L_nodes = np.zeros(n_T)
for el in elements:
    for q in range(Z[el] + 1):
        ion = el if q == 0 else f'{el}{q}+'
        Lam = (1.0 - fn) * lut(ion, j) + fn * lut(ion, j + 1)
        L_nodes += ne * n_el[el] / (Z[el] + 1) * Lam

C = n_tot * KB / (gamma - 1.0)  # energy density per unit temperature


def L_at(T):
    k = min(max(int((np.log10(T) - np.log10(T_min)) / dlogT), 0), n_T - 2)
    f = (T - T_tab[k]) / (T_tab[k + 1] - T_tab[k])
    return (1.0 - f) * L_nodes[k] + f * L_nodes[k + 1]


def cool_time(T_lo, T_hi):
    """Time to cool from T_hi to T_lo, both in the same grid interval."""
    L1, L2 = L_at(T_lo), L_at(T_hi)
    b = (L2 - L1) / (T_hi - T_lo)
    if abs(b) * (T_hi - T_lo) < 1e-10 * L1:
        return C * (T_hi - T_lo) * 2.0 / (L1 + L2)
    return C * np.log(L2 / L1) / b


# Cooling time from T_start down to each grid node below it.
nodes_below = T_tab[T_tab < T_start][::-1]
t_nodes = np.cumsum([
    cool_time(lo, hi)
    for lo, hi in zip(nodes_below, np.r_[T_start, nodes_below[:-1]])
])


def time_to(T):
    """Time to cool from T_start to T."""
    i = np.searchsorted(-nodes_below, -T)  # first node <= T
    T_above = T_start if i == 0 else nodes_below[i - 1]
    t_above = 0.0 if i == 0 else t_nodes[i - 1]
    return t_above + (cool_time(T, T_above) if T < T_above else 0.0)


temps = np.geomspace(T_start, 1.0e4, 50)
times = np.array([time_to(T) for T in temps])
dTdt = np.array([-L_at(T) / C for T in temps])

header = f'''Isochoric radiative cooling: reference data for NEMO's isochoric-cooling test.
Written by tests/isochoric-cooling/make_reference.py from the CHIANTI cooling
tables in source/nemo/chianti_cooling.cpp, using NEMO's interpolation scheme
({n_T} x {n_ne} lookup grid, n_e = {ne:.6f} cm^-3),
with exact integration.

Column 1: elapsed time t [s] at which the gas reaches temperature T.
Column 2: temperature T [K].
Column 3: cooling rate dT/dt [K/s] at temperature T.

        elapsed_time [s]          temperature [K]              dT/dt [K/s]'''
np.savetxt(out_file,
           np.column_stack((times, temps, dTdt)),
           fmt='%24.17e',
           header=header)
print(f'n_e = {ne:.12g}, n_tot = {n_tot:.12g}; wrote {out_file}')
