#!/usr/bin/env python3
"""Plot NEMO's final shell temperatures against radius."""

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "input",
        type=Path,
        default=Path(__file__).resolve().with_name("output_states_pieq.txt"),
        help="radius in cm [0] and temperature in K [1]",
    )
    parser.add_argument(
        "reference",
        type=Path,
        default=Path(__file__).resolve().with_name("pieq_test_refdata.txt"),
        help="radius in cm [0] and temperature in K [1]",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("temperature_profile.png"),
        help="Output figure path (default: temperature_profile.png)",
    )
    args = parser.parse_args()

    data = np.loadtxt(args.input, skiprows=1, ndmin=2)
    refs = np.loadtxt(args.reference, skiprows=1, ndmin=2)

    r1 = np.array(data[:, 0]) / 3.086e18
    r2 = np.array(refs[:, 0]) / 3.086e18
    t1 = np.array(data[:, 1])
    t2 = np.array(refs[:, 1])
    yh1 = np.array(data[:, 9]) / np.array(data[:, 2])
    yh2 = np.array(refs[:, 9]) / np.array(refs[:, 2])
    yhe1 = np.array(data[:, 10]) / np.array(data[:, 3])
    yhe2 = np.array(refs[:, 10]) / np.array(refs[:, 3])
    yhep1 = np.array(data[:, 11]) / np.array(data[:, 3])
    yhep2 = np.array(refs[:, 11]) / np.array(refs[:, 3])

    fig, ax = plt.subplots(figsize=(7, 4.5), constrained_layout=True)
    ax.plot(r1, t1 / 1.0e4, "o-", color="r", markersize=4, linewidth=1.5)
    ax.plot(r2, t2 / 1.0e4, "o-", color="b", markersize=4, linewidth=1.5)
    ax.plot(r1, yh1, "--", color="r", linewidth=1)
    ax.plot(r2, yh2, "--", color="b", linewidth=1)
    ax.plot(r1, yhe1, ":", color="r", linewidth=1)
    ax.plot(r2, yhe2, ":", color="b", linewidth=1)
    ax.plot(r1, yhep1, "-.", color="r", linewidth=1)
    ax.plot(r2, yhep2, "-.", color="b", linewidth=1)
    ax.set_xlabel("Radius (pc)")
    ax.set_ylabel(r"$T\ (10^3\,\mathrm{K})$")
    ax.grid(True, alpha=0.25)
    ax.ticklabel_format(axis="y", style="plain", useOffset=False)
    fig.savefig(args.output, dpi=200)
    plt.close(fig)
    print(f"Saved {args.output.resolve()}")


if __name__ == "__main__":
    main()
