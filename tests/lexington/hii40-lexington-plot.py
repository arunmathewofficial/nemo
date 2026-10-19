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
        nargs="?",
        type=Path,
        default=Path(__file__).resolve().with_name("temperature_profile.txt"),
        help="Two-column profile: radius in pc and temperature in K",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("temperature_profile.png"),
        help="Output figure path (default: temperature_profile.png)",
    )
    args = parser.parse_args()

    data = np.loadtxt(args.input, comments="#", ndmin=2)
    radius, temperature = data.T

    fig, ax = plt.subplots(figsize=(7, 4.5), constrained_layout=True)
    ax.plot(radius,
            temperature / 1.0e3,
            "o-",
            color="#176b9b",
            markersize=4,
            linewidth=1.5)
    ax.set_xlabel("Radius (pc)")
    ax.set_ylabel(r"$T\ (10^3\,\mathrm{K})$")
    ax.grid(True, alpha=0.25)
    ax.ticklabel_format(axis="y", style="plain", useOffset=False)
    fig.savefig(args.output, dpi=200)
    plt.close(fig)
    print(f"Saved {args.output.resolve()}")


if __name__ == "__main__":
    main()
