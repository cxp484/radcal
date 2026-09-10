#!/usr/bin/env python3
"""Plot FDS-style lookup kappa against direct RadCal kappa.

Input is fds_lookup_scatter.csv, written by generate_pl10cm_lookup.f90. The
corresponding temperature-interpolated scatter file is read automatically when
it is beside the FDS-floor scatter file. Two separate linear-scale plots are
saved, one for each lookup-temperature treatment.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("scatter_csv", type=Path, help="fds_lookup_scatter.csv")
    parser.add_argument(
        "--temperature-interpolated-csv",
        type=Path,
        help="fds_lookup_scatter_temperature_interpolated.csv (default: beside scatter_csv)",
    )
    parser.add_argument("--output-directory", type=Path, help="PNG directory (default: beside input)")
    args = parser.parse_args()

    output_directory = args.output_directory or args.scatter_csv.parent
    output_directory.mkdir(parents=True, exist_ok=True)
    interpolated_csv = args.temperature_interpolated_csv or args.scatter_csv.with_name(
        "fds_lookup_scatter_temperature_interpolated.csv"
    )
    direct, fds_lookup, temperature = read_scatter(args.scatter_csv)
    interpolated_direct, interpolated_lookup, interpolated_temperature = read_scatter(interpolated_csv)
    if len(direct) != len(interpolated_direct) or not np.allclose(direct, interpolated_direct):
        raise SystemExit("The two scatter files do not contain the same direct RadCal samples.")
    if not np.allclose(temperature, interpolated_temperature):
        raise SystemExit("The two scatter files do not contain the same temperature samples.")

    fds_valid = np.isfinite(direct) & np.isfinite(fds_lookup) & np.isfinite(temperature)
    interpolated_valid = np.isfinite(direct) & np.isfinite(interpolated_lookup) & np.isfinite(temperature)

    plot_comparison(
        direct[fds_valid],
        fds_lookup[fds_valid],
        temperature[fds_valid],
        output_directory / "fds_lookup_vs_direct_fds_linear.png",
        "FDS Temperature-Bin Lookup vs. Direct RadCal (linear scale)",
    )
    plot_comparison(
        direct[interpolated_valid],
        interpolated_lookup[interpolated_valid],
        temperature[interpolated_valid],
        output_directory / "fds_lookup_vs_direct_temperature_interpolated_linear.png",
        "Temperature-Interpolated Lookup vs. Direct RadCal (linear scale)",
    )


def read_scatter(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Read direct K10, lookup K10, and temperature from one sampled scatter file."""
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    return (
        np.array([float(row["direct_cm-1"]) for row in rows]),
        np.array([float(row["fds_lookup_cm-1"]) for row in rows]),
        np.array([float(row["temperature_k"]) for row in rows]),
    )


def plot_comparison(
    direct: np.ndarray,
    lookup: np.ndarray,
    temperature: np.ndarray,
    output: Path,
    title: str,
) -> None:
    """Save one direct-versus-lookup scatter plot."""
    upper = max(direct.max(), lookup.max()) * 1.02

    figure, axis = plt.subplots(figsize=(7.2, 6.2), constrained_layout=True)
    points = axis.scatter(direct, lookup, c=temperature, s=5, alpha=0.35, linewidths=0, cmap="viridis")
    axis.plot([0.0, upper], [0.0, upper], color="black", linewidth=1.2, label="exact agreement")
    axis.set_xlim(0.0, upper)
    axis.set_ylim(0.0, upper)
    axis.set_title(title)
    axis.set_aspect("equal", adjustable="box")
    axis.set_xlabel(r"Direct RadCal $\kappa_{10\,cm}$ (cm$^{-1}$)")
    axis.set_ylabel(r"FDS-style lookup $\kappa_{10\,cm}$ (cm$^{-1}$)")
    axis.grid(True, which="both", alpha=0.22)
    axis.legend(loc="upper left")
    colorbar = figure.colorbar(points, ax=axis)
    colorbar.set_label("Temperature (K)")
    figure.savefig(output, dpi=220)
    plt.close(figure)
    print(output)


if __name__ == "__main__":
    main()
