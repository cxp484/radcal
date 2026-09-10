#!/usr/bin/env python3
"""Identify the temperature, soot, and gas mixtures behind lookup-table errors.

Run this after generate_pl10cm_lookup.f90 has written an enriched
fds_lookup_scatter.csv.  It saves the largest-relative-error sampled states to
fds_lookup_worst_cases.csv and plots relative error against direct K10, colored
by soot volume fraction.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


STATE_COLUMNS = (
    "temperature_k",
    "x_co2",
    "x_h2o",
    "x_co",
    "x_c2h4",
    "soot_volume_fraction",
    "direct_cm-1",
    "fds_lookup_cm-1",
    "abs_percent_error",
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("scatter_csv", type=Path)
    parser.add_argument("--top", type=int, default=100, help="number of largest-error cases to export")
    parser.add_argument("--output-directory", type=Path)
    args = parser.parse_args()

    out_dir = args.output_directory or args.scatter_csv.parent
    out_dir.mkdir(parents=True, exist_ok=True)
    with args.scatter_csv.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    missing = set(STATE_COLUMNS) - set(rows[0] if rows else ())
    if missing:
        raise SystemExit(
            "The scatter CSV lacks state columns. Re-run generate_pl10cm_lookup.f90 first. "
            f"Missing: {', '.join(sorted(missing))}"
        )

    values = {column: np.array([float(row[column]) for row in rows]) for column in STATE_COLUMNS}
    valid = np.isfinite(values["abs_percent_error"]) & (values["direct_cm-1"] > 1.0e-12)
    order = np.argsort(values["abs_percent_error"][valid])[::-1]
    selected = np.flatnonzero(valid)[order[: args.top]]

    worst_csv = out_dir / "fds_lookup_worst_cases.csv"
    with worst_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=STATE_COLUMNS)
        writer.writeheader()
        for index in selected:
            writer.writerow({column: values[column][index] for column in STATE_COLUMNS})

    print(f"Exported {len(selected)} worst sampled cases: {worst_csv}")
    print("Worst-case ranges:")
    for column in ("temperature_k", "soot_volume_fraction", "x_co2", "x_h2o", "x_co", "x_c2h4", "abs_percent_error"):
        subset = values[column][selected]
        print(f"  {column}: {subset.min():.6e} to {subset.max():.6e}; median {np.median(subset):.6e}")

    print("\nTop 10 worst sampled points by absolute percent error:")
    print(
        " rank  temperature_K          x_CO2          x_H2O           x_CO"
        "         x_C2H4           soot         direct         lookup       error_%"
    )
    for rank, index in enumerate(selected[:10], start=1):
        print(
            f" {rank:>4d}  {values['temperature_k'][index]:13.6e}"
            f" {values['x_co2'][index]:13.6e} {values['x_h2o'][index]:13.6e}"
            f" {values['x_co'][index]:13.6e} {values['x_c2h4'][index]:13.6e}"
            f" {values['soot_volume_fraction'][index]:13.6e}"
            f" {values['direct_cm-1'][index]:13.6e} {values['fds_lookup_cm-1'][index]:13.6e}"
            f" {values['abs_percent_error'][index]:10.4f}"
        )

    figure, axis = plt.subplots(figsize=(7.4, 5.6), constrained_layout=True)
    direct = values["direct_cm-1"][valid]
    error = values["abs_percent_error"][valid]
    soot = values["soot_volume_fraction"][valid]
    positive_soot = soot > 0.0
    color = np.log10(np.maximum(soot, np.min(soot[positive_soot]) if np.any(positive_soot) else 1.0e-20))
    points = axis.scatter(direct, error, c=color, s=5, alpha=0.35, linewidths=0, cmap="plasma")
    axis.set_xscale("log")
    axis.set_yscale("log")
    axis.set_xlabel(r"Direct RadCal $\kappa_{10\,cm}$ (cm$^{-1}$)")
    axis.set_ylabel("Absolute percent error")
    axis.set_title("FDS Lookup Error by Soot Volume Fraction")
    axis.grid(True, which="both", alpha=0.22)
    colorbar = figure.colorbar(points, ax=axis)
    colorbar.set_label(r"$\log_{10}$(soot volume fraction)")
    figure.savefig(out_dir / "fds_lookup_error_by_soot.png", dpi=220)
    plt.close(figure)


if __name__ == "__main__":
    main()
