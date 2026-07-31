import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from read_nn_database import RadcalNNDataset

HERE = Path(__file__).resolve().parent
COLORS = {
    "CO2": "tab:blue",
    "H2O": "tab:orange",
    "CO": "tab:green",
    "C2H4": "tab:purple",
    "SOOT": "tab:red",
}

plt.rcParams.update({
    "font.size": 20,
    "axes.labelsize": 20,
    "legend.fontsize": 14,
    "xtick.labelsize": 18,
    "ytick.labelsize": 18,
})


def read_file_curve(species, target_temperature, tmin, tmax, divisions):
    values = np.loadtxt(HERE / f"kappa_{species}_radcal_pl10cm.dat", comments="#")
    k_index = values[:, 1].astype(int)
    temperature = tmin + k_index * (tmax - tmin) / divisions
    unique_k = np.unique(k_index)
    selected_k = unique_k[np.argmin(np.abs(
        tmin + unique_k * (tmax - tmin) / divisions - target_temperature
    ))]
    mask = (k_index == selected_k) & (values[:, 3] > 0)
    fraction, kappa = values[mask, 3], values[mask, 4]
    ratio = kappa / fraction
    if species == "SOOT":
        ratio *= 1e-5
    order = np.argsort(fraction)
    return fraction[order], ratio[order], temperature[mask][0]


def single_species_indices(data, species_column, chunk_size=1_000_000):
    result = np.full(len(data.gas_fraction), -1, dtype=np.int64)
    for start in range(0, data.n_gas, chunk_size):
        stop = min(start + chunk_size, data.n_gas)
        indices = np.asarray(data.gas_indices[start:stop])
        mask = np.all(np.delete(indices, species_column, axis=1) == 0, axis=1)
        rows = np.flatnonzero(mask)
        result[indices[rows, species_column]] = start + rows
    return result


def zero_gas_index(data, chunk_size=1_000_000):
    for start in range(0, data.n_gas, chunk_size):
        indices = np.asarray(data.gas_indices[start:min(start + chunk_size, data.n_gas)])
        rows = np.flatnonzero(np.all(indices == 0, axis=1))
        if len(rows):
            return start + int(rows[0])
    raise ValueError("The database has no zero-gas composition")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("header")
    parser.add_argument("-o", "--output", default="compare_kappa_10cm_1500K.png")
    parser.add_argument("--temperature", type=float, default=1500.0)
    parser.add_argument("--reference-tmin", type=float, default=270.0)
    parser.add_argument("--reference-tmax", type=float, default=2470.0)
    parser.add_argument("--reference-divisions", type=int, default=50)
    args = parser.parse_args()

    data = RadcalNNDataset(args.header)
    temperature_index = int(np.argmin(np.abs(data.temperature - args.temperature)))
    database_temperature = float(data.temperature[temperature_index])
    zero_soot = int(np.argmin(np.abs(data.soot_fraction)))
    zero_gas = zero_gas_index(data)

    plt.figure(figsize=(9, 7))
    for column, species in enumerate(("CO2", "H2O", "CO", "C2H4")):
        gas_records = single_species_indices(data, column)
        valid = gas_records >= 0
        fraction = data.gas_fraction[valid]
        kappa = np.asarray(data.kappa[temperature_index, gas_records[valid], zero_soot, 1]) * 100.0
        positive = fraction > 0
        x_file, y_file, file_temperature = read_file_curve(
            species, args.temperature, args.reference_tmin, args.reference_tmax,
            args.reference_divisions
        )
        plt.plot(x_file, y_file, "--", linewidth=2.5, color=COLORS[species],
                 label=f"{species} file ({file_temperature:.0f} K)")
        plt.plot(fraction[positive], kappa[positive] / fraction[positive], "s-",
                 markevery=5, linewidth=2, color=COLORS[species],
                 label=f"{species} database ({database_temperature:.0f} K)")

    fraction = data.soot_fraction
    kappa = np.asarray(data.kappa[temperature_index, zero_gas, :, 1]) * 100.0
    positive = fraction > 0
    x_file, y_file, file_temperature = read_file_curve(
        "SOOT", args.temperature, args.reference_tmin, args.reference_tmax,
        args.reference_divisions
    )
    plt.plot(x_file, y_file, "--", linewidth=2.5, color=COLORS["SOOT"],
             label=f"SOOT file ({file_temperature:.0f} K)")
    plt.plot(fraction[positive], kappa[positive] / fraction[positive] * 1e-5, "s-",
             markevery=5, linewidth=2, color=COLORS["SOOT"],
             label=f"SOOT database ({database_temperature:.0f} K)")

    plt.xscale("log")
    plt.xlabel("Mole fraction / soot volume fraction")
    plt.ylabel(r"Effective $\kappa$ (m$^{-1}$)")
    plt.grid(True, which="both", alpha=0.3)
    plt.legend(ncol=2, frameon=False)
    plt.tight_layout()
    plt.savefig(args.output, dpi=200)
    print(f"Database temperature: {database_temperature:.6g} K")


if __name__ == "__main__":
    main()
