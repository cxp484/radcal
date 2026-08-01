import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from read_nn_database import RadcalNNDataset

HERE = Path(__file__).resolve().parent
SPECIES = ("CO2", "H2O", "CO", "C2H4")
TARGET_T = 1527.0
TARGET_CO2 = 8.1854673070690193e-2
TARGET_H2O = 4.9619476030029037e-2
RTOL = 1e-6
ATOL = 1e-8


def gas_lookup(data):
    return {tuple(row): i for i, row in enumerate(data.gas_indices)}


def check(name, x, direct, database, xlabel, ylabel, plot_directory):
    direct = np.asarray(direct, dtype=np.float32)
    database = np.asarray(database, dtype=np.float32)
    filename = name.lower().replace(" ", "_").replace(",", "").replace("/", "_")
    plt.figure(figsize=(8, 6))
    plt.plot(x, direct, "-", linewidth=2.5, label="Direct RADCAL")
    plt.plot(x, database, "o", markersize=7, markerfacecolor="none",
             markeredgewidth=1.5, label="Database")
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(name)
    plt.grid(True, alpha=0.3)
    plt.legend(frameon=False)
    plt.tight_layout()
    plt.savefig(plot_directory / f"{filename}.png", dpi=200)
    plt.close()
    equal = np.allclose(direct, database, rtol=RTOL, atol=ATOL)
    error = 0.0 if equal else float(np.max(np.abs(direct - database)))
    return name, equal, error, len(direct)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("header")
    parser.add_argument("--direct-directory", type=Path, default=HERE)
    parser.add_argument("--plot-directory", type=Path,
                        default=HERE / "direct_comparison_plots")
    args = parser.parse_args()
    args.plot_directory.mkdir(parents=True, exist_ok=True)

    data = RadcalNNDataset(args.header)
    results = []
    lookup = gas_lookup(data)
    n_axis = len(data.gas_fraction)
    zero = (0, 0, 0, 0)
    zero_gas = lookup[zero]
    zero_soot = int(np.argmin(np.abs(data.soot_fraction)))
    temperature_index = int(np.argmin(np.abs(data.temperature - TARGET_T)))
    co2_index = int(np.argmin(np.abs(data.gas_fraction - TARGET_CO2)))
    h2o_index = int(np.argmin(np.abs(data.gas_fraction - TARGET_H2O)))

    for column, species in enumerate(SPECIES):
        direct = np.loadtxt(args.direct_directory / f"direct_{species}_planck.dat")
        indices = [0, 0, 0, 0]
        indices[column] = n_axis - 1
        record = lookup[tuple(indices)]
        results.append(check(
            f"{species} Planck", direct[:, 1], direct[:, 3],
            data.kappa[:, record, zero_soot, 0], "Temperature (K)",
            r"Planck-mean $\kappa$ (cm$^{-1}$)", args.plot_directory
        ))

    direct = np.loadtxt(args.direct_directory / "direct_SOOT_planck.dat")
    soot_index = int(np.argmin(np.abs(data.soot_fraction - direct[0, 2])))
    results.append(check(
        "SOOT Planck", direct[:, 1], direct[:, 3],
        data.kappa[:, zero_gas, soot_index, 0], "Temperature (K)",
        r"Planck-mean $\kappa$ (cm$^{-1}$)", args.plot_directory
    ))

    for column, species in enumerate(SPECIES):
        direct = np.loadtxt(args.direct_directory / f"direct_{species}_kappa10cm.dat")
        records = []
        for axis_index in direct[:, 0].astype(int):
            indices = [0, 0, 0, 0]
            indices[column] = axis_index
            records.append(lookup[tuple(indices)])
        results.append(check(
            f"{species} 10 cm", direct[:, 1], direct[:, 2],
            data.kappa[temperature_index, records, zero_soot, 1],
            "Mole fraction", r"10 cm $\kappa$ (cm$^{-1}$)",
            args.plot_directory
        ))

    direct = np.loadtxt(args.direct_directory / "direct_SOOT_kappa10cm.dat")
    results.append(check(
        "SOOT 10 cm", direct[:, 1], direct[:, 2],
        data.kappa[temperature_index, zero_gas, direct[:, 0].astype(int), 1],
        "Soot volume fraction", r"10 cm $\kappa$ (cm$^{-1}$)",
        args.plot_directory
    ))

    direct = np.loadtxt(args.direct_directory / "direct_CO2_H2O_kappa10cm.dat")
    records = [
        lookup[(co2_index, h2o_axis, 0, 0)]
        for h2o_axis in direct[:, 0].astype(int)
    ]
    results.append(check(
        "fixed CO2, varying H2O 10 cm", direct[:, 1], direct[:, 2],
        data.kappa[temperature_index, records, zero_soot, 1],
        "H2O mole fraction", r"10 cm $\kappa$ (cm$^{-1}$)",
        args.plot_directory
    ))

    direct = np.loadtxt(args.direct_directory / "direct_CO2_H2O_SOOT_kappa10cm.dat")
    record = lookup[(co2_index, h2o_index, 0, 0)]
    results.append(check(
        "fixed CO2/H2O, varying SOOT 10 cm", direct[:, 1], direct[:, 2],
        data.kappa[temperature_index, record, direct[:, 0].astype(int), 1],
        "Soot volume fraction", r"10 cm $\kappa$ (cm$^{-1}$)",
        args.plot_directory
    ))

    failures = []
    for name, equal, error, count in results:
        if equal:
            print(f"PASS {name}: {count} values match within tolerance")
        else:
            failures.append(f"{name}: maximum error={error:.8e}")
            print(f"FAIL {failures[-1]}")
    if failures:
        raise AssertionError(
            f"{len(failures)} of {len(results)} comparisons failed:\n" +
            "\n".join(failures)
        )


if __name__ == "__main__":
    main()
