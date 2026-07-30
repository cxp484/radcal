import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
from read_nn_database import RadcalNNDataset

plt.rcParams.update({
    "font.size": 20,
    "axes.labelsize": 20,
    "legend.fontsize": 14,
    "xtick.labelsize": 20,
    "ytick.labelsize": 20,
})

COLORS = {
    "CO2": "tab:blue",
    "H2O": "tab:orange",
    "CO": "tab:green",
    "C2H4": "tab:purple",
    "SOOT": "tab:red",
}


def nearest_gas(data, target, chunk_size=1_000_000):
    best_index, best_error = 0, np.inf
    for start in range(0, data.n_gas, chunk_size):
        stop = min(start + chunk_size, data.n_gas)
        gas = data.gas_fraction[data.gas_indices[start:stop]]
        error = np.sum((gas - target) ** 2, axis=1)
        local = int(np.argmin(error))
        if error[local] < best_error:
            best_index, best_error = start + local, float(error[local])
    return best_index


def read_radcal_curve(species, rtmp_min, rtmp_max, target_fraction):
    values = np.loadtxt(HERE / f"kappa_{species}_radcal_pm.dat", comments="#")
    k_index, fraction, kappa = values[:, 1], values[:, 3], values[:, 4]
    selected_fraction = np.max(fraction)
    mask = np.isclose(fraction, selected_fraction)
    temperature = rtmp_min + k_index[mask] * (rtmp_max - rtmp_min) / np.max(k_index)
    kappa = kappa[mask] * target_fraction / selected_fraction
    order = np.argsort(temperature)
    return temperature[order], kappa[order], selected_fraction


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("header")
    parser.add_argument("-o", "--output", default="compare_planck_mean.png")
    parser.add_argument("--soot-fv", type=float, default=1e-5)
    parser.add_argument("--reference-tmin", type=float, default=270.0)
    parser.add_argument("--reference-tmax", type=float, default=2470.0)
    args = parser.parse_args()

    data = RadcalNNDataset(args.header)
    targets = {
        "CO2": np.array([1, 0, 0, 0], dtype=np.float32),
        "H2O": np.array([0, 1, 0, 0], dtype=np.float32),
        "CO": np.array([0, 0, 1, 0], dtype=np.float32),
        "C2H4": np.array([0, 0, 0, 1], dtype=np.float32),
    }
    zero_gas = nearest_gas(data, np.zeros(4, dtype=np.float32))
    zero_soot = int(np.argmin(np.abs(data.soot_fraction)))
    soot_index = int(np.argmin(np.abs(data.soot_fraction - args.soot_fv)))

    plt.figure(figsize=(10, 7))
    for species, target in targets.items():
        gas_index = nearest_gas(data, target)
        gas = data.gas_fraction[data.gas_indices[gas_index]]
        fraction = float(gas[np.argmax(target)])
        t_ref, k_ref, source_fraction = read_radcal_curve(
            species, args.reference_tmin, args.reference_tmax, fraction
        )
        k_db = np.asarray(data.kappa[:, gas_index, zero_soot, 0]) * 100.0
        plt.plot(t_ref, k_ref, "--", linewidth=2.5,
                 color=COLORS[species], label=f"{species} FDS")
        plt.plot(data.temperature, k_db, "s-", markevery=5, linewidth=2,
                 color=COLORS[species], label=f"{species} database")
        print(f"{species}: database x={fraction:.7g}, file x={source_fraction:.7g}")

    soot_fraction = float(data.soot_fraction[soot_index])
    t_ref, k_ref, source_fraction = read_radcal_curve(
        "SOOT", args.reference_tmin, args.reference_tmax, soot_fraction
    )
    k_db = np.asarray(data.kappa[:, zero_gas, soot_index, 0]) * 100.0
    plt.plot(t_ref, k_ref, "--", linewidth=2.5,
             color=COLORS["SOOT"], label="SOOT FDS")
    plt.plot(data.temperature, k_db, "s-", markevery=5, linewidth=2,
             color=COLORS["SOOT"], label="SOOT database")
    print(f"SOOT: database fv={soot_fraction:.7g}, file fv={source_fraction:.7g}")

    plt.xlabel("Temperature (K)")
    plt.ylabel(r"Planck-mean $\kappa$ (m$^{-1}$)")
    plt.xlim(300, 2300)
    plt.grid(True, alpha=0.3)
    plt.legend(ncol=2, frameon=False)
    plt.tight_layout()
    plt.savefig(args.output, dpi=200)


if __name__ == "__main__":
    main()
