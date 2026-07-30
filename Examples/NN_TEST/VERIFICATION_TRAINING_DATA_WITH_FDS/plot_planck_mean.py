import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from read_nn_database import RadcalNNDataset


def nearest_gas(data, target, chunk_size=1_000_000):
    best_index, best_error = 0, np.inf
    for start in range(0, data.n_gas, chunk_size):
        stop = min(start+chunk_size, data.n_gas)
        gas = data.gas_fraction[data.gas_indices[start:stop]]
        error = np.sum((gas-target)**2, axis=1)
        local = int(np.argmin(error))
        if error[local] < best_error:
            best_index, best_error = start+local, float(error[local])
    return best_index


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("header")
    parser.add_argument("-o", "--output", default="planck_mean_kappa.png")
    args = parser.parse_args()
    data = RadcalNNDataset(args.header)
    output = Path(args.output)
    targets = {
        "CO2": np.array([1, 0, 0, 0], dtype=np.float32),
        "H2O": np.array([0, 1, 0, 0], dtype=np.float32),
        "CO": np.array([0, 0, 1, 0], dtype=np.float32),
        "C2H4": np.array([0, 0, 0, 1], dtype=np.float32),
    }
    soot_min = 0
    soot_target = int(np.argmin(np.abs(data.soot_fraction-1E-5)))
    zero_gas = nearest_gas(data, np.zeros(4, dtype=np.float32))
    curves = []

    for name, target in targets.items():
        gas_index = nearest_gas(data, target)
        gas = data.gas_fraction[data.gas_indices[gas_index]]
        kappa = np.asarray(data.kappa[:, gas_index, soot_min, 0])
        curves.append((name, gas, data.soot_fraction[soot_min], kappa))

    gas = data.gas_fraction[data.gas_indices[zero_gas]]
    curves.append(("Soot", gas, data.soot_fraction[soot_target],
                   np.asarray(data.kappa[:, zero_gas, soot_target, 0])))

    for name, gas, soot, kappa in curves:
        plt.plot(data.temperature, kappa, marker="o", label=name)
        print(f"{name}: CO2={gas[0]:.7g}, H2O={gas[1]:.7g}, CO={gas[2]:.7g}, "
              f"C2H4={gas[3]:.7g}, soot={soot:.7g}")

    plt.xlabel("Temperature (K)")
    plt.ylabel(r"Planck mean $\kappa$ (cm$^{-1}$)")
    plt.title("Nearest available pure-species compositions")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(output, dpi=200)
    plt.close()

    with output.with_suffix(".csv").open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["temperature_K", *[name for name, *_ in curves]])
        for i, temperature in enumerate(data.temperature):
            writer.writerow([temperature, *[curve[3][i] for curve in curves]])


if __name__ == "__main__":
    main()
