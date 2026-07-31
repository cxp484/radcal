import argparse
from pathlib import Path

import numpy as np

from read_nn_database import RadcalNNDataset

HERE = Path(__file__).resolve().parent
SPECIES = ("CO2", "H2O", "CO", "C2H4")
TARGET_T = 1527.0
TARGET_CO2 = 8.1854673070690193e-2
TARGET_H2O = 4.9619476030029037e-2


def gas_lookup(data):
    return {tuple(row): i for i, row in enumerate(data.gas_indices)}


def check(name, direct, database):
    direct = np.asarray(direct, dtype=np.float32)
    database = np.asarray(database, dtype=np.float32)
    if not np.array_equal(direct, database):
        error = np.max(np.abs(direct - database))
        raise AssertionError(f"{name}: values differ; maximum error={error:.8e}")
    print(f"PASS {name}: {len(direct)} values match exactly")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("header")
    parser.add_argument("--direct-directory", type=Path, default=HERE)
    args = parser.parse_args()

    data = RadcalNNDataset(args.header)
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
        check(f"{species} Planck", direct[:, 3],
              data.kappa[:, record, zero_soot, 0])

    direct = np.loadtxt(args.direct_directory / "direct_SOOT_planck.dat")
    soot_index = int(np.argmin(np.abs(data.soot_fraction - direct[0, 2])))
    check("SOOT Planck", direct[:, 3],
          data.kappa[:, zero_gas, soot_index, 0])

    for column, species in enumerate(SPECIES):
        direct = np.loadtxt(args.direct_directory / f"direct_{species}_kappa10cm.dat")
        records = []
        for axis_index in direct[:, 0].astype(int):
            indices = [0, 0, 0, 0]
            indices[column] = axis_index
            records.append(lookup[tuple(indices)])
        check(f"{species} 10 cm", direct[:, 2],
              data.kappa[temperature_index, records, zero_soot, 1])

    direct = np.loadtxt(args.direct_directory / "direct_SOOT_kappa10cm.dat")
    check("SOOT 10 cm", direct[:, 2],
          data.kappa[temperature_index, zero_gas, direct[:, 0].astype(int), 1])

    direct = np.loadtxt(args.direct_directory / "direct_CO2_H2O_kappa10cm.dat")
    records = [
        lookup[(co2_index, h2o_axis, 0, 0)]
        for h2o_axis in direct[:, 0].astype(int)
    ]
    check("fixed CO2, varying H2O 10 cm", direct[:, 2],
          data.kappa[temperature_index, records, zero_soot, 1])

    direct = np.loadtxt(args.direct_directory / "direct_CO2_H2O_SOOT_kappa10cm.dat")
    record = lookup[(co2_index, h2o_index, 0, 0)]
    check("fixed CO2/H2O, varying SOOT 10 cm", direct[:, 2],
          data.kappa[temperature_index, record, direct[:, 0].astype(int), 1])


if __name__ == "__main__":
    main()
