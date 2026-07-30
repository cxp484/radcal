import argparse
from pathlib import Path

import numpy as np


class RadcalNNDataset:
    input_names = ("temperature_K", "pressure_atm", "x_co2", "x_h2o", "x_co", "soot_fv", "x_c2h4")
    target_names = ("kappa_planck_cm-1", "kappa_10cm_cm-1")

    def __init__(self, header):
        self.header = Path(header)
        values, axes = self._read_header()
        self.n_temperature = int(values["N_TEMPERATURE"])
        self.n_gas = int(values["N_VALID_GAS"])
        self.n_soot = int(values["N_SOOT_AXIS"])
        self.pressure = float(values["PRESSURE_ATM"])
        self.temperature = np.asarray(axes["TEMPERATURE_K"], dtype=np.float32)
        self.gas_fraction = np.asarray(axes["GAS_MOLE_FRACTION"], dtype=np.float32)
        self.soot_fraction = np.asarray(axes["SOOT_VOLUME_FRACTION"], dtype=np.float32)
        gas_file = self._resolve(values["GAS_INDEX_FILE"])
        kappa_file = self._resolve(values["KAPPA_FILE"])
        self.gas_indices = np.memmap(gas_file, dtype=np.int32, mode="r", shape=(self.n_gas, 4))
        self.kappa = np.memmap(
            kappa_file, dtype=np.float32, mode="r",
            shape=(self.n_temperature, self.n_gas, self.n_soot, 2)
        )
        self.size = self.n_temperature*self.n_gas*self.n_soot

    def _resolve(self, name):
        path = Path(name)
        return path if path.exists() else self.header.parent/path.name

    def _read_header(self):
        values = {}
        axes = {}
        section = None
        for raw in self.header.read_text().splitlines():
            line = raw.strip()
            if not line:
                continue
            if line.startswith("["):
                section = line[1:-1]
                axes[section] = []
            elif section:
                axes[section].append(float(line.split()[1]))
            elif "=" in line:
                key, value = line.split("=", 1)
                values[key] = value.strip()
        return values, axes

    def batch(self, start, count):
        stop = min(start+count, self.size)
        record = np.arange(start, stop, dtype=np.int64)
        soot_index = record % self.n_soot
        gas_index = (record//self.n_soot) % self.n_gas
        temperature_index = record//(self.n_soot*self.n_gas)
        gas = self.gas_fraction[self.gas_indices[gas_index]]
        inputs = np.empty((len(record), 7), dtype=np.float32)
        inputs[:, 0] = self.temperature[temperature_index]
        inputs[:, 1] = self.pressure
        inputs[:, 2:5] = gas[:, :3]
        inputs[:, 5] = self.soot_fraction[soot_index]
        inputs[:, 6] = gas[:, 3]
        targets = np.asarray(self.kappa).reshape(-1, 2)[start:stop]
        return inputs, targets

    def batches(self, batch_size):
        for start in range(0, self.size, batch_size):
            yield self.batch(start, batch_size)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("header")
    parser.add_argument("--start", type=int, default=0)
    parser.add_argument("--count", type=int, default=8)
    args = parser.parse_args()
    data = RadcalNNDataset(args.header)
    inputs, targets = data.batch(args.start, args.count)
    print("inputs:", data.input_names)
    print(inputs)
    print("targets:", data.target_names)
    print(targets)


if __name__ == "__main__":
    main()
