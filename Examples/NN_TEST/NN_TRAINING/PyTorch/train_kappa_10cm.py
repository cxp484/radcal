import argparse
import sys
from pathlib import Path

import numpy as np
import torch
from torch import nn

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent.parent / "VERIFICATION_TRAINING_DATA_WITH_FDS"))
from read_nn_database import RadcalNNDataset


class KappaNet(nn.Module):
    def __init__(self):
        super().__init__()
        self.layers = nn.Sequential(
            nn.Linear(6, 64), nn.ReLU(),
            nn.Linear(64, 64), nn.ReLU(),
            nn.Linear(64, 1),
        )

    def forward(self, inputs):
        return self.layers(inputs)


def scales(data):
    gas_positive = data.gas_fraction[data.gas_fraction > 0]
    soot_positive = data.soot_fraction[data.soot_fraction > 0]
    minimum = np.array([
        gas_positive.min(), gas_positive.min(), gas_positive.min(),
        soot_positive.min(), gas_positive.min(),
    ], dtype=np.float32)
    maximum = np.array([
        data.gas_fraction.max(), data.gas_fraction.max(), data.gas_fraction.max(),
        data.soot_fraction.max(), data.gas_fraction.max(),
    ], dtype=np.float32)
    return minimum, maximum


def normalize(inputs, data, minimum, maximum):
    result = np.empty((len(inputs), 6), dtype=np.float32)
    span = data.temperature[-1] - data.temperature[0]
    result[:, 0] = (inputs[:, 0] - data.temperature[0]) / span
    species = inputs[:, [2, 3, 4, 5, 6]]
    result[:, 1:] = np.log1p(species / minimum) / np.log1p(maximum / minimum)
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("header")
    parser.add_argument("--epochs", type=int, default=10)
    parser.add_argument("--batch-size", type=int, default=65536)
    parser.add_argument("--learning-rate", type=float, default=1e-3)
    parser.add_argument("--output", default="kappa_10cm_model.pt")
    parser.add_argument("--max-batches", type=int)
    args = parser.parse_args()

    data = RadcalNNDataset(args.header)
    minimum, maximum = scales(data)
    device = torch.device("cuda" if torch.cuda.is_available() else
                          "mps" if torch.backends.mps.is_available() else "cpu")
    model = KappaNet().to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=args.learning_rate)
    loss_function = nn.MSELoss()

    for epoch in range(args.epochs):
        total_loss = 0.0
        records = 0
        for batch_number, (inputs, targets) in enumerate(data.batches(args.batch_size)):
            if args.max_batches is not None and batch_number >= args.max_batches:
                break
            x = torch.from_numpy(normalize(inputs, data, minimum, maximum)).to(device)
            y = torch.from_numpy(targets[:, 1:2].copy()).to(device)
            order = torch.randperm(len(x), device=device)
            prediction = model(x[order])
            loss = loss_function(prediction, y[order])
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
            total_loss += loss.item() * len(x)
            records += len(x)
        print(f"epoch {epoch + 1}/{args.epochs} loss={total_loss / records:.6e}")

    torch.save({
        "model_state": model.state_dict(),
        "input_names": ("temperature_K", "x_co2", "x_h2o", "x_co", "soot_fv", "x_c2h4"),
        "temperature_min": float(data.temperature[0]),
        "temperature_max": float(data.temperature[-1]),
        "species_min": minimum,
        "species_max": maximum,
        "target_name": "kappa_10cm_cm-1",
    }, args.output)


if __name__ == "__main__":
    main()
