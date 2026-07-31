import argparse
import sys
from pathlib import Path

import numpy as np
import tensorflow as tf

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent.parent / "VERIFICATION_TRAINING_DATA_WITH_FDS"))
from read_nn_database import RadcalNNDataset


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


def create_model(learning_rate):
    model = tf.keras.Sequential([
        tf.keras.layers.Input(shape=(6,)),
        tf.keras.layers.Dense(64, activation="relu"),
        tf.keras.layers.Dense(64, activation="relu"),
        tf.keras.layers.Dense(1),
    ])
    model.compile(
        optimizer=tf.keras.optimizers.Adam(learning_rate),
        loss="mean_squared_error",
    )
    return model


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("header")
    parser.add_argument("--epochs", type=int, default=10)
    parser.add_argument("--batch-size", type=int, default=65536)
    parser.add_argument("--learning-rate", type=float, default=1e-3)
    parser.add_argument("--output", default="kappa_10cm_model.keras")
    parser.add_argument("--max-batches", type=int)
    args = parser.parse_args()

    data = RadcalNNDataset(args.header)
    minimum, maximum = scales(data)
    model = create_model(args.learning_rate)

    for epoch in range(args.epochs):
        total_loss = 0.0
        records = 0
        for batch_number, (inputs, targets) in enumerate(data.batches(args.batch_size)):
            if args.max_batches is not None and batch_number >= args.max_batches:
                break
            x = normalize(inputs, data, minimum, maximum)
            y = targets[:, 1:2].copy()
            order = np.random.permutation(len(x))
            loss = model.train_on_batch(x[order], y[order])
            total_loss += float(loss) * len(x)
            records += len(x)
        print(f"epoch {epoch + 1}/{args.epochs} loss={total_loss / records:.6e}")

    model.save(args.output)
    np.savez(
        Path(args.output).with_suffix(".scales.npz"),
        input_names=np.array(("temperature_K", "x_co2", "x_h2o", "x_co", "soot_fv", "x_c2h4")),
        temperature_min=data.temperature[0],
        temperature_max=data.temperature[-1],
        species_min=minimum,
        species_max=maximum,
        target_name="kappa_10cm_cm-1",
    )


if __name__ == "__main__":
    main()
