# NN database test

This reduced grid produces 20 valid gas compositions and 120 Kappa records.

From the repository root, first change to the test directory:

```bash
cd Examples/NN_TEST
```

On macOS:

```bash
../../Build/intel_osx_64/make_radcal_nn.sh
```

On Linux:

```bash
../../Build/intel_linux_64/make_radcal_nn.sh
```

Generate the test database on macOS:

```bash
python3 generate_nn_database.py \
  ../../Build/intel_osx_64/radcal_nn \
  --input RADCAL_NN_TEST.IN \
  --output-prefix radcal_nn_test
```

On Linux, replace `intel_osx_64` with `intel_linux_64`.

The generated files are:

- `radcal_nn_test_header.txt`
- `radcal_nn_test_valid_gas_indices.bin`
- `radcal_nn_test_kappa.bin`

Read a training batch with:

```bash
python3 read_nn_database.py radcal_nn_test_header.txt --start 0 --count 8
```

Plot nearest-available pure-species Planck mean curves with:

```bash
python3 plot_planck_mean.py \
  radcal_nn_test_header.txt \
  --output planck_mean_kappa.png
```
