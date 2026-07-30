# NN database test

This reduced grid produces 20 valid gas compositions and 120 Kappa records.

On macOS:

```bash
cd Build/intel_osx_64
source ../Scripts/set_compilers.sh
make -f ../makefile clean
make -f ../makefile intel_osx_64_nn VPATH=../../Source
cd ../..
```

On Linux:

```bash
cd Build/intel_linux_64
source ../Scripts/set_compilers.sh
make -f ../makefile clean
make -f ../makefile intel_linux_64_nn VPATH=../../Source
cd ../..
```

Generate the test database:

```bash
python3 Examples/NN_TEST/generate_nn_database.py \
  Build/intel_osx_64/radcal_nn \
  --input Examples/NN_TEST/RADCAL_NN_TEST.IN \
  --output-prefix Examples/NN_TEST/radcal_nn_test
```

The generated files are:

- `radcal_nn_test_header.txt`
- `radcal_nn_test_valid_gas_indices.bin`
- `radcal_nn_test_kappa.bin`

Read a training batch with:

```bash
python3 Examples/NN_TEST/read_nn_database.py \
  Examples/NN_TEST/radcal_nn_test_header.txt \
  --start 0 --count 8
```

Plot nearest-available pure-species Planck mean curves with:

```bash
python3 Examples/NN_TEST/plot_planck_mean.py \
  Examples/NN_TEST/radcal_nn_test_header.txt \
  --output Examples/NN_TEST/planck_mean_kappa.png
```
