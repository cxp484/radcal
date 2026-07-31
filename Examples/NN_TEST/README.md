# NN database test

`RADCAL_TRAINING_DATA_GENERATION/RADCAL_NN_TEST.IN` contains the full
44-temperature configuration.
`RADCAL_NN_PROGRESS_TEST.IN` is the retained reduced fixture used for fast
serial and MPI validation; it produces 20 valid gas compositions and 120
Kappa records.

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
cd RADCAL_TRAINING_DATA_GENERATION
python3 generate_nn_database.py \
  ../../../Build/intel_osx_64/radcal_nn \
  --input RADCAL_NN_TEST.IN \
  --output-prefix radcal_nn_test \
  --processes 4
```

On Linux, replace `intel_osx_64` with `intel_linux_64`.
Temperatures are assigned to MPI ranks. If `--processes` is greater than the
`N_RTMP` value in the input namelist, it is automatically reduced to
`N_RTMP`.
Add `--oversubscribe` when requesting more processes than the local machine
provides as CPU slots.

The generated files are:

- `radcal_nn_test_header.txt`
- `radcal_nn_test_valid_gas_indices.bin`
- `radcal_nn_test_kappa.bin`

Return to the NN test directory:

```bash
cd ..
```

Read a training batch with:

```bash
python3 VERIFICATION_TRAINING_DATA_WITH_FDS/read_nn_database.py \
  RADCAL_TRAINING_DATA_GENERATION/radcal_nn_test_header.txt \
  --start 0 --count 8
```

Plot nearest-available pure-species Planck mean curves with:

```bash
python3 VERIFICATION_TRAINING_DATA_WITH_FDS/compare_planck_mean.py \
  RADCAL_TRAINING_DATA_GENERATION/radcal_nn_test_header.txt \
  --output planck_mean_kappa.png
```

Compare 10 cm path-length Kappa at 1500 K with:

```bash
python3 VERIFICATION_TRAINING_DATA_WITH_FDS/compare_kappa_10cm.py \
  RADCAL_TRAINING_DATA_GENERATION/radcal_nn_test_header.txt \
  --output VERIFICATION_TRAINING_DATA_WITH_FDS/compare_kappa_10cm_1500K.png
```
