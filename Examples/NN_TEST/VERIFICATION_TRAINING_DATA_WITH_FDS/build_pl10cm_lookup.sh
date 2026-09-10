#!/bin/sh
set -eu

cd "$(dirname "$0")"
FC_CMD=${1:-mpiifx}
"$FC_CMD" -O2 -o generate_pl10cm_lookup \
  ../../../Source/rmod.f90 ../../../Source/rcal.f90 \
  generate_pl10cm_lookup.f90
