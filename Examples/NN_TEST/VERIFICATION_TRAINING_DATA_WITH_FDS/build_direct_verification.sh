#!/bin/sh
set -eu

cd "$(dirname "$0")"
FC_CMD=mpiifx
if [ "$#" -gt 0 ]; then FC_CMD=$1; fi
"$FC_CMD" -O2 -cpp -o generate_direct_verification \
  ../../../Source/rmod.f90 ../../../Source/rcal.f90 \
  generate_direct_verification.f90
