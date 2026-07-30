#!/bin/bash
platform=intel64
dir=`pwd`
target=intel_osx_64_nn

source ../Scripts/set_compilers.sh

echo Building $target
make -j4 DEBUG=1 VPATH="../../Source" -f ../makefile $target