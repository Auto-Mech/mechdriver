#!/usr/bin/env bash

# Note that install location assumes one is in a Conda environment

# Chorus was here! Test Test!

mkdir -p build
cd build

cmake -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX -DCMAKE_Fortran_COMPILER=$FC -DCMAKE_Fortran_FLAGS="${FFLAGS}" ..

make VERBOSE=1
make install
