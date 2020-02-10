#!/usr/bin/env bash

# Set library path
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/ajasper/lib/

# Run several auto1dmin.x instances
cd ${run_path}/run1
time ${src_path} < input.dat > output.dat &
% for i in range(nprocs-1):
cd ../run${i+2}
time ${src_path} < input.dat > output.dat &
% endfor
wait
