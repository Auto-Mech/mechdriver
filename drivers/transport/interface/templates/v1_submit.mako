#!/bin/bash
#SBATCH -p bdwall
#SBATCH -N 1
#SBATCH --cpus-per-task=${nprocs}
#SBATCH -t 12:00:00
#SBATCH -J 1dmin
#SBATCH -A CMRP
#SBATCH -o job_%j.log
#SBATCH -e job_%j.err

# Set the libraries
source ~/.bash_profile
## module load molpro/2015.1_170920

# Set library path
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/ajasper/lib/

# Set variables for running Molpro
export TMPDIR=${scratch}
## export MOLPRO_LIB=/soft/molpro/2015.1_170920/bebop/molprop_2015_1_linux_x86_64_i8/lib/
## export MOLPRO_OPTIONS="--nouse-logfile --no-xml-output -L $MOLPRO_LIB -d $TMPDIR -I $TMPDIR -W $TMPDIR -o qc.out -s qc.in"

# Run several auto1dmin.x instances
cd run1
time ./auto1dmin.x < input.dat > output.dat &
% for i in range(nprocs-1):
cd ../run${i+2}
## echo "molpro.exe -c 1 -n 1 --exclusive $MOLPRO_OPTIONS molpro.exe" >> elstruct.x 
time ./auto1dmin.x < input.dat > output.dat &
% endfor
