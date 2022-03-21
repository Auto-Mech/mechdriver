# Mechdriver
This repository houses the main AutoMech driver (bin/automech.py) and also the high-level electronic structure and thermochemical routines and drivers.

## Installation

You can install all dependencies with: 
```bash
conda env create -f environment.yml
```
This will create a new environment `mechdriver-env` to isolate the current installation.
## Usage

You can install `mechdriver` after  activating this environment
```
conda activate mechdriver-env
python setup.py install
```
You need to create input files and place them in `inp` directory. The list of these files are given below.
You can also find examples in the `tests/inp` directory. You can start the calculations by running `automech.py`
in the `bin` directory with the command line option for the path to the directory that contains `inp` directory.
```
python -u $PATH_TO_MECHDRIVER/mechdriver/bin/automech.py $PATH_TO_INP >& output2.txt
```

To deactivate an active environment, use
```
conda deactivate
```
### Input Files
1. run.dat
2. mechanism.dat
3. species.csv
4. theory.dat
5. models.dat
6. species.dat
7. active space, geom.xyz (aux)

## Code Structure
### Drivers
1. esdriver: carries out electronic structure tasks for requested species and theories
2. thermodriver: carries out partition function calculations for requested species and models
3. ktpdriver: carries out rate calculations for requested mechanism and models
4. transdriver: carries out energy transfer calculations for a given mechanism
5. printdriver: creates txt and csv output for requested species, modules, and theories
6. sordriver:

### Mechroutines
1. es
2. pf
3. output

### Mechlib
###
