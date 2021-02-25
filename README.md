# Mechdriver
This repository houses the main AutoMech executable (bin/automech.py) and well as the high-level electronic structure and thermochemical routines and drivers.

## Installation
conda ...

## Usage
### Input Files
1. run.dat
2. mechanism.dat
3. species.csv
4. theory.dat
5. models.dat
6. species.dat

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
