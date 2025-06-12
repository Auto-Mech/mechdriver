#!/usr/env python
""" Modifies the species.csv file in ways requested by the user:

    (0) adds inchis to smiles
    (1) adds required heat-of-formation basis species not present in csv file
    (2) adds stereochemistry to species in the file csv
"""

import os
import sys
import time
import argparse
import ioformat
import mechanalyzer


# Set useful global variables
CWD = os.getcwd()

# Parse the command line
PAR = argparse.ArgumentParser()
PAR.add_argument('-s', '--stereo', default=False, type=bool,
                 help='add stereochemistry to species (False)')
PAR.add_argument('-c', '--canonical', default=False, type=bool,
                 help='add canonical enantiomer (False)')
PAR.add_argument('-a', '--amchi', default=False, type=bool,
                 help='turn bad inchis into amchis')
PAR.add_argument('-b', '--hof-basis', default=False, type=bool,
                 help='add heat-of-formation species (False)')
PAR.add_argument('-u', '--instability', default=False, type=bool,
                 help='add instability product species (False)')
PAR.add_argument('-g', '--sort', default=False, type=bool,
                 help='sort the species in the CSV file by atom counts')
PAR.add_argument('-n', '--nprocs', default=1, type=int,
                 help='number of processors to use for tasks (1)')
PAR.add_argument('-i', '--input', default='species.csv',
                 help='name of input species mechanism file (species.csv)')
PAR.add_argument('-o', '--output', default='mod_species.csv',
                 help='name of outpt species mechanism file (mod_species.csv)')
OPTS = vars(PAR.parse_args())

# Initialize the start time for script execution
t0 = time.time()

# Check if any runtime options
if not any([
        OPTS['hof_basis'], OPTS['stereo'],
        OPTS['instability'], OPTS['canonical'],
        OPTS['amchi']]):
    print('Neither stereo, basis, nor instability job specified.')
    print('Add one of [-b, -s, -u] flags to command.')
    print('Exiting...')
    sys.exit()

# Read input species file into a species dictionary
SPC_STR = ioformat.pathtools.read_file(CWD, OPTS['input'])
mech_spc_dct = mechanalyzer.parser.new_spc.parse_mech_spc_dct(SPC_STR)
# mech_spc_dct = mechanalyzer.parser.spc.build_spc_dct(SPC_STR, 'csv')

# Add species relating to unstable because of nearby radicals
if OPTS['instability']:
    mech_spc_dct = mechanalyzer.parser.spc.add_instability_products(
        mech_spc_dct, nprocs=OPTS['nprocs'], stereo=True)

if OPTS['amchi']:
    mech_spc_dct = mechanalyzer.parser.new_spc.mech_inchi_to_amchi(
        mech_spc_dct)

# Add the stereochemical labels to the species
if OPTS['stereo']:
    mech_spc_dct = mechanalyzer.parser.spc.stereochemical_spc_dct(
        mech_spc_dct, nprocs=OPTS['nprocs'], all_stereo=False)

# Sort the species dictionary, if requested
if OPTS['sort']:
    mech_spc_dct = mechanalyzer.parser.spc.reorder_by_atomcount(
        mech_spc_dct)

HEADERS = ('smiles', 'inchi', 'inchikey', 'mult', 'charge')
if OPTS['canonical']:
    mech_spc_dct = mechanalyzer.parser.new_spc.add_canonical_enantiomer(
        mech_spc_dct)
    HEADERS += ('canon_enant_ich',)


# Add the thermochemical species to the species dictionary
if OPTS['hof_basis']:
    mech_spc_dct = mechanalyzer.parser.new_spc.add_canonical_enantiomer(
        mech_spc_dct, dummy=True)
    # mech_spc_dct = mechanalyzer.parser.spc.add_heat_of_formation_basis(
    #    mech_spc_dct, ref_schemes=('cbh0', 'cbh1'),
    mech_spc_dct = mechanalyzer.parser.spc.add_heat_of_formation_basis(
        mech_spc_dct, ref_schemes=('cbh0', 'cbh1', 'cbh2'),
        nprocs=OPTS['nprocs'])

# Write the new species dictionary to a string
csv_str = mechanalyzer.parser.spc.csv_string(mech_spc_dct, HEADERS)

# Write the string to a file
ioformat.pathtools.write_file(csv_str, CWD, OPTS['output'])

# Compute script run time and print to screen
tf = time.time()
print('\n\nScript executed successfully.')
print(f'Time to complete: {tf-t0:.2f}')
print('Exiting...')
