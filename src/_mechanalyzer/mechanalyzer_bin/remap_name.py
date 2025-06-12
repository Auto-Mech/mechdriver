""" Script to generate files with reactions

    Hashes cannot start with - or a number

"""

import os
import time
import ioformat
import chemkin_io
import mechanalyzer


# Initialize the start time for script execution
t0 = time.time()

# Set useful global variables
CWD = os.getcwd()

# READ AND PARSE INPUT
print('\n---- Parsing the species and mechanism files ---\n')

# Read input species and mechanism files into dictionary
INP_SPC_STR = ioformat.pathtools.read_file(
    CWD, 'ini_species.csv',
    remove_comments='!', remove_whitespace=True)
INP_MECH_STR = ioformat.pathtools.read_file(
    CWD, 'ini_mechanism.dat',
    remove_comments='#', remove_whitespace=True)

# Build the initial dictionaries
mech_spc_dct = mechanalyzer.parser.spc.build_spc_dct(INP_SPC_STR, 'csv')
rxn_param_dct, _, _ = mechanalyzer.parser.mech.parse_mechanism(
    INP_MECH_STR, 'chemkin', mech_spc_dct)

# Build the name mapping dictionary
map_dct = mechanalyzer.builder.functional_group_name_dct(
    mech_spc_dct, rename_rule_dct=None)
re_mech_spc_dct, re_rxn_param_dct = mechanalyzer.builder.remap_mechanism_names(
    mech_spc_dct, rxn_param_dct, map_dct)

# # WRITE OUTPUT AND EXIT

# Write the dictionaries to ordered strings
header_lst = ('smiles', 'inchi', 'mult', 'charge')
csv_str = mechanalyzer.parser.spc.csv_string(
    re_mech_spc_dct, header_lst)
mech_str = chemkin_io.writer.mechanism.write_chemkin_file(
    elem_tuple=(),
    mech_spc_dct=re_mech_spc_dct,
    spc_nasa7_dct=None,
    rxn_param_dct=re_rxn_param_dct,
    rxn_cmts_dct={})

# Write the species and mechanism files
ioformat.pathtools.write_file(csv_str, CWD, 'species.csv')
ioformat.pathtools.write_file(mech_str, CWD, 'mechanism.dat')

# Compute script run time and print to screen
tf = time.time()
print('\n\nScript executed successfully.')
print(f'Time to complete: {tf-t0:.2f}')
print('Exiting...')
