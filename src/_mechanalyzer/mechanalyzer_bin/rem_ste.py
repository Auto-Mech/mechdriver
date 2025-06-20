""" Removal script
"""

import os
import time
import ioformat
import chemkin_io
import mechanalyzer
from mechanalyzer.builder import sorter


# Initialize the start time for script execution
t0 = time.time()

# Set useful global variables
CWD = os.getcwd()

# READ AND PARSE INPUT
print('\n---- Parsing the species and mechanism files ---\n')

# Read the build input file and set up info
BLD_STR = ioformat.pathtools.read_file(
    CWD, 'build.dat',
    remove_comments='#', remove_whitespace=True)

FILE_DCT, _ = mechanalyzer.parser.build_input_file(BLD_STR)

# Read input species and mechanism files into dictionary
INP_SPC_STR = ioformat.pathtools.read_file(
    CWD, FILE_DCT['inp_spc'],
    remove_comments='!', remove_whitespace=True)
INP_MECH_STR = ioformat.pathtools.read_file(
    CWD, FILE_DCT['inp_mech'],
    remove_comments='#', remove_whitespace=True)
SORT_STR = ioformat.pathtools.read_file(
    CWD, FILE_DCT['sort'],
    remove_comments='#', remove_whitespace=True)

# Build the initial dictionaries
inp_mech_spc_dct = mechanalyzer.parser.spc.build_spc_dct(INP_SPC_STR, 'csv')
inp_mech_rxn_dct, _, _ = mechanalyzer.parser.mech.parse_mechanism(
    INP_MECH_STR, 'chemkin', inp_mech_spc_dct)
isolate_spc, sort_lst = mechanalyzer.parser.mech.parse_sort(SORT_STR)

# Remove stereochemistry from the mechanism
print('\n Removing stereo improper reactions')
noste_rxn_dct, noste_spc_dct = mechanalyzer.builder.remove_stereochemistry(
    inp_mech_rxn_dct, inp_mech_spc_dct)

# OBTAIN SORTED SPECIES AND MECHANISMS

# Write the new (sorted) species dictionary to a string
HEADERS = ('smiles', 'inchi', 'inchikey', 'mult', 'charge')
noste_spc_sort_dct = mechanalyzer.parser.spc.reorder_by_atomcount(
    noste_spc_dct)
csv_str = mechanalyzer.parser.spc.csv_string(noste_spc_sort_dct, HEADERS)

# Write initial string to call the sorter
mech_str = chemkin_io.writer.mechanism.write_chemkin_file(
    elem_tuple=None,
    mech_spc_dct=noste_spc_sort_dct,
    spc_nasa7_dct=None,
    rxn_param_dct=noste_rxn_dct,
    rxn_cmts_dct=None)

# Use strings to generate ordered objects
param_dct_sort, _, ste_spc_dct_sort, cmts_dct, _ = sorter.sorted_mech(
    csv_str, mech_str, isolate_spc, sort_lst)
rxn_cmts_dct = chemkin_io.writer.comments.get_rxn_cmts_dct(
    rxn_sort_dct=cmts_dct)

# WRITE OUTPUT AND EXIT

# Write the dictionaries to ordered strings
csv_str = mechanalyzer.parser.spc.csv_string(
    ste_spc_dct_sort, FILE_DCT['headers'])
mech_str = chemkin_io.writer.mechanism.write_chemkin_file(
    elem_tuple=elems,
    mech_spc_dct=ste_spc_dct_sort,
    spc_nasa7_dct=None,
    rxn_param_dct=param_dct_sort,
    rxn_cmts_dct=rxn_cmts_dct)

# Write the species and mechanism files
ioformat.pathtools.write_file(csv_str, CWD, FILE_DCT['out_spc'])
ioformat.pathtools.write_file(mech_str, CWD, FILE_DCT['out_mech'])

# Compute script run time and print to screen
tf = time.time()
print('\n\nScript executed successfully.')
print(f'Time to complete: {tf-t0:.2f}')
print('Exiting...')
