""" Script to generate files with reactions
"""

import sys
import os
import time
import ioformat
import chemkin_io
import mechanalyzer
from mechanalyzer.builder import sorter

# Initialize the start time for script execution
t0 = time.time()

# Set up the paths
CWD = os.getcwd()

# Read the build input file and set up info
BLD_STR = ioformat.pathtools.read_file(
    CWD, 'build.dat',
    remove_comments='#', remove_whitespace=True)

FILE_DCT, RSERIES = mechanalyzer.parser.build_input_file(BLD_STR)
# Read the spc_dct and rxn_param_dct
INP_SPC_STR = ioformat.pathtools.read_file(
    CWD, FILE_DCT['inp_spc'],
    remove_comments='!', remove_whitespace=True)
INP_MECH_STR = ioformat.pathtools.read_file(
    CWD, FILE_DCT['inp_mech'],
    remove_comments='#', remove_whitespace=True)
SORT_STR = ioformat.pathtools.read_file(
    CWD, FILE_DCT['sort'],
    remove_comments='#', remove_whitespace=True)

# Check if the input strings exist
if any(string is None for string in (INP_SPC_STR, INP_MECH_STR, SORT_STR)):
    print('ERROR: Input file(s) species.csv, mechanism.dat, sort.dat missing')
    sys.exit()

mech_spc_dct = mechanalyzer.parser.new_spc.parse_mech_spc_dct(
    INP_SPC_STR)
rxn_param_dct = mechanalyzer.parser.mech.parse_mechanism(
    INP_MECH_STR, 'chemkin')
isolate_spc, sort_lst, _ = mechanalyzer.parser.mech.parse_sort(SORT_STR)
# Generate the requested reactions
mech_spc_dct, rxn_param_dct = mechanalyzer.builder.rxn.build_mechanism(
    mech_spc_dct, rxn_param_dct,
    rxn_series=RSERIES,
    stereo=FILE_DCT.get('stereo', False),
    nprocs=FILE_DCT.get('nprocs', 1))

# Write the dictionaries to original strings
csv_str = mechanalyzer.parser.spc.csv_string(
    mech_spc_dct, FILE_DCT['headers'])
#csv_str = chemkin_io.writer.spc.write_species(mech_spc_dct)
mech_str = chemkin_io.writer.mechanism.write_chemkin_file(
    elem_tuple=None,
    mech_spc_dct=mech_spc_dct,
    spc_nasa7_dct=None,
    rxn_param_dct=rxn_param_dct,
    rxn_cmts_dct=None)

# Write the species and mechanism files
ioformat.pathtools.write_file(csv_str, CWD, FILE_DCT['out_spc'])
ioformat.pathtools.write_file(mech_str, CWD, FILE_DCT['out_mech'])

# Use strings to generate ordered objects
# param_dct_sort, _, mech_spc_dct, cmts_dct, _ = sorter.sorted_mech(
#     csv_str, mech_str, isolate_spc, sort_lst)
# rxn_cmts_dct = chemkin_io.writer.comments.get_rxn_cmts_dct(
#     rxn_sort_dct=cmts_dct)

# Write the dictionaries to ordered strings
#sortd_csv_str = mechanalyzer.parser.spc.csv_string(
#    mech_spc_dct, FILE_DCT['headers'])
# sorted_csv_str = chemkin_io.writer.spc.write_species(mech_spc_dct)
# sortd_mech_str = chemkin_io.writer.mechanism.write_chemkin_file(
#     elem_tuple=None,
#     mech_spc_dct=mech_spc_dct,
#     spc_nasa7_dct=None,
#     rxn_param_dct=param_dct_sort,
#     rxn_cmts_dct=rxn_cmts_dct)
# Write the species and mechanism files
# ioformat.pathtools.write_file(sortd_csv_str, CWD, FILE_DCT['out_spc'])
# ioformat.pathtools.write_file(sortd_mech_str, CWD, FILE_DCT['out_mech'])

# Compute script run time and print to screen
tf = time.time()
print('\n\nScript executed successfully.')
print(f'Time to complete: {tf-t0:.2f}')
print('Exiting...')
