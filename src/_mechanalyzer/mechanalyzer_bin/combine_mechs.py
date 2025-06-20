""" Combines two mechanisms into a single mechanism
"""

import sys
from mechanalyzer.calculator import combine
from mechanalyzer import parser
from mechanalyzer.parser import new_spc as spc_parser
from mechanalyzer.parser import spc as old_spc
from chemkin_io.writer import mechanism as mech_writer
from ioformat import pathtools

# INPUTS
# Filenames
MECH_FILES = [
    'mechanism.dat',
    'healy_2010.ckin',
]
THERM_FILES = [
    'therm.dat',
    'healy_2010.therm',
]
SPC_FILES = [
    'species_no_ste.csv',
    'healy_2010.csv',
]
OUT_PREFIX = 'comb_mech'  # prefix for output files (e.g., 'OUT_PREFIX.ckin')
# Options for combining mechanisms:
# 1.) Whether there is stereo ONLY in mech1 that should be included in mech2
STE_MECH1_ONLY = False
# 2.) Whether to strip stereo when comparing inchis (overwritten internally to
# True if STE_MECH1_ONLY=True)
STRIP_STE = False

# Options for parsing species.csv file:
# 1.) Whether to check inchis for stereo completeness
CHK_STE = False
# 2.) Whether to check for matching smiles and inchis
CHK_MATCH = False


# DON'T CHANGE ANYTHING BELOW THIS LINE

# Load objects
JOB_PATH = sys.argv[1]
RXN_PARAM_DCTS = parser.ckin_.load_rxn_param_dcts(MECH_FILES, JOB_PATH)
SPC_NASA7_DCTS = parser.ckin_.load_spc_nasa7_dcts(THERM_FILES, JOB_PATH)
MECH_SPC_DCTS = spc_parser.load_mech_spc_dcts(
    SPC_FILES, JOB_PATH, chk_ste=CHK_STE, chk_match=CHK_MATCH)

# Combine objects
comb_rxn_param_dct, comb_spc_nasa7_dct, comb_mech_spc_dct = combine.comb_mechs(
    RXN_PARAM_DCTS[0], RXN_PARAM_DCTS[1], SPC_NASA7_DCTS[0], SPC_NASA7_DCTS[1],
    MECH_SPC_DCTS[0], MECH_SPC_DCTS[1], ste_mech1_only=STE_MECH1_ONLY,
    strip_ste=STRIP_STE)

# Create Chemkin-formatted strings
mech_str = mech_writer.write_chemkin_file(
    rxn_param_dct=comb_rxn_param_dct,
    mech_spc_dct=comb_mech_spc_dct)  # mech_spc_dct: for species and elements
therm_str = mech_writer.write_chemkin_file(
    spc_nasa7_dct=comb_spc_nasa7_dct)
headers = ['smiles', 'inchi', 'mult', 'charge', 'exc_flag']
csv_str = old_spc.csv_string(comb_mech_spc_dct, headers)

# Write strings to files
pathtools.write_file(mech_str, JOB_PATH, f'{OUT_PREFIX}.ckin')
pathtools.write_file(therm_str, JOB_PATH, f'{OUT_PREFIX}.therm')
pathtools.write_file(csv_str, JOB_PATH, f'{OUT_PREFIX}.csv')
