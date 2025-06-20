""" Removes species from a csv file that are not in the mechanism

    WARNING: Does not account for species that are listed in the "SPECIES"
    list but do not show up in the reactions (e.g., bath gases). Need to fix
    so that these are retained in the spc.csv file
"""

import sys
from mechanalyzer import parser
from mechanalyzer.parser import new_spc as spc_parser
from mechanalyzer.builder import checker
from ioformat import pathtools

# INPUTS
# Filenames
MECH_FILES = [
    'healy_2010.ckin',
]
THERM_FILES = [
]
SPC_FILES = [
    'aramco3.csv',
]

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


# DON'T CHANGE ANYTHING BELOW THIS LINE ###

# Load objects
JOB_PATH = sys.argv[1]
RXN_PARAM_DCTS = parser.ckin_.load_rxn_param_dcts(MECH_FILES, JOB_PATH)
MECH_SPC_DCTS = spc_parser.load_mech_spc_dcts(
    SPC_FILES, JOB_PATH, chk_ste=CHK_STE, chk_match=CHK_MATCH)

# Get the ones I care about
mech_spc_dct = MECH_SPC_DCTS[0]
rxn_param_dct = RXN_PARAM_DCTS[0]

# Go through the mech_spc_dct and remove any species that don't show up in 
# the mechanism
_, missing_from_mech = checker.get_missing_spcs(rxn_param_dct, mech_spc_dct)
for missing_spc in missing_from_mech:
    mech_spc_dct.pop(missing_spc)

# Write to a new species.csv file
headers = ['smiles', 'inchi', 'mult', 'charge', 'exc_flag', ]
csv_str = parser.spc.csv_string(mech_spc_dct, headers)
pathtools.write_file(csv_str, JOB_PATH, 'clean.csv')
