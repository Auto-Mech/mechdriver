""" test unstable reaction finder
"""

import os
import ioformat
import mechanalyzer


# Set useful global variables
CWD = os.getcwd()

# Read input species and mechanism files into dictionary
INP_SPC_STR = ioformat.pathtools.read_file(
    CWD, 'species.csv',
    remove_comments='!', remove_whitespace=True)
INP_MECH_STR = ioformat.pathtools.read_file(
    CWD, 'mechanism.dat',
    remove_comments='#', remove_whitespace=True)

# Build the initial dictionaries
mech_spc_dct = mechanalyzer.parser.spc.build_spc_dct(INP_SPC_STR, 'csv')
rxn_param_dct, _, _ = mechanalyzer.parser.mech.parse_mechanism(
    INP_MECH_STR, 'chemkin', mech_spc_dct)

# Remove reactions that should not be there
print('\n Removing reactions')
rxn_param_dct = mechanalyzer.builder.remove_unstable_reactions(
    rxn_param_dct, mech_spc_dct)
