""" This script tests the merging of a stereo-specific submechanism (e.g., one
    calculated with AutoMech) with a non-stereo mechanism from the literature 
"""

import os
import tempfile
import numpy
from mechanalyzer.parser import new_spc as spc_parser
from mechanalyzer.parser import ckin_ as ckin_parser
from mechanalyzer.builder import merge_ste
from mechanalyzer.calculator.rates import check_p_t
from chemkin_io.writer import mechanism
from ioformat import pathtools
from autoreact import params

# Set Paths to test/data directory and output directory
DAT_PATH = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'data')
TMP_OUT = tempfile.mkdtemp()

# Filenames
CALC_SPC_CSV = 'merge_ste.csv'
CALC_CKIN = 'merge_ste.ckin'
CALC_THERM = 'merge_ste.therm'
NOSTE_SPC_CSV = 'nuig_no_ste.csv'
NOSTE_CKIN = 'nuig_no_ste.ckin'
NOSTE_THERM = 'nuig_no_ste.therm'

# Load things
CALC_MECH_SPC_DCT = spc_parser.load_mech_spc_dct(
    CALC_SPC_CSV, DAT_PATH, chk_ste=False, chk_match=False)
CALC_RXN_PARAM_DCT = ckin_parser.load_rxn_param_dct(CALC_CKIN, DAT_PATH)
CALC_SPC_NASA7_DCT = ckin_parser.load_spc_nasa7_dct(CALC_THERM, DAT_PATH)
NOSTE_MECH_SPC_DCT = spc_parser.load_mech_spc_dct(
    NOSTE_SPC_CSV, DAT_PATH, chk_ste=False, chk_match=False)
NOSTE_RXN_PARAM_DCT = ckin_parser.load_rxn_param_dct(NOSTE_CKIN, DAT_PATH)
NOSTE_SPC_NASA7_DCT = ckin_parser.load_spc_nasa7_dct(NOSTE_THERM, DAT_PATH)

# Expand the dict
new_rxn_param_dct, only_ste_rxn_param_dct = merge_ste.expand_all_rxns(
    CALC_MECH_SPC_DCT, NOSTE_MECH_SPC_DCT, NOSTE_RXN_PARAM_DCT)

# Rename the spcs in the mech_spc_dct and the thermo
new_mech_spc_dct, new_spc_nasa7_dct = merge_ste.rename_spc(
    CALC_MECH_SPC_DCT, NOSTE_MECH_SPC_DCT, CALC_SPC_NASA7_DCT, 
    NOSTE_SPC_NASA7_DCT)

# Write the mechanism to a Chemkin file
mech_str = mechanism.write_chemkin_file(
    rxn_param_dct=new_rxn_param_dct,
    mech_spc_dct=new_mech_spc_dct,
    spc_nasa7_dct=new_spc_nasa7_dct)
only_ste_mech_str = mechanism.write_chemkin_file(
    rxn_param_dct=only_ste_rxn_param_dct,
    mech_spc_dct=new_mech_spc_dct,
    spc_nasa7_dct=new_spc_nasa7_dct)
pathtools.write_file(mech_str, DAT_PATH, 'merge_ste.out')
pathtools.write_file(only_ste_mech_str, DAT_PATH, 'merge_ste_only_ste.out')

