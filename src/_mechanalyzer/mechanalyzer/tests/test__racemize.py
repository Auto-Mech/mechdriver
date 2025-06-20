""" This script tests the merging of a stereo-specific submechanism (e.g., one
    calculated with AutoMech) with a non-stereo mechanism from the literature 
"""

import os
import tempfile
import numpy
from mechanalyzer.parser import new_spc as spc_parser
from mechanalyzer.parser import spc as old_spc_parser
from mechanalyzer.parser import ckin_ as ckin_parser
from mechanalyzer.builder import merge_ste
from mechanalyzer.builder import racemize
from mechanalyzer.calculator.rates import check_p_t
from chemkin_io.writer import mechanism
from ioformat import pathtools
from autoreact import params
from automol import amchi

# Set Paths to test/data directory and output directory
DAT_PATH = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'data')
TMP_OUT = tempfile.mkdtemp()

PRESSURES = [0.03, 0.1, 0.3, 1, 3, 10, 30, 100]
TEMPS = [numpy.arange(360, 1500, 60)]

# Filenames
SPC_CSV = 'strip_ste.csv'
#SPC_CSV = '230810.csv'
#CKIN = 'amech.ckin'
#CKIN = '230810.ckin'
CKIN = 'strip_ste.ckin'
#CKIN = 'amech_140423.ckin'
#THERM = '230810.therm'
THERM = 'strip_ste.therm'

# Load things
MECH_SPC_DCT = spc_parser.load_mech_spc_dct(SPC_CSV, DAT_PATH, canon_ent=True)
RXN_PARAM_DCT = ckin_parser.load_rxn_param_dct(CKIN, DAT_PATH)
SPC_NASA7_DCT = ckin_parser.load_spc_nasa7_dct(THERM, DAT_PATH)

lump_rxn_param_dct, rac_spc_nasa7_dct, rac_mech_spc_dct = racemize.main(
    RXN_PARAM_DCT, SPC_NASA7_DCT, MECH_SPC_DCT, TEMPS, PRESSURES)

## Run the check for unbalanced reactions. Probably dumb, but whatever
##racemize.check_bal_rxns(rac_rxn_param_dct, rac_mech_spc_dct)
## commenting out so I can do my stupid test reactions
#
## Get isomer sets, racemic_sets, and racemic rxn_param_dct
#iso_sets = racemize.find_iso_sets(MECH_SPC_DCT)
#rac_sets, rac_names, rac_mech_spc_dct = racemize.get_rac_sets(
#    iso_sets, MECH_SPC_DCT)
#rac_rxn_param_dct = racemize.get_rac_rxn_param_dct(
#    rac_sets, rac_names, RXN_PARAM_DCT)
#
## Get the lumped rxn parameter dictionary
#lump_rxn_param_dct = racemize.lump(rac_rxn_param_dct, TEMPS, PRESSURES)
#
## This just checks for missing species in the thermo file
#for spc in MECH_SPC_DCT:
#    if spc not in SPC_NASA7_DCT:
#        print(f'species {spc} not in thermo')
#
## Get the thermo dct with only racemized species names
#rac_spc_nasa7_dct = racemize.get_rac_spc_nasa7_dct(rac_names, SPC_NASA7_DCT)

# Generate strings
mech_str = mechanism.write_chemkin_file(
    mech_spc_dct=rac_mech_spc_dct, rxn_param_dct=lump_rxn_param_dct)
headers = ('smiles', 'inchi', 'mult', 'charge', 'exc_flag')
csv_str = old_spc_parser.csv_string(rac_mech_spc_dct, headers)
therm_str = mechanism.thermo_block(rac_spc_nasa7_dct)

# Write to files
pathtools.write_file(mech_str, DAT_PATH, 'racemize.ckin')
pathtools.write_file(csv_str, DAT_PATH, 'racemize.csv')
pathtools.write_file(therm_str, DAT_PATH, 'racemize.therm')

