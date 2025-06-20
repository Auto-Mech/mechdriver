import sys
import numpy
from mechanalyzer.builder import checker
import mechanalyzer.parser.spc as old_spc_parser
import mechanalyzer.parser.new_spc as spc_parser
import mechanalyzer.parser.ckin_ as ckin_parser
from ioformat import pathtools


# INPUTS
RXNS_OR_THERMO = 'rxns'  # 'rxns' or 'thermo'
FILENAMES = [
    'fullckin.ckin',
]
SPC_FILENAMES = [
    'species.csv',
]
OUTPUT_FILENAME = 'spc_check.txt'
CHK_STER = False
CHK_MATCH = False

# Load dcts
assert len(FILENAMES) == len(SPC_FILENAMES)
assert RXNS_OR_THERMO in ('rxns', 'thermo')
JOB_PATH = sys.argv[1]
MECH_SPC_DCTS = spc_parser.load_mech_spc_dcts(SPC_FILENAMES, JOB_PATH, 
                                              chk_ste=CHK_STER,
                                              chk_match=CHK_MATCH)
if RXNS_OR_THERMO == 'rxns':
    RXN_PARAM_DCTS = ckin_parser.load_rxn_param_dcts(FILENAMES, JOB_PATH)
    # Run checker
    OUTPUT_STR = ''
    for IDX, RXN_PARAM_DCT in enumerate(RXN_PARAM_DCTS):
        MECH_SPC_DCT = MECH_SPC_DCTS[IDX]
        CSV, MECH = checker.get_missing_spcs(RXN_PARAM_DCT, MECH_SPC_DCT)
        OUTPUT_STR += 'Filename:\n' + FILENAMES[IDX] + '\n\n'
        OUTPUT_STR += checker.write_missing_spcs(CSV, MECH)

else:  # 'thermo'
    DUMMY_TEMPS = numpy.array([500, 600])
    SPC_THERMO_DCTS = ckin_parser.load_spc_therm_dcts(FILENAMES, JOB_PATH,
                                                      DUMMY_TEMPS)
    OUTPUT_STR = ''
    for IDX, SPC_THERMO_DCT in enumerate(SPC_THERMO_DCTS):
        MECH_SPC_DCT = MECH_SPC_DCTS[IDX]
        THERMO_SPCS = set(SPC_THERMO_DCT.keys())
        CSV_SPCS = set(MECH_SPC_DCT.keys())
        CSV = list(THERMO_SPCS - CSV_SPCS)  # missing from csv
        THERMO = list(CSV_SPCS - THERMO_SPCS)  # missing from thermo
        OUTPUT_STR += 'Filename:\n' + FILENAMES[IDX] + '\n\n'
        OUTPUT_STR += checker.write_missing_spcs(CSV, THERMO)

pathtools.write_file(OUTPUT_STR, JOB_PATH, OUTPUT_FILENAME)

# Write each mech_spc_dct to a new .csv
headers = ('smiles', 'inchi', 'mult', 'exc_flag', 'charge')
for MECH_IDX, MECH_SPC_DCT in enumerate(MECH_SPC_DCTS):
    FNAME = f'reprocessed{MECH_IDX}.csv'
    CSV_STR = old_spc_parser.csv_string(MECH_SPC_DCT, headers) 
    pathtools.write_file(CSV_STR, JOB_PATH, FNAME)
