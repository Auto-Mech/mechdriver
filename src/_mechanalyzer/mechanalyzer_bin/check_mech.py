""" Runs the mechanism checker on a Chemkin-formatted mechanism file

    Run with "python check_mech.py <path/to/folder>"
"""

import numpy as np
import sys
import ioformat.pathtools as fileio
import mechanalyzer.parser.ckin_ as ckin_parser
from mechanalyzer.builder import checker

# INPUTS
MECH_FILENAME = '209_highT.ckin'
TEMPS = [np.array([360, 400, 600, 800, 1000, 1200, 1400, 1500, 1700, 1900, 2100, 2200, 2400, 2600, 2800, 3000, 3200, 3400])]
PRESSURES = [1, 10, 100]
K_THRESHOLDS = [6e12, 1e15, 1e22]
RXN_NUM_THRESHOLD = 2
OUT_FILENAME = 'mech_check.txt'

# Load dcts
JOB_PATH = sys.argv[1]
RXN_PARAM_DCT = ckin_parser.load_rxn_param_dct(MECH_FILENAME, JOB_PATH)
RXN_KTP_DCT = ckin_parser.load_rxn_ktp_dct(MECH_FILENAME, JOB_PATH,
                                           TEMPS, PRESSURES)

output_str = checker.run_all_checks(
    RXN_PARAM_DCT, RXN_KTP_DCT, K_THRESHOLDS, RXN_NUM_THRESHOLD)
fileio.write_file(output_str, JOB_PATH, OUT_FILENAME)
