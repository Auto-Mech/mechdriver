#!/usr/bin/env python
"""
  Script to fit rates for a standalone
"""

import os
import argparse
import ioformat.pathtools
import mess_io.reader
import chemkin_io
import ratefit


# Set useful global variables
CWD = os.getcwd()

# Parse the command line
PAR = argparse.ArgumentParser()
PAR.add_argument('-m', '--mess', default='rate.out',
                 help='MESS ouput name (rate.out)')
PAR.add_argument('-l', '--label', default='label.inp',
                 help='label dct name (label.inp)')
PAR.add_argument('-c', '--chemkin', default='rate.ckin',
                 help='Chemkin ouput name (rate.ckin)')
PAR.add_argument('-f', '--fit-method', default='plog',
                 help='method to fit the rates (plog, chebyshev)')
OPTS = vars(PAR.parse_args())

# Read label dct
label_dct = ioformat.pathtools.read_json_file(CWD, OPTS['label'])

# Read MESS file and get rate constants
mess_str = ioformat.pathtools.read_file(CWD, OPTS['mess'])
rxn_ktp_dct = mess_io.reader.rates.get_rxn_ktp_dct(
    mess_str, label_dct=label_dct, filter_kts=True
)

# Fit rates
rxn_param_dct, rxn_err_dct = ratefit.fit.fit_rxn_ktp_dct(
    rxn_ktp_dct, OPTS['fit_method'],
)

# Get the comments dct and write the Chemkin string
rxn_cmts_dct = chemkin_io.writer.comments.get_rxn_cmts_dct(
    rxn_err_dct=rxn_err_dct)
ckin_str = chemkin_io.writer.mechanism.write_chemkin_file(
    rxn_param_dct=rxn_param_dct, rxn_cmts_dct=rxn_cmts_dct)

# Write the fitted rate parameter file
ioformat.pathtools.write_file(ckin_str, CWD, OPTS['chemkin'])
