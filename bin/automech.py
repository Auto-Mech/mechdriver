""" Central Execution script to launch a MechDriver process which will
    parse all of the user-supplied input files in a specified directory, then
    launches all of the requested electronic structure, transport,
    thermochemistry and kinetics calculations via their associated
    sub-drivers.
"""

import sys
# import argparse
from mechlib.filesys import prefix_fs
from mechlib.amech_io import parser as ioparser
from mechlib.amech_io import printer as ioprinter
from drivers import esdriver, thermodriver, ktpdriver, transdriver, procdriver
import autofile


# Set runtime options based on user input
JOB_PATH = sys.argv[1]  # Add a check to see if [1] exists; path exits
if len(sys.argv) > 2:
    if sys.argv[2] == 'safemode_off':
        autofile.turn_off_safemode()
        ioprinter.info_message('Running with safemode turned OFF...')

# Print the header message and host name
ioprinter.program_header('amech')
ioprinter.random_cute_animal()
ioprinter.host_name()

# Parse all of the input
ioprinter.program_header('inp')

ioprinter.info_message('\nReading files provided in the inp directory...')
INPUT = ioparser.read_amech_input(JOB_PATH)

ioprinter.info_message('\nParsing input files for runtime parameters...')
THY_DCT = ioparser.thy.theory_dictionary(INPUT['thy'])
KMOD_DCT, SMOD_DCT = ioparser.models.models_dictionary(
    INPUT['mod'], THY_DCT)
INP_KEY_DCT = ioparser.run.input_dictionary(INPUT['run'])
PES_IDX_DCT, SPC_IDX_DCT = ioparser.run.chem_idxs(INPUT['run'])
TSK_LST_DCT = ioparser.run.tasks(INPUT['run'], THY_DCT)
SPC_DCT, GLOB_DCT = ioparser.spc.species_dictionary(
    INPUT['spc'], INPUT['dat'], INPUT['geo'], 'csv')
PES_DCT = ioparser.mech.pes_dictionary(
    INPUT['mech'], 'chemkin', SPC_DCT)

PES_RLST, SPC_RLST = ioparser.rlst.run_lst(
    PES_DCT, SPC_DCT, PES_IDX_DCT, SPC_IDX_DCT)

# Do a check
ioprinter.info_message('\nFinal check if all required input provided...')
ioparser.run.check_inputs(
    TSK_LST_DCT, PES_DCT, KMOD_DCT, SMOD_DCT)

# Build the Run-Save Filesystem Directories
prefix_fs(INP_KEY_DCT['run_prefix'], INP_KEY_DCT['save_prefix'])

# Run Drivers Requested by User
ES_TSKS = TSK_LST_DCT.get('es')
if ES_TSKS is not None:
    ioprinter.program_header('es')
    esdriver.run(
        PES_RLST, SPC_RLST,
        ES_TSKS,
        SPC_DCT, GLOB_DCT, THY_DCT,
        INP_KEY_DCT['run_prefix'], INP_KEY_DCT['save_prefix'],
        print_debug=INP_KEY_DCT['print_debug']
    )
    ioprinter.program_exit('es')

THERM_TSKS = TSK_LST_DCT.get('thermo')
if THERM_TSKS is not None:
    ioprinter.program_header('thermo')
    thermodriver.run(
        PES_RLST, SPC_RLST,
        THERM_TSKS,
        KMOD_DCT, SMOD_DCT,
        SPC_DCT, THY_DCT,
        INP_KEY_DCT['run_prefix'], INP_KEY_DCT['save_prefix'], JOB_PATH
    )
    ioprinter.program_exit('thermo')

TRANS_TSKS = TSK_LST_DCT.get('trans')
if TRANS_TSKS is not None:
    ioprinter.program_header('trans')
    if PES_DCT:
        transdriver.run(
            PES_RLST, SPC_RLST,
            TRANS_TSKS,
            SMOD_DCT,
            SPC_DCT, THY_DCT,
            INP_KEY_DCT['run_prefix'], INP_KEY_DCT['save_prefix'],
        )
    ioprinter.program_exit('trans')

KTP_TSKS = TSK_LST_DCT.get('ktp')
if KTP_TSKS is not None:
    ioprinter.program_header('ktp')
    ktpdriver.run(
        PES_RLST, INPUT['pesgrp'],
        KTP_TSKS,
        SPC_DCT, GLOB_DCT,
        THY_DCT, KMOD_DCT, SMOD_DCT,
        INP_KEY_DCT['run_prefix'], INP_KEY_DCT['save_prefix'], JOB_PATH,
    )
    ioprinter.program_exit('ktp')

PROC_TSKS = TSK_LST_DCT.get('proc')
if PROC_TSKS is not None:
    ioprinter.program_header('proc')
    procdriver.run(
        PES_RLST, SPC_RLST,
        PROC_TSKS,
        SPC_DCT, THY_DCT,
        KMOD_DCT, SMOD_DCT,
        INP_KEY_DCT['run_prefix'], INP_KEY_DCT['save_prefix'], JOB_PATH
    )
    ioprinter.program_exit('proc')

# Check if any drivers were requested to be run
if all(tsks is None
       for tsks in (ES_TSKS, THERM_TSKS, TRANS_TSKS, KTP_TSKS, PROC_TSKS)):
    ioprinter.warning_message(
        'User did not provide (uncommented) driver tasks lists in run.dat')

# Exit Program
ioprinter.obj('vspace')
ioprinter.program_exit('amech')
