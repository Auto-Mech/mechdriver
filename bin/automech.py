"""
   Main Driver to parse and sort the mechanism input files and
   launch the desired drivers
"""

import sys
from drivers import esdriver
from drivers import thermodriver
from drivers import ktpdriver
from drivers import transdriver
from drivers import printdriver
from mechlib.amech_io import parser
from mechlib.amech_io import printer as ioprinter
# from mechlib.amech_io import print_host_name
from mechlib.reaction import rxnid
from mechlib.filesys import prefix_fs


# Set runtime options based on user input
JOB_PATH = sys.argv[1]

# Print the header message and host name
ioprinter.program_header('amech')
ioprinter.random_cute_animal()
ioprinter.host_name()
ioprinter.program_header('inp')

# Parse the run input
ioprinter.reading('run.dat...', newline=1)
RUN_INP_DCT = parser.run.build_run_inp_dct(JOB_PATH)
RUN_OBJ_DCT = parser.run.objects_dct(JOB_PATH)
RUN_JOBS_LST = parser.run.build_run_jobs_lst(JOB_PATH)
ES_TSK_STR = parser.run.read_es_tsks(JOB_PATH)
PRINT_TSK_STR = parser.run.read_print_tsks(JOB_PATH)
TRANS_TSK_STR = parser.run.read_trans_tsks(JOB_PATH)

# Parse the theory input
ioprinter.reading('theory.dat...', newline=1)
THY_DCT = parser.theory.build_thy_dct(JOB_PATH)

# Parse the model input
ioprinter.reading('model.dat...', newline=1)
PES_MODEL_DCT, SPC_MODEL_DCT = parser.model.read_models_sections(JOB_PATH)

# Parse the species input to get a dct with ALL species in mechanism
ioprinter.reading('species.csv...', newline=1)
SPC_DCT = parser.species.build_spc_dct(JOB_PATH, 'csv')

# Parse mechanism input and get a dct with info on PESs user request to run
if RUN_OBJ_DCT['pes']:
    ioprinter.running(
        'Calculations for PESs. Need input for mechanism.', newline=1)
    CLA_DCT = rxnid.parse_rxn_class_file(JOB_PATH)
    ioprinter.reading('mechanism.dat...', newline=1)
    RUN_PES_DCT = parser.mechanism.build_pes_dct(
        JOB_PATH,
        RUN_INP_DCT['mech'],
        SPC_DCT,
        RUN_OBJ_DCT['pes'],
        sort_rxns=True
    )
elif RUN_OBJ_DCT['spc']:
    RUN_PES_DCT = {}
    RUN_SPC_LST_DCT = parser.species.build_run_spc_dct(SPC_DCT, RUN_OBJ_DCT)
    CLA_DCT = {}
else:
    ioprinter.error_message('No Proper Run object specified')
    sys.exit()

# Build a dictionary of submission scripts (to finish)
# SUB_SCRIPT_DCT = build_sub_script_dct(JOB_PATH)

# Kill run if just printing mechanism is wanted
if RUN_INP_DCT['print_mech']:
    ioprinter.obj('vspace')
    ioprinter.program_exit('amech')
    sys.exit()

# Initialize the filesystem
ioprinter.info_message('Building the base Run-Save filesystems at', newline=1)
prefix_fs(RUN_INP_DCT['run_prefix'])
ioprinter.info_message('{}'.format(RUN_INP_DCT['run_prefix']), indent=1)
prefix_fs(RUN_INP_DCT['save_prefix'])
ioprinter.info_message('{}'.format(RUN_INP_DCT['save_prefix']), indent=1)

# Print messages describing drivers and tasks running
ioprinter.info_message(
    'Drivers and tasks user has requested to be run...', newline=1)
RUN_ES = bool('es' in RUN_JOBS_LST)
RUN_PRINT = bool('print' in RUN_JOBS_LST)
WRITE_MESSPF, RUN_MESSPF, RUN_NASA = parser.run.set_thermodriver(RUN_JOBS_LST)
WRITE_MESSRATE, RUN_MESSRATE, RUN_FITS = parser.run.set_ktpdriver(RUN_JOBS_LST)
RUN_TRANS = bool('transport' in RUN_JOBS_LST)
ioprinter.driver_tasks(
    RUN_ES, WRITE_MESSPF, RUN_MESSPF, RUN_NASA,
    WRITE_MESSRATE, RUN_MESSRATE, RUN_FITS, RUN_TRANS)

ioprinter.program_exit('inp')

# ESDriver
if RUN_ES:

    ioprinter.program_header('es')

    # Build the elec struct tsk lst
    ES_TSK_LST = parser.tsks.es_tsk_lst(ES_TSK_STR, THY_DCT)

    # Call ESDriver for spc in each PES or SPC
    if RUN_OBJ_DCT['pes']:
        for (formula, pes_idx, sub_pes_idx), rxn_lst in RUN_PES_DCT.items():

            # Print PES form and SUB PES Channels
            ioprinter.pes(pes_idx, formula, sub_pes_idx)
            for rxn in rxn_lst:
                ioprinter.channel(rxn['chn_idx'], rxn['reacs'], rxn['prods'])

            esdriver.run(
                pes_idx,
                rxn_lst,
                SPC_DCT,
                CLA_DCT,
                ES_TSK_LST,
                THY_DCT,
                RUN_INP_DCT
            )
    else:
        PES_IDX = 0
        esdriver.run(
            PES_IDX,
            RUN_SPC_LST_DCT,
            SPC_DCT,
            CLA_DCT,
            ES_TSK_LST,
            THY_DCT,
            RUN_INP_DCT
        )

    ioprinter.program_exit('es')


# ThermoDriver
if WRITE_MESSPF or RUN_MESSPF or RUN_NASA:

    ioprinter.program_header('thermo')

    # Call ThermoDriver for spc in PES
    if RUN_OBJ_DCT['pes']:
        for _, rxn_lst in RUN_PES_DCT.items():
            thermodriver.run(
                SPC_DCT,
                PES_MODEL_DCT, SPC_MODEL_DCT,
                THY_DCT,
                rxn_lst,
                RUN_INP_DCT,
                write_messpf=WRITE_MESSPF,
                run_messpf=RUN_MESSPF,
                run_nasa=RUN_NASA,
            )
    else:
        for spc in RUN_SPC_LST_DCT:
            ioprinter.info_message(
                'Calculating Thermochem for species: {}'.format(spc),
                newline=1)
        thermodriver.run(
            SPC_DCT,
            PES_MODEL_DCT, SPC_MODEL_DCT,
            THY_DCT,
            RUN_SPC_LST_DCT,
            RUN_INP_DCT,
            write_messpf=WRITE_MESSPF,
            run_messpf=RUN_MESSPF,
            run_nasa=RUN_NASA,
        )

    ioprinter.program_exit('thermo')


# TransportDriver
if RUN_TRANS:

    ioprinter.program_header('trans')

    # Build the elec struct tsk lst
    TRANS_TSK_LST = parser.tsks.trans_tsk_lst(TRANS_TSK_STR, THY_DCT)

    # Call ThermoDriver for spc in PES
    if RUN_OBJ_DCT['pes']:
        for _, rxn_lst in RUN_PES_DCT.items():
            transdriver.run(
                SPC_DCT,
                THY_DCT,
                rxn_lst,
                TRANS_TSK_LST,
                RUN_INP_DCT
            )
    else:
        for spc in RUN_SPC_LST_DCT:
            ioprinter.message('Calculating Transport for species: {}'.format(spc), newline=1)
        transdriver.run(
            SPC_DCT,
            THY_DCT,
            RUN_SPC_LST_DCT,
            TRANS_TSK_LST,
            RUN_INP_DCT
        )

# kTPDriver
if WRITE_MESSRATE or RUN_MESSRATE or RUN_FITS:

    ioprinter.program_header('ktp')

    # Call kTPDriver for each SUB PES
    if RUN_OBJ_DCT['pes']:
        for (formula, pes_idx, sub_pes_idx), rxn_lst in RUN_PES_DCT.items():

            # Print PES form and SUB PES Channels
            ioprinter.pes(pes_idx, formula, sub_pes_idx)
            for chn_idx, rxn in enumerate(rxn_lst):
                ioprinter.channel(
                    rxn['chn_idx'],
                    '+'.join(rxn['reacs']),
                    '+'.join(rxn['prods']))

            ktpdriver.run(
                formula, pes_idx, sub_pes_idx,
                SPC_DCT,
                CLA_DCT,
                THY_DCT,
                rxn_lst,
                PES_MODEL_DCT, SPC_MODEL_DCT,
                RUN_INP_DCT,
                write_messrate=WRITE_MESSRATE,
                run_messrate=RUN_MESSRATE,
                run_fits=RUN_FITS
            )
    else:
        ioprinter.error_message("Can't run kTPDriver without a PES being specified")

    ioprinter.program_exit('ktp')

# Print Driver
if RUN_PRINT:

    ioprinter.program_header('print')

    # Build the elec struct tsk lst
    PRNT_TSK_LST = tsks.prnt_tsk_lst(PRNT_TSK_STR, THY_DCT)

    PES_IDX = 0
    printdriver.run(
        PES_IDX,
        RUN_SPC_LST_DCT,
        SPC_DCT,
        CLA_DCT,
        PRINT_TSK_LST,
        THY_DCT,
        RUN_INP_DCT,
        SPC_MODEL_DCT
    )


# Exit Program
ioprinter.obj('vspace')
ioprinter.program_exit('amech')
