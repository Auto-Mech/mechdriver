"""
   Driver to parse and sort the mechanism input files and
   launch the desired drivers
"""

import sys
from drivers import esdriver
from drivers import thermodriver
from drivers import ktpdriver
from lib.load import run as lrun
from lib.load import theory as lthy
from lib.load import model as lmodel
from lib.load import mechanism as lmech
from lib.load import species as lspc
from lib.filesystem import build as fbuild
from lib import printmsg


# Print the header message for the driver
printmsg.program_header('amech')
printmsg.random_cute_animal()

# Set runtime options based on user input
JOB_PATH = sys.argv[1]

# Parse the run input
print('Parsing the input files...')
print('\nReading run.dat...')
RUN_INP_DCT = lrun.build_run_inp_dct(JOB_PATH)
RUN_OBJ_DCT = lrun.objects_dct(JOB_PATH)
RUN_JOBS_LST = lrun.build_run_jobs_lst(JOB_PATH)
ES_TSK_STR = lrun.read_es_tsks(JOB_PATH)

# Parse the theory input
print('\nReading theory.dat...')
THY_DCT = lthy.build_thy_dct(JOB_PATH)

# Parse the model input
print('\nReading model.dat...')
PES_MODEL_DCT, SPC_MODEL_DCT = lmodel.read_models_sections(JOB_PATH)

# Parse the species input to get a dct with ALL species in mechanism
print('\nReading species.csv...')
SPC_DCT = lspc.build_spc_dct(JOB_PATH, 'csv', check_stereo=False)

# Parse the mechanism input and get a dct with info on PESs user request to run
if RUN_OBJ_DCT['pes']:
    print('\nReaction Channels Needed. Reading mechanism.dat...')
    RUN_PES_DCT = lmech.parse_mechanism_file(
        JOB_PATH,
        RUN_INP_DCT['mech'],
        SPC_DCT,
        RUN_OBJ_DCT['pes'],
        sort_rxns=True
    )
    CLA_DCT = lspc.parse_rxn_class_file(JOB_PATH)
elif RUN_OBJ_DCT['spc']:
    RUN_PES_DCT = {}
    RUN_SPC_LST_DCT = lspc.build_run_spc_dct(SPC_DCT, RUN_OBJ_DCT)
    CLA_DCT = {}
else:
    print('No Proper Run object specified')
    sys.exit()

# Initialize the filesystem
print('\nBuilding the base Run-Save filesystems at')
fbuild.prefix_fs(RUN_INP_DCT['run_prefix'])
print('{}'.format(RUN_INP_DCT['run_prefix']))
fbuild.prefix_fs(RUN_INP_DCT['save_prefix'])
print('{}'.format(RUN_INP_DCT['save_prefix']))

# Run the requested drivers: es, thermo, ktp
print('\n\nRunning the requested drivers...')
if 'es' in RUN_JOBS_LST:
    # Print the header message for the driver
    printmsg.program_header('es')
    # Build the elec struct tsk lst
    ES_TSK_LST = lrun.build_run_es_tsks_lst(
        ES_TSK_STR, SPC_MODEL_DCT, THY_DCT)
    if RUN_OBJ_DCT['pes']:
        # Call ESDriver for spc in each PES
        for (_, pes_idx, sub_pes_idx), rxn_lst in RUN_PES_DCT.items():
            # Do some extra work to prepare the info to pass to the drivers
            # ES_TSK_LST = lrun.build_run_es_tsks_lst(
            #     ES_TSK_STR, SPC_MODEL_DCT, THY_DCT)
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
        # Call ESDriver for all of the species
        pes_idx = 0
        esdriver.run(
            pes_idx,
            RUN_SPC_LST_DCT,
            SPC_DCT,
            CLA_DCT,
            ES_TSK_LST,
            THY_DCT,
            RUN_INP_DCT
        )

WRITE_MESSPF, RUN_MESSPF, RUN_NASA = lrun.set_thermodriver_run(RUN_JOBS_LST)
if WRITE_MESSPF or RUN_MESSPF or RUN_NASA:
    # Print the header message for the driver
    printmsg.program_header('thermo')
    if RUN_OBJ_DCT['pes']:
        # Call ThermoDriver for spc in each PES
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
        # Call ThermoDriver for all of the species
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

WRITE_MESSRATE, RUN_MESSRATE, RUN_FITS = lrun.set_ktpdriver_run(RUN_JOBS_LST)
if WRITE_MESSRATE or RUN_MESSRATE or RUN_FITS:
    # Print the header message for the driver
    printmsg.program_header('ktp')
    if RUN_OBJ_DCT['pes']:
        # Call kTPDriver for spc in each PES
        for (pes_formula, pes_idx, sub_pes_idx), rxn_lst in RUN_PES_DCT.items():
            ktpdriver.run(
                pes_formula, pes_idx, sub_pes_idx,
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
        print("Can't run kTPDriver without a PES being specified")

# Print the program exist message
print('\n\nAutoMech has completed.')
