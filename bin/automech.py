"""
   Main Driver to parse and sort the mechanism input files and
   launch the desired drivers
"""

import sys
from drivers import esdriver
from drivers import thermodriver
from drivers import ktpdriver
from lib.amech_io import parser
from lib.amech_io import printer
from lib.reaction import direction as rxndirn
from lib.filesys.build import prefix_fs


# Set runtime options based on user input
JOB_PATH = sys.argv[1]

# Print the header message and host name
printer.program_header('amech')
printer.random_cute_animal()
printer.host_name()
printer.program_header('inp')

# Parse the run input
print('\nReading run.dat...')
RUN_INP_DCT = parser.run.build_run_inp_dct(JOB_PATH)
RUN_OBJ_DCT = parser.run.objects_dct(JOB_PATH)
RUN_JOBS_LST = parser.run.build_run_jobs_lst(JOB_PATH)
ES_TSK_STR = parser.run.read_es_tsks(JOB_PATH)

# Parse the theory input
print('\nReading theory.dat...')
THY_DCT = parser.theory.build_thy_dct(JOB_PATH)

# Parse the model input
print('\nReading model.dat...')
PES_MODEL_DCT, SPC_MODEL_DCT = parser.model.read_models_sections(JOB_PATH)

# Parse the species input to get a dct with ALL species in mechanism
print('\nReading species.csv...')
SPC_DCT = parser.species.build_spc_dct(JOB_PATH, 'csv')

# Parse mechanism input and get a dct with info on PESs user request to run
if RUN_OBJ_DCT['pes']:
    print('\nRunning Calculations for PESs. Need input for mechanism.')
    CLA_DCT = rxndirn.parse_rxn_class_file(JOB_PATH)
    print('  Reading mechanism.dat...')
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
    print('No Proper Run object specified')
    sys.exit()

# Build a dictionary of submission scripts (to finish)
# SUB_SCRIPT_DCT = build_sub_script_dct(JOB_PATH)

# Kill run if just printing mechanism is wanted
if RUN_INP_DCT['print_mech']:
    print('\n\n')
    printer.program_exit('amech')
    sys.exit()

# Initialize the filesystem
print('\nBuilding the base Run-Save filesystems at')
prefix_fs(RUN_INP_DCT['run_prefix'])
print('{}'.format(RUN_INP_DCT['run_prefix']))
prefix_fs(RUN_INP_DCT['save_prefix'])
print('{}'.format(RUN_INP_DCT['save_prefix']))

# Print messages describing drivers and tasks running
print('\nDrivers and tasks user has requested to be run...')
RUN_ES = bool('es' in RUN_JOBS_LST)
WRITE_MESSPF, RUN_MESSPF, RUN_NASA = parser.run.set_thermodriver(RUN_JOBS_LST)
WRITE_MESSRATE, RUN_MESSRATE, RUN_FITS = parser.run.set_ktpdriver(RUN_JOBS_LST)
if RUN_ES:
    print('  - ESDriver')
if WRITE_MESSPF or RUN_MESSPF or RUN_NASA:
    print('  - ThermoDriver')
    if WRITE_MESSPF:
        print('    - write_messpf')
    if RUN_MESSPF:
        print('    - run_messpf')
    if RUN_MESSPF:
        print('    - run_nasa')
if WRITE_MESSRATE or RUN_MESSRATE or RUN_FITS:
    print('  - kTPDriver')
    if WRITE_MESSRATE:
        print('    - write_messrate')
    if RUN_MESSRATE:
        print('    - run_messrate')
    if RUN_FITS:
        print('    - run_fits')

printer.program_exit('inp')

# ESDriver
if RUN_ES:

    printer.program_header('es')

    # Build the elec struct tsk lst
    ES_TSK_LST = parser.run.build_run_es_tsks_lst(
        ES_TSK_STR, SPC_MODEL_DCT, THY_DCT)

    # Call ESDriver for spc in each PES or SPC
    if RUN_OBJ_DCT['pes']:
        for (formula, pes_idx, sub_pes_idx), rxn_lst in RUN_PES_DCT.items():

            # Print PES form and SUB PES Channels
            print('\nRunning PES {}: {}, SUB PES {}'.format(
                pes_idx, formula, sub_pes_idx))
            for rxn in rxn_lst:
                print('  Running Channel {}: {} = {}'.format(
                    rxn['chn_idx'],
                    '+'.join(rxn['reacs']),
                    '+'.join(rxn['prods'])))

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

    printer.program_exit('es')


# ThermoDriver
if WRITE_MESSPF or RUN_MESSPF or RUN_NASA:

    printer.program_header('thermo')

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
            print('\nCalculating Thermochem for species: {}'.format(spc))
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

    printer.program_exit('thermo')

# kTPDriver
if WRITE_MESSRATE or RUN_MESSRATE or RUN_FITS:

    printer.program_header('ktp')

    # Call kTPDriver for each SUB PES
    if RUN_OBJ_DCT['pes']:
        for (formula, pes_idx, sub_pes_idx), rxn_lst in RUN_PES_DCT.items():

            # Print PES form and SUB PES Channels
            print('\nCalculating Rates for PES {}: {}, SUB PES {}'.format(
                pes_idx, formula, sub_pes_idx))
            for chn_idx, rxn in enumerate(rxn_lst):
                print('  Including Channel {}: {} = {}'.format(
                    rxn['chn_idx'],
                    '+'.join(rxn['reacs']),
                    '+'.join(rxn['prods'])))

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
        print("Can't run kTPDriver without a PES being specified")

    printer.program_exit('ktp')

# Exit Program
print('\n\n')
printer.program_exit('amech')
