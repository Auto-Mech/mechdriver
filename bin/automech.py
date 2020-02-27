"""
   Driver to parse and sort the mechanism input files and
   launch the desired drivers
"""

import sys
from drivers import esdriver
from drivers import thermodriver
from drivers import ktpdriver
from lib.load import run as loadrun
from lib.load import theory as loadthy
from lib.load import model as loadmodel
from lib.load import mechanism as loadmech
from lib.load import species as loadspc
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
RUN_INP_DCT = loadrun.build_run_inp_dct(JOB_PATH)
RUN_OBJ_DCT = loadrun.objects_dct(JOB_PATH)
RUN_OPTIONS_DCT = loadrun.build_run_glob_opts_dct(JOB_PATH)
RUN_JOBS_LST = loadrun.build_run_jobs_lst(JOB_PATH)
ES_TSK_STR = loadrun.read_es_tsks(JOB_PATH)

# Parse the theory input
print('\nReading theory.dat...')
THY_DCT = loadthy.build_thy_dct(JOB_PATH)

# Parse the model input
print('\nReading model.dat...')
MODEL_DCT = loadmodel.read_models_sections(JOB_PATH)

# Parse the species input to get a dct with ALL species in mechanism
print('\nReading species.csv...')
SPC_DCT = loadspc.build_spc_dct(JOB_PATH, RUN_INP_DCT['spc'])

# Parse the mechanism input and get a dct with info on PESs user request to run
if RUN_OBJ_DCT['pes'] or RUN_OBJ_DCT['pspc']:
    print('\nReaction Channels Needed. Reading mechanism.dat...')
    # Prob move this into the fxn below cuz I need the model
    RUN_PES_DCT = loadmech.parse_mechanism_file(
        JOB_PATH,
        RUN_INP_DCT['mech'],
        SPC_DCT,
        RUN_OBJ_DCT['pes'],
        sort_rxns=RUN_INP_DCT['sort_rxns']
    )
elif RUN_OBJ_DCT['spc']:
    RUN_PES_DCT = {}
    RUN_SPC_LST_DCT = loadspc.build_run_spc_dct(SPC_DCT, RUN_OBJ_DCT)
else:
    print('No Proper Run object specified')
    sys.exit()

# Print stuff for test
# print('\n\nEchoing the user input:\n\n')
# print('\nrun inp dct')
# print(RUN_INP_DCT)
# print('\nrun options dct')
# print(RUN_OPTIONS_DCT)
# print('\nrun jobs lst')
# print(RUN_JOBS_LST)
# print('\ntheory dct')
# print(THY_DCT)
# print('\nmodel dct')
# print(MODEL_DCT)
# print('\nspc dct')
# print(SPC_DCT)
# print('\npes dct')
# print(PES_DCT)

# Initialize the filesystem
print('\nBuilding the base Run-Save filesystems at')
fbuild.prefix_filesystem(
    RUN_INP_DCT['run_prefix'],
    RUN_INP_DCT['save_prefix']
)
print('{}'.format(RUN_INP_DCT['run_prefix']))
print('{}'.format(RUN_INP_DCT['save_prefix']))

# Run the requested drivers: es, thermo, ktp
print('\n\nRunning the requested drivers...')
if 'es' in RUN_JOBS_LST:
    if RUN_OBJ_DCT['pes'] or RUN_OBJ_DCT['pspc']:
        # Call ESDriver for spc in each PES
        for pes, rxn_lst in RUN_PES_DCT.items():
            esdriver.run(
                rxn_lst,
                SPC_DCT,
                ES_TSK_STR,
                MODEL_DCT,
                THY_DCT,
                RUN_OPTIONS_DCT,
                RUN_INP_DCT
            )
    else:
        # Call ESDriver for all of the species
        esdriver.run(
            RUN_SPC_LST_DCT,
            SPC_DCT,
            ES_TSK_STR,
            MODEL_DCT,
            THY_DCT,
            RUN_OPTIONS_DCT,
            RUN_INP_DCT
        )

write_messpf=False
run_messpf=False
run_thermo=False
if 'thermochem' in RUN_JOBS_LST:
    write_messpf=True
    run_messpf=True
    run_thermo=True
else:
    if 'write_messpf' in RUN_JOBS_LST:
        write_messpf=True
    if 'run_messpf' in RUN_JOBS_LST:
        run_messpf=True
    if 'run_thermo' in RUN_JOBS_LST:
        run_thermo=True
if write_messpf or run_messpf or run_thermo:
#    if bool('no_write_messpf' in RUN_JOBS_LST) write_messpf=False
#    if bool('no_run_messpf' in RUN_JOBS_LST) run_messpf=False
#    if bool('no_run_thermo' in RUN_JOBS_LST) run_thermo=False
    if RUN_OBJ_DCT['pes'] or RUN_OBJ_DCT['pspc']:
        # Call ThermoDriver for spc in each PES
        for pes, rxn_lst in RUN_PES_DCT.items():
            thermodriver.run(
                SPC_DCT,
                MODEL_DCT,
                THY_DCT,
                rxn_lst,
                RUN_INP_DCT,
                ref_scheme='basic',
                write_messpf=write_messpf,
                run_messpf=run_messpf,
                run_thermo=run_thermo,
            )
    else:
        # Call ThermoDriver for all of the species
        thermodriver.run(
            SPC_DCT,
            MODEL_DCT,
            THY_DCT,
            RUN_SPC_LST_DCT,
            RUN_INP_DCT,
            ref_scheme='basic',
            write_messpf=write_messpf,
            run_messpf=run_messpf,
            run_thermo=run_thermo,
        )

if 'rates' in RUN_JOBS_LST or 'fits' in RUN_JOBS_LST:
    if RUN_OBJ_DCT['pes']:
        # Call kTPDriver for spc in each PES
        for pes_formula, rxn_lst in RUN_PES_DCT.items():
            # Get info for the transition states
            ts_dct = loadspc.build_sadpt_dct(
                rxn_lst, MODEL_DCT, THY_DCT, ES_TSK_STR,
                RUN_INP_DCT, RUN_OPTIONS_DCT, SPC_DCT, {})
            SPC_DCT.update(ts_dct)
            # Run the driver
            ktpdriver.run(
                pes_formula,
                SPC_DCT,
                THY_DCT,
                rxn_lst,
                MODEL_DCT,
                RUN_INP_DCT,
                run_rates=bool('rates' in RUN_JOBS_LST),
                run_fits=bool('fits' in RUN_JOBS_LST)
            )
    else:
        print("Can't run kTPDriver without a PES being specified")

# Print the program exist message
print('\n\nAutoMech has completed.')
