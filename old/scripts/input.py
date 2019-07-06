""" reaction list test
"""

import sys
import os
import pandas
import autofile
from autofile import SFS
import automol
import chemkin_io
import moldr

PROG = 'g09'
SCRIPT_STR = ("#!/usr/bin/env bash\n"
              "g09 run.inp run.out")

METHOD = 'wb97xd'
BASIS = '6-31g*'
RUN_PREFIX = 'run'
SAVE_PREFIX = 'save'

if not os.path.exists(RUN_PREFIX):
    os.mkdir(RUN_PREFIX)

if not os.path.exists(SAVE_PREFIX):
    os.mkdir(SAVE_PREFIX)

def _read_file(file_name):
      with open(file_name, encoding='utf8', errors='ignore') as file_obj:
          file_str = file_obj.read()
      return file_str

PATH = os.path.dirname(os.path.realpath(__file__))

SYNGAS_PATH = os.path.join(PATH, 'data/syngas')
SYNGAS_MECH_STR = _read_file(os.path.join(SYNGAS_PATH, 'mechanism.txt'))
SYNGAS_TAB = pandas.read_csv(os.path.join(SYNGAS_PATH, 'smiles.csv'))
SYNGAS_TAB['inchi'] = list(map(automol.smiles.inchi, SYNGAS_TAB['smiles']))

#NATGAS_PATH = os.path.join(PATH, 'data/natgas')
#NATGAS_MECH_STR = _read_file(os.path.join(NATGAS_PATH, 'mechanism.txt'))
#NATGAS_TAB = pandas.read_csv(os.path.join(SYNGAS_PATH, 'smiles.csv'))
#NATGAS_TAB['inchi'] = list(map(automol.smiles.inchi, SYNGAS_TAB['smiles']))

#HEPTANE_PATH = os.path.join(PATH, 'data/heptane')
#HEPTANE_MECH_STR = _read_file(os.path.join(HEPTANE_PATH, 'mechanism.txt'))
#HEPTANE_TAB = pandas.read_csv(os.path.join(HEPTANE_PATH, 'smiles.csv'))
#HEPTANE_TAB['inchi'] = list(map(automol.smiles.inchi, HEPTANE_TAB['smiles']))

# PROG = 'psi4'
# SCRIPT_STR = ("#!/usr/bin/env bash\n"
#               "psi4 -i run.inp -o run.out >> stdout.log &> stderr.log")

#mech_str = SYNGAS_MECH_STR
#block_str = chemkin_io.reaction_block(mech_str)
#rxn_strs = chemkin_io.reaction.data_strings(block_str)
#assert len(rxn_strs) == 1678

#rct_names_lst = list(
#    map(chemkin_io.reaction.DataString.reactant_names, rxn_strs))
#prd_names_lst = list(
#    map(chemkin_io.reaction.DataString.product_names, rxn_strs))
#coeffs_lst = list(
#    map(chemkin_io.reaction.DataString.high_p_coefficients, rxn_strs))

# make sure we don't have any None's
#assert all(rct_names_lst)
#assert all(prd_names_lst)
#assert all(coeffs_lst)

# generating inchis for testing elsewhere
MECH_STR = SYNGAS_MECH_STR
BLOCK_STR = chemkin_io.reaction_block(MECH_STR)
RXN_STRS = chemkin_io.reaction.data_strings(BLOCK_STR)
RCT_NAMES_LST = list(
    map(chemkin_io.reaction.DataString.reactant_names, RXN_STRS))
PRD_NAMES_LST = list(
    map(chemkin_io.reaction.DataString.product_names, RXN_STRS))

ICH_DCT = dict(zip(SYNGAS_TAB['name'], SYNGAS_TAB['inchi']))
MULT_DCT = dict(zip(SYNGAS_TAB['name'], SYNGAS_TAB['mult']))
ICH_DCT['OHV'] = None
MULT_DCT['OHV'] = None
ICH_DCT['CHV'] = None
MULT_DCT['CHV'] = None

RCT_ICHS_LST = list(
    list(map(ICH_DCT.__getitem__, names)) for names in RCT_NAMES_LST)
RCT_MULTS_LST = list(
    list(map(MULT_DCT.__getitem__, names)) for names in RCT_NAMES_LST)
PRD_ICHS_LST = list(
    list(map(ICH_DCT.__getitem__, names)) for names in PRD_NAMES_LST)
PRD_MULTS_LST = list(
    list(map(MULT_DCT.__getitem__, names)) for names in PRD_NAMES_LST)

RCT_CHGS_LST = [(0,)*len(rct_ichs) for rct_ichs in RCT_ICHS_LST]
PRD_CHGS_LST = [(0,)*len(prd_ichs) for prd_ichs in PRD_ICHS_LST]
ICHS_LST = list(zip(RCT_ICHS_LST, PRD_ICHS_LST))
CHGS_LST = list(zip(RCT_CHGS_LST, PRD_CHGS_LST))
MULTS_LST = list(zip(RCT_MULTS_LST, PRD_MULTS_LST))
REACTION_LIST = list(zip(ICHS_LST, CHGS_LST, MULTS_LST))
#print(REACTION_LIST)

#REACTION_LIST = (
#    ((('InChI=1S/C2H5O/c1-2-3/h2H2,1H3',),
#      ('InChI=1S/CH2O/c1-2/h1H2', 'InChI=1S/CH3/h1H3')),
#     ((0,), (0, 0)),
#     ((2,), (1, 2))),
#    # ((('InChI=1S/CH2O/c1-2/h1H2', 'InChI=1S/CH3/h1H3'),
#    #   ('InChI=1S/C2H5O/c1-2-3/h2H2,1H3',)),
#    #  ((0, 0), (0,)),
#    #  ((1, 2), (2,))),
#    # ((('InChI=1/CH4O/c1-2/h2H,1H3', 'InChI=1/CHO/c1-2/h1H'),
#    #   ('InChI=1/CH3O/c1-2/h2H,1H2', 'InChI=1/CH2O/c1-2/h1H2')),
#    #  ((0, 0), (0, 0)),
#    #  ((1, 2), (2, 1))),
#)
print(REACTION_LIST)
for rxn_inchis, rxn_charges, rxn_mults in REACTION_LIST:
#   ts_mult = sum(rxn_mults[0]) - len(rxn_mults[0]) + 1
    print(rxn_inchis)
    ts_mult0 = rxn_mults[0][0]
    if len(rxn_mults[0]) == 2:
        ts_mult0 = abs(ts_mult0 - rxn_mults[0][1])+1
    ts_mult1 = rxn_mults[0][1]
    if len(rxn_mults[1]) == 2:
        ts_mult1 = abs(ts_mult1 - rxn_mults[1][1])+1
    ts_mult = max(ts_mult0, ts_mult1)
    print(ts_mult0)
    print(ts_mult1)
    print(ts_mult)
    orb_restricted = (ts_mult == 1)
    direction = autofile.system.reaction_direction(
        rxn_inchis, rxn_charges, rxn_mults)
    rxn_inchis, rxn_charges, rxn_mults = autofile.system.sort_together(
        rxn_inchis, rxn_charges, rxn_mults)
    moldr.driver.run_gridopt(
        rxn_inchis=rxn_inchis,
        rxn_charges=rxn_charges,
        rxn_mults=rxn_mults,
        ts_mult=ts_mult,
        method=METHOD,
        basis=BASIS,
        orb_restricted=orb_restricted,
        run_prefix=RUN_PREFIX,
        save_prefix=SAVE_PREFIX,
        script_str=SCRIPT_STR,
        prog=PROG,
    )
    moldr.driver.save_gridopt(
        rxn_inchis=rxn_inchis,
        rxn_charges=rxn_charges,
        rxn_mults=rxn_mults,
        ts_mult=ts_mult,
        method=METHOD,
        basis=BASIS,
        orb_restricted=orb_restricted,
        run_prefix=RUN_PREFIX,
        save_prefix=SAVE_PREFIX,
    )

sys.exit()

# cycle over input species
for smi, mult in SMILES_MULT_LST:
    ich = automol.smiles.inchi(smi)

    moldr.driver.save_conformers(
        ich=ich,
        charge=0,
        mult=mult,
        method=METHOD,
        basis=BASIS,
        orb_restricted=(mult == 1),
        run_prefix=RUN_PREFIX,
        save_prefix=SAVE_PREFIX,)

    # running monte carlo sampling
    moldr.driver.run_conformers(
        ich=ich,
        charge=0,
        mult=mult,
        method=METHOD,
        basis=BASIS,
        orb_restricted=(mult == 1),
        # run arguments
        nsamp=NSAMP,
        run_prefix=RUN_PREFIX,
        save_prefix=SAVE_PREFIX,
        script_str=SCRIPT_STR,
        prog=PROG,)

    moldr.driver.save_conformers(
        ich=ich,
        charge=0,
        mult=mult,
        method=METHOD,
        basis=BASIS,
        orb_restricted=(mult == 1),
        run_prefix=RUN_PREFIX,
        save_prefix=SAVE_PREFIX,)

    ROOT_SPECS = (ich, 0, mult, METHOD, BASIS, (mult == 1))
    specs_lst = SFS.conf.dir.existing(
        prefix=SAVE_PREFIX,
        root_specs=ROOT_SPECS)

    moldr.driver.run_conformer_job(
        ich=ich,
        charge=0,
        mult=mult,
        method=METHOD,
        basis=BASIS,
        orb_restricted=(mult == 1),
        job='gradient',
        run_prefix=RUN_PREFIX,
        save_prefix=SAVE_PREFIX,
        script_str=SCRIPT_STR,
        prog=PROG,)

    # moldr.driver.run_conformer_job(
    #     ich=ich,
    #     charge=0,
    #     mult=mult,
    #     method=METHOD,
    #     basis=BASIS,
    #     orb_restricted=(mult == 1),
    #     job='hessian',
    #     run_prefix=RUN_PREFIX,
    #     save_prefix=SAVE_PREFIX,
    #     script_str=SCRIPT_STR,
    #     prog=PROG,)

    # running hindered rotor scans for each conformer
    for specs in specs_lst:
        cid, = specs
        moldr.driver.run_scan(
            ich=ich,
            charge=0,
            mult=mult,
            method=METHOD,
            basis=BASIS,
            orb_restricted=(mult == 1),
            cid=cid,
            # run arguments
            run_prefix=RUN_PREFIX,
            save_prefix=SAVE_PREFIX,
            script_str=SCRIPT_STR,
            prog=PROG,
            scan_incr=30.)

        moldr.driver.save_scan(
            ich=ich,
            charge=0,
            mult=mult,
            method=METHOD,
            basis=BASIS,
            orb_restricted=(mult == 1),
            cid=cid,
            # run arguments
            run_prefix=RUN_PREFIX,
            save_prefix=SAVE_PREFIX,)

