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

# PROG = 'psi4'
# SCRIPT_STR = ("#!/usr/bin/env bash\n"
#               "psi4 -i run.inp -o run.out >> stdout.log &> stderr.log")

PROG = 'g09'
SCRIPT_STR = ("#!/usr/bin/env bash\n"
              "g09 run.inp run.out")
MACHINE_OPTIONS = ['%NProcShared=10']
MEMORY = 10

NSAMP = 10
METHOD = 'wb97xd'
BASIS = '6-31g*'
RUN_PREFIX = 'run'
SAVE_PREFIX = 'save'
RUN_GRADIENT = False
RUN_HESSIAN = False
RUN_CONFORMER_SCAN = False


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


MECH_STR = SYNGAS_MECH_STR
BLOCK_STR = chemkin_io.species_block(MECH_STR)
SPC_NAMES = chemkin_io.species.names(BLOCK_STR)

# generating inchis, smiles, and multiplicities

ICH_DCT = dict(zip(SYNGAS_TAB['name'], SYNGAS_TAB['inchi']))
MULT_DCT = dict(zip(SYNGAS_TAB['name'], SYNGAS_TAB['mult']))
ICH_DCT['OHV'] = None
MULT_DCT['OHV'] = None
ICH_DCT['CHV'] = None
MULT_DCT['CHV'] = None
ICH_LST = list(map(ICH_DCT.__getitem__, SPC_NAMES))
MULT_LST = list(map(MULT_DCT.__getitem__, SPC_NAMES))
ICH_MULT_LST = list(zip(ICH_LST, MULT_LST))
#ICH_TUP = tuple(map(ICH_DCT.__getitem__, SPC_NAMES))
#MULT_TUP = tuple(map(MULT_DCT.__getitem__, SPC_NAMES))
#ICH_MULT_TUP = tuple((zip(ICH_LST, MULT_LST))

# cycle over species
for ich, mult in ICH_MULT_LST:
#for ich, mult in ICH_MULT_TUP:

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
        prog=PROG,
        memory=MEMORY,
        machine_options=MACHINE_OPTIONS,)

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
    ene_lst = [
        SFS.conf.file.energy.read(
            prefix=SAVE_PREFIX,
            specs=ROOT_SPECS+specs)
        for specs in specs_lst]
    min_idx = ene_lst.index(min(ene_lst))
    min_specs = specs_lst[min_idx]

    if RUN_GRADIENT:
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
            prog=PROG,
            memory=MEMORY,
            machine_options=MACHINE_OPTIONS,)

#       moldr.driver.save_conformer_job(
#           ich=ich,
#           charge=0,
#           mult=mult,
#           method=METHOD,
#           basis=BASIS,
#           orb_restricted=(mult == 1),
#           job='gradient',
#           run_prefix=RUN_PREFIX,
#           save_prefix=SAVE_PREFIX,)

    if RUN_HESSIAN:
        moldr.driver.run_conformer_job(
            ich=ich,
            charge=0,
            mult=mult,
            method=METHOD,
            basis=BASIS,
            orb_restricted=(mult == 1),
            job='hessian',
            run_prefix=RUN_PREFIX,
            save_prefix=SAVE_PREFIX,
            script_str=SCRIPT_STR,
            prog=PROG,
            memory=MEMORY,
            machine_options=MACHINE_OPTIONS,)

#       moldr.driver.save_conformer_job(
#           ich=ich,
#           charge=0,
#           mult=mult,
#           method=METHOD,
#           basis=BASIS,
#           orb_restricted=(mult == 1),
#           job='hessian',
#           run_prefix=RUN_PREFIX,
#           save_prefix=SAVE_PREFIX,)

    # running hindered rotor scans for each conformer
#    for specs in specs_lst:
    min_cid, = min_specs
    moldr.driver.run_scan(
        ich=ich,
        charge=0,
        mult=mult,
        method=METHOD,
        basis=BASIS,
        orb_restricted=(mult == 1),
        cid=min_cid,
        # run arguments
        run_prefix=RUN_PREFIX,
        save_prefix=SAVE_PREFIX,
        script_str=SCRIPT_STR,
        prog=PROG,
        memory=MEMORY,
        machine_options=MACHINE_OPTIONS,
        scan_incr=30.)

    moldr.driver.save_scan(
        ich=ich,
        charge=0,
        mult=mult,
        method=METHOD,
        basis=BASIS,
        orb_restricted=(mult == 1),
        cid=min_cid,
        # run arguments
        run_prefix=RUN_PREFIX,
        save_prefix=SAVE_PREFIX,)

# now cycle over reactions
BLOCK_STR = chemkin_io.reaction_block(MECH_STR)
RXN_STRS = chemkin_io.reaction.data_strings(BLOCK_STR)
RCT_NAMES_LST = list(
    map(chemkin_io.reaction.DataString.reactant_names, RXN_STRS))
PRD_NAMES_LST = list(
    map(chemkin_io.reaction.DataString.product_names, RXN_STRS))

COEFFS_LST = list(
    map(chemkin_io.reaction.DataString.high_p_coefficients, RXN_STRS))

# make sure we don't have any None's
assert all(RCT_NAMES_LST)
assert all(PRD_NAMES_LST)
assert all(COEFFS_LST)

# generating inchis, smiles, and multiplicities

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

for rxn_inchis, rxn_charges, rxn_mults in REACTION_LIST:

    ts_mult_hs = sum(rxn_mults[0]) - len(rxn_mults[0]) + 1
    ts_mult0 = rxn_mults[0][0]
    if len(rxn_mults[0]) == 2:
        ts_mult0 = abs(ts_mult0 - rxn_mults[0][1])+1
    ts_mult1 = rxn_mults[0][1]
    if len(rxn_mults[1]) == 2:
        ts_mult1 = abs(ts_mult1 - rxn_mults[1][1])+1
    ts_mult = max(ts_mult0, ts_mult1)
    ts_mult_ls = max(ts_mult0, ts_mult1)
    ts_mult = ts_mult_ls
    orb_restricted = (ts_mult == 1)
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
        memory=MEMORY,
        machine_options=MACHINE_OPTIONS,
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

#sys.exit()
