""" moldr script
"""
import os
import moldr
import automol

# PROG = 'g09'
# SCRIPT_STR = ("#!/usr/bin/env bash\n"
#               "g09 run.inp run.out")

PROG = 'psi4'
SCRIPT_STR = ("#!/usr/bin/env bash\n"
              "psi4 -i run.inp -o run.out >> stdout.log &> stderr.log")

NSAMP = 2
SMILES_MULT_LST = [
    ('O', 1),
    ('CO[O]', 2),
    ('OCO', 1),
]
METHOD = 'hf'
BASIS = 'sto-3g'
RUN_PREFIX = 'run'
SAVE_PREFIX = 'save'

if not os.path.exists(RUN_PREFIX):
    os.mkdir(RUN_PREFIX)

if not os.path.exists(SAVE_PREFIX):
    os.mkdir(SAVE_PREFIX)

for smi, mult in SMILES_MULT_LST:
    ich = automol.smiles.inchi(smi)

    moldr.driver.run_conformers(
        ich=ich,
        mult=mult,
        method=METHOD,
        basis=BASIS,
        orb_restricted=(mult == 1),
        # run arguments
        run_prefix=RUN_PREFIX,
        save_prefix=SAVE_PREFIX,
        nsamp=NSAMP,
        script_str=SCRIPT_STR,
        prog=PROG,)

    moldr.driver.save_conformers(
        ich=ich,
        mult=mult,
        method=METHOD,
        basis=BASIS,
        orb_restricted=(mult == 1),
        run_prefix=RUN_PREFIX,
        save_prefix=SAVE_PREFIX,)
