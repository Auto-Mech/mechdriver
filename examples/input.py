""" moldr script
"""
import os
# import moldr
import automol
from autofile import fs

# PROG = 'g09'
# SCRIPT_STR = ("#!/usr/bin/env bash\n"
#               "g09 run.inp run.out")

PROG = 'psi4'
SCRIPT_STR = ("#!/usr/bin/env bash\n"
              "psi4 -i run.inp -o run.out >> stdout.log &> stderr.log")

NSAMP = 2
SMILES_MULT_LST = [
    ('CO[O]', 2),
    ('OCO', 1),
]
METHOD = 'hf'
BASIS = 'sto-3g'
CHARGE = 0
RUN_PREFIX = 'run'
SAVE_PREFIX = 'save'

if not os.path.exists(RUN_PREFIX):
    os.mkdir(RUN_PREFIX)

if not os.path.exists(SAVE_PREFIX):
    os.mkdir(SAVE_PREFIX)

for smi, mult in SMILES_MULT_LST:
    ich = automol.smiles.inchi(smi)
    orb_restricted = (mult == 1)
    specs = (ich, mult, METHOD, BASIS, orb_restricted)

    print('specs: ', specs)
    fs.theory.dir.create(RUN_PREFIX, specs)
    fs.theory.dir.create(SAVE_PREFIX, specs)
