""" moldr input
"""
import moldr
import automol

# PROG = 'g09'
# SCRIPT_STR = ("#!/usr/bin/env bash\n"
#               "g09 run.inp run.out")

NSAMP = 3
PROG = 'psi4'
SCRIPT_STR = ("#!/usr/bin/env bash\n"
              "psi4 -i run.inp -o run.out >> stdout.log &> stderr.log")
METHOD = 'hf'
BASIS = 'sto-3g'
RUN_PREFIX = '.'
SAVE_PREFIX = 'save'
GEO = automol.inchi.geometry(automol.smiles.inchi('OCO'))
MULT = 1
CHARGE = 0

moldr.driver.conformers(
    nsamp=NSAMP,
    script_str=SCRIPT_STR,
    run_prefix=RUN_PREFIX,
    save_prefix=SAVE_PREFIX,
    prog=PROG,
    method=METHOD,
    basis=BASIS,
    geo=GEO,
    mult=MULT,
    charge=CHARGE
)
