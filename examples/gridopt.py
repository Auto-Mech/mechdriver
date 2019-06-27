""" moldr script
"""
import os
import moldr

PROG = 'psi4'
SCRIPT_STR = ("#!/usr/bin/env bash\n"
              "psi4 -i run.inp -o run.out >> stdout.log &> stderr.log")

METHOD = 'hf'
BASIS = 'sto-3g'
ORB_RESTRICTED = False
RUN_PREFIX = 'run'
SAVE_PREFIX = 'save'

if not os.path.exists(RUN_PREFIX):
    os.mkdir(RUN_PREFIX)

if not os.path.exists(SAVE_PREFIX):
    os.mkdir(SAVE_PREFIX)

REACTION_LIST = (
    ((('InChI=1S/C2H5O/c1-2-3/h2H2,1H3',),
      ('InChI=1S/CH2O/c1-2/h1H2', 'InChI=1S/CH3/h1H3')),
     ((0,), (0, 0)),
     ((2,), (1, 2))),
    # ((('InChI=1S/CH2O/c1-2/h1H2', 'InChI=1S/CH3/h1H3'),
    #   ('InChI=1S/C2H5O/c1-2-3/h2H2,1H3',)),
    #  ((0, 0), (0,)),
    #  ((1, 2), (2,))),
)

for inchis_pair, charges_pair, mults_pair in REACTION_LIST:
    moldr.driver.run_gridopt(
        inchis_pair=inchis_pair,
        charges_pair=charges_pair,
        mults_pair=mults_pair,
        method=METHOD,
        basis=BASIS,
        orb_restricted=ORB_RESTRICTED,
        run_prefix=RUN_PREFIX,
        save_prefix=SAVE_PREFIX,
        script_str=SCRIPT_STR,
        prog=PROG,
        ts_mult=sum(mults_pair[0]) - len(mults_pair[0]) + 1,
    )
