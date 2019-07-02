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
    # ((('InChI=1/CH4O/c1-2/h2H,1H3', 'InChI=1/CHO/c1-2/h1H'),
    #   ('InChI=1/CH3O/c1-2/h2H,1H2', 'InChI=1/CH2O/c1-2/h1H2')),
    #  ((0, 0), (0, 0)),
    #  ((1, 2), (2, 1))),
)

for rxn_inchis, rxn_charges, rxn_mults in REACTION_LIST:
    ts_mult = sum(rxn_mults[0]) - len(rxn_mults[0]) + 1
    # moldr.driver.run_gridopt(
    #     rxn_inchis=rxn_inchis,
    #     rxn_charges=rxn_charges,
    #     rxn_mults=rxn_mults,
    #     ts_mult=ts_mult,
    #     method=METHOD,
    #     basis=BASIS,
    #     orb_restricted=ORB_RESTRICTED,
    #     run_prefix=RUN_PREFIX,
    #     save_prefix=SAVE_PREFIX,
    #     script_str=SCRIPT_STR,
    #     prog=PROG,
    # )
    moldr.driver.save_gridopt(
        rxn_inchis=rxn_inchis,
        rxn_charges=rxn_charges,
        rxn_mults=rxn_mults,
        ts_mult=ts_mult,
        method=METHOD,
        basis=BASIS,
        orb_restricted=ORB_RESTRICTED,
        run_prefix=RUN_PREFIX,
        save_prefix=SAVE_PREFIX,
    )
