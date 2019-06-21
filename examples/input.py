""" moldr script
"""
import os
import automol
import moldr
from autofile import SFS

# PROG = 'g09'
# SCRIPT_STR = ("#!/usr/bin/env bash\n"
#               "g09 run.inp run.out")

PROG = 'psi4'
SCRIPT_STR = ("#!/usr/bin/env bash\n"
              "psi4 -i run.inp -o run.out >> stdout.log &> stderr.log")

NSAMP = 4
SMILES_MULT_LST = [
    ('CC', 1),
]
#    ('CO[O]', 2),
#    ('OCO', 1),
#    ('O', 1),
METHOD = 'hf'
BASIS = 'sto-3g'
RUN_PREFIX = 'run'
SAVE_PREFIX = 'save'

if not os.path.exists(RUN_PREFIX):
    os.mkdir(RUN_PREFIX)

if not os.path.exists(SAVE_PREFIX):
    os.mkdir(SAVE_PREFIX)

# cycle over input species
for smi, mult in SMILES_MULT_LST:
    ich = automol.smiles.inchi(smi)

# test conformer sampling

# start with save to save any old information in run that wasn't saved
    moldr.driver.save_conformers(
        ich=ich,
        charge=0,
        mult=mult,
        method=METHOD,
        basis=BASIS,
        orb_restricted=(mult == 1),
        run_prefix=RUN_PREFIX,
        save_prefix=SAVE_PREFIX,)

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
#    cid_lst = SFS.conf.dir.existing(RUN_PREFIX,ROOT_SPECS)
    cid_lst = SFS.conf.dir.existing(
        prefix=SAVE_PREFIX,
        root_specs=ROOT_SPECS)
    print(cid_lst)
    for cid in cid_lst:
        print(cid)
        moldr.driver.run_scan(
            ich=ich,
            charge=0,
            mult=mult,
            method=METHOD,
            basis=BASIS,
            orb_restricted=(mult == 1),
            cid=cid[0],
            # run arguments
            run_prefix=RUN_PREFIX,
            save_prefix=SAVE_PREFIX,
            script_str=SCRIPT_STR,
            prog=PROG,
            scan_incr=60.)

# test constrained tau optimization

# start with save to save any old information in run that wasn't saved
#    moldr.driver.save_tau(
#        ich=ich,
#        charge=0,
#        mult=mult,
#        method=METHOD,
#        basis=BASIS,
#        orb_restricted=(mult == 1),
#        run_prefix=RUN_PREFIX,
#        save_prefix=SAVE_PREFIX,)

    moldr.driver.run_tau(
        ich=ich,
        charge=0,
        mult=mult,
        method=METHOD,
        basis=BASIS,
        orb_restricted=(mult == 1),
        nsamp=NSAMP,
        run_prefix=RUN_PREFIX,
        save_prefix=SAVE_PREFIX,
        script_str=SCRIPT_STR,
        prog=PROG,)

    moldr.driver.save_tau(
        ich=ich,
        charge=0,
        mult=mult,
        method=METHOD,
        basis=BASIS,
        orb_restricted=(mult == 1),
        run_prefix=RUN_PREFIX,
        save_prefix=SAVE_PREFIX,)

    print(NSAMP)
    moldr.driver.run_tau_hessian(
        ich=ich,
        charge=0,
        mult=mult,
        method=METHOD,
        basis=BASIS,
        orb_restricted=(mult == 1),
        run_prefix=RUN_PREFIX,
        save_prefix=SAVE_PREFIX,
        script_str=SCRIPT_STR,
        prog=PROG,)
