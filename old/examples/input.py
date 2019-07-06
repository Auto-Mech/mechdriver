""" moldr script
"""
import os
import automol
import moldr
from autofile import SFS

PROG = 'psi4'
SCRIPT_STR = ("#!/usr/bin/env bash\n"
              "psi4 -i run.inp -o run.out >> stdout.log &> stderr.log")

NSAMP = 4
SMILES_MULT_LST = [
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
