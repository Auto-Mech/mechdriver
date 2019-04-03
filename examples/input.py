""" moldr script
"""
import os
import moldr
import automol
import autodir

# PROG = 'g09'
# SCRIPT_STR = ("#!/usr/bin/env bash\n"
#               "g09 run.inp run.out")

PROG = 'psi4'
SCRIPT_STR = ("#!/usr/bin/env bash\n"
              "psi4 -i run.inp -o run.out >> stdout.log &> stderr.log")

NSAMP = 2
SMILES_MULTS = [
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

# loop over species
for smi, mult in SMILES_MULTS:
    open_shell = (mult == 1)
    orb_restricted = (not open_shell)

    # get the geometry
    ich = automol.smiles.inchi(smi)
    geo = automol.inchi.geometry(ich)

    # create the species directories
    run_prefix = autodir.spth.directory_path(
        RUN_PREFIX, ich, mult, METHOD, BASIS, open_shell, orb_restricted)
    print(run_prefix)
    autodir.spth.create(
        RUN_PREFIX, ich, mult, METHOD, BASIS, open_shell, orb_restricted)

    save_prefix = autodir.spth.directory_path(
        SAVE_PREFIX, ich, mult, METHOD, BASIS, open_shell, orb_restricted)
    print(save_prefix)
    autodir.spth.create(
        SAVE_PREFIX, ich, mult, METHOD, BASIS, open_shell, orb_restricted)

    # run the driver for determining conformers via torsional sampling
    moldr.driver.conformers(
        nsamp=NSAMP,
        script_str=SCRIPT_STR,
        run_prefix=run_prefix,
        save_prefix=save_prefix,
        prog=PROG,
        method=METHOD,
        basis=BASIS,
        geo=geo,
        mult=mult,
        charge=CHARGE
    )

    # add gradients
    moldr.driver.add_conformer_gradients(
        script_str=SCRIPT_STR,
        run_prefix=run_prefix,
        save_prefix=save_prefix,
        prog=PROG,
        method=METHOD,
        basis=BASIS,
        mult=mult,
        charge=CHARGE
    )

    # add hessians
    moldr.driver.add_conformer_hessians(
        script_str=SCRIPT_STR,
        run_prefix=run_prefix,
        save_prefix=save_prefix,
        prog=PROG,
        method=METHOD,
        basis=BASIS,
        mult=mult,
        charge=CHARGE
    )
