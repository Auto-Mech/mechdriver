""" moldr script
"""
import tempfile
import moldr
import automol
import autodir

PROG = 'g09'
SCRIPT_STR = ("#!/usr/bin/env bash\n"
              "g09 run.inp run.out")

# PROG = 'psi4'
# SCRIPT_STR = ("#!/usr/bin/env bash\n"
#               "psi4 -i run.inp -o run.out >> stdout.log &> stderr.log")

NSAMP = 2
SMILES = ['OO', 'CO[O]']
METHOD = 'hf'
BASIS = 'sto-3g'
RUN_PREFIX = tempfile.mkdtemp()
SAVE_PREFIX = '.'
CHARGE = 0

print(RUN_PREFIX)

# loop over species
for smi in SMILES:
    # determine the multiplicity
    ich = automol.smiles.inchi(smi)
    gra = automol.inchi.connectivity_graph(ich)
    mult = min(automol.graph.possible_spin_multiplicities(gra))
    open_shell = (mult == 1)
    orb_restricted = (not open_shell)

    # get the geometry
    geo = automol.inchi.geometry(ich)

    # create the species directories
    spc_run_prefix = autodir.species.directory_path(RUN_PREFIX, ich, mult)
    autodir.species.create(RUN_PREFIX, ich, mult)

    spc_save_prefix = autodir.species.directory_path(SAVE_PREFIX, ich, mult)
    autodir.species.create(SAVE_PREFIX, ich, mult)

    # create the theory directories
    run_prefix = autodir.theory.directory_path(
        spc_run_prefix, METHOD, BASIS, open_shell, orb_restricted)
    autodir.theory.create(
        spc_run_prefix, METHOD, BASIS, open_shell, orb_restricted)

    save_prefix = autodir.theory.directory_path(
        spc_save_prefix, METHOD, BASIS, open_shell, orb_restricted)
    autodir.theory.create(
        spc_save_prefix, METHOD, BASIS, open_shell, orb_restricted)

    # run the driver
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
