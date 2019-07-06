""" reaction list test
"""
import os
import sys
import pandas
import chemkin_io
import automol
import autofile
import moldr

# 0. choose which mechanism to run
MECHANISM_NAME = 'syngas'  # options: syngas, natgas, heptane, etc.

# 1. script control parameters
PROG = 'psi4'
SCRIPT_STR = ("#!/usr/bin/env bash\n"
              "psi4 -i run.inp -o run.out >> stdout.log &> stderr.log")
METHOD = 'wb97xd'
BASIS = '6-31g*'
RESTRICT_OPEN_SHELL = False
NSAMP = 5

# 2. create run and save directories
RUN_PREFIX = 'run'
if not os.path.exists(RUN_PREFIX):
    os.mkdir(RUN_PREFIX)

SAVE_PREFIX = 'save'
if not os.path.exists(SAVE_PREFIX):
    os.mkdir(SAVE_PREFIX)

# 3. read in data from the mechanism directory
DATA_PATH = os.path.dirname(os.path.realpath(__file__))
MECH_PATH = os.path.join(DATA_PATH, 'data', MECHANISM_NAME)
MECH_STR = open(os.path.join(MECH_PATH, 'mechanism.txt')).read()
SPC_TAB = pandas.read_csv(os.path.join(MECH_PATH, 'smiles.csv'))

# 4. process species data from the mechanism directory
SPC_TAB['inchi'] = list(map(automol.smiles.inchi, SPC_TAB['smiles']))
SPC_TAB['charge'] = 0
ICH_DCT = dict(zip(SPC_TAB['name'], SPC_TAB['inchi']))
CHG_DCT = dict(zip(SPC_TAB['name'], SPC_TAB['charge']))
MUL_DCT = dict(zip(SPC_TAB['name'], SPC_TAB['mult']))
SPC_BLK_STR = chemkin_io.species_block(MECH_STR)
SPC_NAMES = chemkin_io.species.names(SPC_BLK_STR)

for name in SPC_NAMES:
    # species
    ich = ICH_DCT[name]
    charge = CHG_DCT[name]
    mult = MUL_DCT[name]
    print("Species: {:s}".format(name))

    # theory
    method = METHOD
    basis = BASIS
    if RESTRICT_OPEN_SHELL:
        orb_restr = True
    else:
        orb_restr = (mult == 1)

    # set up the filesystem
    spc_alocs = [ich, charge, mult]         # aloc = absolute locator
    thy_rlocs = [method, basis, orb_restr]  # rloc = relative locator
    thy_alocs = spc_alocs + thy_rlocs

    spc_afs = autofile.fs.species()
    thy_afs = autofile.fs.theory(spc_afs, 'species')

    thy_afs.theory.dir.create(RUN_PREFIX, thy_alocs)
    thy_afs.theory.dir.create(SAVE_PREFIX, thy_alocs)

    run_path = thy_afs.theory.dir.path(RUN_PREFIX, thy_alocs)
    save_path = thy_afs.theory.dir.path(SAVE_PREFIX, thy_alocs)

    # generates the z-matrix and sampling ranges
    geo = automol.inchi.geometry(ich)
    zma = automol.geom.zmatrix(geo)
    tors_names = automol.geom.zmatrix_torsion_coordinate_names(geo)
    tors_ranges = automol.zmatrix.torsional_sampling_ranges(zma, tors_names)
    tors_range_dct = dict(zip(tors_names, tors_ranges))

    moldr.driver.run_conformers(
        zma=zma,
        charge=charge,
        mult=mult,
        method=method,
        basis=basis,
        orb_restr=orb_restr,
        nsamp=NSAMP,
        tors_range_dct=tors_range_dct,
        run_prefix=run_path,
        save_prefix=save_path,
        script_str=SCRIPT_STR,
        prog=PROG,
    )

    print(save_path)
    moldr.driver.save_conformers(
        run_prefix=run_path,
        save_prefix=save_path,
    )

sys.exit()
