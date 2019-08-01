""" reaction list test
"""
import os
import sys
import numpy
import pandas
from qcelemental import constants as qcc
import chemkin_io
import projrot_io
import thermo
import automol
import elstruct
import autofile
import moldr
import mess_io.writer

ANG2BOHR = qcc.conversion_factor('angstrom', 'bohr')
WAVEN2KCAL = qcc.conversion_factor('wavenumber', 'kcal/mol')
EH2KCAL = qcc.conversion_factor('hartree', 'kcal/mol')

# 0. choose which mechanism to run

# MECHANISM_NAME = 'ch4+nh2'  # options: syngas, natgas, heptane
MECHANISM_NAME = 'test'  # options: syngas, natgas, heptane
# MECHANISM_NAME = 'estoktp/add30'  # options: syngas, natgas
# MECHANISM_NAME = 'estoktp/habs65'  # options: syngas, natgas

# 1. script control parameters

# a. Strings to launch executable
# script_strings for electronic structure are obtained from run_qchem_par since
# they vary with method

PROJROT_SCRIPT_STR = ("#!/usr/bin/env bash\n"
                      "RPHt.exe")
PF_SCRIPT_STR = ("#!/usr/bin/env bash\n"
                 "messpf pf.inp build.out >> stdout.log &> stderr.log")
RATE_SCRIPT_STR = ("#!/usr/bin/env bash\n"
                   "mess build.inp build.out >> stdout.log &> stderr.log")
NASA_SCRIPT_STR = ("#!/usr/bin/env bash\n"
                   "cp ../PF/build.out pf.dat\n"
                   "cp /tcghome/sjklipp/PACC/nasa/new.groups .\n"
                   "python /tcghome/sjklipp/PACC/nasa/makepoly.py"
                   " >> stdout.log &> stderr.log")

# b. Electronic structure parameters; code, method, basis, convergence control

# use reference to determine starting geometries, which are used to define
# z-matrices also used for reference geometries in HL calculations
REF_LEVEL = ['wb97xd', '6-31g*']

# use indices to easily choose between a set of standard opts to call
IDX_OPTS = [0, 3]
# set up standard prog, method, basis for opts

OPT_LEVELS = []
OPT_LEVELS.append(['psi4', 'wb97xd', '6-31g*'])
# OPT_LEVELS.append(['psi4', 'wb97xd', 'cc-pVTZ'])
# OPT_LEVELS.append(['psi4', 'b2plypd3', 'cc-pVTZ'])
# OPT_LEVELS.append(['psi4', 'wb97xd', 'cc-pVTZ'])
# OPT_LEVELS.append(['psi4', 'm062x', 'cc-pVTZ'])
# OPT_LEVELS.append(['psi4', 'b3lyp', '6-31g*'])
RUN_OPT_LEVELS = [OPT_LEVELS[0]]

# set up a set of standard hl methods
HIGH_LEVEL_REF = ['wb97xd', '6-31g*']
HIGH_LEVELS = []
HIGH_LEVELS.append(['psi4', 'mp2', 'cc-pVTZ'])
# HIGH_LEVELS.append(['psi4', 'CCSD(T)', 'cc-pVTZ'])
# HIGH_LEVELS.append(['psi4', 'CCSD(T)', 'cc-pVQZ'])
# HIGH_LEVELS.append(['psi4', 'CCSD(T)', 'cc-pV5Z'])
# HIGH_LEVELS.append(['psi4', 'CCSD(T)-F12', 'cc-pVDZ-F12'])
# HIGH_LEVELS.append(['psi4', 'CCSD(T)-F12', 'cc-pVTZ-F12'])
# HIGH_LEVELS.append(['psi4', 'CCSD(T)-F12', 'cc-pVQZ-F12'])
# HIGH_LEVELS.append(['psi4', 'CCSD(T)', 'aug-cc-pVDZ'])
# HIGH_LEVELS.append(['psi4', 'CCSD(T)', 'aug-cc-pVTZ'])
# HIGH_LEVELS.append(['psi4', 'CCSD(T)', 'aug-cc-pVQZ'])
# HIGH_LEVELS.append(['psi4', 'CCSD(T)', 'aug-cc-pV5Z'])
RUN_HIGH_LEVELS = [HIGH_LEVELS[0]]

# c. What type of electronic structure calculations to run
RUN_SPECIES_QCHEM = True
RUN_REACTIONS_QCHEM = True
RUN_TS_KICKS_QCHEM = True
RUN_VDW_QCHEM = True

RUN_INI_GEOM = True
RUN_REMOVE_IMAG = True

KICKOFF_SADDLE = True

RUN_CONF_SAMP = True
RUN_CONF_OPT = True
RUN_MIN_GRAD = True
RUN_MIN_HESS = True
RUN_CONF_GRAD = True
RUN_CONF_HESS = True

RUN_CONF_SCAN = True
RUN_CONF_SCAN_GRAD = False
RUN_CONF_SCAN_HESS = False

RUN_TAU_SAMP = True
RUN_TAU_GRAD = True
RUN_TAU_HESS = True

RUN_TS_CONF_OPT = True
RUN_TS_CONF_SCAN = True
RUN_TS_TAU_SAMP = True

RUN_HL_MIN_ENE = True

# setting these to true turns on corresponding run for min, conf, conf_scan,
# and tau
RUN_GRAD = True
RUN_HESS = True
if RUN_GRAD:
    RUN_MIN_GRAD = True
    RUN_CONF_GRAD = True
    RUN_CONF_SCAN_GRAD = True
    RUN_TAU_GRAD = True

if RUN_HESS:
    RUN_MIN_HESS = True
    RUN_CONF_HESS = True
    RUN_CONF_SCAN_HESS = True
    RUN_TAU_HESS = True
RUN_GRAD_PF = False
RUN_HESS_PF = False

# d. Parameters for number of torsional samplings
NSAMP_CONF = 5
NSAMP_CONF_EXPR = True
NSAMP_CONF_A = 3
NSAMP_CONF_B = 1
NSAMP_CONF_C = 3
NSAMP_CONF_D = 100
NSAMP_CONF_PAR = [NSAMP_CONF_EXPR, NSAMP_CONF_A, NSAMP_CONF_B, NSAMP_CONF_C,
                  NSAMP_CONF_D, NSAMP_CONF]

# NSAMP_TAU = 100
NSAMP_TAU = 10
NSAMP_TAU_EXPR = False
NSAMP_TAU_A = 3
NSAMP_TAU_B = 1
NSAMP_TAU_C = 3
NSAMP_TAU_D = 15
NSAMP_TAU_PAR = [NSAMP_TAU_EXPR, NSAMP_TAU_A, NSAMP_TAU_B, NSAMP_TAU_C,
                 NSAMP_TAU_D, NSAMP_TAU]

NSAMP_VDW = 10
NSAMP_VDW_EXPR = False
NSAMP_VDW_A = 3
NSAMP_VDW_B = 1
NSAMP_VDW_C = 3
NSAMP_VDW_D = 15
NSAMP_VDW_PAR = [NSAMP_VDW_EXPR, NSAMP_VDW_A, NSAMP_VDW_B, NSAMP_VDW_C,
                 NSAMP_VDW_D, NSAMP_VDW]

# e. What to run for thermochemical kinetics
RUN_SPECIES_PF = False
RUN_SPECIES_THERMO = False
RUN_REACTIONS_RATES = False
RUN_VDW_RCT_RATES = False
RUN_VDW_PRD_RATES = False

# f. Partition function parameters
TAU_PF_WRITE = True
SPECIES_RRHO_STR = {}
SPECIES_HR_STR = {}

# Defaults
SCAN_INCREMENT = 30. * qcc.conversion_factor('degree', 'radian')
KICKOFF_SIZE = 0.1
KICKOFF_BACKWARD = False
RESTRICT_OPEN_SHELL = False
OVERWRITE = False
# Temperatue and pressure
TEMPS = [300., 500., 750., 1000., 1250., 1500., 1750., 2000.]
TEMP_STEP = 100.
NTEMPS = 30
PRESS = [0.1, 1., 10., 100.]
# Collisional parameters
EXP_FACTOR = 150.0
EXP_POWER = 50.0
EXP_CUTOFF = 80.0
EPS1 = 100.0
EPS2 = 200.0
SIG1 = 10.0
SIG2 = 20.0
MASS1 = 15.0
MASS2 = 25.0

ELC_DEG_DCT = {
    ('InChI=1S/B', 2): [[0., 2], [16., 4]],
    ('InChI=1S/C', 3): [[0., 1], [16.4, 3], [43.5, 5]],
    ('InChI=1S/N', 2): [[0., 6], [8., 4]],
    ('InChI=1S/O', 3): [[0., 5], [158.5, 3], [226.5, 1]],
    ('InChI=1S/F', 2): [[0., 4], [404.1, 2]],
    ('InChI=1S/Cl', 2): [[0., 4], [883.4, 2]],
    ('InChI=1S/Br', 2): [[0., 4], [685.2, 2]],
    ('InChI=1S/HO/h1H', 2): [[0., 2], [138.9, 2]],
    ('InChI=1S/NO/c1-2', 2): [[0., 2], [123.1, 2]],
    ('InChI=1S/O2/c1-2', 1): [[0., 2]]
}
REF_REF = ['[H],[H]', 'C', 'O', 'N']
CHG_REF = [0, 0, 0, 0]
MULT_REF = [1, 1, 1, 1]

ELC_SIG_LST = {'InChI=1S/CN/c1-2'}
# add CCH to this list

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

# 4. process species data from the mechanism file

SPC_TAB['charge'] = 0
SMI_DCT = dict(zip(SPC_TAB['name'], SPC_TAB['smiles']))
# SMI_DCT['REF_H2'] = '[H][H]'
SMI_DCT['REF_CH4'] = 'C'
SMI_DCT['REF_H2O'] = 'O'
SMI_DCT['REF_NH3'] = 'N'
CHG_DCT = dict(zip(SPC_TAB['name'], SPC_TAB['charge']))
# CHG_DCT['REF_H2'] = 0
CHG_DCT['REF_CH4'] = 0
CHG_DCT['REF_H2O'] = 0
CHG_DCT['REF_NH3'] = 0
MUL_DCT = dict(zip(SPC_TAB['name'], SPC_TAB['mult']))
# MUL_DCT['REF_H2'] = 1
MUL_DCT['REF_CH4'] = 1
MUL_DCT['REF_H2O'] = 1
MUL_DCT['REF_NH3'] = 1
SPC_BLK_STR = chemkin_io.species_block(MECH_STR)
SPC_NAMES = chemkin_io.species.names(SPC_BLK_STR)
# SPC_NAMES += ('REF_H2', 'REF_CH4', 'REF_H2O', 'REF_NH3')
# SPC_NAMES += ('REF_CH4', 'REF_H2O', 'REF_NH3')
print('SPC_NAMES')
print(SPC_NAMES)

GEOM_PATH = os.path.join(DATA_PATH, 'data', 'geoms')
print(GEOM_PATH)
GEOM_DCT = moldr.util.geometry_dictionary(GEOM_PATH)

# take starting geometry from saved directory if possible, otherwise get it
# from inchi via rdkit

if RUN_SPECIES_QCHEM:
    for name in SPC_NAMES:
        # species
        print("Species: {}".format(name))
        smi = SMI_DCT[name]
        ich = automol.smiles.inchi(smi)
        print("smiles: {}".format(smi), "inchi: {}".format(ich))
        # ich = ICH_DCT[name]
        chg = CHG_DCT[name]
        mul = MUL_DCT[name]
        orb_restr = moldr.util.orbital_restriction(mul, RESTRICT_OPEN_SHELL)
        method_ref = REF_LEVEL[0]
        basis_ref = REF_LEVEL[1]
        for prog, method, basis in RUN_OPT_LEVELS:
            # theory
            SCRIPT_STR, OPT_SCRIPT_STR, KWARGS, OPT_KWARGS = (
                moldr.util.run_qchem_par(prog))

            # a. conformer sampling
            spc_run_fs = autofile.fs.species(RUN_PREFIX)
            spc_run_fs.leaf.create([ich, chg, mul])
            spc_run_path = spc_run_fs.leaf.path([ich, chg, mul])

            spc_save_fs = autofile.fs.species(SAVE_PREFIX)
            spc_save_fs.leaf.create([ich, chg, mul])
            spc_save_path = spc_save_fs.leaf.path([ich, chg, mul])

            thy_run_fs = autofile.fs.theory(spc_run_path)
            thy_run_fs.leaf.create([method, basis, orb_restr])
            thy_run_path = thy_run_fs.leaf.path([method, basis, orb_restr])

            thy_save_fs = autofile.fs.theory(spc_save_path)
            thy_save_fs.leaf.create([method, basis, orb_restr])
            thy_save_path = thy_save_fs.leaf.path([method, basis, orb_restr])

            geo_init = moldr.util.reference_geometry(
                ich=ich,
                chg=chg,
                mul=mul,
                method=method_ref,
                basis=basis_ref,
                orb_restr=orb_restr,
                prefix=SAVE_PREFIX,
                geom_dct=GEOM_DCT)

            if RUN_INI_GEOM:
                geo = moldr.driver.run_initial_geometry_opt(
                    ich=ich,
                    chg=chg,
                    mul=mul,
                    method=method,
                    basis=basis,
                    orb_restr=orb_restr,
                    run_prefix=spc_run_path,
                    save_prefix=spc_save_path,
                    script_str=SCRIPT_STR,
                    prog=prog,
                    overwrite=OVERWRITE,
                    geo_init=geo_init,
                    **OPT_KWARGS,
                )

                if RUN_REMOVE_IMAG:
                    moldr.driver.run_remove_imaginary(
                        ich=ich,
                        chg=chg,
                        mul=mul,
                        method=method,
                        basis=basis,
                        orb_restr=orb_restr,
                        run_prefix=spc_run_path,
                        script_str=SCRIPT_STR,
                        prog=prog,
                        overwrite=OVERWRITE,
                        kickoff_backward=KICKOFF_BACKWARD,
                        kickoff_size=KICKOFF_SIZE,
                        **KWARGS
                    )

                moldr.driver.save_initial_geometry(
                    ich=ich,
                    chg=chg,
                    mul=mul,
                    method=method,
                    basis=basis,
                    orb_restr=orb_restr,
                    run_prefix=spc_run_path,
                    save_prefix=spc_save_path,
                    prog=prog,
                )

            if RUN_CONF_SAMP:
                moldr.driver.conformer_sampling(
                    ich=ich,
                    chg=chg,
                    mul=mul,
                    method=method,
                    basis=basis,
                    orb_restr=orb_restr,
                    run_prefix=spc_run_path,
                    save_prefix=spc_save_path,
                    script_str=SCRIPT_STR,
                    prog=prog,
                    overwrite=OVERWRITE,
                    nsamp_par=NSAMP_CONF_PAR,
                    **OPT_KWARGS,
                )

                if RUN_MIN_GRAD:
                    moldr.driver.run_minimum_energy_gradient(
                        ich=ich,
                        chg=chg,
                        mul=mul,
                        method=method,
                        basis=basis, orb_restr=orb_restr,
                        run_prefix=thy_run_path,
                        save_prefix=thy_save_path,
                        script_str=SCRIPT_STR,
                        prog=prog,
                        overwrite=OVERWRITE,
                        **KWARGS,
                    )

                if RUN_MIN_HESS:
                    moldr.driver.run_minimum_energy_hessian(
                        ich=ich,
                        chg=chg,
                        mul=mul,
                        method=method,
                        basis=basis,
                        orb_restr=orb_restr,
                        run_prefix=thy_run_path,
                        save_prefix=thy_save_path,
                        script_str=SCRIPT_STR,
                        prog=prog,
                        overwrite=OVERWRITE,
                        **KWARGS,
                    )

                if RUN_CONF_SCAN:
                    moldr.driver.hindered_rotor_scans(
                        ich=ich,
                        chg=chg,
                        mul=mul,
                        method=method,
                        basis=basis,
                        orb_restr=orb_restr,
                        run_prefix=thy_run_path,
                        save_prefix=thy_save_path,
                        script_str=SCRIPT_STR,
                        prog=prog,
                        overwrite=OVERWRITE,
                        scan_increment=SCAN_INCREMENT,
                        **OPT_KWARGS,
                    )

                if RUN_CONF_GRAD:
                    moldr.driver.run_conformer_gradients(
                        ich=ich,
                        chg=chg,
                        mul=mul,
                        method=method,
                        basis=basis,
                        orb_restr=orb_restr,
                        run_prefix=thy_run_path,
                        save_prefix=thy_save_path,
                        script_str=SCRIPT_STR,
                        prog=prog,
                        overwrite=OVERWRITE,
                        **KWARGS,
                    )

                if RUN_CONF_HESS:
                    moldr.driver.run_conformer_hessians(
                        ich=ich,
                        chg=chg,
                        mul=mul,
                        method=method,
                        basis=basis,
                        orb_restr=orb_restr,
                        run_prefix=thy_run_path,
                        save_prefix=thy_save_path,
                        script_str=SCRIPT_STR,
                        prog=prog,
                        overwrite=OVERWRITE,
                        **KWARGS,
                    )

            if RUN_TAU_SAMP:
                moldr.driver.tau_sampling(
                    ich=ich,
                    chg=chg,
                    mul=mul,
                    method=method,
                    basis=basis,
                    orb_restr=orb_restr,
                    run_prefix=spc_run_path,
                    save_prefix=spc_save_path,
                    script_str=SCRIPT_STR,
                    prog=prog,
                    overwrite=OVERWRITE,
                    nsamp_par=NSAMP_TAU_PAR,
                    **OPT_KWARGS,
                )

                if RUN_TAU_GRAD:
                    moldr.driver.run_tau_gradients(
                        ich=ich,
                        chg=chg,
                        mul=mul,
                        method=method,
                        basis=basis,
                        orb_restr=orb_restr,
                        run_prefix=thy_run_path,
                        save_prefix=thy_save_path,
                        script_str=SCRIPT_STR,
                        prog=prog,
                        overwrite=OVERWRITE,
                        **KWARGS,
                    )

                if RUN_TAU_HESS:
                    moldr.driver.run_tau_hessians(
                        ich=ich,
                        chg=chg,
                        mul=mul,
                        method=method,
                        basis=basis,
                        orb_restr=orb_restr,
                        run_prefix=thy_run_path,
                        save_prefix=thy_save_path,
                        script_str=SCRIPT_STR,
                        prog=prog,
                        overwrite=OVERWRITE,
                        **KWARGS,
                    )

            if TAU_PF_WRITE:
                moldr.driver.tau_pf_write(
                    name=name,
                    ich=ich,
                    chg=chg,
                    mul=mul,
                    method=method,
                    basis=basis,
                    orb_restr=orb_restr,
                    save_prefix=thy_save_path,
                    run_grad=RUN_GRAD_PF,
                    run_hess=RUN_HESS_PF,
                    **KWARGS,
                )

        for prog, method, basis, in RUN_HIGH_LEVELS:
            SCRIPT_STR, OPT_SCRIPT_STR, KWARGS, OPT_KWARGS = (
                moldr.util.run_qchem_par(prog))
            method_ref = HIGH_LEVEL_REF[0]
            basis_ref = HIGH_LEVEL_REF[1]

            min_cnf_locs = moldr.util.min_energy_conformer_locators(
                thy_save_path)

            cnf_run_fs = autofile.fs.conformer(thy_run_path)
            cnf_run_path = cnf_run_fs.leaf.path(min_cnf_locs)

            cnf_save_fs = autofile.fs.conformer(thy_save_path)
            cnf_save_path = cnf_save_fs.leaf.path(min_cnf_locs)
            min_cnf_geo = cnf_save_fs.leaf.file.geometry.read(min_cnf_locs)
            print('HL test')
            if RUN_HL_MIN_ENE:
                moldr.driver.run_single_point_energy(
                    geo=min_cnf_geo,
                    chg=chg,
                    mul=mul,
                    method=method,
                    basis=basis,
                    orb_restr=orb_restr,
                    run_prefix=thy_run_path,
                    save_prefix=thy_save_path,
                    script_str=SCRIPT_STR,
                    prog=prog,
                    overwrite=OVERWRITE,
                    **KWARGS,
                )

# everything through here is converted over

if RUN_SPECIES_PF:
    for name in SPC_NAMES:
        for prog, method, basis in RUN_OPT_LEVELS:
            # set up species information
            smi = SMI_DCT[name]
            ich = automol.smiles.inchi(smi)
            print("smiles: {}".format(smi), "inchi: {}".format(ich))
            chg = CHG_DCT[name]
            mult = MUL_DCT[name]
            orb_restr = moldr.util.orbital_restriction(
                mult, RESTRICT_OPEN_SHELL)
            # read in geometry, hessian and hindered rotor potentials for
            # minimum energy conformer
            spc_save_fs = autofile.fs.species(SAVE_PREFIX)
            spc_save_path = spc_save_fs.leaf.path([ich, chg, mul])

            thy_save_fs = autofile.fs.theory(spc_save_path)
            thy_save_path = thy_save_fs.leaf.path([method, basis, orb_restr])

            min_cnf_locs = moldr.util.min_energy_conformer_locators(
                thy_save_path)
            cnf_save_fs = autofile.fs.conformer(thy_save_path)
            cnf_save_path = cnf_save_fs.leaf.path(min_cnf_locs)

            # I think we need something for if it is none
            if min_cnf_locs is not None:
                geo = cnf_save_fs.leaf.file.geometry.read(min_cnf_locs)
                grad = cnf_save_fs.leaf.file.gradient.read(min_cnf_locs)
                hess = cnf_save_fs.leaf.file.hessian.read(min_cnf_locs)
                freqs = elstruct.util.harmonic_frequencies(geo, hess)
                zpe = sum(freqs)*WAVEN2KCAL/2.
                zma = automol.geom.zmatrix(geo)
                gra = automol.zmatrix.graph(zma, remove_stereo=True)
                min_ene = cnf_save_fs.leaf.file.energy.read(min_cnf_locs)
                cnf_save_path = cnf_save_fs.leaf.path(min_cnf_locs)
                tors_names = automol.geom.zmatrix_torsion_coordinate_names(geo)
                coo_dct = automol.zmatrix.coordinates(zma, multi=False)
                hind_rot_str = ""
                rotors_str = ""

                # prepare axis, group, and projection info
                scn_save_fs = autofile.fs.scan(cnf_save_path)
                for tors_name in tors_names:
                    enes = [scn_save_fs.leaf.file.energy.read(locs) for locs
                            in scn_save_fs.leaf.existing([[tors_name]])]
                    enes = numpy.subtract(enes, min_ene)
                    pot = list(enes*EH2KCAL)
                    axis = coo_dct[tors_name][1:3]
                    group = list(
                        automol.graph.branch_atom_keys(gra, axis[0], axis) -
                        set(axis))
                    group = list(numpy.add(group, 1))
                    axis = list(numpy.add(axis, 1))
                    sym = 1
                    hind_rot_str += mess_io.writer.write_rotor_hindered(
                        group, axis, sym, pot)
                    rotors_str += projrot_io._write.write_rotors_str(
                        axis, group)

                print('hind_rot_str test')
                print(hind_rot_str)

                # Write the string for the ProjRot input
                COORD_PROJ = 'cartesian'
                print('grad')
                print(grad)
                print('hess')
                print(hess)
                projrot_inp_str = projrot_io._write.write_rpht_input(
                    geo, grad, hess, rotors_str=rotors_str,
                    coord_proj=COORD_PROJ)

                bld_locs = ['PROJROT', 0]
                bld_save_fs = autofile.fs.build(thy_save_path)
                bld_save_fs.leaf.create(bld_locs)
                path = bld_save_fs.leaf.path(bld_locs)
                print('Build Path for Partition Functions')
                print(path)
                proj_file_path = os.path.join(path, 'RPHt_input_data.dat')
                with open(proj_file_path, 'w') as proj_file:
                    proj_file.write(projrot_inp_str)

                moldr.util.run_script(PROJROT_SCRIPT_STR, path)

                rtproj_freqs, _ = projrot_io._read.read_rpht_output(
                    path+'/RTproj_freq.dat')
                rthrproj_freqs, _ = projrot_io._read.read_rpht_output(
                    path+'/hrproj_freq.dat')
                # the second variable above is the imaginary frequency list
                print('Projection test')
                print(rtproj_freqs)
                print(rthrproj_freqs)

                # set up messpf input
                elec_levels = [[0., mult]]
                if (ich, mult) in ELC_DEG_DCT:
                    elec_levels = ELC_DEG_DCT[(ich, mult)]
                # to be generalized
                symfactor = 1.

                # create a messpf input file
                temp_step = TEMP_STEP
                ntemps = NTEMPS
                global_pf_str = mess_io.writer.write_global_pf(
                    [], temp_step, ntemps, rel_temp_inc=0.001,
                    atom_dist_min=0.6)
                print(global_pf_str)
                species_head_str = 'Species ' + name
                print(species_head_str)
                if automol.geom.is_atom(geo):
                    print('This is an atom')
                    SPECIES_RRHO_STR[name] = mess_io.writer.write_atom(
                        name, elec_levels)
                else:
                    core = mess_io.writer.write_core_rigidrotor(geo, symfactor)
                    SPECIES_RRHO_STR[name] = mess_io.writer.write_molecule(
                        core, rtproj_freqs, zpe, elec_levels,
                        hind_rot='',
                        )
                    pf_rrho_inp_str = '\n'.join(
                        [global_pf_str, species_head_str,
                         SPECIES_RRHO_STR[name]])
                    print(SPECIES_RRHO_STR[name])

                    if pot is not None:
                        SPECIES_HR_STR[name] = mess_io.writer.write_molecule(
                            core, rthrproj_freqs, zpe, elec_levels,
                            hind_rot=hind_rot_str,
                        )
                    pf_hr_inp_str = '\n'.join(
                        [global_pf_str, species_head_str,
                         SPECIES_HR_STR[name]])
                    print(SPECIES_HR_STR[name])

                bld_locs = ['PF', 0]
                bld_save_fs = autofile.fs.build(thy_save_path)
                bld_save_fs.leaf.create(bld_locs)
                path = bld_save_fs.leaf.path(bld_locs)
                print('Build Path for Partition Functions')
                print(path)

                with open(os.path.join(path, 'pf_rrho.inp'), 'w') as pf_file:
                    pf_file.write(pf_rrho_inp_str)
                pf_script_str = PF_SCRIPT_STR.replace('pf.inp', 'pf_rrho.inp')
                moldr.util.run_script(pf_script_str, path)

                if hind_rot_str != '':
                    pf_1dhr_inp_str = '\n'.join(
                        [global_pf_str, species_head_str,
                         SPECIES_HR_STR[name]])
                    pf_file_path = os.path.join(path, 'pf_1dhr.inp')
                    with open(pf_file_path, 'w') as pf_file:
                        pf_file.write(pf_1dhr_inp_str)

                    pf_script_str = PF_SCRIPT_STR.replace(
                        'pf.inp', 'pf_1dhr.inp')
                    moldr.util.run_script(pf_script_str, path)

# AVC: didn't fix this section
if RUN_SPECIES_THERMO:
    for name in SPC_NAMES:
        # set up species information
        smi = SMI_DCT[name]
        ich = automol.smiles.inchi(smi)
        print("smiles: {}".format(smi), "inchi: {}".format(ich))
        chg = CHG_DCT[name]
        for prog, method, basis in RUN_OPT_LEVELS:
            mult = MUL_DCT[name]
            orb_restr = moldr.util.orbital_restriction(mult, RESTRICT_OPEN_SHELL)
            # read in geometry, hess and hindered rotor potentials for minimum energy conformer
            spc_save_path = moldr.util.species_path(ich, chg, mult, SAVE_PREFIX)
            thy_save_path = moldr.util.theory_path(method, basis, orb_restr, spc_save_path)
            spc_save_path = moldr.util.species_path(ich, chg, mult, SAVE_PREFIX)
            thy_save_path = moldr.util.theory_path(method, basis, orb_restr, spc_save_path)


            bld_locs = ['PF', 0]
            bld_afs.build.dir.create(thy_save_path, bld_locs)
            pf_path = bld_afs.build.dir.path(thy_save_path, bld_locs)
            print('pf build path')
            print(pf_path)

            nasa_inp_str = ('nasa')
            bld_afs = autofile.fs.build()
            bld_locs = ['NASA_POLY', 0]
            bld_afs.build.dir.create(thy_save_path, bld_locs)
            nasa_path = bld_afs.build.dir.path(thy_save_path, bld_locs)
            print('NASA build path')
            print(path)

            cnf_afs = autofile.fs.conformer()
            min_cnf_locs = moldr.util.min_energy_conformer_locators(thy_save_path)
    # I think we need something for if it is none
            if min_cnf_locs is not None:
                min_ene = cnf_afs.conf.file.energy.read(thy_save_path, min_cnf_locs)

            formula = thermo.util.inchi_formula(ich)
            print('\nformula:')
            print(formula)

            # Get atom count dictionary
            atom_dict = thermo.util.get_atom_counts_dict(formula)
            print('\natom dict:')
            print(atom_dict)

            # Get the list of the basis
            basis = thermo.heatform.select_basis(atom_dict)
            print('\nbasis:')
            print(basis)

            # Get the coefficients for the balanced heat-of-formation eqn
            coeff = thermo.heatform.calc_coefficients(basis, atom_dict)
            print('\ncoeff:')
            print(coeff)

            # Get the energy for each of the basis species
            e_basis = thermo.heatform.get_basis_energy(basis)
            print('\ne_basis:')
            print(e_basis)

            # Get the 0 K heat of formation
            h0form = thermo.heatform.calc_hform_0k(min_ene, e_basis, coeff)

            os.chdir(path)

            # Write thermp input file
            ENTHALPYT = 0.
            BREAKT = 1000.
            thermo.runner.write_thermp_input(
                    formula=formula,
                    deltaH=h0form,
                    enthalpyT=ENTHALPYT,
                    breakT=BREAKT,
                    thermp_file_name='thermp.dat')

            PF_TYPES = ['pf_rrho.dat', 'pf_1dhr.dat']
            for pf_type in PF_TYPES:

                # Run thermp
                thermo.runner.run_thermp(
                        pf_path=pf_path,
                        thermp_path=path,
                        thermp_file_name='thermp.dat',
                        pf_file_name=pf_type
                        )

                # Run pac99
                print('formula test')
                print(formula)
                print(path)
                FORMULA = 'CH4O'
                thermo.runner.run_pac99(path, FORMULA)
#                thermo.runner.run_pac99(path, formula)

                with open(os.path.join(path, 'thermp.out'), 'r') as thermp_outfile:
                    thermp_out_str = thermp_outfile.read()

                # Get the 298 K heat of formation
                h298form = thermo.heatform.get_hform_298k_thermp(thermp_out_str)
                print('\nhform(298 K):')
                print(h298form)

                with open(os.path.join(path, formula+'.o97'), 'r') as pac99_file:
                    pac99_str = pac99_file.read()

                # Get the pac99 polynomial
                pac99_poly_str = thermo.nasapoly.get_pac99_polynomial(pac99_str)
                print('\nPAC99 Polynomial:')
                print(pac99_poly_str)

                # Convert the pac99 polynomial to chemkin polynomial
                chemkin_poly_str = thermo.nasapoly.convert_pac_to_chemkin(pac99_poly_str)
                print('\nCHEMKIN Polynomial:')
                print(chemkin_poly_str)

                bld_afs.build.file.input.write(nasa_inp_str, thy_save_path, bld_locs)
                moldr.util.run_script(NASA_SCRIPT_STR, path)

# 5. process reaction data from the mechanism file
RXN_BLOCK_STR = chemkin_io.reaction_block(MECH_STR)
RXN_STRS = chemkin_io.reaction.data_strings(RXN_BLOCK_STR)
RCT_NAMES_LST = list(
    map(chemkin_io.reaction.DataString.reactant_names, RXN_STRS))
PRD_NAMES_LST = list(
    map(chemkin_io.reaction.DataString.product_names, RXN_STRS))

if RUN_REACTIONS_QCHEM:
    for rct_names, prd_names in zip(RCT_NAMES_LST, PRD_NAMES_LST):
        # print the CHEMKIN reaction name for reference
        rxn_name = '='.join(['+'.join(rct_names), '+'.join(prd_names)])
        print()
        print("Reaction: {}".format(rxn_name))

        # determine inchis, charges, and multiplicities

        rct_smis = list(map(SMI_DCT.__getitem__, rct_names))
        prd_smis = list(map(SMI_DCT.__getitem__, prd_names))
        rct_ichs = list(map(automol.smiles.inchi, rct_smis))
        prd_ichs = list(map(automol.smiles.inchi, prd_smis))
        rct_chgs = list(map(CHG_DCT.__getitem__, rct_names))
        prd_chgs = list(map(CHG_DCT.__getitem__, prd_names))
        rct_muls = list(map(MUL_DCT.__getitem__, rct_names))
        prd_muls = list(map(MUL_DCT.__getitem__, prd_names))

        # determine the transition state multiplicity
        ts_mul = automol.mult.ts.low(rct_muls, prd_muls)

        # theory
        for prog, method, basis in RUN_OPT_LEVELS:
            ts_orb_restr = moldr.util.orbital_restriction(
                ts_mul, RESTRICT_OPEN_SHELL)

            # check direction of reaction
            rxn_ichs = [rct_ichs, prd_ichs]
            rxn_chgs = [rct_chgs, prd_chgs]
            rxn_muls = [rct_muls, prd_muls]
            rxn_exo = moldr.util.reaction_energy(
                SAVE_PREFIX, rxn_ichs, rxn_chgs, rxn_muls, method, basis,
                RESTRICT_OPEN_SHELL)
            print(rxn_exo)
            if rxn_exo > 0:
                rct_ichs, prd_ichs = prd_ichs, rct_ichs
                rct_chgs, prd_chgs = prd_chgs, rct_chgs
                rct_muls, prd_muls = prd_muls, rct_muls
                print('ts search will be performed in reverse direction')

            # obtain geometries from a hierachy of (i) data directory and (ii)
            # previous species calculation
            rct_geos = []
            for ich, chg, mul in zip(rct_ichs, rct_chgs, rct_muls):
                orb_restr = moldr.util.orbital_restriction(
                    mul, RESTRICT_OPEN_SHELL)
                geo = moldr.util.reference_geometry(
                    ich, chg, mul, method, basis, orb_restr, SAVE_PREFIX,
                    GEOM_DCT)
                rct_geos.append(geo)

            prd_geos = []
            for ich, chg, mul in zip(prd_ichs, prd_chgs, prd_muls):
                orb_restr = moldr.util.orbital_restriction(
                    mul, RESTRICT_OPEN_SHELL)
                geo = moldr.util.reference_geometry(
                    ich, chg, mul, method, basis, orb_restr, SAVE_PREFIX,
                    GEOM_DCT)
                prd_geos.append(geo)

            # determine the transition state z-matrix
            # replace this with save values if they are available
            rct_zmas = list(map(automol.geom.zmatrix, rct_geos))
            prd_zmas = list(map(automol.geom.zmatrix, prd_geos))

            typ = None

            # # (migrations are not yet implemented)
            # ret = automol.zmatrix.ts.hydrogen_migration(rct_zmas, prd_zmas)
            # if ret and typ is None:
            #     typ = 'hydrogen migration'

            ret = automol.zmatrix.ts.beta_scission(rct_zmas, prd_zmas)
            if ret and typ is None:
                typ = 'beta scission'
                ts_zma, dist_name, tors_names = ret

            ret = automol.zmatrix.ts.addition(rct_zmas, prd_zmas)
            if ret and typ is None:
                typ = 'addition'
                ts_zma, dist_name, tors_names = ret

            # fix this later
            # ret = automol.zmatrix.ts.hydrogen_abstraction(rct_zmas, prd_zmas,
            #                                               sigma=True)
            ret = automol.zmatrix.ts.hydrogen_abstraction(rct_zmas, prd_zmas,
                                                          sigma=False)
            if ret and typ is None:
                typ = 'hydrogen abstraction'
                ts_zma, dist_name, tors_names = ret

            if typ is None:
                print("Failed to classify reaction.")
            else:
                print("Type: {}".format(typ))

                # determine the grid
                dist_coo, = automol.zmatrix.coordinates(ts_zma)[dist_name]
                syms = automol.zmatrix.symbols(ts_zma)
                bnd_len_key = tuple(sorted(map(syms.__getitem__, dist_coo)))

                bnd_len_dct = {
                    ('C', 'C'): 1.54 * ANG2BOHR,
                    ('C', 'H'): 1.09 * ANG2BOHR,
                    ('H', 'H'): 0.74 * ANG2BOHR,
                    ('N', 'N'): 1.45 * ANG2BOHR,
                    ('O', 'O'): 1.48 * ANG2BOHR,
                    ('C', 'N'): 1.47 * ANG2BOHR,
                    ('C', 'O'): 1.43 * ANG2BOHR,
                    ('H', 'O'): 1.20 * ANG2BOHR,
                    ('H', 'N'): 0.99 * ANG2BOHR,
                }

                if typ in ('beta scission', 'addition'):
                    rmin = 1.4 * ANG2BOHR
                    rmin = 2.8 * ANG2BOHR
                    if bnd_len_key in bnd_len_dct:
                        bnd_len = bnd_len_dct[bnd_len_key]
                        rmin = bnd_len + 0.2 * ANG2BOHR
                        rmax = bnd_len + 1.6 * ANG2BOHR
                elif typ == 'hydrogen abstraction':
                    rmin = 0.7 * ANG2BOHR
                    rmax = 2.2 * ANG2BOHR
                    if bnd_len_key in bnd_len_dct:
                        bnd_len = bnd_len_dct[bnd_len_key]
                        rmin = bnd_len
                        rmax = bnd_len + 1.0 * ANG2BOHR

                npoints = 8
                grid = numpy.linspace(rmin, rmax, npoints)

                # construct the filesystem
                rxn_ichs = [rct_ichs, prd_ichs]
                rxn_chgs = [rct_chgs, prd_chgs]
                rxn_muls = [rct_muls, prd_muls]

                # set up the filesystem
                is_rev = autofile.system.reaction_is_reversed(
                    rxn_ichs, rxn_chgs, rxn_muls)
                rxn_ichs, rxn_chgs, rxn_muls = autofile.system.sort_together(
                    rxn_ichs, rxn_chgs, rxn_muls)
                print(" - The reaction direction is {}"
                      .format('backward' if is_rev else 'forward'))

                rxn_run_fs = autofile.fs.reaction(RUN_PREFIX)
                rxn_run_fs.leaf.create([rxn_ichs, rxn_chgs, rxn_muls, ts_mul])
                rxn_run_path = rxn_run_fs.leaf.path(
                    [rxn_ichs, rxn_chgs, rxn_muls, ts_mul])

                rxn_save_fs = autofile.fs.reaction(RUN_PREFIX)
                rxn_save_fs.leaf.create([rxn_ichs, rxn_chgs, rxn_muls, ts_mul])
                rxn_save_path = rxn_save_fs.leaf.path(
                    [rxn_ichs, rxn_chgs, rxn_muls, ts_mul])

                thy_run_fs = autofile.fs.theory(rxn_run_path)
                thy_run_fs.leaf.create([method, basis, orb_restr])
                thy_run_path = thy_run_fs.leaf.path(
                    [method, basis, orb_restr])

                thy_save_fs = autofile.fs.theory(rxn_save_path)
                thy_save_fs.leaf.create([method, basis, orb_restr])
                thy_save_path = thy_save_fs.leaf.path(
                    [method, basis, orb_restr])

                moldr.driver.run_scan(
                    zma=ts_zma,
                    chg=0,
                    mul=ts_mul,
                    method=method,
                    basis=basis,
                    orb_restr=ts_orb_restr,
                    grid_dct={dist_name: grid},
                    run_prefix=thy_run_path,
                    save_prefix=thy_save_path,
                    script_str=SCRIPT_STR,
                    prog=prog,
                    overwrite=OVERWRITE,
                    update_guess=False,
                    reverse_sweep=False,
                    **OPT_KWARGS
                )

                moldr.driver.save_scan(
                    run_prefix=thy_run_path,
                    save_prefix=thy_save_path,
                    coo_names=[dist_name],
                )

                scn_save_fs = autofile.fs.scan(thy_save_path)
                locs_lst = [
                    locs for locs in scn_save_fs.leaf.existing([[dist_name]])
                    if scn_save_fs.leaf.file.energy.exists(locs)]
                print(locs_lst)
                enes = [scn_save_fs.leaf.file.energy.read(locs)
                        for locs in locs_lst]
                max_locs = locs_lst[enes.index(max(enes))]
                max_ene = max(enes)
                max_zma = scn_save_fs.leaf.file.zmatrix.read(max_locs)

                print('optimizing ts')
                # find saddlepoint from maximum on the grid opt scan
                ts_run_fs = autofile.fs.ts(thy_run_path)
                ts_run_fs.trunk.create()
                ts_run_path = ts_run_fs.trunk.path()

                ts_save_fs = autofile.fs.ts(thy_save_path)
                ts_save_fs.trunk.create()
                ts_save_path = ts_save_fs.trunk.path()

                moldr.driver.run_job(
                    job='optimization',
                    script_str=SCRIPT_STR,
                    prefix=ts_run_path,
                    geom=max_zma,
                    chg=0,
                    mult=ts_mul,
                    method=method,
                    basis=basis,
                    orb_restr=ts_orb_restr,
                    prog=prog,
                    saddle=True,
                    overwrite=OVERWRITE,
                    **OPT_KWARGS,
                )
                opt_ret = moldr.driver.read_job(
                    job='optimization',
                    prefix=ts_run_path,
                )
                if opt_ret is not None:
                    inf_obj, inp_str, out_str = opt_ret
                    prog = inf_obj.prog
                    method = inf_obj.method
                    ene = elstruct.reader.energy(prog, method, out_str)
                    geo = elstruct.reader.opt_geometry(prog, out_str)
                    zma = elstruct.reader.opt_zmatrix(prog, out_str)

                    print(" - Saving...")
                    print(" - Save path: {}".format(ts_save_path))

                    ts_save_fs.leaf.file.geometry_info.write(inf_obj)
                    ts_save_fs.leaf.file.geometry_input.write(inp_str)
                    ts_save_fs.leaf.file.energy.write(ene)
                    ts_save_fs.leaf.file.geometry.write(geo)
                    ts_save_fs.leaf.file.zmatrix.write(zma)

                    moldr.driver.run_job(
                        job='hessian',
                        script_str=SCRIPT_STR,
                        prefix=ts_run_path,
                        geom=geo,
                        chg=0,
                        mult=ts_mul,
                        method=method,
                        basis=basis,
                        orb_restr=orb_restr,
                        prog=prog,
                        overwrite=OVERWRITE,
                        **KWARGS,
                    )
                    hess_ret = moldr.driver.read_job(
                        job='hessian',
                        prefix=ts_run_path,
                    )
                    if hess_ret is not None:
                        inf_obj, inp_str, out_str = hess_ret
                        prog = inf_obj.prog
                        method = inf_obj.method
                        hess = elstruct.reader.hessian(prog, out_str)
                        freqs = elstruct.util.harmonic_frequencies(geo, hess)

                        print(" - Saving hessian...")
                        print(" - Save path: {}".format(ts_save_path))

                        ts_save_fs.leaf.file.hessian_info.write(inf_obj)
                        ts_save_fs.leaf.file.hessian_input.write(inp_str)
                        ts_save_fs.leaf.file.hessian.write(hess)
                        ts_save_fs.leaf.file.harmonic_frequencies.write(freqs)

                if RUN_TS_TAU_SAMP:

                    moldr.driver.save_tau(
                        run_prefix=thy_run_path,
                        save_prefix=thy_save_path,
                    )

                    zma = ts_save_fs.leaf.file.zmatrix.read()
                    tors_ranges = automol.zmatrix.torsional_sampling_ranges(
                        zma, tors_names)
                    tors_range_dct = dict(zip(tors_names, tors_ranges))

                    moldr.driver.run_tau(
                        zma=zma,
                        chg=chg,
                        mult=ts_mul,
                        method=method,
                        basis=basis,
                        orb_restr=orb_restr,
                        nsamp=nsamp,
                        tors_range_dct=tors_range_dct,
                        run_prefix=thy_run_path,
                        save_prefix=thy_save_path,
                        script_str=SCRIPT_STR,
                        prog=prog,
                        saddle=True,
                        overwrite=OVERWRITE,
                        **OPT_KWARGS,
                    )

                    moldr.driver.save_tau(
                        run_prefix=thy_run_path,
                        save_prefix=thy_save_path,
                    )

                if RUN_TS_CONF_SCAN:
                    zma = ts_save_fs.leaf.file.zmatrix.read(thy_save_path)
                    val_dct = automol.zmatrix.values(zma)
                    tors_linspaces = automol.zmatrix.torsional_scan_linspaces(
                        zma, tors_names, SCAN_INCREMENT)
                    tors_grids = [
                        numpy.linspace(*linspace) + val_dct[name]
                        for name, linspace in zip(tors_name, tors_linspaces)]
                    for tors_name, tors_grid in zip(tors_names, tors_grids):
                        moldr.driver.run_scan(
                            zma=zma,
                            chg=chg,
                            mult=ts_mul,
                            method=method,
                            basis=basis,
                            orb_restr=orb_restr,
                            grid_dct={tors_name: tors_grid},
                            run_prefix=ts_run_path,
                            save_prefix=ts_save_path,
                            script_str=SCRIPT_STR,
                            prog=prog,
                            saddle=True,
                            overwrite=OVERWRITE,
                            **OPT_KWARGS,
                        )

                        moldr.driver.save_scan(
                            run_prefix=ts_run_path,
                            save_prefix=ts_save_path,
                            coo_names=[tors_name],
                        )
# AVC: debugged to here
                    hind_rot_dct = {}
                    scn_run_fs = autofile.fs.scan(ts_run_path)
                    scn_save_fs = autofile.fs.scan(ts_save_path)

                    min_ene = ts_save_fs.leaf.file.energy.read()
                    for tors_name in tors_names:
                        enes = [scn_save_fs.leaf.file.energy.read(locs)
                                for locs in scn_save_fs.leaf.existing()]
                        enes = numpy.subtract(enes, min_ene)
                        hind_rot_dct[tors_name] = enes*EH2KCAL

                    print('ts hindered rotor potential')
                    print(hind_rot_dct)

                if RUN_TS_KICKS_QCHEM:
                    ret = moldr.driver.read_job(
                        job=elstruct.Job.HESSIAN, prefix=ts_run_path)
                    if ret:
                        inf_obj, _, out_str = ret
                        prog = inf_obj.prog
                        hess = elstruct.reader.hessian(prog, out_str)
                        freqs = elstruct.util.harmonic_frequencies(geo, hess, project=True)
                        norm_coos = elstruct.util.normal_coordinates(geo, hess, project=True)
                        assert freqs[0] < -100

                        print('Kicking off from saddle in forward direction')
                        im_norm_coo = numpy.array(norm_coos)[:, 0]
                        disp_xyzs = numpy.reshape(im_norm_coo, (-1, 3))
                        dir_afs = autofile.fs.direction()
                        fwd_run_path = dir_afs.direction.dir.path(thy_run_path, [True])
                        dir_afs.direction.dir.create(thy_run_path, [True])
                        fwd_save_path = dir_afs.direction.dir.path(thy_save_path, [True])
                        dir_afs.direction.dir.create(thy_save_path, [True])
                        print(automol.geom.string(geo))
                        print(disp_xyzs)
                        moldr.driver.run_kickoff_saddle(
                            geo, disp_xyzs, chg, mult, method, basis, orb_restr,
                            fwd_run_path, SCRIPT_STR, prog, OVERWRITE,
                            kickoff_size=KICKOFF_SIZE, kickoff_backward=False,
                            opt_cart=True, **OPT_KWARGS)
                        print('Saving product of kick off from saddle in forward direction')
                        ret = moldr.driver.read_job(job=elstruct.Job.OPTIMIZATION, prefix=fwd_run_path)
                        if ret:
                            inf_obj, inp_str, out_str = ret
                            prog = inf_obj.prog
                            method = inf_obj.method
                            ene = elstruct.reader.energy(prog, method, out_str)
                            geo = elstruct.reader.opt_geometry(prog, out_str)
                            fwd_save_path = dir_afs.direction.dir.path(thy_save_path, [True])
                            print('save path', fwd_save_path)
                            dir_afs.direction.file.geometry_info.write(inf_obj, thy_save_path, [True])
                            dir_afs.direction.file.geometry_input.write(inp_str, thy_save_path, [True])
                            dir_afs.direction.file.geometry.write(geo, thy_save_path, [True])
                            dir_afs.direction.file.energy.write(ene, thy_save_path, [True])

                        print('Kicking off from saddle in backward direction')
                        bwd_run_path = dir_afs.direction.dir.path(thy_run_path, [False])
                        dir_afs.direction.dir.create(thy_run_path, [False])
                        bwd_save_path = dir_afs.direction.dir.path(thy_save_path, [False])
                        dir_afs.direction.dir.create(thy_save_path, [False])
                        moldr.driver.run_kickoff_saddle(
                            geo, disp_xyzs, chg, mult, method, basis,
                            orb_restr, bwd_run_path, SCRIPT_STR, prog,
                            OVERWRITE, kickoff_size=KICKOFF_SIZE,
                            kickoff_backward=True, **OPT_KWARGS)
                        print('Saving product of kick off from saddle in backward direction')
                        ret = moldr.driver.read_job(job=elstruct.Job.OPTIMIZATION, prefix=bwd_run_path)
                        if ret:
                            inf_obj, inp_str, out_str = ret
                            prog = inf_obj.prog
                            method = inf_obj.method
                            ene = elstruct.reader.energy(prog, method, out_str)
                            geo = elstruct.reader.opt_geometry(prog, out_str)
                            bwd_save_path = dir_afs.direction.dir.path(thy_save_path, [False])
                            print('save path', bwd_save_path)
                            dir_afs.direction.file.geometry_info.write(inf_obj, thy_save_path, [False])
                            dir_afs.direction.file.geometry_input.write(inp_str, thy_save_path, [False])
                            dir_afs.direction.file.geometry.write(geo, thy_save_path, [False])
                            dir_afs.direction.file.energy.write(ene, thy_save_path, [False])

if RUN_VDW_QCHEM:
    if NSAMP_VDW_EXPR:
        nsamp = min(NSAMP_VDW_A + NSAMP_VDW_B * NSAMP_VDW_C**ntaudof, NSAMP_VDW_D)
    else:
        nsamp = NSAMP_VDW

    VDW_NAMES_LST = []
    for rct_names, prd_names in zip(RCT_NAMES_LST, PRD_NAMES_LST):
        rct_muls = list(map(MUL_DCT.__getitem__, rct_names))
        prd_muls = list(map(MUL_DCT.__getitem__, prd_names))
        ts_mul = automol.mult.ts.low(rct_muls, prd_muls)
        if len(rct_names) == 2:
            if sorted(rct_names) not in VDW_NAMES_LST:
                VDW_NAMES_LST.append([sorted(rct_names), ts_mul])
        if len(prd_names) == 2:
            if sorted(prd_names) not in VDW_NAMES_LST:
                VDW_NAMES_LST.append([sorted(prd_names), ts_mul])

    for names, ts_mul in VDW_NAMES_LST:
        smis = list(map(SMI_DCT.__getitem__, names))
        ichs = list(map(automol.smiles.inchi, smis))
        chgs = list(map(CHG_DCT.__getitem__, names))
        muls = list(map(MUL_DCT.__getitem__, names))

        method = METHOD
        basis = BASIS
        geos = []
        for ich, chg, mult in zip(ichs, chgs, muls):
            orb_restr = moldr.util.orbital_restriction(mult, RESTRICT_OPEN_SHELL)
            geo = moldr.util.reference_geometry(
                ich, chg, mult, method, basis, orb_restr, SAVE_PREFIX, GEOM_DCT)
            geos.append(geo)
           
        geo1, geo2 = geos
        geo1 = automol.geom.mass_centered(geo1)
        geo2 = automol.geom.mass_centered(geo2)
        for idx in range(nsamp):
            print('Optimizing vdw geometry {}/{}'.format(idx+1,nsamp))
            angs1 = numpy.multiply(
                numpy.random.rand(3), [1*numpy.pi, 2*numpy.pi, 2*numpy.pi])
            angs2 = numpy.multiply(
                numpy.random.rand(3), [1*numpy.pi, 2*numpy.pi, 2*numpy.pi])
            angs12 = numpy.multiply(
                numpy.random.rand(2), [1*numpy.pi, 2*numpy.pi])
            geo1 = automol.geom.euler_rotated(geo1, *angs1)
            geo2 = automol.geom.euler_rotated(geo2, *angs2)
            dist_cutoff = 3.*qcc.conversion_factor('angstrom', 'bohr')

            geo = automol.geom.join(geo1, geo2, dist_cutoff, *angs12)
            print("Species: {}".format('+'.join(names)))
            print('vdw starting geometry')
            print(automol.geom.xyz_string(geo))

    # set up the filesystem
            ich = automol.inchi.recalculate(automol.inchi.join(ichs))
            chg = sum(chgs)
            mul = ts_mul
            orb_restr = moldr.util.orbital_restriction(mul, RESTRICT_OPEN_SHELL)
            spc_run_path = moldr.util.species_path(ich, chg, mul, RUN_PREFIX)
            spc_save_path = moldr.util.species_path(ich, chg, mul, SAVE_PREFIX)
            thy_run_path = moldr.util.theory_path(method, basis, orb_restr, spc_run_path)
            thy_save_path = moldr.util.theory_path(method, basis, orb_restr, spc_save_path)

    # generate reference geometry
    # generate the z-matrix and sampling ranges

            moldr.driver.run_job(
                job=elstruct.Job.OPTIMIZATION,
                geom=geo,
                chg=chg,
                mult=mul,
                method=method,
                basis=basis,
                orb_restr=orb_restr,
                prefix=thy_run_path,
                script_str=SCRIPT_STR,
                prog=PROG,
                overwrite=OVERWRITE,
                **OPT_KWARGS,
            )

    # save info for the initial geometry (from inchi or from save directory)
            ret = moldr.driver.read_job(job=elstruct.Job.OPTIMIZATION, prefix=thy_run_path)
            if ret:
                print('Saving reference geometry')
                print(" - Save path: {}".format(thy_save_path))

                inf_obj, inp_str, out_str = ret
                prog = inf_obj.prog
                method = inf_obj.method
                geo = elstruct.reader.opt_geometry(prog, out_str)
                thy_afs = autofile.fs.theory()
                thy_afs.theory.file.geometry.write(geo, spc_save_path, [method, basis, orb_restr])
                ene = elstruct.reader.energy(prog, method, out_str)
                print('ene test in vdw')
                print(ene)
                thy_afs.theory.file.energy.write(ene, spc_save_path, [method, basis, orb_restr])


if RUN_REACTIONS_RATES:
    for rct_names, prd_names in zip(RCT_NAMES_LST, PRD_NAMES_LST):
        # print the CHEMKIN reaction name for reference
        rxn_name = '='.join(['+'.join(rct_names), '+'.join(prd_names)])
        print()
        print("Reaction: {}".format(rxn_name))
        print('Mess Input for')
        print(rxn_name)
        rct_smis = list(map(SMI_DCT.__getitem__, rct_names))
        prd_smis = list(map(SMI_DCT.__getitem__, prd_names))
        rct_ichs = list(map(automol.smiles.inchi, rct_smis))
        prd_ichs = list(map(automol.smiles.inchi, prd_smis))
        rct_chgs = list(map(CHG_DCT.__getitem__, rct_names))
        prd_chgs = list(map(CHG_DCT.__getitem__, prd_names))
        rct_muls = list(map(MUL_DCT.__getitem__, rct_names))
        prd_muls = list(map(MUL_DCT.__getitem__, prd_names))

# header section
        temperatures = TEMPS
        pressures = PRESS
        header_section_str = mess_io.writer.write_global_reaction(temperatures, pressures)
        print(header_section_str)

# energy transfer section
        exp_factor = EXP_FACTOR 
        exp_ppower = EXP_POWER 
        exp_cutoff = EXP_CUTOFF 
        eps1 = EPS1 
        eps2 = EPS2 
        sig1 = SIG1 
        sig2 = SIG2 
        mass1 = MASS1 
        mass2 = MASS2 
        energy_trans_section_str = mess_io.writer.write_energy_transfer(
            exp_factor, exp_power, exp_cutoff, eps1, eps2, sig1, sig2, mass1, mass2)
        print(energy_trans_section_str)
# cycle over reactant and product species
# check if unimolecular or bimolecular species
        print('Unimolecular or bimolecular test')
        indxw = 0
        indxp = 0
        if rct_ichs[1]:
            print(rct_ichs)
            indxw += 1
# write W_indxw 
        else:
            indxp += 1
# write P_indxp
        if prd_ichs[1]:
            print(prd_ichs)
            indxw += 1
# write W_indxw 
        else:
            indxp += 1
# write P_indxp

# cycle over vdw Species
# cycle over transition states
# for now - have only one TS

#        if reaction_typ = addition:
# reactants
#            print(species_str(rct1))
#            print(species_str(rct2))
# well
#            print(species_str(prod1))
#        if reaction_typ = abstraction:
# reactants
#            print(species_str(rct1))
#            print(species_str(rct2))
# vdw
#           for vdw_species in ...
#                   print(species_str(vdwi))
# products
#            print(species_str(prod1))
#            print(species_str(prod2))

#        if reaction_typ = abstraction
# reactants
#            print(species_str(rct1))
#            print(species_str(rct2))
# vdw
#            print(species_str(vdw1))
#            print(species_str(vdw2))
# products
#            print(species_str(prod1))
#            print(species_str(prod2))
# ts
        core = mess_io.writer.write_core_rigidrotor(geo, symfactor)
        molecule_section_str1 = mess_io.writer.write_molecule(
            core, freqs, zpe, elec_levels, hind_rot='',
        )
        print(molecule_section_str1)

sys.exit()
