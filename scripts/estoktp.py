""" reaction list test
"""
import os
import sys
import numpy
import pandas
from qcelemental import constants as qcc
import chemkin_io
import automol
import elstruct
import autofile
import moldr
import mess_io.writer

ANG2BOHR = qcc.conversion_factor('angstrom', 'bohr')
WAVEN2KCAL = qcc.conversion_factor('wavenumber', 'kcal/mol')
EH2KCAL = qcc.conversion_factor('hartree', 'kcal/mol')

# 0. choose which mechanism to run
MECHANISM_NAME = 'ch4+nh2'  # options: syngas, natgas, heptane, test, estoktp, ...
#MECHANISM_NAME = 'onereac'  # options: syngas, natgas, heptane, test, estoktp, ...
#MECHANISM_NAME = 'estoktp/add30'  # options: syngas, natgas, heptane, test, estoktp, ...
#MECHANISM_NAME = 'estoktp/habs65'  # options: syngas, natgas, heptane, test, estoktp, ...

# 1. script control parameters
# a. Electronic structure parameters; code, Program, method, basis, convergence control

#PROG = 'psi4'
#SCRIPT_STR = ("#!/usr/bin/env bash\n"
#              "psi4 -i run.inp -o run.out >> stdout.log &> stderr.log")
#KWARGS = {}

PROG = 'g09'
SCRIPT_STR = ("#!/usr/bin/env bash\n"
              "g09 run.inp run.out >> stdout.log &> stderr.log")
PF_SCRIPT_STR = ("#!/usr/bin/env bash\n"
              "messpf build.inp build.out >> stdout.log &> stderr.log")
NASA_SCRIPT_STR = ("#!/usr/bin/env bash\n"
              "cp ../PF/build.out pf.dat"
              "cp /tcghome/sjklipp/PACC/nasa/new.groups ."
              "python /tcghome/sjklipp/PACC/nasa/makepoly.py >> stdout.log &> stderr.log")
METHOD = 'wb97xd'
BASIS = '6-31g*'
OPT_KWARGS = {
    'memory': 10,
    'machine_options': ['%NProcShared=10'],
    'gen_lines': ['# int=ultrafine'],
    'feedback': True,
    'errors': [
        elstruct.Error.OPT_NOCONV
    ],
    'options_mat': [
        [{},
         {},
         {},
         {'job_options': ['calcfc']},
         {'job_options': ['calcfc']},
         {'job_options': ['calcall']}]
    ],
}
KWARGS = {
    'memory': 10,
    'machine_options': ['%NProcShared=10'],
    'gen_lines': ['# int=ultrafine'],
}

#    'job_options': ['verytight'],

# b. What to run for electronic structure calculations
RUN_SPECIES_QCHEM = True
RUN_REACTIONS_QCHEM = True
RUN_TS_KICKS_QCHEM = True
RUN_VDW_QCHEM = True

KICKOFF_SADDLE = True
KICKOFF_BACKWARD = False
KICKOFF_SIZE = 0.1

RUN_CONF_OPT = True
RUN_MIN_GRAD = True
RUN_MIN_HESS = True
RUN_CONF_GRAD = False
RUN_CONF_HESS = False

RUN_CONF_SCAN = True
RUN_CONF_SCAN_GRAD = False
RUN_CONF_SCAN_HESS = False

RUN_TAU_SAMP = False
RUN_TAU_GRAD = False
RUN_TAU_HESS = False

RUN_TS_CONF_OPT = False
RUN_TS_CONF_SCAN = False
RUN_TS_TAU_SAMP = False

RUN_GRAD = False
if RUN_GRAD:
    RUN_MIN_GRAD = True
    RUN_CONF_GRAD = True
    RUN_CONF_SCAN_GRAD = True
    RUN_TAU_GRAD = True

RUN_HESS = False
if RUN_HESS:
    RUN_MIN_HESS = True
    RUN_CONF_HESS = True
    RUN_CONF_SCAN_HESS = True
    RUN_TAU_HESS = True

# c. What to run for thermochemical kinetics
RUN_SPECIES_PF = False
RUN_SPECIES_THERMO = False
RUN_REACTIONS_RATES = False
RUN_VDW_RCT_RATES = False
RUN_VDW_PRD_RATES = False

# d. Parameters for number of torsional samplings
NSAMP_CONF = 5
NSAMP_CONF_EXPR = False
NSAMP_CONF_A = 3
NSAMP_CONF_B = 1
NSAMP_CONF_C = 3
NSAMP_CONF_D = 15

NSAMP_TAU = 10
NSAMP_TAU_EXPR = False
NSAMP_TAU_A = 3
NSAMP_TAU_B = 1
NSAMP_TAU_C = 3
NSAMP_TAU_D = 15

NSAMP_VDW = 10
NSAMP_VDW_EXPR = False
NSAMP_VDW_A = 3
NSAMP_VDW_B = 1
NSAMP_VDW_C = 3
NSAMP_VDW_D = 15

# Parameters for partition function and kinetics
SPECIES_STR={}
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

# Defaults
SCAN_INCREMENT = 30. * qcc.conversion_factor('degree', 'radian')
RESTRICT_OPEN_SHELL = False
OVERWRITE = False
#SMILES_LST = ['[O]', '[OH]', '[N]=O', '[CH2]', '[C]', '[B]', '[N]', '[F]', '[Cl]', '[Br]',
#              '[BH]', '[BH2]', 'C[C]', '[O][O]']
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
#SMILES_LST = ['[C]#N']
#for smiles in SMILES_LST:
#    ich = automol.convert.smiles.inchi(smiles)
#    print(ich)

ELC_SIG_LST = {'InChI=1S/CN/c1-2'}

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

#SMILES_LST = ['[H]', '[OH]', 'O[O]', '[CH3]', '[O]', 'C', 'CC', 'C[CH2]', 'C=C', 'C=[CH]',
#              'C#C', 'C#[C]', 'CO', '[CH2]=O', 'C[O]', 'OC=O', 'OC=[O]', 'O[C]O', 'COC',
#              'CO[CH2]', 'C=O', 'O=[CH]', 'CCl', '[CH2]Cl', 'S', '[SH]', 'N', '[NH2]']

SPC_TAB['charge'] = 0
SMI_DCT = dict(zip(SPC_TAB['name'], SPC_TAB['smiles']))
CHG_DCT = dict(zip(SPC_TAB['name'], SPC_TAB['charge']))
MUL_DCT = dict(zip(SPC_TAB['name'], SPC_TAB['mult']))
SPC_BLK_STR = chemkin_io.species_block(MECH_STR)
SPC_NAMES = chemkin_io.species.names(SPC_BLK_STR)

GEOM_PATH = os.path.join(DATA_PATH, 'data', 'geoms')
print(GEOM_PATH)
GEOM_DCT = moldr.util.geometry_dictionary(GEOM_PATH)

# take starting geometry from saved directory if possible, otherwise get it from inchi via rdkit

if RUN_SPECIES_QCHEM:
    for name in SPC_NAMES:
# species
        print("Species: {}".format(name))
        smi = SMI_DCT[name]
        ich = automol.smiles.inchi(smi)
        print("smiles: {}".format(smi), "inchi: {}".format(ich))
#        ich = ICH_DCT[name]
        charge = CHG_DCT[name]
        mult = MUL_DCT[name]

# theory
        method = METHOD
        basis = BASIS
        orb_restr = moldr.util.orbital_restriction(mult, RESTRICT_OPEN_SHELL)

# set up the filesystem
        spc_run_path = moldr.util.species_path(ich, charge, mult, RUN_PREFIX)
        spc_save_path = moldr.util.species_path(ich, charge, mult, SAVE_PREFIX)
        thy_run_path = moldr.util.theory_path(method, basis, orb_restr, spc_run_path)
        thy_save_path = moldr.util.theory_path(method, basis, orb_restr, spc_save_path)

# generate reference geometry
# generate the z-matrix and sampling ranges

        geo = moldr.util.reference_geometry(
            ich, charge, mult, method, basis, orb_restr, SAVE_PREFIX, GEOM_DCT)
        zma = automol.geom.zmatrix(geo)
        moldr.driver.run_job(
            job=elstruct.Job.OPTIMIZATION,
            geom=zma,
            charge=charge,
            mult=mult,
            method=method,
            basis=basis,
            orb_restr=orb_restr,
            prefix=thy_run_path,
            script_str=SCRIPT_STR,
            prog=PROG,
            overwrite=OVERWRITE,
            **OPT_KWARGS,
        )
        if KICKOFF_SADDLE:
# check if optimized geometry has negative frequencies
# if it does then kick in direction of imaginary mode
            ret = moldr.driver.read_job(job=elstruct.Job.OPTIMIZATION, prefix=thy_run_path)
            if ret:
                print('Checking for saddle')
                inf_obj, _, out_str = ret
                prog = inf_obj.prog
                geo = elstruct.reader.opt_geometry(prog, out_str)
                moldr.driver.run_job(
                    job=elstruct.Job.HESSIAN,
                    geom=geo,
                    charge=charge,
                    mult=mult,
                    method=method,
                    basis=basis,
                    orb_restr=orb_restr,
                    prefix=thy_run_path,
                    script_str=SCRIPT_STR,
                    prog=PROG,
                    overwrite=OVERWRITE,
                    **KWARGS,
                )
                ret = moldr.driver.read_job(job=elstruct.Job.HESSIAN, prefix=thy_run_path)
                if ret:
                    inf_obj, _, out_str = ret
                    prog = inf_obj.prog
                    hess = elstruct.reader.hessian(prog, out_str)
                    print('hess test')
                    print(automol.geom.string(geo))
                    print(numpy.array(hess))
                    freqs = elstruct.util.harmonic_frequencies(geo, hess, project=True)
                    print('projected freqs')
                    print(freqs)
                    norm_coos = elstruct.util.normal_coordinates(geo, hess, project=True)

# if there's an imaginary frequency, optimize again after displacing along the mode
                    if freqs[0] < -100:
                        print('Imaginary mode found: Attempting to kickoff from saddle')
                        im_norm_coo = numpy.array(norm_coos)[:, 0]
                        disp_xyzs = numpy.reshape(im_norm_coo, (-1, 3))
                        moldr.driver.run_kickoff_saddle(
                            geo, disp_xyzs, charge, mult, method, basis,
                            orb_restr, thy_run_path, SCRIPT_STR, prog,
                            OVERWRITE, kickoff_size=KICKOFF_SIZE,
                            kickoff_backward=KICKOFF_BACKWARD, opt_cart=False,
                            **OPT_KWARGS)
                        print('removing saddlepoint hessian')
                        run_afs = autofile.fs.run()
                        run_afs.run.dir.remove(thy_run_path, ['hessian'])
#                        print('writing corrected geometry to data directory')
#                        geom_dfile = autofile.system.file_.geometry(name)
#                        geom_dfile.write(geo, GEOM_PATH)

# save info for the initial geometry (from inchi or from save directory)
        ret = moldr.driver.read_job(job=elstruct.Job.OPTIMIZATION, prefix=thy_run_path)
        if ret:
            print('Saving reference geometry')
            print(" - Save path: {}".format(thy_save_path))

            inf_obj, inp_str, out_str = ret
            prog = inf_obj.prog
            method = inf_obj.method
            ene = elstruct.reader.energy(prog, method, out_str)
            geo = elstruct.reader.opt_geometry(prog, out_str)
            zma = automol.geom.zmatrix(geo)
            thy_afs = autofile.fs.theory()
            thy_afs.theory.file.geometry.write(geo, spc_save_path, [method, basis, orb_restr])
            thy_afs.theory.file.zmatrix.write(zma, spc_save_path, [method, basis, orb_restr])

# a. conformer sampling
        tors_names = automol.geom.zmatrix_torsion_coordinate_names(geo)
        tors_ranges = automol.zmatrix.torsional_sampling_ranges(
            zma, tors_names)
        tors_range_dct = dict(zip(tors_names, tors_ranges))

        gra = automol.inchi.graph(ich)
        ntaudof = len(automol.graph.rotational_bond_keys(gra, with_h_rotors=False))
        if NSAMP_CONF_EXPR:
            nsamp = min(NSAMP_CONF_A + NSAMP_CONF_B * NSAMP_CONF_C**ntaudof, NSAMP_CONF_D)
        else:
            nsamp = NSAMP_CONF

        moldr.driver.save_conformers(
            run_prefix=thy_run_path,
            save_prefix=thy_save_path,
        )
        print('zma test in conformer sampling')
        print(zma)

        moldr.driver.run_conformers(
            zma=zma,
            charge=charge,
            mult=mult,
            method=method,
            basis=basis,
            orb_restr=orb_restr,
            nsamp=nsamp,
            tors_range_dct=tors_range_dct,
            run_prefix=thy_run_path,
            save_prefix=thy_save_path,
            script_str=SCRIPT_STR,
            prog=PROG,
            overwrite=OVERWRITE,
            **OPT_KWARGS,
        )

        moldr.driver.save_conformers(
            run_prefix=thy_run_path,
            save_prefix=thy_save_path,
        )

# b. save information about the minimum energy conformer in top directory
        thy_afs = autofile.fs.theory()
        cnf_afs = autofile.fs.conformer()
        cnf_alocs_lst = cnf_afs.conf.dir.existing(thy_save_path)
        min_cnf_alocs = moldr.util.min_energy_conformer_locators(thy_save_path)
        if min_cnf_alocs:
            min_cnf_run_path = cnf_afs.conf.dir.path(thy_run_path, min_cnf_alocs)
            min_cnf_save_path = cnf_afs.conf.dir.path(thy_save_path, min_cnf_alocs)
            geo = cnf_afs.conf.file.geometry.read(thy_save_path, min_cnf_alocs)
            zma = cnf_afs.conf.file.zmatrix.read(thy_save_path, min_cnf_alocs)
            assert automol.zmatrix.almost_equal(zma, automol.geom.zmatrix(geo))
            thy_afs.theory.file.geometry.write(geo, spc_save_path, [method, basis, orb_restr])
            thy_afs.theory.file.zmatrix.write(zma, spc_save_path, [method, basis, orb_restr])


# c. gradients and hessians for minimum energy conformer
            if RUN_MIN_GRAD:
                print('Minimum energy conformer gradient')
                moldr.driver.run_job(
                    job='gradient',
                    script_str=SCRIPT_STR,
                    prefix=min_cnf_run_path,
                    geom=geo,
                    charge=charge,
                    mult=mult,
                    method=method,
                    basis=basis,
                    orb_restr=orb_restr,
                    prog=PROG,
                    overwrite=OVERWRITE,
                    **KWARGS,
                )

                ret = moldr.driver.read_job(
                    job='gradient',
                    prefix=min_cnf_run_path,
                )

                if ret is not None:
                    inf_obj, inp_str, out_str = ret

                    print(" - Reading gradient from output...")
                    grad = elstruct.reader.gradient(inf_obj.prog, out_str)

                    print(" - Saving gradient...")
                    print(" - Save path: {}".format(min_cnf_save_path))
                    cnf_afs.conf.file.gradient_info.write(inf_obj, thy_save_path, min_cnf_alocs)
                    cnf_afs.conf.file.gradient_input.write(inp_str, thy_save_path, min_cnf_alocs)
                    cnf_afs.conf.file.gradient.write(grad, thy_save_path, min_cnf_alocs)

            if RUN_MIN_HESS:
                print('Minimum energy conformer hessian')
                moldr.driver.run_job(
                    job='hessian',
                    script_str=SCRIPT_STR,
                    prefix=min_cnf_run_path,
                    geom=geo,
                    charge=charge,
                    mult=mult,
                    method=method,
                    basis=basis,
                    orb_restr=orb_restr,
                    prog=PROG,
                    overwrite=OVERWRITE,
                )

                ret = moldr.driver.read_job(
                    job='hessian',
                    prefix=min_cnf_run_path,
                )

                if ret is not None:
                    inf_obj, inp_str, out_str = ret

                    print(" - Reading hessian from output...")
                    hess = elstruct.reader.hessian(inf_obj.prog, out_str)
                    freqs = elstruct.util.harmonic_frequencies(geo, hess, project=False)
                    print('Freqs test')
                    print(freqs)
                    freqs = elstruct.util.harmonic_frequencies(geo, hess, project=True)
                    print('Projected freqs test')
                    print(freqs)

                    print(" - Saving hessian...")
                    print(" - Save path: {}".format(min_cnf_save_path))
                    cnf_afs.conf.file.hessian_info.write(inf_obj, thy_save_path, min_cnf_alocs)
                    cnf_afs.conf.file.hessian_input.write(inp_str, thy_save_path, min_cnf_alocs)
                    cnf_afs.conf.file.hessian.write(hess, thy_save_path, min_cnf_alocs)
                    cnf_afs.conf.file.harmonic_frequencies.write(freqs, thy_save_path, min_cnf_alocs)

# d. hindered rotor scans

            val_dct = automol.zmatrix.values(zma)
            tors_names = automol.geom.zmatrix_torsion_coordinate_names(geo)
            tors_linspaces = automol.zmatrix.torsional_scan_linspaces(
                zma, tors_names, SCAN_INCREMENT)
            tors_grids = [
                numpy.linspace(*linspace) +val_dct[name]
                for name, linspace in zip(tors_names, tors_linspaces)]

            if RUN_CONF_SCAN:
                for tors_name, tors_grid in zip(tors_names, tors_grids):
                    moldr.driver.run_scan(
                        zma=zma,
                        charge=charge,
                        mult=mult,
                        method=method,
                        basis=basis,
                        orb_restr=orb_restr,
                        grid_dct={tors_name: tors_grid},
                        run_prefix=min_cnf_run_path,
                        save_prefix=min_cnf_save_path,
                        script_str=SCRIPT_STR,
                        prog=PROG,
                        overwrite=OVERWRITE,
                        **OPT_KWARGS,
                    )

                    moldr.driver.save_scan(
                        run_prefix=min_cnf_run_path,
                        save_prefix=min_cnf_save_path,
                        coo_names=[tors_name],
                    )
# need to add grad and hessian for each point on the conf_scan
#                    if RUN_CONF_SCAN_GRAD:

#                    if RUN_CONF_SCAN_HESS:

# determine gradient and hessian for each of the conformers
            if RUN_CONF_GRAD:
                for alocs in cnf_locs_lst:
                    cnf_run_path = cnf_afs.conf.dir.path(thy_run_path, alocs)
                    cnf_save_path = cnf_afs.conf.dir.path(thy_save_path, alocs)
                    geo = cnf_afs.conf.file.geometry.read(thy_save_path, alocs)

                    print('Running conformer gradient')
                    moldr.driver.run_job(
                        job='gradient',
                        script_str=SCRIPT_STR,
                        prefix=cnf_run_path,
                        geom=geo,
                        charge=charge,
                        mult=mult,
                        method=method,
                        basis=basis,
                        orb_restr=orb_restr,
                        prog=PROG,
                        overwrite=OVERWRITE,
                        **KWARGS,
                    )

                    ret = moldr.driver.read_job(
                        job='gradient',
                        prefix=cnf_run_path,
                    )

                    if ret is not None:
                        inf_obj, inp_str, out_str = ret

                        print(" - Reading gradient from output...")
                        grad = elstruct.reader.gradient(inf_obj.prog, out_str)

                        print(" - Saving gradient...")
                        print(" - Save path: {}".format(cnf_save_path))
                        cnf_afs.conf.file.gradient_info.write(inf_obj, thy_save_path, alocs)
                        cnf_afs.conf.file.gradient_input.write(inp_str, thy_save_path, alocs)
                        cnf_afs.conf.file.gradient.write(grad, thy_save_path, alocs)

            if RUN_CONF_HESS:
                for alocs in cnf_alocs_lst:
                    cnf_run_path = cnf_afs.conf.dir.path(thy_run_path, alocs)
                    cnf_save_path = cnf_afs.conf.dir.path(thy_save_path, alocs)
                    geo = cnf_afs.conf.file.geometry.read(thy_save_path, alocs)

                    print('Running conformer hessian')
                    moldr.driver.run_job(
                        job='hessian',
                        script_str=SCRIPT_STR,
                        prefix=cnf_run_path,
                        geom=geo,
                        charge=charge,
                        mult=mult,
                        method=method,
                        basis=basis,
                        orb_restr=orb_restr,
                        prog=PROG,
                        overwrite=OVERWRITE,
                    )

                    ret = moldr.driver.read_job(
                        job='hessian',
                        prefix=cnf_run_path,
                    )

                    if ret is not None:
                        inf_obj, inp_str, out_str = ret

                        print(" - Reading hessian from output...")
                        hess = elstruct.reader.hessian(inf_obj.prog, out_str)
                        freqs = elstruct.util.harmonic_frequencies(geo, hess, project=False)
                        print('Conformer Freqs test')
                        print(freqs)

                        print(" - Saving hessian...")
                        print(" - Save path: {}".format(cnf_save_path))
                        cnf_afs.conf.file.hessian_info.write(inf_obj, thy_save_path, alocs)
                        cnf_afs.conf.file.hessian_input.write(inp_str, thy_save_path, alocs)
                        cnf_afs.conf.file.hessian.write(hess, thy_save_path, alocs)
                        cnf_afs.conf.file.harmonic_frequencies.write(freqs, thy_save_path, alocs)

        if RUN_TAU_SAMP:
            if NSAMP_TAU_EXPR:
                nsamp = min(NSAMP_TAU_A + NSAMP_TAU_B * NSAMP_TAU_C**ntaudof, NSAMP_TAU_D)
            else:
                nsamp = NSAMP_TAU
            moldr.driver.save_tau(
                run_prefix=thy_run_path,
                save_prefix=thy_save_path,
            )

            tors_names = automol.geom.zmatrix_torsion_coordinate_names(geo)
            tors_ranges = automol.zmatrix.torsional_sampling_ranges(
                zma, tors_names)
            tors_range_dct = dict(zip(tors_names, tors_ranges))
            moldr.driver.run_tau(
                zma=zma,
                charge=charge,
                mult=mult,
                method=method,
                basis=basis,
                orb_restr=orb_restr,
                nsamp=nsamp,
                tors_range_dct=tors_range_dct,
                run_prefix=thy_run_path,
                save_prefix=thy_save_path,
                script_str=SCRIPT_STR,
                prog=PROG,
                overwrite=OVERWRITE,
                **OPT_KWARGS,
            )

            moldr.driver.save_tau(
                run_prefix=thy_run_path,
                save_prefix=thy_save_path,
            )

            if RUN_TAU_GRAD:
# cycle through saved tau geometries
                afs = autofile.fs.tau()
                for alocs in afs.tau.dir.existing(thy_save_path):
#                    tau_run_path = afs.tau.dir.path(thy_run_path, alocs)
#                    print('path test')
#                    print(thy_save_path)
#                    print(alocs)
#                    print(tau_run_path)
#                    print(tau_save_path)
#                    tau_save_path = afs.tau.dir.path(thy_save_path, alocs)
#                    geo = afs.tau.file.geometry.read(tau_save_path)
                    geo = afs.tau.file.geometry.read(thy_save_path, alocs)
                    print('Running tau gradient')
                    moldr.driver.run_job(
                        job='gradient',
                        script_str=SCRIPT_STR,
                        prefix=thy_run_path,
                        geom=geo,
                        charge=charge,
                        mult=mult,
                        method=method,
                        basis=basis,
                        orb_restr=orb_restr,
                        prog=PROG,
                        overwrite=OVERWRITE,
                        **KWARGS,
                    )

                    ret = moldr.driver.read_job(
                        job='gradient',
                        prefix=thy_run_path,
                    )

                    if ret is not None:
                        inf_obj, inp_str, out_str = ret

                        print(" - Reading gradient from output...")
                        grad = elstruct.reader.gradient(inf_obj.prog, out_str)

                        print(" - Saving gradient...")
                        print(" - Save path: {}".format(thy_save_path))
                        afs.tau.file.gradient_info.write(inf_obj, thy_save_path, alocs)
                        afs.tau.file.gradient_input.write(inp_str, thy_save_path, alocs)
                        afs.tau.file.gradient.write(grad, thy_save_path, alocs)

# cycle through saved tau geometries
            if RUN_TAU_HESS:
# cycle through saved tau geometries
                afs = autofile.fs.tau()
                for alocs in afs.tau.dir.existing(thy_save_path):
#                    tau_run_path = afs.tau.dir.path(thy_run_path, alocs)
#                    tau_save_path = afs.tau.dir.path(thy_save_path, alocs)
                    geo = afs.tau.file.geometry.read(thy_save_path, alocs)
                    print('Running tau Hessian')
                    moldr.driver.run_job(
                        job='hessian',
                        script_str=SCRIPT_STR,
                        prefix=thy_run_path,
                        geom=geo,
                        charge=charge,
                        mult=mult,
                        method=method,
                        basis=basis,
                        orb_restr=orb_restr,
                        prog=PROG,
                        overwrite=OVERWRITE,
                        **KWARGS,
                    )

                    ret = moldr.driver.read_job(
                        job='hessian',
                        prefix=thy_run_path,
                    )

                    if ret is not None:
                        inf_obj, inp_str, out_str = ret

                        print(" - Reading hessian from output...")
                        hess = elstruct.reader.hessian(inf_obj.prog, out_str)

                        print(" - Saving hessian...")
                        print(" - Save path: {}".format(thy_save_path))
                        afs.tau.file.hessian_info.write(inf_obj, thy_save_path, alocs)
                        afs.tau.file.hessian_input.write(inp_str, thy_save_path, alocs)
                        afs.tau.file.hessian.write(hess, thy_save_path, alocs)

if RUN_SPECIES_PF:
    for name in SPC_NAMES:
# set up species information
        smi = SMI_DCT[name]
        ich = automol.smiles.inchi(smi)
        print("smiles: {}".format(smi), "inchi: {}".format(ich))
        charge = CHG_DCT[name]
        mult = MUL_DCT[name]
# specify electronic structure method used
        method = METHOD
        basis = BASIS
        orb_restr = moldr.util.orbital_restriction(mult, RESTRICT_OPEN_SHELL)
# read in geometry, hessian and hindered rotor potentials for minimum energy conformer
        spc_save_path = moldr.util.species_path(ich, charge, mult, SAVE_PREFIX)
        thy_save_path = moldr.util.theory_path(method, basis, orb_restr, spc_save_path)
        cnf_afs = autofile.fs.conformer()
        min_cnf_alocs = moldr.util.min_energy_conformer_locators(thy_save_path)
# I think we need something for if it is none
        if min_cnf_alocs is not None:
            geo = cnf_afs.conf.file.geometry.read(thy_save_path, min_cnf_alocs)
            hess = cnf_afs.conf.file.hessian.read(thy_save_path, min_cnf_alocs)
            freqs = elstruct.util.harmonic_frequencies(geo, hess)
            zpe = sum(freqs)*WAVEN2KCAL/2.
            zma = automol.geom.zmatrix(geo)
            gra = automol.zmatrix.graph(zma, remove_stereo=True)
            scan_afs = autofile.fs.scan()
            min_ene = cnf_afs.conf.file.energy.read(thy_save_path, min_cnf_alocs)
            cnf_save_path = cnf_afs.conf.dir.path(thy_save_path, min_cnf_alocs)
            tors_names = automol.geom.zmatrix_torsion_coordinate_names(geo)
            coo_dct = automol.zmatrix.coordinates(zma, multi=False)
            hind_rot_str = ""

            for tors_name in tors_names:
                enes = [scan_afs.scan.file.energy.read(cnf_save_path, [[tors_name]] + rlocs)
                        for rlocs in scan_afs.scan.dir.existing(cnf_save_path, [[tors_name]])]
                enes = numpy.subtract(enes, min_ene)
                pot = list(enes*EH2KCAL)
                axis = coo_dct[tors_name][1:3]
                group = list(automol.graph.branch_atom_keys(gra, axis[0], axis) - set(axis))
                group = list(numpy.add(group, 1))
                axis = list(numpy.add(axis, 1))
                sym = 1
                hind_rot_str += mess_io.writer.write_rotor_hindered(
                    group, axis, sym, pot)

            print('hina_rot_str test')
            print(hind_rot_str)

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
            [], temp_step, ntemps, rel_temp_inc=0.001, atom_dist_min=0.6)
        print(global_pf_str)
        species_head_str = 'Species ' + name
        print(species_head_str)
        if automol.geom.is_atom(geo):
            print('This is an atom')
            SPECIES_STR[name] = mess_io.writer.write_atom(name, elec_levels)
        else:
            if len(automol.geom.symbols(geo)) == 2:
                freq_offset = 5
            else:
                freq_offset = 6
            core = mess_io.writer.write_core_rigidrotor(geo, symfactor)
            if pot is not None:
                SPECIES_STR[name] = mess_io.writer.write_molecule(
                    core, freqs[freq_offset:], zpe, elec_levels,
                    hind_rot=hind_rot_str,
                )
                freq_offset += 1
        print(SPECIES_STR[name])
        pf_inp_str = '\n'.join(
            [global_pf_str, species_head_str, SPECIES_STR[name]])
        bld_afs = autofile.fs.build()
        bld_alocs = ['PF', 0]
        bld_afs.build.dir.create(thy_save_path, bld_alocs)
        path = bld_afs.build.dir.path(thy_save_path, bld_alocs)
        print('Build Path for Partition Functions')
        print(path)
        bld_afs.build.file.input.write(pf_inp_str, thy_save_path, bld_alocs)
        moldr.util.run_script(PF_SCRIPT_STR, path)

if RUN_SPECIES_THERMO:
    for name in SPC_NAMES:
# set up species information
        smi = SMI_DCT[name]
        ich = automol.smiles.inchi(smi)
        print("smiles: {}".format(smi), "inchi: {}".format(ich))
        charge = CHG_DCT[name]
        mult = MUL_DCT[name]
# specify electronic structure method used
        method = METHOD
        basis = BASIS
        orb_restr = moldr.util.orbital_restriction(mult, RESTRICT_OPEN_SHELL)
# read in geometry, hess and hindered rotor potentials for minimum energy conformer
        spc_save_path = moldr.util.species_path(ich, charge, mult, SAVE_PREFIX)
        thy_save_path = moldr.util.theory_path(method, basis, orb_restr, spc_save_path)
        nasa_inp_str=('nasa')
        # above needs to fixed
        bld_afs = autofile.fs.build()
        bld_alocs = ['NASA_POLY', 0]
        bld_afs.build.dir.create(thy_save_path, bld_alocs)
        path = bld_afs.build.dir.path(thy_save_path, bld_alocs)
        print('Build Path for NASA Polynomials')
        print(path)
        bld_afs.build.file.input.write(nasa_inp_str, thy_save_path, bld_alocs)
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
        method = METHOD
        basis = BASIS
        ts_orb_restr = moldr.util.orbital_restriction(ts_mul, RESTRICT_OPEN_SHELL)

        # check direction of reaction
        rxn_ichs = [rct_ichs, prd_ichs]
        rxn_chgs = [rct_chgs, prd_chgs]
        rxn_muls = [rct_muls, prd_muls]
        rxn_exo = moldr.util.reaction_energy(
            SAVE_PREFIX, rxn_ichs, rxn_chgs, rxn_muls, method, basis, RESTRICT_OPEN_SHELL
            )
        print(rxn_exo)
        if rxn_exo > 0:
            rct_ichs, prd_ichs = prd_ichs, rct_ichs
            rct_chgs, prd_chgs = prd_chgs, rct_chgs
            rct_muls, prd_muls = prd_muls, rct_muls
            print('ts search will be performed in reverse direction')

# obtain geometries from a hierachy of (i) data directory and (ii) previous species calculation
        rct_geos = []
        for ich, chg, mult in zip(rct_ichs, rct_chgs, rct_muls):
            orb_restr = moldr.util.orbital_restriction(mult, RESTRICT_OPEN_SHELL)
            geo = moldr.util.reference_geometry(
                ich, chg, mult, method, basis, orb_restr, SAVE_PREFIX, GEOM_DCT)
            rct_geos.append(geo)

        prd_geos = []
        for ich, chg, mult in zip(prd_ichs, prd_chgs, prd_muls):
            orb_restr = moldr.util.orbital_restriction(mult, RESTRICT_OPEN_SHELL)
            geo = moldr.util.reference_geometry(
                ich, chg, mult, method, basis, orb_restr, SAVE_PREFIX, GEOM_DCT)
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

#        ret = automol.zmatrix.ts.hydrogen_abstraction(rct_zmas, prd_zmas, sigma=True)
# fix this later
        ret = automol.zmatrix.ts.hydrogen_abstraction(rct_zmas, prd_zmas, sigma=False)
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

            rxn_run_path = moldr.util.reaction_path(
                rxn_ichs, rxn_chgs, rxn_muls, ts_mul, RUN_PREFIX)
            rxn_save_path = moldr.util.reaction_path(
                rxn_ichs, rxn_chgs, rxn_muls, ts_mul, SAVE_PREFIX)
            thy_run_path = moldr.util.theory_path(method, basis, orb_restr, rxn_run_path)
            thy_save_path = moldr.util.theory_path(method, basis, orb_restr, rxn_save_path)
            moldr.driver.run_scan(
                zma=ts_zma,
                charge=0,
                mult=ts_mul,
                method=method,
                basis=basis,
                orb_restr=ts_orb_restr,
                grid_dct={dist_name: grid},
                run_prefix=thy_run_path,
                save_prefix=thy_save_path,
                script_str=SCRIPT_STR,
                prog=PROG,
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

            scan_afs = autofile.fs.scan()
            rlocs_lst = scan_afs.scan.dir.existing(thy_save_path, [[dist_name]])
            enes = [scan_afs.scan.file.energy.read(thy_save_path, [[dist_name]] + rlocs)
                    for rlocs in rlocs_lst]
            max_ene = max(enes)
            max_rlocs = rlocs_lst[enes.index(max(enes))]
            geos = [scan_afs.scan.file.energy.read(thy_save_path, [[dist_name]] + rlocs)
                    for rlocs in rlocs_lst]
            max_ene = max(enes)
            max_zma = scan_afs.scan.file.zmatrix.read(thy_save_path, [[dist_name]] + max_rlocs)
            print(max_ene)
            print(max_rlocs)
            print('optimizing ts')
# find saddlepoint from maximum on the grid opt scan
            ts_afs = autofile.fs.ts()
            ts_afs.ts.dir.create(thy_run_path)
            ts_afs.ts.dir.create(thy_save_path)
            ts_run_path = ts_afs.ts.dir.path(thy_run_path)
            ts_save_path = ts_afs.ts.dir.path(thy_save_path)
            moldr.driver.run_job(
                job='optimization',
                script_str=SCRIPT_STR,
                prefix=ts_run_path,
                geom=max_zma,
                charge=0,
                mult=ts_mul,
                method=method,
                basis=basis,
                orb_restr=ts_orb_restr,
                prog=PROG,
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

                ts_afs.ts.file.geometry_info.write(inf_obj, thy_save_path)
                ts_afs.ts.file.geometry_input.write(inp_str, thy_save_path)
                ts_afs.ts.file.energy.write(ene, thy_save_path)
                ts_afs.ts.file.geometry.write(geo, thy_save_path)
                ts_afs.ts.file.zmatrix.write(zma, thy_save_path)

                moldr.driver.run_job(
                    job='hessian',
                    script_str=SCRIPT_STR,
                    prefix=ts_run_path,
                    geom=geo,
                    charge=0,
                    mult=ts_mul,
                    method=method,
                    basis=basis,
                    orb_restr=orb_restr,
                    prog=PROG,
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

                    ts_afs.ts.file.hessian_info.write(inf_obj, thy_save_path)
                    ts_afs.ts.file.hessian_input.write(inp_str, thy_save_path)
                    ts_afs.ts.file.hessian.write(hess, thy_save_path)
                    ts_afs.ts.file.harmonic_frequencies.write(freqs, thy_save_path)

            if RUN_TS_TAU_SAMP:

                moldr.driver.save_tau(
                    run_prefix=thy_run_path,
                    save_prefix=thy_save_path,
                )

                zma = ts_afs.ts.file.zmatrix.read(thy_save_path)
                tors_ranges = automol.zmatrix.torsional_sampling_ranges(
                    zma, tors_names)
                tors_range_dct = dict(zip(tors_names, tors_ranges))

                moldr.driver.run_tau(
                    zma=zma,
                    charge=charge,
                    mult=ts_mul,
                    method=method,
                    basis=basis,
                    orb_restr=orb_restr,
                    nsamp=nsamp,
                    tors_range_dct=tors_range_dct,
                    run_prefix=thy_run_path,
                    save_prefix=thy_save_path,
                    script_str=SCRIPT_STR,
                    prog=PROG,
                    saddle=True,
                    overwrite=OVERWRITE,
                    **OPT_KWARGS,
                )

                moldr.driver.save_tau(
                    run_prefix=thy_run_path,
                    save_prefix=thy_save_path,
                )

            if RUN_TS_CONF_SCAN:
                zma = ts_afs.ts.file.zmatrix.read(thy_save_path)
                val_dct = automol.zmatrix.values(zma)
                tors_linspaces = automol.zmatrix.torsional_scan_linspaces(
                    zma, tors_names, SCAN_INCREMENT)
                tors_grids = [
                    numpy.linspace(*linspace) +val_dct[name]
                    for name, linspace in zip(tors_name, tors_linspaces)]
                for tors_name, tors_grid in zip(tors_names, tors_grids):
                    moldr.driver.run_scan(
                        zma=zma,
                        charge=charge,
                        mult=ts_mul,
                        method=method,
                        basis=basis,
                        orb_restr=orb_restr,
                        grid_dct={tors_name: tors_grid},
                        run_prefix=thy_run_path,
                        save_prefix=thy_save_path,
                        script_str=SCRIPT_STR,
                        prog=PROG,
                        saddle=True,
                        overwrite=OVERWRITE,
                        **OPT_KWARGS,
                    )

                    moldr.driver.save_scan(
                        run_prefix=thy_run_path,
                        save_prefix=thy_save_path,
                        coo_names=[tors_name],
                    )
                hind_rot_dct = {}
                scan_afs = autofile.fs.scan()
                min_ene = ts_afs.ts.file.energy.read(thy_save_path)
#                min_ene = cnf_afs.conf.file.energy.read(thy_save_path, min_cnf_alocs)
#                cnf_afs.conf_trunk.file.energy.write(min_ene, thy_save_path)
                for tors_name in tors_names:
                    enes = [scan_afs.scan.file.energy.read(thy_save_path, [[tors_name]] + rlocs)
                            for rlocs in scan_afs.scan.dir.existing(thy_save_path, [[tors_name]])]
                    enes = numpy.subtract(enes, min_ene)
                    hind_rot_dct[tors_name] = enes*EH2KCAL

                print('ts hindered rotor potential')
                print(hind_rot_dct)

            if RUN_TS_KICKS_QCHEM:
                ret = moldr.driver.read_job(job=elstruct.Job.HESSIAN, prefix=ts_run_path)
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
                        geo, disp_xyzs, charge, mult, method, basis, orb_restr,
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
                        geo, disp_xyzs, charge, mult, method, basis,
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
                charge=chg,
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
