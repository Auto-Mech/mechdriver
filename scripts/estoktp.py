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
print('EH2KCAL')
print(EH2KCAL)
#WAVEN2KCAL = 1./349.7

# 0. choose which mechanism to run
MECHANISM_NAME = 'onereac'  # options: syngas, natgas, heptane, test, estoktp, ...
#MECHANISM_NAME = 'estoktp/add30'  # options: syngas, natgas, heptane, test, estoktp, ...
#MECHANISM_NAME = 'estoktp/habs65'  # options: syngas, natgas, heptane, test, estoktp, ...

# 1. script control parameters
# Program, method, basis, convergence control

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
KWARGS = {
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

# What to run for electronic structure calculations
RUN_SPECIES_QCHEM = True
RUN_REACTIONS_QCHEM = False
RUN_VDW_QCHEM = False
KICKOFF_SADDLE = False
KICKOFF_BACKWARD = False
KICKOFF_SIZE = 0.1

RUN_CONFORMER_OPT = True
RUN_CONFORMER_SCAN = True
RUN_TAU_SAMPLING = False

RUN_TS_CONFORMER_OPT = True
RUN_TS_CONFORMER_SCAN = True
RUN_TS_TAU_SAMPLING = False

RUN_GRADIENT = True
RUN_HESSIAN = True

# What to run for thermochemical kinetics
RUN_SPECIES_PF = True
RUN_SPECIES_THERMO = False
RUN_REACTIONS_RATES = False
RUN_VDW_RCT_RATES = False
RUN_VDW_PRD_RATES = False

# Parameters for number of torsional samplings
NSAMP_A = 3
NSAMP_B = 1
NSAMP_C = 3
NSAMP_D = 15

# Parameters for partition function and kinetics
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
SMILES_LST = ['[O]', '[OH]', '[N]=O', '[CH2]', '[C]', '[B]', '[N]', '[F]', '[Cl]', '[Br]',
              '[BH]', '[BH2]', 'C[C]', '[O][O]']
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
#for smiles in SMILES_LST:
#    ich = automol.convert.smiles.inchi(smiles)
#    ich2 = automol.smiles.inchi(smiles)
#    print(ich)

SPC_TAB['charge'] = 0
SMI_DCT = dict(zip(SPC_TAB['name'], SPC_TAB['smiles']))
CHG_DCT = dict(zip(SPC_TAB['name'], SPC_TAB['charge']))
MUL_DCT = dict(zip(SPC_TAB['name'], SPC_TAB['mult']))
SPC_BLK_STR = chemkin_io.species_block(MECH_STR)
SPC_NAMES = chemkin_io.species.names(SPC_BLK_STR)

GEOM_PATH = os.path.join(DATA_PATH, 'data', 'geoms')
print(GEOM_PATH)
GEOM_DCT = {}
for dir_path, _, file_names in os.walk(GEOM_PATH):
    for file_name in file_names:
        file_path = os.path.join(dir_path, file_name)
        if file_path.endswith('.xyz'):
            xyz_str = autofile.file.read_file(file_path)
            print(file_path)
            print(xyz_str)
            geo = automol.geom.from_xyz_string(xyz_str)
            ich = automol.geom.inchi(geo)
            if ich in GEOM_DCT:
                print('Warning: Dupilicate xyz geometry for ', ich)
            GEOM_DCT[ich] = geo
            print(geo)

# take starting geometry from saved directory if possible, otherwise get it from inchi via rdkit

def inchi_to_geometry(ich):
    if ich in GEOM_DCT:
        geo = GEOM_DCT[ich]
    else:
        geo = automol.inchi.geometry(ich)
    return geo


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
        if RESTRICT_OPEN_SHELL:
            orb_restr = True
        else:
            orb_restr = (mult == 1)

        # set up the filesystem
        thy_run_path = moldr.util.species_theory_path(ich, charge, mult, method, basis, orb_restr, RUN_PREFIX)
        thy_save_path = moldr.util.species_theory_path(ich, charge, mult, method, basis, orb_restr, SAVE_PREFIX)

        # generate reference geometry
        # generate the z-matrix and sampling ranges

        geo = inchi_to_geometry(ich)
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
            **KWARGS,
        )

        if KICKOFF_SADDLE:
            print('Checking for saddle')
            ret = moldr.driver.read_job(job=elstruct.Job.OPTIMIZATION, prefix=thy_run_path)
            if ret:
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
                    freqs = elstruct.util.harmonic_frequencies(geo, hess, project=True)
                    norm_coos = elstruct.util.normal_coordinates(geo, hess, project=True)

                    # if there's an imaginary frequency, try again after displacing along
                    # the mode
                    if freqs[0] < -10:
                        print('Imaginary mode found: Attempting to kickoff from saddle')
                        im_norm_coo = numpy.array(norm_coos)[:, 0]
                        disp_len = KICKOFF_SIZE * qcc.conversion_factor('angstrom', 'bohr')
                        if KICKOFF_BACKWARD:
                            disp_len *= -1
                        disp_xyzs = numpy.reshape(im_norm_coo, (-1, 3)) * disp_len
                        geo = automol.geom.displaced(geo, disp_xyzs)
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
                            **KWARGS,
                        )

        ret = moldr.driver.read_job(job=elstruct.Job.OPTIMIZATION, prefix=thy_run_path)
        if ret:
            inf_obj, _, out_str = ret
            prog = inf_obj.prog
            geo = elstruct.reader.opt_geometry(prog, out_str)
            zma = automol.geom.zmatrix(geo)
#            save_path = afs.conf.dir.path(save_prefix, alocs)
#            print(" - Geometry is unique. Saving...")
#            print(" - Save path: {}".format(save_path))

#            afs.dir.create(save_prefix, alocs)
#            afs.file.geometry_info.write(
#                inf_obj, save_prefix, alocs)
#            afs.file.geometry_input.write(
#                inp_str, save_prefix, alocs)
#            afs.file.energy.write(ene, save_prefix, alocs)
#            afs.file.geometry.write(geo, save_prefix, alocs)
#            afs.file.zmatrix.write(zma, save_prefix, alocs)


        # a. conformer sampling
        tors_names = automol.geom.zmatrix_torsion_coordinate_names(geo)
        tors_ranges = automol.zmatrix.torsional_sampling_ranges(
            zma, tors_names)
        tors_range_dct = dict(zip(tors_names, tors_ranges))

        gra = automol.inchi.graph(ich)
        ntaudof = len(automol.graph.rotational_bond_keys(gra, with_h_rotors=False))
        nsamp = min(NSAMP_A + NSAMP_B * NSAMP_C**ntaudof, NSAMP_D)

        moldr.driver.save_conformers(
            run_prefix=thy_run_path,
            save_prefix=thy_save_path,
        )

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
            **KWARGS,
        )

        moldr.driver.save_conformers(
            run_prefix=thy_run_path,
            save_prefix=thy_save_path,
        )

        # get the list of saved conformers, and their filesystem
        cnf_afs = autofile.fs.conformer()
        cnf_alocs_lst = cnf_afs.conf.dir.existing(thy_save_path)

        # b. conformer gradients and hessians
        for cnf_alocs in cnf_alocs_lst:
            geo = cnf_afs.conf.file.geometry.read(thy_save_path, cnf_alocs)

            cnf_run_path = cnf_afs.conf.dir.path(thy_run_path, cnf_alocs)
            cnf_save_path = cnf_afs.conf.dir.path(thy_save_path, cnf_alocs)

            if RUN_GRADIENT:
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

                    save_path = cnf_afs.conf.dir.path(thy_save_path, cnf_alocs)
                    print(" - Saving gradient...")
                    print(" - Save path: {}".format(save_path))
                    print(cnf_afs.conf.dir.exists(thy_save_path, cnf_alocs))
                    cnf_afs.conf.file.gradient_info.write(
                        inf_obj, thy_save_path, cnf_alocs)
                    cnf_afs.conf.file.gradient_input.write(
                        inp_str, thy_save_path, cnf_alocs)
                    cnf_afs.conf.file.gradient.write(
                        grad, thy_save_path, cnf_alocs)

            if RUN_HESSIAN:
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
                    grad = elstruct.reader.hessian(inf_obj.prog, out_str)

                    save_path = cnf_afs.conf.dir.path(thy_save_path, cnf_alocs)
                    print(" - Saving hessian...")
                    print(" - Save path: {}".format(save_path))
                    print(cnf_afs.conf.dir.exists(thy_save_path, cnf_alocs))
                    cnf_afs.conf.file.hessian_info.write(
                        inf_obj, thy_save_path, cnf_alocs)
                    cnf_afs.conf.file.hessian_input.write(
                        inp_str, thy_save_path, cnf_alocs)
                    cnf_afs.conf.file.hessian.write(
                        grad, thy_save_path, cnf_alocs)

        # d. hindered rotor scans
        # determine the lowest energy conformer to get the correct path
#        cnf_enes = [cnf_afs.conf.file.energy.read(thy_save_path, alocs)
#                    for alocs in cnf_alocs_lst]
#        min_cnf_alocs = cnf_alocs_lst[cnf_enes.index(min(cnf_enes))]
        min_cnf_alocs = moldr.util.min_energy_conformer_locators(thy_save_path)
        if min_cnf_alocs is not None:
            cnf_run_path = cnf_afs.conf.dir.path(thy_run_path, min_cnf_alocs)
            cnf_save_path = cnf_afs.conf.dir.path(thy_save_path, min_cnf_alocs)

            # generate the z-matrix and sampling grids (grids)
            geo = cnf_afs.conf.file.geometry.read(thy_save_path, min_cnf_alocs)
            zma = automol.geom.zmatrix(geo)
            val_dct = automol.zmatrix.values(zma)
            tors_names = automol.geom.zmatrix_torsion_coordinate_names(geo)
            tors_linspaces = automol.zmatrix.torsional_scan_linspaces(
                zma, tors_names, SCAN_INCREMENT)
            print('SCAN_INCREMENT test')
            print(SCAN_INCREMENT)
            tors_grids = [
                numpy.linspace(*linspace) +val_dct[name]
                for name, linspace in zip(tors_names, tors_linspaces)]

            print(tors_grids)
            print(2*numpy.pi/12.)

            # run one-dimensional scans for each torsional coordinate
            if RUN_CONFORMER_SCAN:
                for tors_name, tors_grid in zip(tors_names, tors_grids):
                    moldr.driver.run_scan(
                        zma=zma,
                        charge=charge,
                        mult=mult,
                        method=method,
                        basis=basis,
                        orb_restr=orb_restr,
                        grid_dct={tors_name: tors_grid},
                        run_prefix=cnf_run_path,
                        save_prefix=cnf_save_path,
                        script_str=SCRIPT_STR,
                        prog=PROG,
                        overwrite=OVERWRITE,
                        **KWARGS,
                    )

                    moldr.driver.save_scan(
                        run_prefix=cnf_run_path,
                        save_prefix=cnf_save_path,
                        coo_names=[tors_name],
                    )

        if RUN_TAU_SAMPLING:
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
                **KWARGS,
            )

            moldr.driver.save_tau(
                run_prefix=thy_run_path,
                save_prefix=thy_save_path,
            )

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
        if RESTRICT_OPEN_SHELL:
            orb_restr = True
        else:
            orb_restr = (mult == 1)
# read in geometry, hess and hindered rotor potentials for minimum energy conformer
        thy_save_path = moldr.util.species_theory_path(
            ich, charge, mult, method, basis, orb_restr, SAVE_PREFIX)
        cnf_afs = autofile.fs.conformer()
        min_cnf_alocs = moldr.util.min_energy_conformer_locators(thy_save_path)
# I think we need something for if it is none
        if min_cnf_alocs is not None:
            geo = cnf_afs.conf.file.geometry.read(thy_save_path, min_cnf_alocs)
            hess = cnf_afs.conf.file.hessian.read(thy_save_path, min_cnf_alocs)
            zpe = sum(freqs)*WAVEN2KCAL/2.
            zma = automol.geom.zmatrix(geo)
            gra = automol.zmatrix.graph(zma, remove_stereo=True)
            scan_afs = autofile.fs.scan()
            min_ene = cnf_afs.conf.file.energy.read(thy_save_path, min_cnf_alocs)
            cnf_afs.conf_trunk.file.energy.write(min_ene, thy_save_path)
            cnf_save_path = cnf_afs.conf.dir.path(thy_save_path, min_cnf_alocs)
            tors_names = automol.geom.zmatrix_torsion_coordinate_names(geo)
            coo_dct = automol.zmatrix.coordinates(zma, multi=False)
            group_dct = {}
            axis_dct = {}
            sym_dct = {}
            pot_dct = {}
            for tors_name in tors_names:
                enes = [scan_afs.scan.file.energy.read(cnf_save_path, [[tors_name]] + rlocs)
                        for rlocs in scan_afs.scan.dir.existing(cnf_save_path, [[tors_name]])]
                enes = numpy.subtract(enes, min_ene)
                pot_dct[tors_name] = enes*EH2KCAL
                axis = coo_dct[tors_name][1:3]
                axis_dct[tors_name] = numpy.add(axis, 1)
                group = list(automol.graph.branch_atom_keys(gra, axis[0], axis) - set(axis))
                group_dct[tors_name] = numpy.add(group, 1)
                sym_dct[tors_name] = 1
                hind_rot_dct = {}
                hind_rot_dct[tors_name] = mess.io.writer.write_hind_rot(
                    group_dct[tors_name], axis_dct[tors_name],
                    sym_dct[tors_name], pot_dct[tors_name])

            print('Species hindered rotor potential')
            print(pot_dct)
            print(group_dct)
            print(axis_dct)

#            hind_rot_key_str = '/n'.join([potential_dct,group_dct,axis_dct])
#            print(hind_rot_key_str)


# set up messpf input
        elec_levels = [[0., mult]]
        if (ich, mult) in ELC_DEG_DCT:
            elec_levels = ELC_DEG_DCT[(ich, mult)]
        freqs = elstruct.util.harmonic_frequencies(geo, hess)
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
            species_str = mess_io.writer.write_atom(name, elec_levels)
        else:
            if len(automol.geom.symbols(geo)) == 2:
                freq_offset = 5
            else:
                freq_offset = 6
            core = mess_io.writer.write_core_rigidrotor(geo, symfactor)
            if pot_dct is not None:
                species_str = mess_io.writer.write_molecule(
                    core, freqs[freq_offset:], zpe, elec_levels,
                    hind_rot='hind_rot_key_str',
                )
        print(species_str)
        pf_inp_str = '\n'.join([global_pf_str, species_head_str, species_str])
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
        if RESTRICT_OPEN_SHELL:
            orb_restr = True
        else:
            orb_restr = (mult == 1)
# read in geometry, hess and hindered rotor potentials for minimum energy conformer
        thy_save_path = moldr.util.species_theory_path(
            ich, charge, mult, method, basis, orb_restr, SAVE_PREFIX)
        bld_afs = autofile.fs.build()
        bld_alocs = ['NASA_POLY', 0]
        bld_afs.build.dir.create(thy_save_path, bld_alocs)
        path = bld_afs.build.dir.path(thy_save_path, bld_alocs)
        print('Build Path for NASA Polynomials')
        print(path)
        bld_afs.build.file.input.write(nasa_inp_str, thy_save_path, bld_alocs)
        moldr.util.run_script(NASA_SCRIPT_STR, path)

sys.exit()


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

        # determine the transition state z-matrix
        rct_zmas = list(
            map(automol.geom.zmatrix, map(inchi_to_geometry, rct_ichs)))
        prd_zmas = list(
            map(automol.geom.zmatrix, map(inchi_to_geometry, prd_ichs)))

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

        ret = automol.zmatrix.ts.hydrogen_abstraction(rct_zmas, prd_zmas)
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

            # determine the transition state multiplicity
            ts_mul = automol.mult.ts.low(rct_muls, prd_muls)

            # theory
            method = METHOD
            basis = BASIS
            if RESTRICT_OPEN_SHELL:
                orb_restr = True
            else:
                orb_restr = (ts_mul == 1)

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

            thy_run_path = moldr.util.reaction_theory_path(
                rxn_ichs, rxn_chgs, rxn_muls, ts_mul, method, basis, orb_restr,
                RUN_PREFIX)
            thy_save_path = moldr.util.reaction_theory_path(
                rxn_ichs, rxn_chgs, rxn_muls, ts_mul, method, basis, orb_restr,
                SAVE_PREFIX)
            moldr.driver.run_scan(
                zma=ts_zma,
                charge=0,
                mult=ts_mul,
                method=method,
                basis=basis,
                orb_restr=orb_restr,
                grid_dct={dist_name: grid},
                run_prefix=thy_run_path,
                save_prefix=thy_save_path,
                script_str=SCRIPT_STR,
                prog=PROG,
                overwrite=OVERWRITE,
                update_guess=False,
                reverse_sweep=False,
                **KWARGS
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
                orb_restr=orb_restr,
                prog=PROG,
                saddle=True,
                overwrite=OVERWRITE,
                **KWARGS,
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

                print(" - Running hessian")
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

            if RUN_TS_TAU_SAMPLING:

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
                    **KWARGS,
                )

                moldr.driver.save_tau(
                    run_prefix=thy_run_path,
                    save_prefix=thy_save_path,
                )

            if RUN_TS_CONFORMER_SCAN:
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
                        **KWARGS,
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

        sys.exit()

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
