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

# 0. choose which mechanism to run
MECHANISM_NAME = 'test'  # options: syngas, natgas, heptane, etc.

# 1. script control parameters
#PROG = 'psi4'
#SCRIPT_STR = ("#!/usr/bin/env bash\n"
#              "psi4 -i run.inp -o run.out >> stdout.log &> stderr.log")
#KWARGS = {}

PROG = 'g09'
SCRIPT_STR = ("#!/usr/bin/env bash\n"
              "g09 run.inp run.out >> stdout.log &> stderr.log")
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

METHOD = 'wb97xd'
BASIS = '6-31g*'
RESTRICT_OPEN_SHELL = False
NSAMP = 5
RUN_SPECIES = True
RUN_REACTIONS = True

RUN_GRADIENT = True
RUN_HESSIAN = True
RUN_CONFORMER_SCAN = True
SCAN_INCREMENT = 30. * qcc.conversion_factor('degree', 'radian')

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
SPC_TAB['inchi'] = list(map(automol.smiles.inchi, SPC_TAB['smiles']))
SPC_TAB['charge'] = 0
ICH_DCT = dict(zip(SPC_TAB['name'], SPC_TAB['inchi']))
CHG_DCT = dict(zip(SPC_TAB['name'], SPC_TAB['charge']))
MUL_DCT = dict(zip(SPC_TAB['name'], SPC_TAB['mult']))
SPC_BLK_STR = chemkin_io.species_block(MECH_STR)
SPC_NAMES = chemkin_io.species.names(SPC_BLK_STR)

GEOM_PATH = os.path.join(DATA_PATH, 'data', 'geoms')
GEOM_DCT = {}
for dir_path, _, file_names in os.walk(GEOM_PATH):
    for file_name in file_names:
        file_path = os.path.join(dir_path, file_name)
        if file_path.endswith('.xyz'):
            xyz_str = autofile.file.read_file(file_path)
            geo = automol.geom.from_xyz_string(xyz_str)
            ich = automol.geom.inchi(geo)
            if ich in GEOM_DCT:
                print('Warning: Dupilicate xyz geometry for ', ich)
            GEOM_DCT[ich] = geo


def inchi_to_geometry(ich):
    if ich in GEOM_DCT:
        geo = GEOM_DCT[ich]
    else:
        geo = automol.inchi.geometry(ich)
    return geo


if RUN_SPECIES:
    for name in SPC_NAMES:
        # species
        ich = ICH_DCT[name]
        charge = CHG_DCT[name]
        mult = MUL_DCT[name]
        print("Species: {}".format(name))

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

        thy_run_path = thy_afs.theory.dir.path(RUN_PREFIX, thy_alocs)
        thy_save_path = thy_afs.theory.dir.path(SAVE_PREFIX, thy_alocs)

        # a. conformer sampling
        # generate the z-matrix and sampling ranges

        geo = inchi_to_geometry(ich)

        zma = automol.geom.zmatrix(geo)
        tors_names = automol.geom.zmatrix_torsion_coordinate_names(geo)
        tors_ranges = automol.zmatrix.torsional_sampling_ranges(zma,
                                                                tors_names)
        tors_range_dct = dict(zip(tors_names, tors_ranges))

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
            nsamp=NSAMP,
            tors_range_dct=tors_range_dct,
            run_prefix=thy_run_path,
            save_prefix=thy_save_path,
            script_str=SCRIPT_STR,
            prog=PROG,
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
        cnf_enes = [cnf_afs.conf.file.energy.read(thy_save_path, alocs)
                    for alocs in cnf_alocs_lst]
        min_cnf_alocs = cnf_alocs_lst[cnf_enes.index(min(cnf_enes))]
        cnf_run_path = cnf_afs.conf.dir.path(thy_run_path, min_cnf_alocs)
        cnf_save_path = cnf_afs.conf.dir.path(thy_save_path, min_cnf_alocs)

        # generate the z-matrix and sampling grids (grids)
        geo = cnf_afs.conf.file.geometry.read(thy_save_path, min_cnf_alocs)
        zma = automol.geom.zmatrix(geo)
        tors_names = automol.geom.zmatrix_torsion_coordinate_names(geo)
        tors_linspaces = automol.zmatrix.torsional_scan_linspaces(
            zma, tors_names, SCAN_INCREMENT)
        tors_grids = [
            numpy.linspace(*linspace) for linspace in tors_linspaces]

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
                    **KWARGS,
                )

                moldr.driver.save_scan(
                    run_prefix=cnf_run_path,
                    save_prefix=cnf_save_path,
                    coo_names=[tors_name],
                )
            hind_rot_dct = {}
            scan_afs = autofile.fs.scan()
            min_ene = cnf_afs.conf.file.energy.read(thy_save_path, min_cnf_alocs)
            for tors_name in tors_names:
                enes = [scan_afs.scan.file.energy.read(cnf_save_path, [[tors_name]] + rlocs)
                    for rlocs in scan_afs.scan.dir.existing(cnf_save_path, [[tors_name]])]
                enes = numpy.subtract(enes, min_ene)
                hind_rot_dct[tors_name] = enes

            print(hind_rot_dct)

        hess = cnf_afs.conf.file.hessian.read(thy_save_path, min_cnf_alocs)
        freqs = elstruct.util.harmonic_frequencies(geo, hess)
        zpe = sum(freqs)/2.

#       hind_rot = 
# to be generalized
        symfactor = 1.
# in pyx2z but not in automol yet.
        elec_levels = [[mult, 0.]]
        
        core = mess_io.writer.write_core_rigidrotor(geo, symfactor)
        molecule_section_str1 = mess_io.writer.write_molecule(
            core, freqs, zpe, elec_levels,
            hind_rot='',
        ) 
        print (molecule_section_str1)

        sys.exit()

# 5. process reaction data from the mechanism file
RXN_BLOCK_STR = chemkin_io.reaction_block(MECH_STR)
RXN_STRS = chemkin_io.reaction.data_strings(RXN_BLOCK_STR)
RCT_NAMES_LST = list(
    map(chemkin_io.reaction.DataString.reactant_names, RXN_STRS))
PRD_NAMES_LST = list(
    map(chemkin_io.reaction.DataString.product_names, RXN_STRS))

if RUN_REACTIONS:
    for rct_names, prd_names in zip(RCT_NAMES_LST, PRD_NAMES_LST):
        # print the CHEMKIN reaction name for reference
        rxn_name = '='.join(['+'.join(rct_names), '+'.join(prd_names)])
        print()
        print("Reaction: {}".format(rxn_name))

        # determine inchis, charges, and multiplicities
        rct_ichs = list(map(ICH_DCT.__getitem__, rct_names))
        prd_ichs = list(map(ICH_DCT.__getitem__, prd_names))
        rct_chgs = list(map(CHG_DCT.__getitem__, rct_names))
        prd_chgs = list(map(CHG_DCT.__getitem__, prd_names))
        rct_muls = list(map(MUL_DCT.__getitem__, rct_names))
        prd_muls = list(map(MUL_DCT.__getitem__, prd_names))

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
            ts_zma, dist_name = ret

        ret = automol.zmatrix.ts.addition(rct_zmas, prd_zmas)
        if ret and typ is None:
            typ = 'addition'
            ts_zma, dist_name = ret

        ret = automol.zmatrix.ts.hydrogen_abstraction(rct_zmas, prd_zmas)
        if ret and typ is None:
            typ = 'hydrogen abstraction'
            ts_zma, dist_name = ret

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
            direction = autofile.system.reaction_direction(
                rxn_ichs, rxn_chgs, rxn_muls)
            rxn_ichs, rxn_chgs, rxn_muls = autofile.system.sort_together(
                rxn_ichs, rxn_chgs, rxn_muls)
            print(" - The reaction direction is {}"
                  .format('forward' if direction else 'backward'))
            rxn_alocs = [rxn_ichs, rxn_chgs, rxn_muls, ts_mul]
            thy_rlocs = [method, basis, orb_restr]
            thy_alocs = rxn_alocs + thy_rlocs

            rxn_afs = autofile.fs.reaction()
            thy_afs = autofile.fs.theory(rxn_afs, 'reaction')

            thy_afs.theory.dir.create(RUN_PREFIX, thy_alocs)
            thy_afs.theory.dir.create(SAVE_PREFIX, thy_alocs)

            thy_run_path = thy_afs.theory.dir.path(RUN_PREFIX, thy_alocs)
            thy_save_path = thy_afs.theory.dir.path(SAVE_PREFIX, thy_alocs)
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
                update_guess=False,
                reverse_sweep=False,
                **KWARGS
            )

            moldr.driver.save_scan(
                run_prefix=thy_run_path,
                save_prefix=thy_save_path,
                coo_names=[dist_name],
            )

sys.exit()
