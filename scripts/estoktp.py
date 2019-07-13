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
#WAVEN2KCAL = 1./349.7

# 0. choose which mechanism to run
MECHANISM_NAME = 'estoktp/add30'  # options: syngas, natgas, heptane, test, estoktp, ...
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

# What to run
RUN_SPECIES = True
RUN_REACTIONS = True
RUN_GRADIENT = True
RUN_HESSIAN = True
RUN_CONFORMER_SCAN = True
RUN_TAU_SAMPLING = True

# Parameters for number of torsional samplings
NSAMP_A = 3
NSAMP_B = 1
NSAMP_C = 3
NSAMP_D = 15

# Defaults
SCAN_INCREMENT = 30. * qcc.conversion_factor('degree', 'radian')
RESTRICT_OPEN_SHELL = False
OVERWRITE = False
RUN_OVER = False

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
#    print(ich)
#SPC_TAB['inchi'] = list(map(automol.smiles.inchi, SPC_TAB['smiles']))

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


SPC_AFS = autofile.fs.species()
THY_AFS = autofile.fs.theory(SPC_AFS, 'species')

if RUN_SPECIES:
    for name in SPC_NAMES:
        # species
        print("Species: {}".format(name))
        smi = SMI_DCT[name]
        ich = automol.smiles.inchi(smi)
        print("smiles: {}".format(smi),"inchi: {}".format(ich))
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
        spc_alocs = [ich, charge, mult]         # aloc = absolute locator
        thy_rlocs = [method, basis, orb_restr]  # rloc = relative locator
        thy_alocs = spc_alocs + thy_rlocs
        print(RUN_PREFIX)
        print(thy_alocs)
        THY_AFS.theory.dir.create(RUN_PREFIX, thy_alocs)
        THY_AFS.theory.dir.create(SAVE_PREFIX, thy_alocs)

        thy_run_path = THY_AFS.theory.dir.path(RUN_PREFIX, thy_alocs)
        thy_save_path = THY_AFS.theory.dir.path(SAVE_PREFIX, thy_alocs)

        # a. conformer sampling
        # generate the z-matrix and sampling ranges

        geo = inchi_to_geometry(ich)

        zma = automol.geom.zmatrix(geo)
        tors_names = automol.geom.zmatrix_torsion_coordinate_names(geo)
        tors_ranges = automol.zmatrix.torsional_sampling_ranges(
            zma, tors_names)
        tors_range_dct = dict(zip(tors_names, tors_ranges))

        gra = automol.inchi.graph(ich)
        ntaudof = len(automol.graph.rotational_bond_keys(gra, with_h_rotors=False))
        print ('nsamp generation')
        print(ntaudof)
        nsamp = min(NSAMP_A + NSAMP_B * NSAMP_C**ntaudof, NSAMP_D)
        print(nsamp)

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
            run_over=RUN_OVER,
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
                    run_over=RUN_OVER,
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
                    run_over=RUN_OVER,
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
                    overwrite=OVERWRITE,
                    run_over=RUN_OVER,
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
            cnf_afs.conf_trunk.file.energy.write(min_ene, thy_save_path)
            for tors_name in tors_names:
                enes = [scan_afs.scan.file.energy.read(cnf_save_path, [[tors_name]] + rlocs)
                    for rlocs in scan_afs.scan.dir.existing(cnf_save_path, [[tors_name]])]
                enes = numpy.subtract(enes, min_ene)
                hind_rot_dct[tors_name] = enes

            print(hind_rot_dct)

        if RUN_TAU_SAMPLING:
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
                run_over=RUN_OVER,
                **KWARGS,
            )

            moldr.driver.save_tau(
                run_prefix=thy_run_path,
                save_prefix=thy_save_path,
            )

        hess = cnf_afs.conf.file.hessian.read(thy_save_path, min_cnf_alocs)
        freqs = elstruct.util.harmonic_frequencies(geo, hess)
        zpe = sum(freqs)*WAVEN2KCAL/2.

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

        rct_smis = list(map(SMI_DCT.__getitem__, rct_names))
        prd_smis = list(map(SMI_DCT.__getitem__, prd_names))
        rct_ichs = list(map(automol.smiles.inchi,rct_smis))
        prd_ichs = list(map(automol.smiles.inchi,prd_smis))
        rct_chgs = list(map(CHG_DCT.__getitem__, rct_names))
        prd_chgs = list(map(CHG_DCT.__getitem__, prd_names))
        rct_muls = list(map(MUL_DCT.__getitem__, rct_names))
        prd_muls = list(map(MUL_DCT.__getitem__, prd_names))

        # check direction of reaction
        rxn_ichs = [rct_ichs, prd_ichs]
        rxn_chgs = [rct_chgs, prd_chgs]
        rxn_muls = [rct_muls, prd_muls]
        rxn_exo = moldr.util.reaction_energy(SAVE_PREFIX, rxn_ichs, rxn_chgs, rxn_muls, method, basis, RESTRICT_OPEN_SHELL)
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
            THY_AFS = autofile.fs.theory(rxn_afs, 'reaction')

            THY_AFS.theory.dir.create(RUN_PREFIX, thy_alocs)
            THY_AFS.theory.dir.create(SAVE_PREFIX, thy_alocs)

            thy_run_path = THY_AFS.theory.dir.path(RUN_PREFIX, thy_alocs)
            thy_save_path = THY_AFS.theory.dir.path(SAVE_PREFIX, thy_alocs)
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
                run_over=RUN_OVER,
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
                run_over=RUN_OVER,
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
                    run_over=RUN_OVER,
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


sys.exit()
