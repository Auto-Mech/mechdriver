""" Generate the information necessary to product the vrctst input files
"""

import os
import autofile
import automol
import elstruct
import varecof_io
from routines.es._routines import _scan as scan
from routines.es._routines import _wfn as wfn
from routines.es import runner as es_runner
from lib import filesys
from lib.submission import run_script
from lib.submission import DEFAULT_SCRIPT_DCT
from lib.phydat import phycon


# CENTRAL FUNCTION TO WRITE THE VARECOF INPUT FILES AND RUN THE PROGRAM
def calc_vrctst_flux(ts_zma, ts_formula, ts_info, ts_dct, spc_dct,
                     high_mul, grid1, grid2, dist_name,
                     multi_level, multi_info, multi_sp_info,
                     ini_thy_info, thy_info,
                     thy_run_path, thy_save_path,
                     overwrite, update_guess,
                     run_prefix, save_prefix,
                     vrc_dct,
                     corr_pot=True):
    """ Set up n VRC-TST calculations to get the flux file
    """

    print('TS DCT')
    print(ts_dct.keys())

    # Set the active space
    num_act_orb, num_act_elc = wfn.active_space(
        ts_dct, spc_dct, ts_dct['high_mul'])

    # Input stuff
    # sp_thy_level = ['molpro2015', 'caspt2c', 'cc-pvdz', 'RR']
    # fortran_compiler = 'gfortran'
    bnd_frm_idxs = automol.zmatrix.coord_idxs(ts_zma, dist_name)
    min_idx, max_idx = min(bnd_frm_idxs), max(bnd_frm_idxs)
    bnd_frm_idxs = (bnd_frm_idxs[0]+1, bnd_frm_idxs[1]+1)
    # spc_name = 'mol'
    # memory = 4.0
    # basis = 'cc-pvdz'
    # method = '{rs2c, shift=0.25}'
    # r1dists_lr = [8., 6., 5., 4.5, 4.]
    # r1dists_sr = [4., 3.8, 3.6, 3.4, 3.2, 3., 2.8, 2.6, 2.4, 2.2]
    # r2dists_sr = [4., 3.8, 3.6, 3.4, 3.2, 3., 2.8, 2.6, 2.4, 2.2]
    # d1dists = [0.01, 0.5, 1.]
    # d2dists = [0.01, 0.5, 1.]
    # conditions = {}
    # nsamp_max = 2000
    # nsamp_min = 50
    # flux_err = 10
    # pes_size = 2
    # exe_path = '/blues/gpfs/home/sjklipp/bin/molpro'
    # base_name = 'mol'

    # Set the run directory
    ts_run_fs = autofile.fs.ts(thy_run_path)
    ts_run_fs[0].create()
    ts_run_path = ts_run_fs[0].path()

    # Set the scan directories
    scn_run_fs = autofile.fs.scan(thy_run_path)
    scn_save_fs = autofile.fs.scan(thy_save_path)
    cscn_run_fs = autofile.fs.cscan(thy_run_path)
    cscn_save_fs = autofile.fs.cscan(thy_save_path)

    # Build the VRC-TST run directory, including the needed scr dir
    bld_locs = ['VARECOF', 0]
    bld_save_fs = autofile.fs.build(ts_run_path)
    bld_save_fs[-1].create(bld_locs)
    vrc_path = bld_save_fs[-1].path(bld_locs)
    os.makedirs(os.path.join(vrc_path, 'scratch'), exist_ok=True)
    print('Build Path for VaReCoF calculations')
    print(vrc_path)

    # Set the orb type
    if multi_sp_info is not None:
        orb_restr = filesys.inf.orbital_restriction(ts_info, multi_sp_info)
        sp_level = multi_sp_info[0:3]
        sp_level.append(orb_restr)

    # Set up the spc infor for the reactants
    rcts = ts_dct['reacs']
    spc_1_info = [spc_dct[rcts[0]]['inchi'],
                  spc_dct[rcts[0]]['chg'],
                  spc_dct[rcts[0]]['mul']]
    spc_2_info = [spc_dct[rcts[1]]['inchi'],
                  spc_dct[rcts[1]]['chg'],
                  spc_dct[rcts[1]]['mul']]

    max_grid = max([max(grid1), max(grid2)])
    locs = [[dist_name], [max_grid]]
    inf_sep_ene = scan.infinite_separation_energy(
        spc_1_info, spc_2_info, ts_info, high_mul, ts_zma,
        ini_thy_info, thy_info, multi_info,
        run_prefix, save_prefix, scn_run_fs, scn_save_fs, locs,
        overwrite=False, num_act_elc=None, num_act_orb=None)

    # Calculate the correction potential
    if corr_pot:

        # Build the constraint dictionary
        constraint_dct = set_alt_constraints(ts_zma, ts_dct['rct_zmas'])

        # Run the potentials
        run_potentials(
            ts_zma,
            ts_formula, ts_info, high_mul,
            multi_level,
            dist_name, grid1, grid2,
            scn_run_fs, scn_save_fs,
            cscn_run_fs, cscn_save_fs,
            overwrite, update_guess,
            num_act_elc, num_act_orb,
            constraint_dct,
            sp_thy_level=sp_level)

        # Combine and sort the grids for organization
        full_grid = list(grid1) + list(grid2)
        full_grid.sort()

        # Read the values for the correction potential from filesystem
        potentials, pot_labels = read_potentials(
            scn_save_fs, cscn_save_fs,
            sp_level, dist_name, full_grid,
            constraint_dct)

        # Build correction potential .so file used by VaReCoF
        build_correction_potential(
            full_grid, potentials,
            bnd_frm_idxs, vrc_dct['fortran_compiler'], vrc_path,
            dist_restrict_idxs=(),
            pot_labels=pot_labels,
            pot_file_names=[vrc_dct['spc_name']],
            spc_name=vrc_dct['spc_name'])

    # Write the electronic structure template file
    memory, basis, method = '', '', ''
    tml_inp_str = write_molpro_template_str(
        ts_zma, ts_info, ts_dct['high_mul'],
        memory, basis, method, inf_sep_ene)

    # Write the remaining VaReCoF input file strings
    input_strs = input_prep(
        ts_zma, ts_dct['rct_zmas'], len(potentials),
        min_idx, max_idx,
        vrc_dct['r1dists_lr'], vrc_dct['r1dists_sr'], vrc_dct['r2dists_sr'],
        vrc_dct['d1dists'], vrc_dct['d2dists'],
        vrc_dct['conditions'],
        vrc_dct['nsamp_max'], vrc_dct['nsamp_min'],
        vrc_dct['flux_err'], vrc_dct['pes_size'],
        vrc_dct['base_name'], vrc_dct['exe_path'], vrc_path)
    [struct_inp_str, lr_divsur_inp_str, tst_inp_str,
     els_inp_str, mc_flux_inp_str, conv_inp_str] = input_strs

    # Write machines file to set compute nodes
    machine_file_str = write_machinefile_str()

    # Write the files for VaReCoF in the build directory
    with open(os.path.join(vrc_path, 'structure.inp'), 'w') as inp_file:
        inp_file.write(struct_inp_str)
    with open(os.path.join(vrc_path, 'lr_divsur.inp'), 'w') as inp_file:
        inp_file.write(lr_divsur_inp_str)
    with open(os.path.join(vrc_path, 'tst.inp'), 'w') as inp_file:
        inp_file.write(tst_inp_str)
    with open(os.path.join(vrc_path, 'molpro.inp'), 'w') as inp_file:
        inp_file.write(els_inp_str)
    with open(os.path.join(vrc_path, 'mc_flux.inp'), 'w') as inp_file:
        inp_file.write(mc_flux_inp_str)
    with open(os.path.join(vrc_path, 'convert.inp'), 'w') as inp_file:
        inp_file.write(conv_inp_str)
    with open(os.path.join(vrc_path, 'mol.tml'), 'w') as inp_file:
        inp_file.write(tml_inp_str)
    with open(os.path.join(vrc_path, 'machines'), 'w') as inp_file:
        inp_file.write(machine_file_str)

    # Run VaReCoF
    run_script(DEFAULT_SCRIPT_DCT['varecof'], vrc_path)

    # Calculate the flux file from the output
    print('Generating flux file with TS N(E) from VaReCoF output...')
    run_script(DEFAULT_SCRIPT_DCT['mcflux'], vrc_path)

    # Check for success of the VaReCoF run and flux file generation
    # mcflux_file = os.path.join(vrc_path, 'mc_flux.out')
    # if os.path.exists(mcflux_file):
    #     with open(mcflux_file, 'r') as fluxfile:
    #         flux_str = mcflux_file.read()
    #     if flux_str != '':
    #         ts_found = True

    #         # Set the run directory
    #         ts_save_fs = autofile.fs.ts(thy_save_path)
    #         ts_save_fs[0].create()
    #         ts_save_path = ts_save_fs[0].path()

    #         # Write the VRC-Flux and info files
    #         ts_save_fs[0].file.vrc_flux.write(flux_str)
    #         ts_save_fs[0].file.vrc_flux_info.write(flux_str)

    # return ts_found


# FUNCTIONS TO WRITE THE STRINGS FOR ALL OF THE VARECOF INPUT FILE
def input_prep(ts_zma, rct_zmas, npot, min_idx, max_idx,
               r1dists_lr, r1dists_sr, r2dists_sr,
               d1dists, d2dists,
               conditions,
               nsamp_max, nsamp_min, flux_err, pes_size,
               base_name, exe_path, vrc_path):
    """ prepare all the input files for a vrc-tst calculation
    """

    # Build geometries needed for the varecof run
    total_geom, frag_geoms, frag_geoms_wdummy = fragment_geometries(
        ts_zma, rct_zmas, min_idx, max_idx)

    # Set information for the pivot points needed in divsur.inp
    frames, npivots = build_pivot_frames(
        min_idx, max_idx, total_geom, frag_geoms, frag_geoms_wdummy)
    pivot_angles = calc_pivot_angles(frag_geoms, frag_geoms_wdummy, frames)
    pivot_xyzs = calc_pivot_xyzs(min_idx, max_idx, total_geom, frag_geoms)

    # Write the long- and short-range divsur input files
    lr_divsur_inp_str = varecof_io.writer.input_file.divsur(
        r1dists_lr, 1, 1, [0.0, 0.0, 0.0], [0.0, 0.0, 0.0])

    # Write the short-range divsur files
    t1angs = [pivot_angles[0]] if pivot_angles[0] is not None else []
    t2angs = [pivot_angles[1]] if pivot_angles[1] is not None else []
    if automol.geom.is_atom(frag_geoms[0]):
        d1dists = []
        t1angs = []
    if automol.geom.is_atom(frag_geoms[1]):
        d2dists = []
        t2angs = []
    if automol.geom.is_linear(frag_geoms[0]):
        d1dists = [0.]
        t1angs = []
    if automol.geom.is_linear(frag_geoms[1]):
        d2dists = [0.]
        t2angs = []
    sr_divsur_inp_str = varecof_io.writer.input_file.divsur(
        r1dists_sr, npivots[0], npivots[1], pivot_xyzs[0], pivot_xyzs[1],
        frame1=frames[0], frame2=frames[1],
        d1dists=d1dists, d2dists=d2dists,
        t1angs=t1angs, t2angs=t2angs,
        r2dists=r2dists_sr,
        **conditions)

    # Write the structure input files
    struct_inp_str = varecof_io.writer.input_file.structure(
        frag_geoms_wdummy[0], frag_geoms_wdummy[1])
    with open(os.path.join(vrc_path, 'structure.inp'), 'w') as inp_file:
        inp_file.write(struct_inp_str)

    # Write the divsur input file with determined frames
    with open(os.path.join(vrc_path, 'divsur.inp'), 'w') as inp_file:
        inp_file.write(sr_divsur_inp_str)

    # Obtain the divsur.out file with divsur-frame fragment geoms
    divsur_out_str = build_divsur_out_file(vrc_path, os.getcwd())

    # Write the tst.inp file
    faces, faces_symm = assess_face_symmetries(divsur_out_str)
    tst_inp_str = varecof_io.writer.input_file.tst(
        nsamp_max, nsamp_min, flux_err, pes_size,
        faces=faces, faces_symm=faces_symm)

    # Write the potential energy surface input string
    els_inp_str = varecof_io.writer.input_file.elec_struct(
        exe_path, vrc_path, base_name, npot,
        dummy_name='dummy_corr_', lib_name='libcorrpot.so',
        geom_ptt='GEOMETRY_HERE', ene_ptt='molpro_energy')

    # Write the mc_flux.inp input string
    mc_flux_inp_str = varecof_io.writer.input_file.mc_flux()

    # Write the convert.inp input string
    conv_inp_str = varecof_io.writer.input_file.convert()

    # Collate the input strings together
    input_strs = [
        struct_inp_str, lr_divsur_inp_str, tst_inp_str,
        els_inp_str, mc_flux_inp_str, conv_inp_str]

    return input_strs


def write_molpro_template_str(ts_zma, ts_info, high_mul,
                              memory, basis, method, inf_sep_ene):
    """ Write the electronic structure template file
    """
    num_act_elc = high_mul
    num_act_orb = num_act_elc
    ts_formula = automol.geom.formula(automol.zmatrix.geometry(ts_zma))
    wfn_str = wfn.wfn_string(
        ts_info, ts_formula, num_act_elc, num_act_orb,
        high_mul, add_two_closed=False)
    tml_inp_str = varecof_io.writer.input_file.tml(
        memory, basis, wfn_str, method, inf_sep_ene)

    return tml_inp_str


def write_machinefile_str():
    """ Take machine list and write the string for the machine file
    """
    machines = ['b450:8', 'b451:8', 'b452:8', 'b453:8']
    machine_file_str = ''
    for machine in machines:
        machine_file_str += machine + '\n'

    return machine_file_str


# FUNCTIONS TO SET UP THE libcorrpot.so FILE USED BY VARECOF
def build_correction_potential(mep_distances, potentials,
                               bnd_frm_idxs, fortran_compiler, vrc_path,
                               dist_restrict_idxs=(),
                               pot_labels=(),
                               pot_file_names=(),
                               spc_name='mol'):
    """  use the MEP potentials to compile the correction potential .so file
    """

    # Change the coordinates of the MEP distances
    mep_distances = [dist * phycon.BOHR2ANG for dist in mep_distances]

    # Build string Fortan src file containing correction potentials
    species_corr_str = varecof_io.writer.corr_potentials.species(
        mep_distances,
        potentials,
        bnd_frm_idxs,
        dist_restrict_idxs=dist_restrict_idxs,
        pot_labels=pot_labels,
        species_name=spc_name)

    # Build string dummy corr file where no correction used
    dummy_corr_str = varecof_io.writer.corr_potentials.dummy()

    # Build string for auxiliary file needed for spline fitting
    pot_aux_str = varecof_io.writer.corr_potentials.auxiliary()

    # Build string for makefile to compile corr pot file into .so file
    makefile_str = varecof_io.writer.corr_potentials.makefile(
        fortran_compiler, pot_file_names=pot_file_names)

    # Write all of the files needed to build the correction potential
    with open(os.path.join(vrc_path, spc_name+'_corr.f'), 'w') as corr_file:
        corr_file.write(species_corr_str)
    with open(os.path.join(vrc_path, 'dummy_corr.f'), 'w') as corr_file:
        corr_file.write(dummy_corr_str)
    with open(os.path.join(vrc_path, 'pot_aux.f'), 'w') as corr_file:
        corr_file.write(pot_aux_str)
    with open(os.path.join(vrc_path, 'makefile'), 'w') as corr_file:
        corr_file.write(makefile_str)

    # Compile the correction potential
    varecof_io.writer.corr_potentials.compile_corr_pot(vrc_path)


def run_potentials(inf_sep_zma,
                   ts_formula, ts_info, high_mul,
                   multi_level,
                   dist_name, grid1, grid2,
                   scn_run_fs, scn_save_fs,
                   cscn_run_fs, cscn_save_fs,
                   overwrite, update_guess,
                   num_act_elc, num_act_orb,
                   constraint_dct,
                   sp_thy_level=None):
    """ Run all of the scans
    """

    # Run and save the scan while constraining only reaction coordinate
    scan.run_multiref_rscan(
        formula=ts_formula,
        high_mul=high_mul,
        zma=inf_sep_zma,
        spc_info=ts_info,
        multi_level=multi_level,
        dist_name=dist_name,
        grid1=grid1,
        grid2=grid2,
        scn_run_fs=scn_run_fs,
        scn_save_fs=scn_save_fs,
        overwrite=overwrite,
        update_guess=update_guess,
        gradient=False,
        hessian=False,
        num_act_elc=num_act_elc,
        num_act_orb=num_act_orb,
        constraint_dct=None
    )

    scan.save_scan(
        scn_run_fs=scn_run_fs,
        scn_save_fs=scn_save_fs,
        coo_names=[dist_name],
        thy_info=multi_level
    )

    # Run and save the scan while constraining all intermolecular coordinates
    scan.run_multiref_rscan(
        formula=ts_formula,
        high_mul=high_mul,
        zma=inf_sep_zma,
        spc_info=ts_info,
        multi_level=multi_level,
        dist_name=dist_name,
        grid1=grid1,
        grid2=grid2,
        scn_run_fs=cscn_run_fs,
        scn_save_fs=cscn_save_fs,
        overwrite=overwrite,
        update_guess=update_guess,
        gradient=False,
        hessian=False,
        num_act_elc=num_act_elc,
        num_act_orb=num_act_orb,
        constraint_dct=constraint_dct
    )

    scan.save_cscan(
        cscn_run_fs=cscn_run_fs,
        cscn_save_fs=cscn_save_fs,
        coo_names=[dist_name],
        thy_info=multi_level)

    # Run the single points on top of the initial scan
    if sp_thy_level is not None:
        scan_sp(inf_sep_zma, ts_info, ts_formula, scn_run_fs, scn_save_fs,
                multi_level, sp_thy_level, [dist_name], overwrite,
                num_act_elc, num_act_orb, high_mul)


def scan_sp(ts_zma, ts_info, ts_formula, scn_run_fs, scn_save_fs,
            multi_level, sp_thy_level, dist_name, overwrite,
            num_act_elc, num_act_orb, high_mul):
    """ get sps for the scan
    """

    # Set up script and kwargs for the irc run
    script_str, _, kwargs, _ = es_runner.par.run_qchem_par(*sp_thy_level[0:2])

    # Build the elstruct CASSCF options list for multiref calcs
    cas_opt = []
    cas_opt.append(
        wfn.cas_options(
            ts_info, ts_formula, num_act_elc, num_act_orb,
            add_two_closed=False))
    cas_opt.append(
        wfn.cas_options(
            ts_info, ts_formula, num_act_elc, num_act_orb,
            add_two_closed=True))

    # Write the lines containing all the calcs for a guess wfn
    guess_str = wfn.multiref_wavefunction_guess(
        high_mul, ts_zma, ts_info, multi_level, cas_opt)
    guess_lines = guess_str.splitlines()

    # Add the above-built objects to the elstruct opt_kwargs dct
    kwargs['casscf_options'] = cas_opt[1]
    kwargs['gen_lines'] = {1: guess_lines}

    # Add option to opt_kwargs dct to turn off symmetry
    kwargs['mol_options'] = ['nosym']

    # Compute the single-point energies along the scan
    for locs in scn_run_fs[-1].existing([dist_name]):

        # Set up single point filesys
        scn_run_fs[-1].create(locs)
        scn_run_path = scn_run_fs[-1].path(locs)
        scn_save_path = scn_save_fs[-1].path(locs)
        print('scn_run_path')
        print(scn_run_path)

        # Set up the run filesys for the job
        sp_run_fs = autofile.fs.single_point(scn_run_path)
        sp_save_fs = autofile.fs.single_point(scn_save_path)
        sp_run_fs[-1].create(sp_thy_level[1:4])
        sp_run_path = sp_run_fs[-1].path(sp_thy_level[1:4])
        run_fs = autofile.fs.run(sp_run_path)

        # Read the geometry from the save filesys
        exists = sp_save_fs[-1].file.energy.exists(sp_thy_level[1:4])
        if not exists or overwrite:
            geo = scn_save_fs[-1].file.geometry.read(locs)
            es_runner.run_job(
                job='energy',
                script_str=script_str,
                run_fs=run_fs,
                geom=geo,
                spc_info=ts_info,
                thy_level=sp_thy_level,
                overwrite=overwrite,
                **kwargs,
            )

        ret = es_runner.read_job(
            job='energy',
            run_fs=run_fs,
        )

        if ret is not None:
            inf_obj, inp_str, out_str = ret

            print(" - Reading energy from output...")
            ene = elstruct.reader.energy(inf_obj.prog, inf_obj.method, out_str)

            print(" - Saving energy...")
            sp_save_fs[-1].create(sp_thy_level[1:4])
            sp_save_fs[-1].file.input.write(inp_str, sp_thy_level[1:4])
            sp_save_fs[-1].file.info.write(inf_obj, sp_thy_level[1:4])
            sp_save_fs[-1].file.energy.write(ene, sp_thy_level[1:4])


def set_alt_constraints(inf_sep_zma, rct_zmas):
    """ Set the additional constraints for the constrained MEP
    """

    frag1_natom = automol.zmatrix.count(rct_zmas[0])
    frag2_natom = automol.zmatrix.count(rct_zmas[1])

    # Build pairs for intermolecular coords to optimize:
    #   (zma_atom_idx, coord_idx in row) (uses 0-indexing)
    # frag1_natom, 0 is the scan coord already accounted for
    no_frz_idxs = []
    no_frz_idxs.append([frag1_natom, 0])
    no_frz_idxs.append([frag1_natom, 1])
    no_frz_idxs.append([frag1_natom, 2])
    if frag2_natom == 2:
        no_frz_idxs.append([frag1_natom+1, 1])
        no_frz_idxs.append([frag1_natom+1, 2])
    elif frag2_natom > 2:
        no_frz_idxs.append([frag1_natom+1, 1])
        no_frz_idxs.append([frag1_natom+1, 2])
        no_frz_idxs.append([frag1_natom+2, 1])

    # Now grab the coordinates NOT in the opt coord idxs
    alt_froz_coords = []
    name_matrix = automol.zmatrix.name_matrix(inf_sep_zma)
    for row_idx, row in enumerate(name_matrix):
        for coord_idx, coord in enumerate(row):
            if [row_idx, coord_idx] not in no_frz_idxs:
                if coord is not None:
                    alt_froz_coords.append(coord)

    # Now build the constraint dictionary
    zma_vals = automol.zmatrix.values(inf_sep_zma)
    constraint_dct = dict(zip(
        alt_froz_coords, (zma_vals[name] for name in alt_froz_coords)))

    return constraint_dct


def read_potentials(scn_save_fs, cscn_save_fs,
                    sp_thy_level, dist_name, full_grid,
                    constraint_dct):
    """ Read values form the filesystem to get the values to
        correct ht MEP
    """

    # Read the energies from the full and constrained opts along MEP
    smp_pot = []
    const_pot = []
    sp_pot = []
    for grid_val in full_grid:

        # Set the locs for the full scan and constrained scan
        locs = [[dist_name], [grid_val]]
        const_locs = [[dist_name], [grid_val], constraint_dct]

        # Read the energies from the scan and constrained scan
        if scn_save_fs[-1].file.energy.exists(locs):
            smp_pot.append(scn_save_fs[-1].file.energy.read(locs))
        else:
            print('No scan energy')
        if cscn_save_fs[-1].file.energy.exists(const_locs):
            const_pot.append(cscn_save_fs[-1].file.energy.read(const_locs))
        else:
            print('No constrained scan energy')

        # Read the single point energy from the potential
        if sp_thy_level is not None:
            scn_save_path = scn_save_fs[-1].path(locs)
            sp_save_fs = autofile.fs.single_point(scn_save_path)
            sp_save_fs[-1].create(sp_thy_level[1:4])
            if sp_save_fs[-1].file.energy.read(sp_thy_level[1:4]):
                sp_pot.append(
                    sp_save_fs[-1].file.energy.read(sp_thy_level[1:4]))

    # Calculate each of the correction potentials
    relax_corr_pot = []
    sp_corr_pot = []
    full_corr_pot = []
    for i, _ in enumerate(smp_pot):
        relax_corr = (smp_pot[i] - const_pot[i]) * phycon.EH2KCAL
        relax_corr_pot.append(relax_corr)
        if sp_thy_level is not None:
            sp_corr = (sp_pot[i] - smp_pot[i]) * phycon.EH2KCAL
            sp_corr_pot.append(sp_corr)
        else:
            sp_corr = 0.0
        full_corr_pot.append(relax_corr + sp_corr)

    # Collate the potentials together in a list
    if sp_thy_level is not None:
        potentials = [relax_corr_pot, sp_corr_pot, full_corr_pot]
        potential_labels = ['relax', 'sp', 'full']
    else:
        potentials = [relax_corr_pot, full_corr_pot]
        potential_labels = ['relax', 'full']

    return potentials, potential_labels


# FUNCTION TO SET UP THE FRAGMENT GEOMETRIES FOR THE STRUCTURE.INP FILE
def fragment_geometries(ts_zma, rct_zmas, min_idx, max_idx):
    """ Generate the fragment geometries from the ts Z-matrix and the
        indices involved in the forming bond
    """

    # Get geometries of fragments from the ts_zma from the MEP
    mep_total_geo = automol.zmatrix.geometry(ts_zma)
    mep_fgeos = [mep_total_geo[:max_idx], mep_total_geo[max_idx:]]

    # Get geometries of isolated fragments at infinite sepearation
    iso_fgeos = [automol.zmatrix.geometry(zma) for zma in rct_zmas]

    # Reorder the iso_fgeos to line up with the mep_frag_geos
    (iso1_symbs, iso2_symbs) = (automol.geom.symbols(geo) for geo in iso_fgeos)
    (mep1_symbs, mep2_symbs) = (automol.geom.symbols(geo) for geo in mep_fgeos)
    if iso1_symbs != mep1_symbs or iso2_symbs != mep2_symbs:
        iso_fgeos[0], iso_fgeos[1] = iso_fgeos[1], iso_fgeos[0]

    # Get the geometries for the structure.inp file
    iso_fgeos_wdummy = []
    mol_data = zip(mep_fgeos, iso_fgeos, (max_idx, min_idx))
    for i, (mep_fgeo, iso_fgeo, idx) in enumerate(mol_data):

        if not automol.geom.is_atom(mep_fgeo):

            # Build MEPFragGeom+X coordinates using MEP geometry
            # x_idx: index for geom to place dummy X atom
            # a1_idx: index corresponding to "bonding" atom in geometry
            x_coord = mep_total_geo[idx][1]
            dummy_row = ('X', x_coord)
            if i == 0:
                mep_geo_wdummy = mep_fgeo + (dummy_row,)
                x_idx = len(mep_geo_wdummy) - 1
                a1_idx = 0
            else:
                mep_geo_wdummy = (dummy_row,) + mep_fgeo
                x_idx = 0
                a1_idx = 1

            # Set a2_idx to a1_idx + 1; should not be any restrictions
            a2_idx = a1_idx + 1

            # Set a3_idx.
            # Need to ensure idx does NOT correspond to atom where x = 0.0
            # The internal xyzp routine dies in this case
            for idx2 in range(a2_idx+1, len(iso_fgeo)):
                if not iso_fgeo[idx2][1][0] == 0.0:
                    a3_idx = idx2
                    break

            # Calculate coords to define X position in IsoFragGeom structure
            xyz1 = iso_fgeo[a1_idx][1]
            xyz2 = iso_fgeo[a2_idx][1]
            xdistance = automol.geom.distance(
                mep_geo_wdummy, x_idx, a1_idx)
            xangle = automol.geom.central_angle(
                mep_geo_wdummy, x_idx, a1_idx, a2_idx)
            if len(mep_fgeo) > 2:
                xyz3 = iso_fgeo[a3_idx][1]
                xdihedral = automol.geom.dihedral_angle(
                    mep_geo_wdummy, x_idx, a1_idx, a2_idx, a3_idx)
            else:
                xyz3 = 0.0
                xdihedral = 0.0

            # Calculate the X Position for the IsoFrag structure
            xyzp = automol.geom.find_xyzp_using_internals(
                xyz1, xyz2, xyz3, xdistance, xangle, xdihedral)

            # Generate the IsoFragGeom+X coordinates for the structure.inp file
            if i == 0:
                iso_geo_wdummy = iso_fgeo + (('X', xyzp),)
            else:
                iso_geo_wdummy = (('X', xyzp),) + iso_fgeo

            # Append to final geoms
            iso_fgeos_wdummy.append(iso_geo_wdummy)

        else:
            # If atom, set IsoFragGeom+X coords equal to mep_geo
            iso_fgeos_wdummy.append(mep_fgeo)

    return mep_total_geo, iso_fgeos, iso_fgeos_wdummy


# FUNCTIONS TO SET UP THE SYMMETRY FOR THE TST.INP FILE
def build_divsur_out_file(vrc_path, work_path):
    """ get the divsur.out string containing divsur-frame geoms
    """

    # Have to to path with divsur.inp to run script (maybe can fix)
    os.chdir(vrc_path)

    # Run the VaReCoF utility script to get the divsur.out file
    # Contains the fragment geometries in the divsur-defined coord sys
    varecof_io.writer.util.divsur_frame_geom_script()

    # Read fragment geoms from divsur.out with coordinates in the divsur frame
    with open(os.path.join(vrc_path, 'divsur.out'), 'r') as divsur_file:
        output_string = divsur_file.read()

    os.chdir(work_path)

    return output_string


def assess_face_symmetries(divsur_out_string):
    """ check the symmetry of the faces for each fragment
    """

    # Read fragment geoms from divsur.out with coordinates in the divsur frame
    fgeo1, fgeo2 = varecof_io.reader.divsur.frag_geoms_divsur_frame(
        divsur_out_string)
    fgeos = [automol.geom.from_string(fgeo1), automol.geom.from_string(fgeo2)]

    # Check facial symmetry if fragments are molecules
    symms = [False, False]
    for i, fgeo in enumerate(fgeos):
        if not automol.geom.is_atom(fgeo):
            # Reflect the dummy atom (pivot position) about the xy plane
            if i == 0:
                dummy_idx = len(fgeo) - 1
            else:
                dummy_idx = 0
            fgeo_reflect = automol.geom.reflect_coordinates(
                fgeo, [dummy_idx], ['x', 'y'])
            # Compute Coloumb spectrum for each geom to its reflected version
            symms[i] = automol.geom.almost_equal_coulomb_spectrum(
                fgeo, fgeo_reflect, rtol=5e-2)

    # Set the face and face_sym keywords based on the above tests
    [symm1, symm2] = symms
    if symm1 and symm2:
        faces = [0, 1]
        face_symm = 4
    elif symm1 and not symm2:
        faces = [0, 1]
        face_symm = 2
    elif not symm1 and symm2:
        faces = [0, 1]
        face_symm = 2
    elif not symm1 and not symm2:
        faces = [0]
        face_symm = 1

    return faces, face_symm


# FUNCTIONS TO SET UP THE DIVIDING SURFACE FRAMES
def build_pivot_frames(min_idx, max_idx,
                       total_geom, frag_geoms, frag_geoms_wdummy):
    """ Use geometries to get pivot info only set up for 1 or 2 pivot points
    """

    frames, npivots, = [], []
    geom_data = zip([min_idx, max_idx], frag_geoms, frag_geoms_wdummy)
    for i, (rxn_idx, geom, _) in enumerate(geom_data):

        # Single pivot point centered on atom
        print('geom test in build_pivot:', geom)
        if automol.geom.is_atom(geom):
            npivot = 1
            frame = [0, 0, 0, 0]
        # For linear species we place the pivot point on radical
        # with no displacment, so no need to find coordinates
        elif automol.geom.is_linear(geom):
            npivot = 2
            frame = [0, 0, 0, 0]
        else:
            # else we build an xy frame to easily place pivot point
            npivot = 2

            # Find the idx in each fragment bonded to the atom at the pivot pt
            for j, coords in enumerate(geom):
                if coords == total_geom[rxn_idx]:
                    coord_idx = j
                    break

            # For each fragment, get indices for a
            # chain (up to three atoms, that terminates at the dummy atom)
            gra = automol.geom.graph(geom)
            gra_neighbor_dct = automol.graph.atom_neighbor_keys(gra)
            bond_neighbors = gra_neighbor_dct[coord_idx]

            # Find idx in each fragment geom that corresponds to the bond index
            for j, idx in enumerate(bond_neighbors):
                if geom[idx][0] != 'H':
                    bond_neighbor_idx = idx
                    break
                if geom[idx][0] == 'H' and j == (len(bond_neighbors) - 1):
                    bond_neighbor_idx = idx

            # Set up the frame indices for the divsur file
            if i == 0:
                pivot_idx = len(geom)
                frame = [coord_idx, bond_neighbor_idx, pivot_idx, coord_idx]
            else:
                pivot_idx = 0
                coord_idx += 1
                bond_neighbor_idx += 1
                frame = [coord_idx, bond_neighbor_idx, pivot_idx, coord_idx]
            frame = [val+1 for val in frame]

        # Append to lists
        frames.append(frame)
        npivots.append(npivot)

    return frames, npivots


def calc_pivot_angles(frag_geoms, frag_geoms_wdummy, frames):
    """ get the angle for the three atoms definining the frame
    """
    angles = []
    for geom, geom_wdummy, frame in zip(frag_geoms, frag_geoms_wdummy, frames):
        if automol.geom.is_atom(geom) or automol.geom.is_linear(geom):
            angle = None
        else:
            frame = [val-1 for val in frame]
            angle = automol.geom.central_angle(
                geom_wdummy, frame[2], frame[0], frame[1])
        angles.append(angle)

    return angles


def calc_pivot_xyzs(min_idx, max_idx, total_geom, frag_geoms):
    """ figure out where pivot point will be centered
        only linear speces need to have non-zero xyz, as they
        will not have a frame set up for them like atoms and
        polyatomics
    """
    xyzs = []
    for rxn_idx, geom in zip([min_idx, max_idx], frag_geoms):
        if automol.geom.is_linear(geom):
            xyz = total_geom[rxn_idx][1]
        else:
            xyz = [0.0, 0.0, 0.0]

        xyzs.append(xyz)

    return xyzs
