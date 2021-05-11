""" Generate the information necessary to product the vrctst input files
"""

import os
import autofile
import automol
import varecof_io
from phydat import phycon
from autorun import run_script
from mechroutines.es.runner import scan
from mechroutines.es._routines import _wfn as wfn
from mechroutines.es._routines import sp
from mechlib.amech_io import printer as ioprinter


# CENTRAL FUNCTION TO WRITE THE VARECOF INPUT FILES AND RUN THE PROGRAM
def calc_vrctst_flux(ini_zma, ts_info, hs_info,
                     ts_formula, high_mul, active_space,
                     rct_info, rct_ichs, rct_zmas, rcts_cnf_fs,
                     grid1, grid2, coord_name,
                     mod_var_scn_thy_info,
                     mod_var_sp1_thy_info, var_sp2_thy_info,
                     hs_var_sp1_thy_info, hs_var_sp2_thy_info,
                     vscnlvl_thy_save_fs,
                     vscnlvl_ts_run_fs,
                     vscnlvl_scn_run_fs, vscnlvl_scn_save_fs,
                     vscnlvl_cscn_run_fs, vscnlvl_cscn_save_fs,
                     overwrite, update_guess):
    """ Set up n VRC-TST calculations to get the flux file
    """

    # Set vrc tst dct
    vrc_dct = _vrc_dct()

    # Set up the casscf options
    ref_zma = automol.zmat.set_values_by_name(ini_zma, {coord_name: grid1[0]})
    cas_kwargs = wfn.build_wfn(ref_zma, ts_info, ts_formula, high_mul,
                               rct_ichs, rct_info,
                               active_space, mod_var_scn_thy_info)
    _, script_str, _, _ = qchem_params(
        mod_var_sp1_thy_info[0], mod_var_sp1_thy_info[1])

    # Get indices for potentials and input
    bnd_frm_idxs = automol.zmat.coord_idxs(ini_zma, coord_name)
    min_idx, max_idx = min(bnd_frm_idxs), max(bnd_frm_idxs)
    bnd_frm_idxs = (bnd_frm_idxs[0]+1, bnd_frm_idxs[1]+1)

    # Build the VRC-TST run directory, including the needed scr dir
    vrc_path = _build_vrctst_fs(vscnlvl_ts_run_fs)

    # Calculate the correction potential along the MEP
    inf_sep_ene, npot, zma_for_inp = _build_correction_potential(
        ts_info, hs_info, ref_zma,
        coord_name, bnd_frm_idxs,
        grid1, grid2,
        rct_info, rcts_cnf_fs, rct_zmas,
        mod_var_scn_thy_info, mod_var_sp1_thy_info,
        hs_var_sp1_thy_info,
        var_sp2_thy_info,
        hs_var_sp2_thy_info,
        vscnlvl_scn_run_fs, vscnlvl_scn_save_fs,
        vscnlvl_cscn_run_fs, vscnlvl_cscn_save_fs,
        vscnlvl_thy_save_fs,
        overwrite, update_guess,
        vrc_dct, vrc_path,
        cas_kwargs)

    # Write remaining VaReCoF input files
    _write_varecof_input(zma_for_inp, ts_info, ts_formula, high_mul,
                         rct_ichs, rct_info, rct_zmas,
                         active_space, mod_var_sp1_thy_info,
                         npot, inf_sep_ene,
                         min_idx, max_idx,
                         vrc_dct, vrc_path, script_str)

    # Run VaReCoF to generate flux file
    _save_flux()


# FUNCTIONS TO SET UP THE libcorrpot.so FILE USED BY VARECOF
def _build_correction_potential(ts_info, hs_info, ref_zma,
                                coord_name, bnd_frm_idxs,
                                grid1, grid2,
                                rct_info, rcts_cnf_fs, rct_zmas,
                                mod_var_scn_thy_info, mod_var_sp1_thy_info,
                                hs_var_sp1_thy_info,
                                var_sp2_thy_info,
                                hs_var_sp2_thy_info,
                                vscnlvl_scn_run_fs, vscnlvl_scn_save_fs,
                                vscnlvl_cscn_run_fs, vscnlvl_cscn_save_fs,
                                vscnlvl_thy_save_fs,
                                overwrite, update_guess,
                                vrc_dct, vrc_path,
                                cas_kwargs):
    """  use the MEP potentials to compile the correction potential .so file
    """

    # Build the constraint dictionary
    constraint_dct = _set_alt_constraints(ref_zma, rct_zmas)

    # Run the potentials
    _run_potentials(ref_zma, ts_info,
                    mod_var_scn_thy_info,
                    coord_name, grid1, grid2,
                    vscnlvl_scn_run_fs, vscnlvl_scn_save_fs,
                    vscnlvl_cscn_run_fs, vscnlvl_cscn_save_fs,
                    vscnlvl_thy_save_fs,
                    overwrite, update_guess,
                    constraint_dct,
                    cas_kwargs,
                    sp_thy_info=mod_var_sp1_thy_info)

    # Calculate and store the infinite separation energy
    # max_grid = max([max(grid1), max(grid2)])
    # locs = [[dist_name], [max_grid]]
    inf_locs = [[coord_name], [grid1[0]]]
    ts_zma = vscnlvl_scn_save_fs[-1].file.zmatrix.read(inf_locs)
    geo = vscnlvl_scn_save_fs[-1].file.geometry.read(inf_locs)

    geo_run_path = vscnlvl_scn_run_fs[-1].path(inf_locs)
    geo_save_path = vscnlvl_scn_save_fs[-1].path(inf_locs)

    inf_sep_ene = scan.radrad_inf_sep_ene(
        hs_info, ts_zma,
        rct_info, rcts_cnf_fs,
        var_sp2_thy_info,
        hs_var_sp1_thy_info, hs_var_sp2_thy_info,
        geo, geo_run_path, geo_save_path,
        vscnlvl_scn_save_fs, inf_locs,
        overwrite, **cas_kwargs)

    # Combine and sort the grids for organization
    full_grid = list(grid1) + list(grid2)
    full_grid.sort()

    # Get grid val for zma used to make structure.inp and divsur.inp
    ioprinter.debug_message('grid1', grid1)
    grid_val_for_zma = grid1[-1]

    # Read the values for the correction potential from filesystem
    potentials, pot_labels, zma_for_inp = _read_potentials(
        vscnlvl_scn_save_fs, vscnlvl_cscn_save_fs,
        mod_var_scn_thy_info, mod_var_sp1_thy_info,
        coord_name, full_grid, inf_sep_ene,
        constraint_dct, grid_val_for_zma)

    # Build correction potential .so file used by VaReCoF
    _compile_potentials(
        full_grid, potentials,
        bnd_frm_idxs, vrc_dct['fortran_compiler'], vrc_path,
        dist_restrict_idxs=(),
        pot_labels=pot_labels,
        pot_file_names=[vrc_dct['spc_name']],
        spc_name=vrc_dct['spc_name'])

    # Set zma if needed
    if zma_for_inp is None:
        zma_for_inp = ref_zma

    return inf_sep_ene, len(potentials), zma_for_inp


def _run_potentials(inf_sep_zma, ts_info,
                    mod_var_scn_thy_info,
                    coord_name, grid1, grid2,
                    vscnlvl_scn_run_fs, vscnlvl_scn_save_fs,
                    vscnlvl_cscn_run_fs, vscnlvl_cscn_save_fs,
                    vscnlvl_thy_save_fs,
                    overwrite, update_guess,
                    constraint_dct,
                    cas_kwargs,
                    sp_thy_info=None):
    """ Run and save the scan along both grids while
          (1) constraining only reaction coordinate, then
          (2) constraining all intermolecular coordinates
    """

    for constraints in (None, constraint_dct):
        if constraints is None:
            scn_run_fs, scn_save_fs = vscnlvl_scn_run_fs, vscnlvl_scn_save_fs
            ioprinter.info_message('Running full scans..', newline=1)
        else:
            scn_run_fs, scn_save_fs = vscnlvl_cscn_run_fs, vscnlvl_cscn_save_fs
            ioprinter.info_message('Running constrained scans..', newline=1)

        scan.multiref_rscan(
            ts_zma=inf_sep_zma,
            ts_info=ts_info,
            grid1=grid1,
            grid2=grid2,
            coord_name=coord_name,
            mod_var_scn_thy_info=mod_var_scn_thy_info,
            vscnlvl_thy_save_fs=vscnlvl_thy_save_fs,
            scn_run_fs=scn_run_fs,
            scn_save_fs=scn_save_fs,
            overwrite=overwrite,
            update_guess=update_guess,
            constraint_dct=constraints,
            **cas_kwargs
        )

    # Run the single points on top of the initial scan
    if sp_thy_info is not None:
        _scan_sp(ts_info, coord_name,
                 vscnlvl_scn_run_fs, vscnlvl_scn_save_fs,
                 sp_thy_info, overwrite,
                 cas_kwargs)


def _scan_sp(ts_info, coord_name,
             vscnlvl_scn_run_fs, vscnlvl_scn_save_fs,
             mod_var_sp1_thy_info, overwrite,
             cas_kwargs):
    """ get sps for the scan; cas options and gen lines should be same
    """

    # Set up script and kwargs for the irc run
    script_str, _, sp_kwargs, _ = qchem_params(
        *mod_var_sp1_thy_info[0:2])
    sp_kwargs.update(cas_kwargs)

    # Compute the single-point energies along the scan
    for locs in vscnlvl_scn_save_fs[-1].existing([[coord_name]]):

        # Set up single point filesys
        vscnlvl_scn_run_fs[-1].create(locs)
        geo_run_path = vscnlvl_scn_run_fs[-1].path(locs)
        geo_save_path = vscnlvl_scn_save_fs[-1].path(locs)
        geo = vscnlvl_scn_save_fs[-1].file.geometry.read(locs)
        zma = vscnlvl_scn_save_fs[-1].file.zmatrix.read(locs)

        # Run the energy
        sp.run_energy(zma, geo, ts_info, mod_var_sp1_thy_info,
                      vscnlvl_scn_save_fs,
                      geo_run_path, geo_save_path, locs,
                      script_str, overwrite,
                      highspin=False, **sp_kwargs)


def _read_potentials(scn_save_fs, cscn_save_fs,
                     mod_var_scn_thy_info, mod_var_sp1_thy_info,
                     dist_name, full_grid,
                     constraint_dct, grid_val_for_zma):
    """ Read values form the filesystem to get the values to
        correct ht MEP
    """

    # Read the energies from the full and constrained opts along MEP
    smp_pot = []
    const_pot = []
    sp_pot = []

    # Put scans info together
    scans = (
        (scn_save_fs, mod_var_scn_thy_info[1:4]),
        (cscn_save_fs, mod_var_scn_thy_info[1:4])
    )
    if mod_var_sp1_thy_info is not None:
        scans += ((scn_save_fs, mod_var_sp1_thy_info[1:4]),)

    smp_pot, const_pot, sp_pot = [], [], []
    for idx, (scn_fs, thy_info) in enumerate(scans):

        for grid_val in full_grid:

            # Set the locs for the full scan and constrained scan
            if idx in (0, 2):
                locs = [[dist_name], [grid_val]]
            else:
                locs = [constraint_dct, [dist_name], [grid_val]]

            # Read the energies from the scan and constrained scan
            scn_path = scn_fs[-1].path(locs)
            sp_fs = autofile.fs.single_point(scn_path)
            if sp_fs[-1].file.energy.exists(thy_info):
                sp_ene = sp_fs[-1].file.energy.read(thy_info)
                ioprinter.debug_message('sp_path', sp_fs[-1].path(thy_info))
                ioprinter.debug_message('sp_ene', sp_ene)
            else:
                ioprinter.info_message('No scan energy')

            # Store the energy in a lst
            if idx == 0:
                smp_pot.append(sp_ene)
            elif idx == 1:
                const_pot.append(sp_ene)
            elif idx == 2:
                sp_pot.append(sp_ene)

        ioprinter.info_message()

    # Calculate each of the correction potentials
    relax_corr_pot = []
    sp_corr_pot = []
    full_corr_pot = []
    for i, _ in enumerate(smp_pot):
        relax_corr = (smp_pot[i] - const_pot[i]) * phycon.EH2KCAL
        relax_corr_pot.append(relax_corr)
        if sp_pot:
            sp_corr = (sp_pot[i] - smp_pot[i]) * phycon.EH2KCAL
            sp_corr_pot.append(sp_corr)
        else:
            sp_corr = 0.0
        full_corr_pot.append(relax_corr + sp_corr)

    # Collate the potentials together in a list
    if sp_pot:
        potentials = [relax_corr_pot, sp_corr_pot, full_corr_pot]
        potential_labels = ['relax', 'sp', 'full']
    else:
        potentials = [relax_corr_pot, full_corr_pot]
        potential_labels = ['relax', 'full']

    # Get zma used to make structure.inp and divsur.inp
    inp_zma_locs = [[dist_name], [grid_val_for_zma]]
    if scn_save_fs[-1].file.zmatrix.exists(inp_zma_locs):
        zma_for_inp = scn_save_fs[-1].file.zmatrix.read(inp_zma_locs)
    else:
        zma_for_inp = None

    return potentials, potential_labels, zma_for_inp


# FUNCTIONS TO WRITE THE STRINGS FOR ALL OF THE VARECOF INPUT FILE
def _write_varecof_input(ref_zma, ts_info, ts_formula, high_mul,
                         rct_ichs, rct_info, rct_zmas,
                         active_space, mod_var_sp1_thy_info,
                         npot, inf_sep_ene,
                         min_idx, max_idx,
                         vrc_dct, vrc_path, script_str):
    """ prepare all the input files for a vrc-tst calculation
    """

    r1dists_lr = vrc_dct['r1dists_lr']
    r1dists_sr = vrc_dct['r1dists_sr']
    r2dists_sr = vrc_dct['r2dists_sr']
    d1dists = vrc_dct['d1dists']
    d2dists = vrc_dct['d2dists']
    conditions = vrc_dct['conditions']
    nsamp_max = vrc_dct['nsamp_max']
    nsamp_min = vrc_dct['nsamp_min']
    flux_err = vrc_dct['flux_err']
    pes_size = vrc_dct['pes_size']
    base_name = vrc_dct['base_name']
    # exe_path = vrc_dct['exe_path']

    # Build geometries needed for the varecof run
    total_geom, frag_geoms, frag_geoms_wdummy = fragment_geometries(
        ref_zma, rct_zmas, min_idx, max_idx)

    # Set information for the pivot points needed in divsur.inp
    frames, npivots = build_pivot_frames(
        min_idx, max_idx, total_geom, frag_geoms, frag_geoms_wdummy)
    pivot_angles = calc_pivot_angles(frag_geoms, frag_geoms_wdummy, frames)
    pivot_xyzs = calc_pivot_xyzs(min_idx, max_idx, total_geom, frag_geoms)

    # Write the long- and short-range divsur input files
    lrdivsur_inp_str = varecof_io.writer.input_file.divsur(
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
    if all(npiv > 1 for npiv in npivots):
        r2dists = r2dists_sr
    else:
        r2dists = []
        ioprinter.warning_message('no r2dist')

    srdivsur_inp_str = varecof_io.writer.input_file.divsur(
        r1dists_sr, npivots[0], npivots[1], pivot_xyzs[0], pivot_xyzs[1],
        frame1=frames[0], frame2=frames[1],
        d1dists=d1dists, d2dists=d2dists,
        t1angs=t1angs, t2angs=t2angs,
        r2dists=r2dists,
        **conditions)

    # Build the structure input file string
    struct_inp_str = varecof_io.writer.input_file.structure(
        frag_geoms_wdummy[0], frag_geoms_wdummy[1])

    # Write the structure and divsur files to get the divsur out file
    inp = ((struct_inp_str, 'structure.inp'), (srdivsur_inp_str, 'divsur.inp'))
    _write_varecof_inp(inp, vrc_path)

    # Obtain the divsur.out file with divsur-frame fragment geoms
    divsur_out_str = build_divsur_out_file(vrc_path, os.getcwd())

    # Write the tst.inp file
    faces, faces_symm = assess_face_symmetries(divsur_out_str)
    tst_inp_str = varecof_io.writer.input_file.tst(
        nsamp_max, nsamp_min, flux_err, pes_size,
        faces=faces, faces_symm=faces_symm)

    # Write the molpro executable and potential energy surface input string
    els_inp_str = varecof_io.writer.input_file.elec_struct(
        vrc_path, base_name, npot,
        dummy_name='dummy_corr_', lib_name='libcorrpot.so',
        exe_name='molpro.sh',
        geom_ptt='GEOMETRY_HERE', ene_ptt='molpro_energy')

    # Write the electronic structure template file
    tml_inp_str = _build_molpro_template_str(
        ref_zma, ts_info, ts_formula, high_mul,
        rct_ichs, rct_info,
        active_space, mod_var_sp1_thy_info,
        inf_sep_ene)

    # Write the mc_flux.inp input string
    mc_flux_inp_str = varecof_io.writer.input_file.mc_flux()

    # Write the convert.inp input string
    conv_inp_str = varecof_io.writer.input_file.convert()

    # Write machines file to set compute nodes
    machine_file_str = build_machinefile_str()

    # Collate the input strings and write the remaining files
    input_strs = (
        lrdivsur_inp_str, tst_inp_str,
        els_inp_str, tml_inp_str,
        mc_flux_inp_str, conv_inp_str,
        machine_file_str, script_str)
    input_names = (
        'lr_divsur.inp', 'tst.inp',
        'molpro.inp', 'mol.tml',
        'mc_flux.inp', 'convert.inp',
        'machines', 'molpro.sh')
    inp = tuple(zip(input_strs, input_names))
    _write_varecof_inp(inp, vrc_path)


def _build_molpro_template_str(ref_zma, ts_info, ts_formula, high_mul,
                               rct_ichs, rct_info,
                               active_space, mod_var_sp1_thy_info,
                               inf_sep_ene):
    """ Write the electronic structure template file
    """

    cas_kwargs = wfn.build_wfn(ref_zma, ts_info, ts_formula, high_mul,
                               rct_ichs, rct_info,
                               active_space, mod_var_sp1_thy_info)

    tml_inp_str = wfn.wfn_string(
        ts_info, mod_var_sp1_thy_info, inf_sep_ene, cas_kwargs)

    return tml_inp_str


# FUNCTIONS TO SET UP THE DIVIDING SURFACE FRAMES
def _save_flux(vrc_ret, ts_run_fs, ts_save_fs, ts_locs=(0,), vrc_locs=(0,)):
    """ Save the VaReCoF flux file
    """
    mcflux_file = os.path.join(vrc_path, 'mc_flux.out')
    if os.path.exists(mcflux_file):
        with open(mcflux_file, 'r') as fluxfile:
            flux_str = mcflux_file.read()
        if flux_str != '':
            ts_save_fs[0].create()
            ts_save_path = ts_save_fs[0].path()

            # Write the VRC-Flux and info files
            ts_save_fs[0].file.vrc_flux.write(flux_str)
            ts_save_fs[0].file.vrc_flux_info.write(flux_str)

    # Unpack the ret

    # Save the files
    ts_save_path = ts_save_fs[-1].path(ts_locs)

    vrc_fs = autofile.fs.vrctst(ts_save_path)
    vrc_fs[-1].create(vrc_locs)
    vrc_fs[-1].file.vrctst_tst.write(ref_tst_str, vrc_locs)
    vrc_fs[-1].file.vrctst_divsur.write(ref_divsur_str, vrc_locs)
    vrc_fs[-1].file.vrctst_molpro.write(ref_molpro_str, vrc_locs)
    vrc_fs[-1].file.vrctst_tml.write(ref_tml_str, vrc_locs)
    vrc_fs[-1].file.vrctst_struct.write(ref_struct_str, vrc_locs)
    vrc_fs[-1].file.vrctst_pot.write(ref_pot_str, vrc_locs)
    vrc_fs[-1].file.vrctst_flux.write(ref_flux_str, vrc_locs)
