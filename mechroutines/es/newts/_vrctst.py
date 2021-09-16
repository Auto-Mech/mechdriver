""" Generate the information necessary to product the vrctst input files
"""

import autofile
import automol
from phydat import phycon
import varecof_io
import autorun
from mechroutines.es.runner import scan
from mechroutines.es.runner import qchem_params
from mechroutines.es._routines import _wfn as wfn
from mechroutines.es._routines import sp
from mechlib.amech_io import printer as ioprinter
from mechlib import filesys


# CENTRAL FUNCTION TO WRITE THE VARECOF INPUT FILES AND RUN THE PROGRAM
def calc_vrctst_flux(ini_zma, ts_info, hs_info,
                     ts_formula, high_mul, active_space,
                     rct_info, rct_ichs, rct_zmas, rcts_cnf_fs,
                     grid1, grid2, coord_name,
                     mod_var_scn_thy_info,
                     mod_var_sp1_thy_info, var_sp2_thy_info,
                     hs_var_sp1_thy_info, hs_var_sp2_thy_info,
                     vscnlvl_thy_save_fs,
                     vscnlvl_ts_run_fs, vscnlvl_ts_save_fs,
                     vscnlvl_scn_run_fs, vscnlvl_scn_save_fs,
                     vscnlvl_cscn_run_fs, vscnlvl_cscn_save_fs,
                     overwrite, update_guess):
    """ Set up n VRC-TST calculations to get the flux file
    """

    # Set vrc tst dct (just using defaults for now)
    vrc_dct = autorun.varecof.VRC_DCT

    # Set up machine dct (how to do this?)
    machine_dct = {}

    # Set up the casscf options
    ref_zma = automol.zmat.set_values_by_name(ini_zma, {coord_name: grid1[0]})
    cas_kwargs = wfn.build_wfn(ref_zma, ts_info, ts_formula, high_mul,
                               rct_ichs, rct_info,
                               active_space, mod_var_scn_thy_info)

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

    # Write the electronic structure template file

    # Write remaining VaReCoF input files
    inp_strs = autorun.varecof.write_varecof_input(
        vrc_path,
        zma_for_inp, rct_zmas,
        npot, min_idx, max_idx,
        machine_dct, vrc_dct)

    inp_str += (('', els_str),)

    # Run VaReCoF to generate flux file
    varecof_script_str = autorun.SCRIPT_DCT['varecof']
    mcflux_script_str = autorun.SCRIPT_DCT['mcflux']
    flux_str = autorun.varecof.flux_file(
       varecof_script_str, mcflux_script_str,
       vrc_path, input_strs)

    # Save the flux file
    if flux_str is not None:
        filesys.save.flux(flux_str, inp_strs, ts_save_fs,
                          ts_locs=(0,), vrc_locs=(0,))


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
    constraint_dct = varecof_io.writer.intramolecular_constraint_dct(
        ref_zma, rct_zmas)

    # Run the constrained and full opt potential scans
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

    inf_sep_ene = rpath.inf_sep_ene(
        ts_dct, thy_inf_dct, savefs_dct, runfs_dct, es_keyword_dct)

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
    autorun.varecof.compile_potentials(
        vrc_path, full_grid, potentials,
        bnd_frm_idxs, vrc_dct['fortran_compiler'],
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

        for grid in (grid1, grid2):
            scan.execute_scan(
                zma=inf_sep_zma,
                spc_info=ts_info,
                mod_thy_info=mod_var_scn_thy_info,
                coord_names=[coord_name],
                coord_grids=[grid],
                scn_run_fs=scn_run_fs,
                scn_save_fs=scn_save_fs,
                scn_typ='relaxed',
                script_str=opt_script_str,
                overwrite=overwrite,
                update_guess=update_guess,
                reverse_sweep=False,
                saddle=False,
                constraint_dct=constraint_dct,
                retryfail=True,
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
