""" Generate the information necessary to product the vrctst input files
"""

import os
import autofile
import automol
import varecof_io
from routines.es._routines import _scan as scan
from routines.es._routines import _wfn as wfn
from routines.es._routines import sp
from routines.es import runner as es_runner
from lib import filesys
from lib.submission import run_script
from lib.submission import DEFAULT_SCRIPT_DCT
from lib.phydat import phycon
from lib import get_host_node


# CENTRAL FUNCTION TO WRITE THE VARECOF INPUT FILES AND RUN THE PROGRAM
def calc_vrctst_flux(ini_zma, ts_info, ts_formula, high_mul, active_space,
                     rct_info, rct_ichs, rct_zmas, rcts_cnf_fs,
                     grid1, grid2, coord_name,
                     mod_var_scn_thy_info,
                     mod_var_sp1_thy_info,
                     mod_var_sp2_thy_info,
                     hs_var_scn_thy_info,
                     hs_var_sp1_thy_info,
                     hs_var_sp2_thy_info,
                     mod_ini_thy_info, mod_thy_info,
                     vscnlvl_thy_save_fs,
                     vscnlvl_ts_save_fs,
                     vscnlvl_ts_run_fs,
                     vscnlvl_scn_run_fs, vscnlvl_scn_save_fs,
                     vscnlvl_cscn_run_fs, vscnlvl_cscn_save_fs,
                     run_prefix, save_prefix,
                     overwrite, update_guess):
    """ Set up n VRC-TST calculations to get the flux file
    """

    # Set vrc tst dct
    vrc_dct = _vrc_dct()

    # Set up the casscf options
    ref_zma = automol.zmatrix.set_values(ini_zma, {coord_name: grid1[0]})
    cas_kwargs = wfn.build_wfn(ref_zma, ts_info, ts_formula, high_mul,
                               rct_ichs, rct_info,
                               active_space, mod_var_scn_thy_info)
    _, script_str, _, _ = es_runner.qchem_params(
        mod_var_sp1_thy_info[0], mod_var_sp1_thy_info[1])

    # Get indices for potentials and input
    bnd_frm_idxs = automol.zmatrix.coord_idxs(ini_zma, coord_name)
    min_idx, max_idx = min(bnd_frm_idxs), max(bnd_frm_idxs)
    bnd_frm_idxs = (bnd_frm_idxs[0]+1, bnd_frm_idxs[1]+1)

    # Build the VRC-TST run directory, including the needed scr dir
    vrc_path = _build_vrctst_fs(vscnlvl_ts_run_fs)

    # Calculate the correction potential along the MEP
    inf_sep_ene, npot, zma_for_inp = _build_correction_potential(
        ts_info, high_mul, ref_zma,
        coord_name, bnd_frm_idxs,
        grid1, grid2,
        rct_info, rcts_cnf_fs, rct_zmas,
        mod_var_scn_thy_info, mod_var_sp1_thy_info,
        hs_var_scn_thy_info, hs_var_sp1_thy_info,
        mod_ini_thy_info, mod_thy_info,
        mod_var_sp2_thy_info,
        hs_var_sp2_thy_info,
        vscnlvl_scn_run_fs, vscnlvl_scn_save_fs,
        vscnlvl_cscn_run_fs, vscnlvl_cscn_save_fs,
        vscnlvl_thy_save_fs,
        overwrite, update_guess,
        vrc_dct, vrc_path,
        **cas_kwargs)

    # Write remaining VaReCoF input files
    _write_varecof_input(zma_for_inp, ts_info, ts_formula, high_mul,
                         rct_ichs, rct_info, rct_zmas,
                         active_space, mod_var_sp1_thy_info,
                         npot, inf_sep_ene,
                         min_idx, max_idx,
                         vrc_dct, vrc_path, script_str)

    # Run VaReCoF to generate flux file
    _run_varecof(vrc_path)

    # Check for success and save the flux file if so
    # if _varecof_success(vrc_path):
    #    _save_flux()


# FUNCTIONS TO SET UP THE libcorrpot.so FILE USED BY VARECOF
def _build_correction_potential(ts_info, high_mul, ref_zma,
                                coord_name, bnd_frm_idxs,
                                grid1, grid2,
                                rct_info, rcts_cnf_fs, rct_zmas,
                                mod_var_scn_thy_info, mod_var_sp1_thy_info,
                                hs_var_scn_thy_info, hs_var_sp1_thy_info,
                                mod_ini_thy_info, mod_thy_info,
                                mod_var_sp2_thy_info,
                                hs_var_sp2_thy_info,
                                vscnlvl_scn_run_fs, vscnlvl_scn_save_fs,
                                vscnlvl_cscn_run_fs, vscnlvl_cscn_save_fs,
                                vscnlvl_thy_save_fs,
                                overwrite, update_guess,
                                vrc_dct, vrc_path,
                                **cas_kwargs):
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
                    sp_thy_info=mod_var_sp1_thy_info,
                    **cas_kwargs)

    # Calculate and store the infinite separation energy
    # max_grid = max([max(grid1), max(grid2)])
    # locs = [[dist_name], [max_grid]]
    inf_locs = [[coord_name], [grid1[0]]]
    ts_zma = vscnlvl_scn_save_fs[-1].file.zmatrix.read(inf_locs)
    geo = vscnlvl_scn_save_fs[-1].file.geometry.read(inf_locs)

    geo_run_path = vscnlvl_scn_run_fs[-1].path(inf_locs)
    geo_save_path = vscnlvl_scn_save_fs[-1].path(inf_locs)

    inf_sep_ene = scan.radrad_inf_sep_ene(
        ts_info, high_mul, ts_zma,
        rct_info, rcts_cnf_fs,
        mod_var_scn_thy_info, mod_var_sp1_thy_info,
        hs_var_scn_thy_info, hs_var_sp1_thy_info,
        mod_ini_thy_info, mod_thy_info,
        mod_var_sp2_thy_info,
        hs_var_sp2_thy_info,
        geo, geo_run_path, geo_save_path,
        vscnlvl_scn_save_fs, inf_locs,
        overwrite, **cas_kwargs)

    # Combine and sort the grids for organization
    full_grid = list(grid1) + list(grid2)
    full_grid.sort()
    
    # Get grid val for zma used to make structure.inp and divsur.inp
    grid_val_for_zma = grid1[-1]
    
    # Read the values for the correction potential from filesystem
    potentials, pot_labels, zma_for_inp = _read_potentials(
        vscnlvl_scn_save_fs, vscnlvl_cscn_save_fs,
        mod_var_sp1_thy_info, coord_name, full_grid, inf_sep_ene,
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


def _set_alt_constraints(inf_sep_zma, rct_zmas):
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


def _run_potentials(inf_sep_zma, ts_info,
                    mod_var_scn_thy_info,
                    coord_name, grid1, grid2,
                    vscnlvl_scn_run_fs, vscnlvl_scn_save_fs,
                    vscnlvl_cscn_run_fs, vscnlvl_cscn_save_fs,
                    vscnlvl_thy_save_fs,
                    overwrite, update_guess,
                    constraint_dct,
                    sp_thy_info=None,
                    **cas_kwargs):
    """ Run and save the scan along both grids while
          (1) constraining only reaction coordinate, then
          (2) constraining all intermolecular coordinates
    """

    for constraints in (None, constraint_dct):
        if constraints is None:
            scn_run_fs, scn_save_fs = vscnlvl_scn_run_fs, vscnlvl_scn_save_fs
            print('\nRunning full scans..')
        else:
            scn_run_fs, scn_save_fs = vscnlvl_cscn_run_fs, vscnlvl_cscn_save_fs
            print('\nRunning constrained scans..')

        print('scn thy test', mod_var_scn_thy_info)
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
    print('spthy', sp_thy_info)
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
    script_str, _, _, _ = es_runner.qchem_params(
        *mod_var_sp1_thy_info[0:2])

    # Compute the single-point energies along the scan
    for locs in vscnlvl_scn_save_fs[-1].existing([[coord_name]]):

        print('splocs', locs)

        # Set up single point filesys
        vscnlvl_scn_run_fs[-1].create(locs)
        geo_run_path = vscnlvl_scn_run_fs[-1].path(locs)
        geo_save_path = vscnlvl_scn_save_fs[-1].path(locs)
        geo = vscnlvl_scn_save_fs[-1].file.geometry.read(locs)
        zma = vscnlvl_scn_save_fs[-1].file.zmatrix.read(locs)
        # print('scn_run_path')
        # print(geo_run_path)

        # Run the energy
        sp.run_energy(zma, geo, ts_info, mod_var_sp1_thy_info,
                      vscnlvl_scn_save_fs,
                      geo_run_path, geo_save_path, locs,
                      script_str, overwrite,
                      highspin=False, **cas_kwargs)


def _read_potentials(scn_save_fs, cscn_save_fs,
                     sp_thy_info, dist_name, full_grid, inf_sep_ene,
                     constraint_dct, grid_val_for_zma):
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
        const_locs = [constraint_dct, [dist_name], [grid_val]]

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
        if sp_thy_info is not None:
            scn_save_path = scn_save_fs[-1].path(locs)
            sp_save_fs = autofile.fs.single_point(scn_save_path)
            sp_save_fs[-1].create(sp_thy_info[1:4])
            # print('sp', sp_save_fs[-1].path(sp_thy_info[1:4]))
            if sp_save_fs[-1].file.energy.exists(sp_thy_info[1:4]):
                sp_pot.append(
                    sp_save_fs[-1].file.energy.read(sp_thy_info[1:4]))

    # Calculate each of the correction potentials
    relax_corr_pot = []
    sp_corr_pot = []
    full_corr_pot = []
    for i, _ in enumerate(smp_pot):
        relax_corr = (smp_pot[i] - const_pot[i]) * phycon.EH2KCAL
        relax_corr_pot.append(relax_corr)
        if sp_thy_info is not None:
            sp_corr = (sp_pot[i] - smp_pot[i]) * phycon.EH2KCAL
            sp_corr_pot.append(sp_corr)
        else:
            sp_corr = 0.0
        full_corr_pot.append(relax_corr + sp_corr)

    # Collate the potentials together in a list
    if sp_thy_info is not None:
        potentials = [relax_corr_pot, sp_corr_pot, full_corr_pot]
        potential_labels = ['relax', 'sp', 'full']
    else:
        potentials = [relax_corr_pot, full_corr_pot]
        potential_labels = ['relax', 'full']

    # Get zma used to make structure.inp and divsur.inp
    inp_zma_locs = [[dist_name], [grid_val]]
    if scn_save_fs[-1].file.zmatrix.exists(inp_zma_locs):
        zma_for_inp = scn_save_fs[-1].file.zmatrix.read(inp_zma_locs)
    else:
        zma_for_inp = None

    return potentials, potential_labels, zma_for_inp


def _compile_potentials(mep_distances, potentials,
                        bnd_frm_idxs, fortran_compiler, vrc_path,
                        dist_restrict_idxs=(),
                        pot_labels=(),
                        pot_file_names=(),
                        spc_name='mol'):
    """  use the MEP potentials to compile the correction potential .so file
    """

    # Change the coordinates of the MEP distances
    # mep_distances = [dist * phycon.BOHR2ANG for dist in mep_distances]

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
    srdivsur_inp_str = varecof_io.writer.input_file.divsur(
        r1dists_sr, npivots[0], npivots[1], pivot_xyzs[0], pivot_xyzs[1],
        frame1=frames[0], frame2=frames[1],
        d1dists=d1dists, d2dists=d2dists,
        t1angs=t1angs, t2angs=t2angs,
        r2dists=r2dists_sr,
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


def build_machinefile_str():
    """ Take machine list and write the string for the machine file
    """

    host_node = get_host_node()
    num_cores = '10'

    machines = ['{}:{}'.format(host_node, num_cores)]
    machine_file_str = ''
    for machine in machines:
        machine_file_str += machine + '\n'

    return machine_file_str


# FUNCTION TO SET UP THE FRAGMENT GEOMETRIES FOR THE STRUCTURE.INP FILE
def fragment_geometries(ts_zma, rct_zmas, min_idx, max_idx):
    """ Generate the fragment geometries from the ts Z-matrix and the
        indices involved in the forming bond
    """

    # Get geometries of fragments from the ts_zma from the MEP
    print('amol zma\n', ts_zma)
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


# VARIOUS JOB RUN FUNCTIONS
def _build_vrctst_fs(ts_run_fs):
    """ build the filesystem and return the path
    """

    ts_fs, _ = ts_run_fs
    ts_run_path = ts_fs[0].path()
    bld_locs = ['VARECOF', 0]
    bld_save_fs = autofile.fs.build(ts_run_path)
    bld_save_fs[-1].create(bld_locs)
    vrc_path = bld_save_fs[-1].path(bld_locs)
    os.makedirs(os.path.join(vrc_path, 'scratch'), exist_ok=True)

    print('Build Path for VaReCoF calculations')
    print(vrc_path)

    return vrc_path


def _write_varecof_inp(varecof_inp, vrc_path):
    """ Write all of the VaReCoF inut files and run the code
    """

    exe_names = ('molpro.sh')

    # Write all of the VaReCoF input files
    for (inp_str, inp_name) in varecof_inp:
        file_name = os.path.join(vrc_path, inp_name)
        # Write the file
        with open(file_name, 'w') as inp_file:
            inp_file.write(inp_str)
        # Make file an executable
        if inp_name in exe_names:
            os.chmod(file_name, mode=os.stat(file_name).st_mode | stat.S_IEXEC)


def _run_varecof(vrc_path):
    """ Write all of the VaReCoF inut files and run the code
    """

    # Run VaReCoF
    run_script(DEFAULT_SCRIPT_DCT['varecof'], vrc_path)

    # Calculate the flux file from the output
    print('Generating flux file with TS N(E) from VaReCoF output...')
    run_script(DEFAULT_SCRIPT_DCT['mcflux'], vrc_path)


# def _varecof_success(vrc_path):
#     """ Check for success of the VaReCoF run and flux file generation
#     """
#     mcflux_file = os.path.join(vrc_path, 'mc_flux.out')
#     if os.path.exists(mcflux_file):
#         with open(mcflux_file, 'r') as fluxfile:
#             flux_str = mcflux_file.read()
#         if flux_str != '':
#             ts_found = True
#
#     return ts_found

# def _save_flux():
#     """ Save the VaReCoF flux file
#     """
#     ts_save_fs[0].create()
#     ts_save_path = ts_save_fs[0].path()
#
#     # Write the VRC-Flux and info files
#     ts_save_fs[0].file.vrc_flux.write(flux_str)
#     ts_save_fs[0].file.vrc_flux_info.write(flux_str)


def _vrc_dct():
    """ Build VRC dict
    """
    return {
        'fortran_compiler': 'gfortran',
        'base_name': 'mol',
        'spc_name': 'mol',
        'memory': 4.0,
        'r1dists_lr': [8., 6., 5., 4.5, 4.],
        'r1dists_sr': [4., 3.8, 3.6, 3.4, 3.2, 3., 2.8, 2.6, 2.4, 2.2],
        'r2dists_sr': [4., 3.8, 3.6, 3.4, 3.2, 3., 2.8, 2.6, 2.4, 2.2],
        'd1dists': [0.01, 0.5, 1.],
        'd2dists': [0.01, 0.5, 1.],
        'conditions': {},
        'nsamp_max': 2000,
        'nsamp_min': 50,
        'flux_err': 10,
        'pes_size': 2,
    }
