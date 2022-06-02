""" Generate the information necessary to product the vrctst input files
"""

import ioformat
import automol
import autofile
from phydat import phycon
import varecof_io
import elstruct
import autorun
from mechanalyzer.inf import thy as tinfo
from mechanalyzer.inf import rxn as rinfo
from mechlib.amech_io.printer import info_message, warning_message, running
from mechlib import filesys
from mechroutines.es.runner import scan
from mechroutines.es.runner import qchem_params
from mechroutines.es._routines import sp
from mechroutines.es.ts import _rpath as rpath


# CENTRAL FUNCTION TO WRITE THE VARECOF INPUT FILES AND RUN THE PROGRAM
def calc_vrctst_flux(ts_dct,
                     thy_inf_dct, thy_method_dct, mref_params,
                     es_keyword_dct,
                     runfs_dct, savefs_dct):
    """ Set up n VRC-TST calculations to get the flux file
    """

    # Build VRC-TST filesys stuff
    vrc_fs = runfs_dct['vrctst']
    vrc_path = vrc_fs[-1].create((0,))
    vrc_path = vrc_fs[-1].path((0,))
    vrc_dct = autorun.varecof.VRC_DCT  # need code to input one

    # Get a bunch of info that describes the grid
    scan_inf_dct = _scan_inf_dct(ts_dct, savefs_dct)

    # print the guess ts zma
    print('guess zma\n')
    print(automol.geom.string(automol.zmat.geometry(ts_dct['zma'])))

    # Run and Read all of the info for the correction potential
    inf_sep_ene, potentials, pot_lbls, zma_for_inp = _correction_pot_data(
        ts_dct, scan_inf_dct,
        thy_inf_dct, thy_method_dct, mref_params,
        es_keyword_dct,
        runfs_dct, savefs_dct)

    # Write the VaReCoF input files
    inp_strs = autorun.varecof.write_input(
        vrc_path,
        zma_for_inp, scan_inf_dct['rct_zmas'],
        len(potentials), scan_inf_dct['rxn_bond_keys'],
        vrc_dct, machine_dct=None)

    rxn_info = ts_dct['rxn_info']
    ts_info = rinfo.ts_info(rxn_info)
    mod_var_scn_thy_info = thy_inf_dct['mod_var_scnlvl']
    cas_kwargs = mref_params['var_scnlvl']
    molp_tmpl_str = varecof_io.writer.molpro_template(
        ts_info, mod_var_scn_thy_info, inf_sep_ene, cas_kwargs)

    inp_strs.update({'run.tml': molp_tmpl_str})

    # Build correction potential .so file used by VaReCoF
    ndummy_added = 0
    for zma in scan_inf_dct['rct_zmas']:
        if automol.zmat.count(zma) > 2:
            ndummy_added += 1
    aidx = scan_inf_dct['rxn_bond_keys'][0]+1
    bidx = scan_inf_dct['rxn_bond_keys'][1]+1+ndummy_added

    autorun.varecof.compile_potentials(
        vrc_path, scan_inf_dct['full_grid'], potentials,
        aidx, bidx, vrc_dct['fortran_compiler'],
        dist_restrict_idxs=(),
        pot_labels=pot_lbls,
        pot_file_names=[vrc_dct['spc_name']],
        spc_name=vrc_dct['spc_name'])

    # Run VaReCoF to generate flux file
    running(f'VaReCoF at {vrc_path}')
    flux_str = autorun.varecof.flux_file(
        autorun.SCRIPT_DCT['varecof_multi'],
        autorun.SCRIPT_DCT['varecof_conv_multi'],
        autorun.SCRIPT_DCT['varecof_mcflux'],
        vrc_path, inp_strs,
        nprocs=es_keyword_dct['varecof_nprocs'])

    # Read the corr pot file to send to save
    pot_corr_str = ioformat.pathtools.read_file(
        vrc_path, 'run_corr.f')

    # Save the flux file
    if flux_str is not None:
        filesys.save.flux(flux_str, pot_corr_str, inp_strs,
                          savefs_dct['vrctst'], vrc_locs=(0,))
        success = True
    else:
        success = False

    return success


def _scan_inf_dct(ts_dct, savefs_dct):
    """ Determine all informationa about the scans and guess information
    """

    # Build initial coord, grid, and other info
    zrxn, ts_zma = ts_dct['zrxn'], ts_dct['zma']

    cls = ts_dct['class']
    scan_inf = automol.reac.build_scan_info(
        zrxn, ts_zma,
        var=(automol.par.is_radrad(cls) and automol.par.is_low_spin(cls)))
    coord_names, _, coord_grids, update_guess = scan_inf

    # Get fol constraint dct
    _rcts_cnf_fs = savefs_dct['rcts_cnf']
    rct_zmas = ()
    for (_, cnf_save_fs, min_locs, _) in _rcts_cnf_fs:
        zma_fs = autofile.fs.zmatrix(cnf_save_fs[-1].path(min_locs))
        rct_zmas += (zma_fs[-1].file.zmatrix.read((0,)),)
    constraint_dct = varecof_io.writer.intramolecular_constraint_dct(
        ts_zma, rct_zmas)

    # Get indices for potentials and input
    frm_bnd_key, = automol.graph.ts.forming_bond_keys(zrxn.forward_ts_graph)

    # Set up zma for the scan
    inf_sep_zma = automol.zmat.set_values_by_name(
        ts_zma, {coord_names[0]: coord_grids[1][-2]}, angstrom=False)

    # set up grid from lowest R to largest R
    full_grid = tuple(sorted(list(coord_grids[0]) + list(coord_grids[1][1:])))

    return {
        'coord_names': coord_names,
        'coord_grids': coord_grids,
        'inf_sep_zma': inf_sep_zma,
        'grid_val_for_zma': coord_grids[0][0],  # first val should be ~ 2 Ang.
        'inf_locs': (coord_names, (coord_grids[1][-1],)),
        'full_grid': full_grid,
        'update_guess': update_guess,
        'constraint_dct': constraint_dct,
        'rxn_bond_keys': (min(frm_bnd_key), max(frm_bnd_key)),
        'rct_zmas': rct_zmas
    }


# FUNCTIONS TO SET UP THE libcorrpot.so FILE USED BY VARECOF
def _correction_pot_data(ts_dct, scan_inf_dct,
                         thy_inf_dct, thy_method_dct, mref_params,
                         es_keyword_dct,
                         runfs_dct, savefs_dct):
    """  use the MEP potentials to compile the correction potential .so file
    """
    rxn_info = ts_dct['rxn_info']
    ts_info = rinfo.ts_info(rxn_info)
    ts_geo = automol.zmat.geometry(ts_dct['zma'])

    # Run the constrained and full opt potential scans
    _run_potentials(ts_info, ts_geo,
                    scan_inf_dct,
                    thy_inf_dct, thy_method_dct, mref_params,
                    es_keyword_dct, runfs_dct, savefs_dct)

    # Obtain the energy at infinite separation
    inf_sep_ene = rpath.inf_sep_ene(
        ts_dct, thy_inf_dct, thy_method_dct, mref_params,
        savefs_dct, runfs_dct, es_keyword_dct)

    # Read the values for the correction potential from filesystem
    potentials, pot_labels, zma_for_inp = _read_potentials(
        scan_inf_dct, thy_inf_dct, savefs_dct)

    # Set zma if needed
    if zma_for_inp is None:
        zma_for_inp = ts_dct['zma']

    return inf_sep_ene, potentials, pot_labels, zma_for_inp


def _run_potentials(ts_info, ts_geo,
                    scan_inf_dct,
                    thy_inf_dct, thy_method_dct, mref_params,
                    es_keyword_dct, runfs_dct, savefs_dct):
    """ Run and save the scan along both grids while
          (1) optimization: constraining only reaction coordinate, then
          (2) optimization: constraining all intermolecular coordinates
          (3) single-point energy on scan (1)
    """

    # Get fs and method objects
    scn_run_fs = runfs_dct['vscnlvl_scn']
    scn_save_fs = savefs_dct['vscnlvl_scn']
    cscn_run_fs = runfs_dct['vscnlvl_cscn']
    cscn_save_fs = savefs_dct['vscnlvl_cscn']
    # sp_scn_save_fs = savefs_dct['vscnlvl_scn']
    scn_thy_info = thy_inf_dct['mod_var_scnlvl']
    sp_thy_info = thy_inf_dct['mod_var_splvl1']

    # Set up the options for the Molpro input strings
    opt_script_str, opt_kwargs = qchem_params(
        thy_method_dct['var_scnlvl'], elstruct.Job.OPTIMIZATION,
        geo=ts_geo, spc_info=ts_info)
    cas_kwargs = mref_params['var_scnlvl']
    opt_kwargs.update(cas_kwargs)

    sp_script_str, sp_kwargs = qchem_params(
        thy_method_dct['var_splvl1'],
        geo=ts_geo, spc_info=ts_info)
    sp_cas_kwargs = mref_params['var_splvl1']
    sp_kwargs.update(sp_cas_kwargs)

    # Run optimization scans
    for constraints in (None, scan_inf_dct['constraint_dct']):
        if constraints is None:
            _run_fs = scn_run_fs
            _save_fs = scn_save_fs
            info_message('Running full scans..', newline=1)
        else:
            _run_fs = cscn_run_fs
            _save_fs = cscn_save_fs
            info_message('Running constrained scans..', newline=1)

        info_message('Method:', tinfo.string(scn_thy_info))

        # Loop over grids (both should start at same point and go in and out)
        for grid in scan_inf_dct['coord_grids']:
            info_message(f'Grid: {grid}')
            scan.execute_scan(
                zma=scan_inf_dct['inf_sep_zma'],
                spc_info=ts_info,
                mod_thy_info=thy_inf_dct['mod_var_scnlvl'],
                coord_names=scan_inf_dct['coord_names'],
                coord_grids=(grid,),
                scn_run_fs=_run_fs,
                scn_save_fs=_save_fs,
                scn_typ='relaxed',
                script_str=opt_script_str,
                overwrite=es_keyword_dct['overwrite'],
                update_guess=scan_inf_dct['update_guess'],
                reverse_sweep=False,
                saddle=False,
                constraint_dct=constraints,
                retryfail=True,
                **opt_kwargs
            )
            info_message('')

    # Run the single points on top of the initial, full scan
    if sp_thy_info is not None:
        info_message('')
        info_message('Running single-point calculations on the full scan...')
        info_message('Method:', tinfo.string(scn_thy_info, sp_thy_info))
        for locs in scn_save_fs[-1].existing((scan_inf_dct['coord_names'],)):
            scn_run_fs[-1].create(locs)
            geo = scn_save_fs[-1].file.geometry.read(locs)
            zma = scn_save_fs[-1].file.zmatrix.read(locs)
            sp.run_energy(zma, geo, ts_info, sp_thy_info,
                          scn_run_fs, scn_save_fs, locs, runfs_dct['prefix'],
                          sp_script_str, es_keyword_dct['overwrite'],
                          highspin=False, **sp_kwargs)


def _read_potentials(scan_inf_dct, thy_inf_dct, savefs_dct):
    """ Read values form the filesystem to get the values to
        correct ht MEP
    # Read the energies from the full and constrained opts along MEP
    """

    scn_save_fs = savefs_dct['vscnlvl_scn']
    cscn_save_fs = savefs_dct['vscnlvl_cscn']
    mod_var_scn_thy_info = thy_inf_dct['mod_var_scnlvl']
    mod_var_sp1_thy_info = thy_inf_dct['mod_var_splvl1']
    coord_name = scan_inf_dct['coord_names'][0]
    full_grid = scan_inf_dct['full_grid']
    constraint_dct = scan_inf_dct['constraint_dct']
    grid_val_for_zma = scan_inf_dct['grid_val_for_zma']

    # build objects for loops
    smp_pot, const_pot, sp_pot = [], [], []
    scans = (
        (scn_save_fs, mod_var_scn_thy_info),
        (cscn_save_fs, mod_var_scn_thy_info)
    )
    if mod_var_sp1_thy_info is not None:
        scans += ((scn_save_fs, mod_var_sp1_thy_info),)

    for idx, (scn_fs, thy_info) in enumerate(scans):
        for grid_val in full_grid:
            if idx in (0, 2):
                locs = [[coord_name], [grid_val]]
            else:
                locs = [constraint_dct, [coord_name], [grid_val]]
            sp_ene = filesys.read.energy(scn_fs, locs, thy_info)
            print('scn_fs: ', scn_fs)
            print('locs: ', locs)
            print('thy_info: ', thy_info)
            print('sp_ene: ', sp_ene)

            # Store the energy in a lst
            if idx == 0:
                smp_pot.append(sp_ene)
            elif idx == 1:
                const_pot.append(sp_ene)
            elif idx == 2:
                sp_pot.append(sp_ene)

    print('SUMMARY OF POTENTIALS:')
    print(' - SAMPLING POT:', smp_pot)
    print(' - CONSTRAINED POT:', const_pot)
    print(' - SP POT:', sp_pot)

    sp_corr_inf = (sp_pot[-1] - smp_pot[-1])

    # Calculate each of the correction potentials
    relax_corr_pot, sp_corr_pot, full_corr_pot = [], [], []
    for i, _ in enumerate(smp_pot):
        # We assume relax_corr_inf = 0
        relax_corr = (smp_pot[i] - const_pot[i]) * phycon.EH2KCAL
        relax_corr_pot.append(relax_corr)
        if all(ene is not None for ene in sp_pot):
            sp_corr = ((sp_pot[i] - smp_pot[i]) - sp_corr_inf) * phycon.EH2KCAL
            sp_corr_pot.append(sp_corr)
        else:
            warning_message("No single point correction applied to potential")
            sp_corr = 0.0
        full_corr_pot.append(relax_corr + sp_corr)

    # Collate the potentials together in a list
    if sp_pot:
        potentials = [full_corr_pot, relax_corr_pot, sp_corr_pot]
        potential_labels = ['full', 'relax', 'sp']
    else:
        potentials = [full_corr_pot, relax_corr_pot]
        potential_labels = ['full', 'relax']

    # Get zma used to make structure.inp and divsur.inp
    inp_zma_locs = [[coord_name], [grid_val_for_zma]]
    if scn_save_fs[-1].file.zmatrix.exists(inp_zma_locs):
        zma_for_inp = scn_save_fs[-1].file.zmatrix.read(inp_zma_locs)
        print('Path for getting Z-Matrix to set dummy atom location'
              'for structure.inp file for VaReCoF:')
        print('  ', scn_save_fs[-1].file.zmatrix.path(inp_zma_locs))
    else:
        zma_for_inp = None

    return potentials, potential_labels, zma_for_inp
