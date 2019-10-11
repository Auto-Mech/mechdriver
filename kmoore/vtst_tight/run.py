"""
Do VTST for tight transition states
"""

import elstruct
import mess_io


TS_ZMA = ''
TS_LABEL = 'B1'
REAC_LABEL = 'R1'
PROD_LABEL = 'P1'


def _run_irc(
        spc_info, thy_level, cnf_run_fs, cnf_save_fs,
        script_str, overwrite, scan_increment=30., saddle=False, tors_names='', new_grid=False, **opt_kwargs):
    """ Run the IRC 
    """
        
    min_cnf_locs = moldr.util.min_energy_conformer_locators(cnf_save_fs)
    if min_cnf_locs:
        min_cnf_run_path = cnf_run_fs.leaf.path(min_cnf_locs)
        min_cnf_save_path = cnf_save_fs.leaf.path(min_cnf_locs)
        scn_run_fs = autofile.fs.scan(min_cnf_run_path)
        scn_save_fs = autofile.fs.scan(min_cnf_save_path)
        zma = cnf_save_fs.leaf.file.zmatrix.read(min_cnf_locs)

    # Run the IRC
    moldr.driver.run_job(
        job='irc',
        script_str=opt_script_str,
        run_fs=run_fs,
        geom=zma,
        spc_info=ts_info,
        thy_level=ref_level,
        saddle=True,
        overwrite=overwrite,
        **opt_kwargs,
        )

    opt_ret = moldr.driver.read_job(
        job='irc',
        run_fs=run_fs,
    )
    if opt_ret is not None:
        inf_obj, _, out_str = opt_ret
        prog = inf_obj.prog
        method = inf_obj.method
        geos, gras, hessians = elstruct.reader.irc_points(prog, out_str)

        print(" - Saving...")
        print(" - Save path: {}".format(ts_save_path))

        for locs in locs:
        for geo, gra, hess in zip(geos, gras, hessians):
            ts_save_fs.trunk.file.energy.write(ene)
            ts_save_fs.trunk.file.geometry.write(geo)
            ts_save_fs.trunk.file.zmatrix.write(zma)
                
        save_scan(
            scn_run_fs=scn_run_fs,
            scn_save_fs=scn_save_fs,
            coo_names=[tors_name],
        )


    return geoms, grads, hessians


def _project_frequencies(geoms, grads, hessians):
    """ fake fxn that gets frequencies (sets up and runs projrot)
    """
    return freqs


def _get_irc_high_energy(geoms, high_e_thy_lvl):
    """ fake fxn to get the high level energies
    """
    return energies


def _build_mess_input(geoms, frequencies, energies):
    """ fake fxn to build the mess input for a variational barrier
    """
    # Determine the the number of points along the irc
    nirc = 21
    
    # Loop over all the points of the irc and build MESS strings
    irc_pt_strings = []
    for i in range(nirc):

        # Iniialize the header of the string
        irc_pt_string = '!-----------------------------------------------'
        irc_pt_string += '! IRC Point {0}\n'.format(str(i+1))

        # Write the molecule section for each irc point
        core = mess_io.writer.mol_data.core_rigidrotor(geom1, sym_factor, interp_emax=''):
        irc_pt_str += mess_io.writer.species.molecule(core, freqs, elec_levels,
             hind_rot='', xmat=None, rovib_coups='', rot_dists=''):

        # Append the zero point energy for the molecule
        irc_pt_str += mess_io.writer.species.molecule(core, freqs, elec_levels,
        irc_pt_str += '    ZeroEnergy[kcal/mol]      {0:}'.format(zero_energy)

        # Append string to list
        irc_pt_strings.append(irc_pt_string)
    
    # Write the MESS string for the variational sections
    varational_str = mess_io.writer.rxnchan.ts_variational(
        ts_label, reac_label, prod_label, irc_pt_strings)

    return variational_str


def _do_sct_projections(geoms, grads, hessians):
    """ fake fxn to do the projections needed for SCT
    """
    return freqs_sct


def _make_mess_sct_files():
    """ writes the files that MESS needs to do SCT
    """
    return mess_sct_string


# Call the job runner and do the variational calculations
geoms, grads, hessians = _run_irc(ts_zma, thy_lvl)
energies = _get_irc_high_energy(geoms, high_e_thy_lvl)
if SCT:
    _do_sct_projections(geoms, grads, hessians)
    _make_mess_sct_files()
_build_mess_input(geoms, frequencies, energies)
