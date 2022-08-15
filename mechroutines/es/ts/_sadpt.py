""" Saddle point finding functions x
"""

import autofile
import autorun
import elstruct
import automol
from mechanalyzer.inf import thy as tinfo
from mechanalyzer.inf import rxn as rinfo
from mechlib.amech_io import printer as ioprinter
from mechlib import filesys
from mechroutines.es import runner as es_runner
from mechroutines.es.runner import qchem_params
from mechroutines.es.ts import _rpath as rpath


# Functions to assess the status of existing saddle point structures in SAVE
def read_existing_saddle_points(spc_dct, tsname, savefs_dct, zma_locs=None):
    """ Searches for and reads out, if present, the Z-matrix for a
        conformer of the saddle-point in the SAVE filesystem for
        the electronic structure method specified in the
        user-defined 'runlvl' and 'inplvl' keywords that are
        stored in the savefs_dct.

        :param savefs_dct: filesystem objects under SAVE prefix
        :type savefs_dct: dict[str: autofile.fs objects]
        :rtype: tuple(automol.zmat object)
    """
    zmas = [[], None]
    for idx, _fs in enumerate(('runlvl_cnf_tuple', 'inplvl_cnf_tuple')):
        lvl = 'runlvl' if idx == 0 else 'inplvl'
        ioprinter.info_message(
            '\nSearching save filesys for Z-Matrix calculated',
            f'at {lvl} level...')
        cnf_fs, cnf_locs = savefs_dct[_fs]
        if any(cnf_locs):
            zma_fs = autofile.fs.zmatrix(cnf_fs[-1].path(cnf_locs))
            for zma_locs in zma_fs[-1].existing():
                geo_path = zma_fs[-1].file.zmatrix.path(zma_locs)
                ioprinter.info_message(
                    ' - Z-Matrix found.')
                ioprinter.info_message(
                    f' - Reading Z-Matrix from path {geo_path}')
                if idx == 1:
                    zmas[idx] = zma_fs[-1].file.zmatrix.read(zma_locs)
                else:
                    zmas[idx].append((
                        cnf_fs, cnf_locs, zma_locs,
                        zma_fs[-1].file.zmatrix.read(zma_locs),))
    return tuple(zmas)


def search_required(runlvl_zma, es_keyword_dct):
    """ Determine if MechDriver needs to search for the saddle
        point of some reaction.

        If a saddle point was found, the function assesses whether
        the overwrite keyword has been set to True, which would
        requiring a new search from scratch.

        :param savefs_dct: filesystem objects under SAVE prefix
        :type savefs_dct: dict[str: autofile.fs objects]
        :rtype: bool
    """

    overwrite = es_keyword_dct['overwrite']
    if runlvl_zma == []:
        ioprinter.info_message(
            '\nSince no transition state found in filesys',
            f'at {es_keyword_dct["runlvl"]} level',
            'proceeding to find it...')
        _run = True
    else:
        if overwrite:
            print('\nUser specified to overwrite transition state search.'
                  'Redoing task...')
            _run = True
        else:
            print('\nSince transition state found and saved previously,',
                  'proceeding to next task.')
            _run = False

    return _run


# Functions to attempt to find, optimize, and save valid saddle point
def search(ini_zma, spc_dct, tsname,
           thy_inf_dct, thy_method_dct, mref_dct,
           es_keyword_dct,
           runfs_dct, savefs_dct):
    """ Attempt to locate and optimize a proper transition state.
    """
    # CHECK
    _ = thy_inf_dct

    # Initialize success variable to False
    success = False

    # Build ts dct for functions
    ts_dct = spc_dct[tsname]

    if ini_zma is not None:
        # Use the zma from the ini level as a guess
        guess_zmas = (ini_zma,)
    else:
        # Generate a guess zma by scanning along rxncoord and finding max
        guess_zmas = rpath.internal_coordinates_scan(
            ts_zma=ts_dct['zma'],
            zrxn=ts_dct['zrxn'],
            ts_info=rinfo.ts_info(ts_dct['canon_rxn_info']),
            rxn_class=ts_dct['class'],
            method_dct=thy_method_dct['runlvl'],
            mref_params=mref_dct['runlvl'],
            scn_run_fs=runfs_dct['runlvl_scn'],
            scn_save_fs=savefs_dct['runlvl_scn'],
            cscn_run_fs=runfs_dct['runlvl_cscn'],
            cscn_save_fs=savefs_dct['runlvl_cscn'],
            es_keyword_dct=es_keyword_dct)

    # Optimize guess and check if saddle point good to save
    if guess_zmas is not None:

        cnf_locs = (autofile.schema.generate_new_ring_id(),
                    autofile.schema.generate_new_conformer_id())

        opt_ret, hess_ret = optimize_saddle_point(
            guess_zmas, ts_dct,
            thy_method_dct['runlvl'], mref_dct['runlvl'],
            runfs_dct, es_keyword_dct,
            cnf_locs)

        status = assess_saddle_point(opt_ret, hess_ret,
                                     runfs_dct, cnf_locs)
        if status == 'save':
            save_saddle_point(opt_ret, hess_ret,
                              ts_dct, thy_method_dct['runlvl'],
                              savefs_dct, cnf_locs)
            success = True
            # print the saddle point

    return success


def save_reversed_zma(
        spc_dct, tsname,
        cnf_fs, cnf_locs, orig_zma_locs, zma):
    """
    """
    ts_dct = spc_dct[tsname]
    # Pull info from the dictionaries to save
    zrxn = ts_dct['zrxn']
    zma_fs = autofile.fs.zmatrix(cnf_fs[-1].path(cnf_locs))
    zma_locs = [zma_fs[-1].existing()[-1][0] + 1]
    inf_obj = zma_fs[-1].file.geometry_info.read(orig_zma_locs)
    inp_str = zma_fs[-1].file.geometry_input.read(orig_zma_locs)

    filesys.save._save_zmatrix_parsed(zma, inf_obj, inp_str, zma_fs, zma_locs)
    filesys.save._conformer_aux_info(zma_fs, zma_locs, zrxn=zrxn)
    # # Save initial saddle point conformer



def save_opt_from_run(spc_dct, tsname,
                      thy_method_dct,
                      runfs_dct, savefs_dct):
    """ Look in the run filesystem for an opt+hess job and try and read it
        to save the conformer
    """

    ts_dct = spc_dct[tsname]
    runlvl_cnf_run_fs = runfs_dct['runlvl_cnf']

    # just grab the first one for now, don't have better idea:
    existing_locs = runlvl_cnf_run_fs[-1].existing()
    if any(existing_locs):
        # just grab the first one for now, don't have better idea:
        cnf_locs = existing_locs[0]
    else:
        cnf_locs = None

    # If locs exist then try and read the existing job
    success = False
    if cnf_locs is not None:
        runlvl_cnf_run_fs = runfs_dct['runlvl_cnf']
        run_fs = autofile.fs.run(runlvl_cnf_run_fs[-1].path(cnf_locs))

        ioprinter.info_message(
            '\nSearching for OPT job to save in the RUN filesys '
            f'at path {run_fs[0].path()}')
        opt_success, opt_ret = es_runner.read_job(
            job='optimization',
            run_fs=run_fs)
        ioprinter.info_message(
            'Searching for HESS job to save in the RUN filesys '
            f'at path {run_fs[0].path()}')
        hess_success, hess_ret = es_runner.read_job(
            job='hessian',
            run_fs=run_fs)

        if opt_success and hess_success:
            ioprinter.info_message(
                'Saving geometry and hessian of TS conformer from RUN filesys '
                f'at path {run_fs[0].path()}')
            status = assess_saddle_point(opt_ret, hess_ret,
                                         runfs_dct, cnf_locs)
            if status == 'save':
                save_saddle_point(opt_ret, hess_ret,
                                  ts_dct, thy_method_dct['runlvl'],
                                  savefs_dct, cnf_locs)
                success = True
    else:
        ioprinter.info_message(
            '\nNo optimization job of TS conformer found in '
            'RUN filesys to save')

    return success


def optimize_saddle_point(guess_zmas, ts_dct,
                          method_dct, mref_kwargs,
                          runfs_dct, es_keyword_dct,
                          cnf_locs):
    """ Optimize the transition state structure obtained from the grid search
    """

    # Get info (move later)
    ts_info = rinfo.ts_info(ts_dct['canon_rxn_info'])
    mod_thy_info = tinfo.modify_orb_label(tinfo.from_dct(method_dct), ts_info)

    overwrite = es_keyword_dct['overwrite']
    ts_info = rinfo.ts_info(ts_dct['canon_rxn_info'])

    # Set the run filesystem for the job
    runlvl_cnf_run_fs = runfs_dct['runlvl_cnf']
    runlvl_cnf_run_fs[-1].create(cnf_locs)
    run_fs = autofile.fs.run(runlvl_cnf_run_fs[-1].path(cnf_locs))
    
    ioprinter.info_message(
        '\nAttempting to get optimized TS from guess Z-Matrices')
    ioprinter.info_message(
        f'There are {len(guess_zmas)} guess Z-Matrices '
        'to attempt to find saddle point.', newline=1)

    # Loop over all the guess zmas to find a TS
    opt_ret, hess_ret = None, None
    for idx, zma in enumerate(guess_zmas):
        ioprinter.info_message(
            f'\nOptimizing guess Z-Matrix {idx+1}...')

        # Run the transition state optimization
        script_str, kwargs = qchem_params(
            method_dct, job=elstruct.Job.OPTIMIZATION,
            geo=automol.zmat.geometry(zma), spc_info=ts_info)
        kwargs.update(mref_kwargs)

        opt_success, opt_ret = es_runner.execute_job(
            job='optimization',
            script_str=script_str,
            run_fs=run_fs,
            geo=zma,
            spc_info=ts_info,
            thy_info=mod_thy_info,
            saddle=True,
            overwrite=overwrite,
            **kwargs,
            )

        if opt_success:
            break

    if opt_success:
        # Obtain geometry from optimization
        opt_inf_obj, _, opt_out_str = opt_ret
        opt_prog = opt_inf_obj.prog
        geo = elstruct.reader.opt_geometry(opt_prog, opt_out_str)

        # Set up the script str
        script_str, kwargs = qchem_params(
            method_dct,
            geo=geo, spc_info=ts_info)

        # Run a Hessian
        _, hess_ret = es_runner.execute_job(
            job='hessian',
            script_str=script_str,
            run_fs=run_fs,
            geo=geo,
            spc_info=ts_info,
            thy_info=mod_thy_info,
            overwrite=overwrite,
            **kwargs,
            )

    return opt_ret, hess_ret


# Checker functions
# def assess_saddle_point(opt_ret, hess_ret, ts_dct, runfs_dct, cnf_locs):
def assess_saddle_point(opt_ret, hess_ret, runfs_dct, cnf_locs):
    """ run things for checking Hessian
        If successful, Read the geom and energy from the optimization
    """

    status = 'failure'
    if hess_ret is not None:

        # Get the physical info used for the checks
        # zrxn = ts_dct['zrxn']

        opt_inf, _, opt_out_str = opt_ret
        hess_inf, _, hess_out_str = hess_ret

        # zma = elstruct.reader.opt_zmatrix(opt_inf.prog, opt_out_str)
        geo = elstruct.reader.opt_geometry(opt_inf.prog, opt_out_str)
        hess = elstruct.reader.hessian(hess_inf.prog, hess_out_str)

        # Set filesys information
        runlvl_cnf_run_fs = runfs_dct['runlvl_cnf']
        run_fs = autofile.fs.run(runlvl_cnf_run_fs[-1].path(cnf_locs))
        freq_run_path = run_fs[-1].path(['hessian'])
        run_fs[-1].create(['hessian'])

        # Freq magnitude check
        script_str = autorun.SCRIPT_DCT['projrot']
        _, _, imags, _ = autorun.projrot.frequencies(
            script_str, freq_run_path, [geo], [[]], [hess])

        freq_success = _check_freqs(imags)

        # ted coordinate check
        # if not automol.zmat.dummy_keys(zma):
        #     script_str = autorun.SCRIPT_DCT['intder']
        #     ted_names = autorun.intder.ted_zmatrix_coordinates(
        #         script_str, freq_run_path, geo, zma, hess, 0)
        #     ted_success = _ted_coordinate_check(ted_names, zrxn, zma)
        # else:
        #     print('Z-Matrix has dummy atoms, cannot do TED check')
        #     ted_success = True
        ted_success = True

        # Set overall success value to return
        if freq_success == 'kick':
            status = 'kick'
        else:
            if freq_success and ted_success:
                status = 'save'
            else:
                status = 'failure'

    return status


def _check_freqs(imags):
    """ Check the magnitude of the imaginary modes.
    """

    big_imag, kick_imag = 0, 0

    ioprinter.checking('the imaginary frequencies of the saddle point...')
    if len(imags) < 1:
        ioprinter.warning_message('No imaginary modes for geometry')
        status = 'fail'
    else:
        if len(imags) > 1:
            ioprinter.warning_message(
                'More than one imaginary mode for geometry')
            status = 'fail'
        for idx, imag in enumerate(imags):
            if imag <= 50.0:
                ioprinter.warning_message(
                    f'Mode {idx+1} {imag} cm-1 is low,')
            elif 50.0 < imag <= 200.0:
                lowstr = f'Mode {idx+1} {imag} cm-1 is low,'
                ioprinter.debug_message(
                    lowstr + ' check mode and see if it should be corrected')
                big_imag += 1
                # Adding to the kick counter kills code for good TSs
                # Some addditions of big species have low mode of this
                # ioprinter.warning_message(
                #     lowstr + 'need a kickoff procedure to remove')
                # kick_imag += 1
            else:
                ioprinter.debug_message(
                    f'Mode {idx+1} {imag} cm-1 is likely fine,')
                big_imag += 1

        if big_imag > 1:
            ioprinter.warning_message(
                'More than one imaginary mode for geometry')
            if kick_imag >= 1:
                ioprinter.debug_message('Will kickoff to get saddle point')
                status = 'kick'
            else:
                status = False
        elif big_imag == 1:
            status = True
        elif big_imag == 0:
            status = False

    return status


def _ted_coordinate_check(ted_names, zrxn, zma):
    """ Assess if the coordinate check can be performed
    """

    if ted_names is not None:

        # Get the zmat names corresponding to frm/brk keys; remove Nones
        rxn_names = automol.reac.zmatrix_coordinate_names(zrxn, zma)
        rxn_names = (tuple(x for x in rxn_names[0] if x is not None) +
                     tuple(x for x in rxn_names[1] if x is not None))

        print('Comparing Z-Matrix Coordinates.')
        tedname_str = ' '.join(ted_names)
        print(f'- TED: {tedname_str}'.format(tedname_str))
        rname_str = ' '.join(rxn_names)
        print(f'- Forming/Breaking Bonds: {rname_str}')

        if set(ted_names) & set(rxn_names):
            print('Overlap of coordinates found, possible success')
            success = True
        else:
            print('No similarity of coords, likely something is wrong')
            success = False
    else:
        print('INTDER had some error, skipping TED check')
        success = True

    return success


def ts_zma_locs(spc_dct, tsname, zma_fs, spc_dct_i=None):
    if spc_dct_i is None:
        spc_dct_i = spc_dct[tsname]
    zma_locs = None
    for zma_locs_i in zma_fs[-1].existing():
        if zma_fs[-1].file.zmatrix.exists(zma_locs_i):
            zma_i = zma_fs[-1].file.zmatrix.read(zma_locs_i)
            if automol.zmat.vmatrix(zma_i) == automol.zmat.vmatrix(spc_dct_i['zma']):
                zma_locs = zma_locs_i
    return zma_locs


# Save the saddle point
def save_saddle_point(opt_ret, hess_ret,
                      ts_dct, method_dct, savefs_dct,
                      cnf_locs):
    """ Given the saddle point guess structure, obtain a
        proper saddle point
    """

    # Pull info from the dictionaries to save
    zrxn = ts_dct['zrxn']
    runlvl_cnf_save_fs, _ = savefs_dct['runlvl_cnf_tuple']
    ts_info = rinfo.ts_info(ts_dct['canon_rxn_info'])
    mod_thy_info = tinfo.modify_orb_label(tinfo.from_dct(method_dct), ts_info)

    # Save initial saddle point conformer
    filesys.save.conformer(
        opt_ret, hess_ret, runlvl_cnf_save_fs, mod_thy_info[1:],
        zrxn=zrxn,
        rng_locs=(cnf_locs[0],),
        tors_locs=(cnf_locs[1],),
        zma_locs=None)
