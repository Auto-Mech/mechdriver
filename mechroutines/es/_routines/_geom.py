""" es_runners for initial geometry optimization
"""

import numpy
import automol
import elstruct
import autorun
from phydat import phycon
from mechanalyzer.inf import thy as tinfo
from mechlib.amech_io import printer as ioprinter
from mechroutines.es import runner as es_runner
from mechroutines.es.runner import qchem_params
# from mechroutines.es._routines.sp import rerun_hessian_and_opt


def remove_imag(geo, ini_ret, spc_info, method_dct, run_fs,
                kickoff_size=0.1, kickoff_backward=False, kickoff_mode=0,
                overwrite=False):
    """ if there is an imaginary frequency displace geometry along the imaginary
    mode and then reoptimize
    """

    thy_info = tinfo.from_dct(method_dct)
    mod_thy_info = tinfo.modify_orb_label(thy_info, spc_info)
    opt_script_str, opt_kwargs = qchem_params(
        method_dct, job=elstruct.Job.OPTIMIZATION)
    script_str, kwargs = qchem_params(
        method_dct)

    ioprinter.info_message(
        'The initial geometries will be checked for imaginary frequencies')

    imag, norm_coords = _check_imaginary(
        spc_info, geo, mod_thy_info, run_fs, script_str,
        overwrite, **kwargs)

    # Make five attempts to remove imag mode if found
    count, kick_ret = 1, None
    while imag and count <= 2:

        # Either attempt a kickoff-reopt or tight-reopt
        # Only attempt kickoff-reopt on first two tries
        if count <= 2:
            ioprinter.info_message(
                f'Attempting kick off along mode, attempt {count}...')

            geo, kick_ret = _kickoff_saddle(
                geo, norm_coords, spc_info, mod_thy_info,
                run_fs, opt_script_str,
                kickoff_size=kickoff_size,
                kickoff_backward=kickoff_backward,
                kickoff_mode=kickoff_mode,
                opt_cart=True, **opt_kwargs)

            ioprinter.info_message(
                'Removing faulty geometry from SAVE filesys. '
                'Rerunning Hessian...')
            run_fs[-1].remove([elstruct.Job.HESSIAN])
        else:
            pass
            # rerun_hessian_and_opt(
            #     zma, spc_info, thy_info,
            #     geo_run_fs, geo_save_fs, locs,
            #     script_str, zrxn=None,
            #     retryfail=True, method_dct=None, attempt=0,
            #     hess_script_str=None, **hess_kwargs)

        # Assess the imaginary mode after the reoptimization
        ioprinter.info_message('Rerunning Hessian...')
        imag, norm_coords = _check_imaginary(
            spc_info, geo, mod_thy_info, run_fs, script_str,
            overwrite, **kwargs)

        # Update the loop counter and kickoff size
        count += 1
        kickoff_size *= 2.0

    if kick_ret is None:
        ret = ini_ret
    else:
        ret = kick_ret

    return geo, ret


def _check_imaginary(spc_info, geo, mod_thy_info, run_fs, script_str,
                     overwrite=False, **kwargs):
    """ check if species has an imaginary frequency
    """

    # Initialize info
    has_imag = False
    norm_coords = []
    hess = ((), ())

    # Run Hessian calculation
    success, ret = es_runner.execute_job(
        job=elstruct.Job.HESSIAN,
        spc_info=spc_info,
        thy_info=mod_thy_info,
        geo=geo,
        run_fs=run_fs,
        script_str=script_str,
        overwrite=overwrite,
        **kwargs,
        )

    # Check for imaginary modes
    if success:
        inf_obj, _, out_str = ret
        prog = inf_obj.prog
        hess = elstruct.reader.hessian(prog, out_str)

        # Calculate vibrational frequencies
        if hess:
            run_path = run_fs[-1].path([elstruct.Job.HESSIAN])
            # run_path = run_fs[-1].path(['VIB'])
            script_str = autorun.SCRIPT_DCT['projrot']
            _, _, imag_freq, _ = autorun.projrot.frequencies(
               script_str, run_path, [geo], [[]], [hess])

            # Mode for now set the imaginary frequency check to -100:
            # Should decrease once freq projector functions properly
            if imag_freq:
                ioprinter.warning_message(f'Imaginary mode found: {imag_freq}')
                norm_coords = elstruct.reader.normal_coordinates(prog, out_str)
                has_imag = True

    return has_imag, norm_coords


def _kickoff_saddle(geo, norm_coords, spc_info, mod_thy_info,
                    run_fs, opt_script_str,
                    kickoff_size=0.1, kickoff_backward=False, kickoff_mode=0,
                    opt_cart=True, **kwargs):
    """ kickoff from saddle to find connected minima
    """

    # Choose the displacement xyzs and set kickoff direction and size
    kickoff = kickoff_size if not kickoff_backward else -1.0*kickoff_size
    disp_xyzs = norm_coords[kickoff_mode]
    disp_xyzs = numpy.multiply(disp_xyzs, kickoff)
    
    disp_geo = automol.geom.translate_along_matrix(geo, disp_xyzs)

    ioprinter.debug_message(
        f'Creating displacement geometry with kickoff size of {kickoff}'
        f'\ninitial geometry:\n{automol.geom.string(geo)}'
        f'\ndisplaced geometry:\n{automol.geom.string(disp_geo)}'
    )

    # Optimize displaced geometry
    geom = disp_geo if opt_cart else automol.geom.zmatrix(disp_geo)
    success, ret = es_runner.execute_job(
        job=elstruct.Job.OPTIMIZATION,
        script_str=opt_script_str,
        run_fs=run_fs,
        geo=geom,
        spc_info=spc_info,
        thy_info=mod_thy_info,
        overwrite=True,
        **kwargs,
    )

    if success:
        inf_obj, _, out_str = ret
        prog = inf_obj.prog
        ret_geo = elstruct.reader.opt_geometry(prog, out_str)
    else:
        ret_geo = None

    return ret_geo, ret
