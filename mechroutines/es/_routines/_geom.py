""" es_runners for initial geometry optimization
"""

import numpy
import automol
import elstruct
import autorun
from mechanalyzer.inf import thy as tinfo
from mechlib.amech_io import printer as ioprinter
from mechroutines.es import runner as es_runner
from mechroutines.es.runner import qchem_params


def remove_imag(geo, ini_ret, spc_info, method_dct, run_fs,
                kickoff_size=0.1, kickoff_backward=False, kickoff_mode=0):
    """ if there is an imaginary frequency displace geometry along the imaginary
    mode and then reoptimize
    """

    thy_info = tinfo.from_dct(method_dct)
    mod_thy_info = tinfo.modify_orb_label(thy_info, spc_info)
    hess_script_str, kwargs = qchem_params(method_dct)

    ioprinter.info_message(
        'Checking the initial geometry for imaginary frequencies...')
    hess, hess_ret = _hess(spc_info, mod_thy_info, geo, run_fs,
                           hess_script_str, **kwargs)
    imags, norm_coords = _check_imaginary(geo, hess, hess_ret, run_fs)

    # Make five attempts to remove imag mode if found
    count, opt_ret = 1, None
    while imags and count <= 4:

        # Either attempt a kickoff-reopt or tight-reopt
        # Only attempt kickoff-reopt on first two tries
        if count <= 2:
            ioprinter.info_message(
                f'Attempting kick off along mode, attempt {count}...')
            opt_script_str, opt_kwargs = qchem_params(
                method_dct, job=elstruct.Job.OPTIMIZATION)
            disp_geo = _kickoff_saddle(
                geo, norm_coords,
                size=kickoff_size, backward=kickoff_backward,
                mode=kickoff_mode)
            geo, opt_ret = _opt(spc_info, mod_thy_info, disp_geo, run_fs,
                                opt_script_str, opt_cart=True, **opt_kwargs)
            hess_script_str, hess_kwargs = qchem_params(
                method_dct)
        else:
            ioprinter.info_message(
                f'Attempting tight opt, attempt {count-2}...')
            opt_script_str, opt_kwargs = qchem_params(
                method_dct, job='tightopt')
            geo, opt_ret = _opt(spc_info, mod_thy_info, geo, run_fs,
                                opt_script_str, opt_cart=True, **opt_kwargs)
            hess_script_str, hess_kwargs = qchem_params(
                method_dct, job='tightfreq')

        # Assess the imaginary mode after the reoptimization
        ioprinter.info_message('Rerunning Hessian...')
        hess, hess_ret = _hess(spc_info, mod_thy_info, geo, run_fs,
                               hess_script_str, **hess_kwargs)
        imags, norm_coords = _check_imaginary(geo, hess, hess_ret, run_fs)

        # Update the loop counter and kickoff size
        count += 1
        kickoff_size *= 2.0

    if opt_ret is None:
        ret = ini_ret
    else:
        ret = opt_ret

    return geo, ret


# Utility functions
def _kickoff_saddle(geo, norm_coords, size=0.1, backward=False, mode=0):
    """ kickoff from saddle to find connected minima
    """

    # Choose the displacement xyzs and set kickoff direction and size
    kickoff = size if not backward else -1.0*size
    disp_xyzs = norm_coords[mode]
    disp_xyzs = numpy.multiply(disp_xyzs, kickoff)

    disp_geo = automol.geom.translate_along_matrix(geo, disp_xyzs)

    ioprinter.debug_message(
        f'Creating displacement geometry with kickoff size of {kickoff}'
        f'\ninitial geometry:\n{automol.geom.string(geo)}'
        f'\ndisplaced geometry:\n{automol.geom.string(disp_geo)}'
    )

    return disp_geo


def _check_imaginary(geo, hess, hess_ret, run_fs):
    """ Assess the imaginary modes to decide whether to kick or not
    """
    _, _, imag_freq, _ = autorun.projrot.frequencies(
       autorun.SCRIPT_DCT['projrot'],
       run_fs[-1].path([elstruct.Job.HESSIAN]),
       [geo], [[]], [hess])

    # Mode for now set the imaginary frequency check to -100:
    # Should decrease once freq projector functions properly
    if imag_freq:
        _, norm_coords = autorun.projrot.displacements(
           autorun.SCRIPT_DCT['projrot'],
           run_fs[-1].path([elstruct.Job.HESSIAN]),
           [geo], [[]], [hess])
    else:
        norm_coords = None

    return imag_freq, norm_coords


# Simple runner functions that should be replaced at some point
def _opt(spc_info, mod_thy_info, geo, run_fs,
         script_str, opt_cart=True, **kwargs):
    """ Run an optimization
    """

    # Optimize displaced geometry
    geom = geo if opt_cart else automol.geom.zmatrix(geo)
    success, ret = es_runner.execute_job(
        job=elstruct.Job.OPTIMIZATION,
        script_str=script_str,
        run_fs=run_fs,
        geo=geom,
        spc_info=spc_info,
        thy_info=mod_thy_info,
        zrxn=None,
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


def _hess(spc_info, mod_thy_info, geo, run_fs,
          script_str, **kwargs):
    """ Calculate the Hessian
    """

    # Run Hessian calculation
    success, ret = es_runner.execute_job(
        job=elstruct.Job.HESSIAN,
        spc_info=spc_info,
        thy_info=mod_thy_info,
        geo=geo,
        zrxn=None,
        run_fs=run_fs,
        script_str=script_str,
        overwrite=True,
        **kwargs,
        )

    # Read the Hessian
    if success:
        inf_obj, _, out_str = ret
        prog = inf_obj.prog
        ret_hess = elstruct.reader.hessian(prog, out_str)
    else:
        ret_hess = None

    return ret_hess, ret
