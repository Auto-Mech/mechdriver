""" es_runners for initial geometry optimization
"""

import numpy
import automol
import elstruct
import autorun
from phydat import phycon
from mechanalyzer.inf import thy as tinfo
from mechroutines.es import runner as es_runner
from mechroutines.es.runner import qchem_params
from mechlib import structure
from mechlib.amech_io import printer as ioprinter


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
    chk_idx = 0
    kick_ret = None
    while imag and chk_idx < 5:
        chk_idx += 1
        ioprinter.info_message(
            'Attempting kick off along mode, attempt {}...'.format(chk_idx))

        geo, kick_ret = _kickoff_saddle(
            geo, norm_coords, spc_info, mod_thy_info,
            run_fs, opt_script_str,
            kickoff_size, kickoff_backward, kickoff_mode,
            opt_cart=True, **opt_kwargs)

        ioprinter.info_message(
            'Removing faulty geometry from filesystem. Rerunning Hessian...')
        run_fs[-1].remove([elstruct.Job.HESSIAN])

        ioprinter.info_message('Rerunning Hessian...')
        imag, norm_coords = _check_imaginary(
            spc_info, geo, mod_thy_info, run_fs, script_str,
            overwrite, **kwargs)

        # Update kickoff size
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
    imag = False
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
            # _, _, imag_freq, _ = structure.vib.projrot_freqs(
            #     [geo], [hess], run_path)

            # Mode for now set the imaginary frequency check to -100:
            # Should decrease once freq projector functions properly
            if imag_freq:
                ioprinter.warning_message('Imaginary mode found:')
                norm_coords = elstruct.reader.normal_coordinates(prog, out_str)

    return bool(imag_freq), norm_coords


def _kickoff_saddle(geo, norm_coords, spc_info, mod_thy_info,
                    run_fs, opt_script_str,
                    kickoff_size=0.1, kickoff_backward=False, kickoff_mode=0,
                    opt_cart=True, **kwargs):
    """ kickoff from saddle to find connected minima
    """

    # Choose the displacement xyzs
    disp_xyzs = norm_coords[kickoff_mode]

    # Set the displacement vectors and displace geometry
    disp_len = kickoff_size * phycon.ANG2BOHR
    if kickoff_backward:
        disp_len *= -1
    disp_xyzs = numpy.multiply(disp_xyzs, disp_len)
    ioprinter.debug_message(
        'geo test in kickoff_saddle:',
        automol.geom.string(geo), disp_xyzs)

    geo = automol.geom.displace(geo, disp_xyzs)

    # Optimize displaced geometry
    if opt_cart:
        geom = geo
    else:
        geom = automol.geom.zmatrix(geo)
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
        geo = elstruct.reader.opt_geometry(prog, out_str)

    return geo, ret
