""" es_runners for initial geometry optimization
"""

import numpy
import automol
import elstruct
import autofile
from phydat import phycon
from mechroutines.es import runner as es_runner
from mechlib import structure
from mechlib.submission import qchem_params
from mechlib.amech_io import printer as ioprinter


def remove_imag(geo, spc_info, mod_thy_info, thy_run_fs, run_fs,
                kickoff_size=0.1, kickoff_backward=False, kickoff_mode=0,
                overwrite=False):
    """ if there is an imaginary frequency displace geometry along the imaginary
    mode and then reoptimize
    """

    ioprinter.info_message(
        'The initial geometries will be checked for imaginary frequencies')
    script_str, opt_script_str, kwargs, opt_kwargs = qchem_params(
        *mod_thy_info[0:2])

    imag, norm_coords = _check_imaginary(
        spc_info, geo, mod_thy_info, thy_run_fs, script_str,
        overwrite, **kwargs)

    # Make var to fix the imaginary mode if needed to pass to other functions
    imag_fix_needed = bool(imag)

    # Make five attempts to remove imag mode if found
    chk_idx = 0
    while imag and chk_idx < 5:
        chk_idx += 1
        ioprinter.info_message(
            'Attempting kick off along mode, attempt {}...'.format(chk_idx))

        geo = _kickoff_saddle(
            geo, norm_coords, spc_info, mod_thy_info,
            run_fs, thy_run_fs, opt_script_str,
            kickoff_size, kickoff_backward, kickoff_mode,
            opt_cart=True, **opt_kwargs)

        ioprinter.info_message(
            'Removing faulty geometry from filesystem. Rerunning Hessian...')
        thy_run_path = thy_run_fs[-1].path(mod_thy_info[1:4])
        run_fs = autofile.fs.run(thy_run_path)
        run_fs[-1].remove([elstruct.Job.HESSIAN])

        ioprinter.info_message('Rerunning Hessian...')
        imag, norm_coords = _check_imaginary(
            spc_info, geo, mod_thy_info, thy_run_fs, script_str,
            overwrite, **kwargs)

        # Update kickoff size
        kickoff_size *= 2.0

    return geo, imag_fix_needed


def _check_imaginary(spc_info, geo, mod_thy_info, thy_run_fs, script_str,
                     overwrite=False, **kwargs):
    """ check if species has an imaginary frequency
    """

    # Handle filesystem
    thy_run_fs[-1].create(mod_thy_info[1:4])
    thy_run_path = thy_run_fs[-1].path(mod_thy_info[1:4])
    run_fs = autofile.fs.run(thy_run_path)

    # Initialize info
    imag = False
    norm_coords = []
    hess = ((), ())

    # Run Hessian calculation
    es_runner.run_job(
        job=elstruct.Job.HESSIAN,
        spc_info=spc_info,
        thy_info=mod_thy_info,
        geom=geo,
        run_fs=run_fs,
        script_str=script_str,
        overwrite=overwrite,
        **kwargs,
        )

    # Check for imaginary modes
    success, ret = es_runner.read_job(job=elstruct.Job.HESSIAN, run_fs=run_fs)
    if success:
        inf_obj, _, out_str = ret
        prog = inf_obj.prog
        hess = elstruct.reader.hessian(prog, out_str)

        # Calculate vibrational frequencies
        if hess:
            imag = False
            _, _, imag_freq, _ = structure.vib.projrot_freqs(
                [geo], [hess], thy_run_path)
            if imag_freq:
                imag = True

            # Mode for now set the imaginary frequency check to -100:
            # Should decrease once freq projector functions properly
            if imag:
                imag = True
                ioprinter.warning_message('Imaginary mode found:')
                norm_coords = elstruct.reader.normal_coordinates(prog, out_str)

    return imag, norm_coords


def _kickoff_saddle(geo, norm_coords, spc_info, mod_thy_info,
                    run_fs, thy_run_fs, opt_script_str,
                    kickoff_size=0.1, kickoff_backward=False, kickoff_mode=0,
                    opt_cart=True, **kwargs):
    """ kickoff from saddle to find connected minima
    """

    # Set the filesys
    thy_run_fs[-1].create(mod_thy_info[1:4])
    thy_run_path = thy_run_fs[-1].path(mod_thy_info[1:4])
    run_fs = autofile.fs.run(thy_run_path)

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
    es_runner.run_job(
        job=elstruct.Job.OPTIMIZATION,
        script_str=opt_script_str,
        run_fs=run_fs,
        geom=geom,
        spc_info=spc_info,
        thy_info=mod_thy_info,
        overwrite=True,
        **kwargs,
    )

    success, ret = es_runner.read_job(
        job=elstruct.Job.OPTIMIZATION, run_fs=run_fs)
    if success:
        inf_obj, _, out_str = ret
        prog = inf_obj.prog
        geo = elstruct.reader.opt_geometry(prog, out_str)

    return geo
