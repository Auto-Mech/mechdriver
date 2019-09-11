""" drivers for initial geometry optimization
"""
import numpy
from qcelemental import constants as qcc
import automol
import elstruct
import autofile
import moldr

WAVEN2KCAL = qcc.conversion_factor('wavenumber', 'kcal/mol')
EH2KCAL = qcc.conversion_factor('hartree', 'kcal/mol')

def run_initial_geometry_opt(
        spc_info, thy_level, run_fs, thy_run_fs, thy_save_fs, 
        script_str, overwrite, geo_init, return_msg = False, **kwargs):
    """ generate initial geometry via optimization from either reference
    geometries or from inchi
    """
    msg = ''
    # set up the filesystem
    print(thy_run_fs)
    thy_run_fs.leaf.create(thy_level[1:4])
    thy_run_path = thy_run_fs.leaf.path(thy_level[1:4])

    thy_save_fs.leaf.create(thy_level[1:4])
    thy_save_path = thy_save_fs.leaf.path(thy_level[1:4])

    # check if geometry has already been saved
    if thy_save_fs.leaf.file.geometry.exists(thy_level[1:4]):
        geo = thy_save_fs.leaf.file.geometry.read(thy_level[1:4])
    # or if geometry has already been run
    elif thy_run_fs.leaf.file.geometry.exists(thy_level[1:4]):
        geo = thy_run_fs.leaf.file.geometry.read(thy_level[1:4])
        zma = automol.geom.zmatrix(geo)
        thy_save_fs.leaf.file.geometry.write(geo, thy_level[1:4])
        thy_save_fs.leaf.file.zmatrix.write(zma, thy_level[1:4])
    else:
        # if not call the electronic structure optimizer
        zma = automol.geom.zmatrix(geo_init)
        run_fs = autofile.fs.run(thy_run_path)
        moldr.driver.run_job(
            job=elstruct.Job.OPTIMIZATION,
            script_str=script_str,
            run_fs=run_fs,
            geom=zma,
            spc_info=spc_info,
            thy_level=thy_level,
            overwrite=overwrite,
            **kwargs,
        )
        ret = moldr.driver.read_job(job=elstruct.Job.OPTIMIZATION, run_fs=run_fs)
        geo = None
        if ret:
            msg = 'Saving reference geometry'
            msg +="\n - Save path: {}".format(thy_save_path)
            inf_obj, _, out_str = ret
            prog = inf_obj.prog
            geo = elstruct.reader.opt_geometry(prog, out_str)
            zma = automol.geom.zmatrix(geo)
            thy_save_fs.leaf.file.geometry.write(geo, thy_level[1:4])
            thy_save_fs.leaf.file.zmatrix.write(zma, thy_level[1:4])
    if return_msg:
        return geo, msg
    print(msg)
    return geo


def run_check_imaginary(
        spc_info, thy_level, thy_run_fs, thy_save_fs, script_str,
        overwrite, **kwargs):
    """ check if species has an imaginary frequency
    """
    thy_run_fs.leaf.create(thy_level[1:4])
    thy_run_path = thy_run_fs.leaf.path(thy_level[1:4])

    thy_save_fs.leaf.create(thy_level[1:4])
    thy_save_path = thy_save_fs.leaf.path(thy_level[1:4])

    run_fs = autofile.fs.run(thy_run_path)

    print('thy_run_path in check imag:', thy_run_path)
    print('thy_save_path in check imag:', thy_save_path)
    print('thy_level:', thy_level[1:4])

    ret = moldr.driver.read_job(job=elstruct.Job.OPTIMIZATION, run_fs=run_fs)
    geo = None
    imag = False
    disp_xyzs = []
    if ret:
        print('Checking for saddle')
        inf_obj, _, out_str = ret
        prog = inf_obj.prog
        geo = elstruct.reader.opt_geometry(prog, out_str)

        if thy_save_fs.leaf.file.hessian.exists(thy_level[1:4]):
            hess = thy_save_fs.leaf.file.hessian.read(thy_level[1:4])
        else:
            hess = None
        ret = moldr.driver.read_job(job=elstruct.Job.HESSIAN, run_fs=run_fs)
        if not hess:
            if not ret:
                moldr.driver.run_job(
                    job=elstruct.Job.HESSIAN,
                    spc_info=spc_info,
                    thy_level=thy_level,
                    geom=geo,
                    run_fs=run_fs,
                    script_str=script_str,
                    overwrite=overwrite,
                    **kwargs,
                    )
                ret = moldr.driver.read_job(job=elstruct.Job.HESSIAN, run_fs=run_fs)
            if ret:
                inf_obj, _, out_str = ret
                prog = inf_obj.prog
                hess = elstruct.reader.hessian(prog, out_str)

        if hess:
            if len(geo) > 1:
                freqs = elstruct.util.harmonic_frequencies(geo, hess, project=False)
#                freqs = elstruct.util.harmonic_frequencies(geo, hess, project=True)
#            print('projected freqs')
#            print(freqs)

# mode for now set the imaginary frequency check to -100:
# Ultimately should decrease once frequency projector is functioning properly
                if freqs[0] < -100:
                    imag = True
                    print('Imaginary mode found:')
                    norm_coos = elstruct.util.normal_coordinates(
                        geo, hess, project=True)
                    im_norm_coo = numpy.array(norm_coos)[:, 0]
                    disp_xyzs = numpy.reshape(im_norm_coo, (-1, 3))
    return imag, geo, disp_xyzs


def run_kickoff_saddle(
        geo, disp_xyzs, spc_info, thy_level, run_fs, thy_run_fs,
        script_str, kickoff_size=0.1, kickoff_backward=False,
        opt_cart=True, **kwargs):
    """ kickoff from saddle to find connected minima
    """
    print('kickoff from saddle')
    thy_run_fs.leaf.create(thy_level[1:4])
    thy_run_path = thy_run_fs.leaf.path(thy_level[1:4])
    run_fs = autofile.fs.run(thy_run_path)
    disp_len = kickoff_size * qcc.conversion_factor('angstrom', 'bohr')
    if kickoff_backward:
        disp_len *= -1
    disp_xyzs = numpy.multiply(disp_xyzs, disp_len)
    geo = automol.geom.displaced(geo, disp_xyzs)
    if opt_cart:
        geom = geo
    else:
        geom = automol.geom.zmatrix(geo)
    moldr.driver.run_job(
        job=elstruct.Job.OPTIMIZATION,
        script_str=script_str,
        run_fs=run_fs,
        geom=geom,
        spc_info=spc_info,
        thy_level=thy_level,
        overwrite=True,
        **kwargs,
    )


def save_initial_geometry(
        spc_info, thy_level, run_fs, thy_run_fs, thy_save_fs):
    """ save the geometry from the initial optimization as a reference geometry
    """
    thy_run_fs.leaf.create(thy_level[1:4])
    thy_run_path = thy_run_fs.leaf.path(thy_level[1:4])
    run_fs = autofile.fs.run(thy_run_path)

    thy_save_fs.leaf.create(thy_level[1:4])
    thy_save_path = thy_save_fs.leaf.path(thy_level[1:4])

    ret = moldr.driver.read_job(job=elstruct.Job.OPTIMIZATION, run_fs=run_fs)
    if ret:
        print('Saving reference geometry')
        print(" - Save path: {}".format(thy_save_path))

        inf_obj, _, out_str = ret
        prog = inf_obj.prog
        geo = elstruct.reader.opt_geometry(prog, out_str)
        zma = automol.geom.zmatrix(geo)
        thy_save_fs.leaf.file.geometry.write(geo, thy_level[1:4])
        thy_save_fs.leaf.file.zmatrix.write(zma, thy_level[1:4])


