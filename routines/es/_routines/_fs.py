"""
  Common functions for reading and writing the filesystem
"""

import automol
import elstruct
import autofile
from routines.es import runner as es_runner


def save_struct(run_fs, save_fs, locs, job, mod_thy_info,
                inzma, zma_locs=(0,), in_zma_fs=False):
    """ Save a geometry and associated information from some
        electronic structure routine into the filesystem.
    """

    ret = es_runner.read_job(job=job, run_fs=run_fs)

    if ret:

        # Get the geo, zma, and ene based on job type
        ene, geo, zma = _read(run_fs, job, inzma)

        # Obtain inf obj and inp str to write in filesys
        inf_obj, inp_str, _ = ret

        # Set and print the save path information
        save_path = save_fs[-1].path(locs)
        print(" - Saving...")
        print(" - Save path: {}".format(save_path))

        # Save the geometry information
        save_fs[-1].create(locs)
        save_fs[-1].file.geometry_info.write(inf_obj, locs)
        save_fs[-1].file.geometry_input.write(inp_str, locs)
        save_fs[-1].file.geometry.write(geo, locs)

        # Save zma information seperately, if required
        if not in_zma_fs:
            zma_save_fs = autofile.fs.manager(save_path, 'ZMATRIX')
            zma_save_fs[-1].create(zma_locs)
            zma_save_fs[-1].file.geometry_info.write(inf_obj, zma_locs)
            zma_save_fs[-1].file.geometry_input.write(inp_str, zma_locs)
            zma_save_fs[-1].file.zmatrix.write(zma, zma_locs)
        else:
            save_fs[-1].file.zmatrix.write(zma, locs)

        # Saving the energy to an SP filesys
        print(" - Saving energy...")
        sp_save_fs = autofile.fs.single_point(save_path)
        sp_save_fs[-1].create(mod_thy_info[1:4])
        sp_save_fs[-1].file.input.write(inp_str, mod_thy_info[1:4])
        sp_save_fs[-1].file.info.write(inf_obj, mod_thy_info[1:4])
        sp_save_fs[-1].file.energy.write(ene, mod_thy_info[1:4])

        saved = True

    else:
        saved = False

    return saved


def _read(run_fs, job, inzma):
    """ Read the output
    """

    assert job in (elstruct.Job.OPTIMIZATION, elstruct.Job.ENERGY), (
        'Job not opt or energy'
    )
    ret = es_runner.read_job(
        job=job, run_fs=run_fs)
    if ret:
        inf_obj, _, out_str = ret
        prog = inf_obj.prog
        method = inf_obj.method
        ene = elstruct.reader.energy(prog, method, out_str)
        if job == elstruct.Job.OPTIMIZATION:
            geo = elstruct.reader.opt_geometry(prog, out_str)
            zma = elstruct.reader.opt_zmatrix(prog, out_str)
        elif job == elstruct.Job.ENERGY:
            zma = inzma
            geo = automol.zmatrix.geometry(zma)
    else:
        ene, zma, geo = None, None, None

    return ene, zma, geo


def _hr_traj(coord_names, scn_save_fs, locs_lst):
    """
    """
    idxs_lst = [locs[-1] for locs in locs_lst]
    enes = [scn_save_fs[-1].file.energy.read(locs)
            for locs in locs_lst]
    geos = [scn_save_fs[-1].file.geometry.read(locs)
            for locs in locs_lst]

    traj = []
    for idxs, ene, geo in zip(idxs_lst, enes, geos):
        comment = (
            'energy: {:>15.10f}, '.format(ene) +
            'grid idxs: {}'.format(idxs)
        )
        traj.append((comment, geo))

    traj_path = scn_save_fs[1].file.trajectory.path([coord_names])
    print("Updating scan trajectory file at {}".format(traj_path))
    scn_save_fs[1].file.trajectory.write(traj, [coord_names])


def _check_vma(zma, save_fs):
    """ Ensure that the ZMA is being used
    """
    vma = automol.zmatrix.var_(zma)
    if save_fs[0].file.vmatrix.exists():
        existing_vma = save_fs[0].file.vmatrix.read()
        assert vma == existing_vma
