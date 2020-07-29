""" es_runners for initial geometry optimization
"""

import numpy
import automol
import elstruct
import autofile
from routines.es import runner as es_runner
from routines.es._routines._fs import save_struct
from routines.es._routines._fs import save_instab
from lib import structure
from lib.phydat import phycon


def reference_geometry(spc_dct_i, spc_info,
                       mod_thy_info,
                       thy_run_fs, thy_save_fs,
                       cnf_run_fs, cnf_save_fs,
                       run_fs,
                       opt_script_str, overwrite,
                       kickoff_size=0.1, kickoff_backward=False,
                       **opt_kwargs):
    """ determine what to use as the reference geometry for all future runs
    If ini_thy_info refers to geometry dictionary then use that,
    otherwise values are from a hierarchy of:
    running level of theory, input level of theory, inchis.
    From the hierarchy an optimization is performed followed by a check for
    an imaginary frequency and then a conformer file system is set up.
    """

    # Initialize empty return
    ret = None

    if run_fs[0].file.info.exists([]):
        inf_obj = run_fs[0].file.info.read([])
        if inf_obj.status == autofile.schema.RunStatus.RUNNING:
            print('Reference geometry already running')
            return ret
    else:
        [prog, method, basis, _] = mod_thy_info
        status = autofile.schema.RunStatus.RUNNING
        inf_obj = autofile.schema.info_objects.run(
            job='', prog=prog, version='version', method=method, basis=basis,
            status=status)
        run_fs[0].file.info.write(inf_obj, [])

    if not thy_save_fs[-1].file.geometry.exists(mod_thy_info[1:4]):
        print('No geometry in filesys. Attempting to initialize new geometry.')
        geo_init = _obtain_ini_geom(spc_dct_i)
        if geo_init is not None:
            zma_init = automol.geom.zmatrix(geo_init)
            if not automol.geom.is_atom(geo_init):
                geo_found = _optimize_molecule(
                    spc_info, zma_init,
                    mod_thy_info, thy_run_fs, thy_save_fs,
                    cnf_save_fs,
                    run_fs,
                    opt_script_str, overwrite,
                    kickoff_size=kickoff_size,
                    kickoff_backward=kickoff_backward,
                    **opt_kwargs)
            else:
                geo_found = _optimize_atom(
                    spc_info, zma_init,
                    mod_thy_info, thy_run_fs,
                    cnf_save_fs, run_fs,
                    overwrite, opt_script_str, **opt_kwargs)
        else:
            geo_found = False
            print('Unable to obtain an initial guess geometry')
    else:
        geo_found = True
        thy_path = thy_save_fs[-1].path(mod_thy_info[1:4])
        print('Initial geometry found and saved previously at {}'.format(
            thy_path))

    # Write the job status into the run filesystem
    if geo_found:
        inf_obj.status = autofile.schema.RunStatus.SUCCESS
        run_fs[0].file.info.write(inf_obj, [])
    else:
        inf_obj.status = autofile.schema.RunStatus.FAILURE
        run_fs[0].file.info.write(inf_obj, [])

    return geo_found


def _obtain_ini_geom(spc_dct_i):
    """ Obtain an initial geometry to be optimized. Checks a hieratchy
        of places to obtain the initial geom.
            (1) Geom dict which is the input from the user
            (2) Geom from inchi
    """

    geo_init = None

    # Obtain geom from user input geometry or inchi
    if 'geo_obj' in spc_dct_i:
        geo_init = spc_dct_i['geo_obj']
        print('Getting geometry from geom dictionary')

    if geo_init is None:
        geo_init = automol.inchi.geometry(spc_dct_i['inchi'])
        print('Getting reference geometry from inchi')

    # Check if the init geometry is connected
    if geo_init is not None:
        if not automol.geom.connected(geo_init):
            geo_init = None

    return geo_init


def _optimize_atom(spc_info, zma_init,
                   mod_thy_info, thy_run_fs,
                   cnf_save_fs, run_fs,
                   overwrite, opt_script_str, **opt_kwargs):
    """ Deal with an atom separately
    """

    geo, zma = _init_geom_opt(zma_init, spc_info, mod_thy_info,
                              run_fs, thy_run_fs,
                              opt_script_str, overwrite, **opt_kwargs)

    if geo is not None and zma is not None:
        locs = [autofile.schema.generate_new_conformer_id()]
        job = elstruct.Job.OPTIMIZATION
        save_struct(run_fs, cnf_save_fs, locs, job, mod_thy_info,
                    zma_locs=(0,), in_zma_fs=False)
        conf_found = True
    else:
        conf_found = False

    return conf_found


def _optimize_molecule(spc_info, zma_init,
                       mod_thy_info, thy_run_fs, thy_save_fs,
                       cnf_save_fs,
                       run_fs,
                       opt_script_str, overwrite,
                       kickoff_size=0.1, kickoff_backward=False,
                       **opt_kwargs):
    """ Optimize a proper geometry
    """

    # Optimize the initial geometry
    geo, zma = _init_geom_opt(
        zma_init, spc_info, mod_thy_info, run_fs, thy_run_fs,
        opt_script_str, overwrite, **opt_kwargs)

    # If connected, check for imaginary modes and fix them if possible
    if automol.geom.connected(geo):

        # Remove the imaginary mode
        geo, imag_fix_needed = _remove_imag(
            spc_info, geo, mod_thy_info, thy_run_fs,
            run_fs, kickoff_size, kickoff_backward,
            overwrite=overwrite)

        # Recheck connectivity for imag-checked geometry
        if geo is not None:

            conf_found = True
            if automol.geom.connected(geo):

                print('\nSaving structure as the first conformer...')
                locs = [autofile.schema.generate_new_conformer_id()]
                job = elstruct.Job.OPTIMIZATION
                save_struct(run_fs, cnf_save_fs, locs, job, mod_thy_info,
                            zma_locs=(0,), in_zma_fs=False,
                            cart_to_zma=imag_fix_needed)

            else:

                print('Saving disconnected species...')
                locs = [autofile.schema.generate_new_conformer_id()]
                job = elstruct.Job.OPTIMIZATION
                save_instab(zma_init, run_fs, cnf_save_fs, locs,
                            mod_thy_info)
                structure.instab.write_instab(
                    zma_init, zma, thy_save_fs, mod_thy_info[1:4])

        else:

            print('\n No geom found...')
            conf_found = False

    else:

        print('Saving disconnected species...')
        conf_found = False
        locs = [autofile.schema.generate_new_conformer_id()]
        job = elstruct.Job.OPTIMIZATION
        save_instab(zma_init, run_fs, cnf_save_fs, locs,
                    mod_thy_info)
        structure.instab.write_instab(
            zma_init, zma, thy_save_fs, mod_thy_info[1:4])

    # Save geom in thy filesys if a good geom is found
    if conf_found:
        thy_save_fs[-1].create(mod_thy_info[1:4])
        thy_save_path = thy_save_fs[-1].path(mod_thy_info[1:4])
        thy_save_fs[-1].file.geometry.write(geo, mod_thy_info[1:4])

        print('Saving reference geometry')
        print(" - Save path: {}".format(thy_save_path))

    return conf_found


def _init_geom_opt(zma_init, spc_info, mod_thy_info,
                   run_fs, thy_run_fs,
                   opt_script_str, overwrite, **opt_kwargs):
    """ Generate initial geometry via optimization from either reference
    geometries or from inchi
    """

    # Set up the filesystem
    thy_run_fs[-1].create(mod_thy_info[1:4])
    thy_run_path = thy_run_fs[-1].path(mod_thy_info[1:4])

    # Call the electronic structure optimizer
    run_fs = autofile.fs.run(thy_run_path)
    es_runner.run_job(
        job=elstruct.Job.OPTIMIZATION,
        script_str=opt_script_str,
        run_fs=run_fs,
        geom=zma_init,
        spc_info=spc_info,
        thy_info=mod_thy_info,
        overwrite=overwrite,
        **opt_kwargs,
    )
    success, ret = es_runner.read_job(
        job=elstruct.Job.OPTIMIZATION, run_fs=run_fs)

    geo, zma = None, None
    if success:
        print('Succesful reference geometry optimization')
        inf_obj, _, out_str = ret
        prog = inf_obj.prog
        geo = elstruct.reader.opt_geometry(prog, out_str)
        zma = elstruct.reader.opt_zmatrix(prog, out_str)

    return geo, zma


# Functions to remove imaginary mode from a structure
def _remove_imag(spc_info, geo, mod_thy_info, thy_run_fs, run_fs,
                 kickoff_size=0.1, kickoff_backward=False,
                 overwrite=False):
    """ if there is an imaginary frequency displace geometry along the imaginary
    mode and then reoptimize
    """

    print('The initial geometries will be checked for imaginary frequencies')
    script_str, opt_script_str, kwargs, opt_kwargs = es_runner.qchem_params(
        *mod_thy_info[0:2])

    imag, disp_xyzs = _check_imaginary(
        spc_info, geo, mod_thy_info, thy_run_fs, script_str,
        overwrite, **kwargs)

    # Make a variable to fix the imaginary mode if needed to pass to other functions
    imag_fix_needed = bool(imag)

    # Make five attempts to remove imag mode if found
    chk_idx = 0
    while imag and chk_idx < 5:
        chk_idx += 1
        print('Attempting kick off along mode, attempt {}...'.format(chk_idx))

        geo = _kickoff_saddle(
            geo, disp_xyzs, spc_info, mod_thy_info, run_fs, thy_run_fs,
            opt_script_str, kickoff_size, kickoff_backward,
            opt_cart=True, **opt_kwargs)

        print('Removing faulty geometry from filesystem. Rerunning Hessian...')
        thy_run_path = thy_run_fs[-1].path(mod_thy_info[1:4])
        run_fs = autofile.fs.run(thy_run_path)
        run_fs[-1].remove([elstruct.Job.HESSIAN])

        print('Rerunning Hessian...')
        imag, disp_xyzs = _check_imaginary(
            spc_info, geo, mod_thy_info, thy_run_fs, script_str,
            overwrite, **kwargs)

        # Update kickoff size
        kickoff_size *= 2

    return geo, imag_fix_needed


def _check_imaginary(
        spc_info, geo, mod_thy_info, thy_run_fs, script_str,
        overwrite=False, **kwargs):
    """ check if species has an imaginary frequency
    """

    # Handle filesystem
    thy_run_fs[-1].create(mod_thy_info[1:4])
    thy_run_path = thy_run_fs[-1].path(mod_thy_info[1:4])
    run_fs = autofile.fs.run(thy_run_path)

    # Initialize info
    imag = False
    disp_xyzs = []
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
                print('Imaginary mode found:')
                # norm_coos = elstruct.util.normal_coordinates(
                #     geo, hess, project=True)
                # im_norm_coo = numpy.array(norm_coos)[:, 0]
                # disp_xyzs = numpy.reshape(im_norm_coo, (-1, 3))
                norm_coos = elstruct.reader.normal_coords(prog, out_str)
                im_norm_coo = norm_coos[0]
                disp_xyzs = im_norm_coo

    return imag, disp_xyzs


def _kickoff_saddle(
        geo, disp_xyzs, spc_info, mod_thy_info, run_fs, thy_run_fs,
        opt_script_str, kickoff_size=0.1, kickoff_backward=False,
        opt_cart=True, **kwargs):
    """ kickoff from saddle to find connected minima
    """

    # Set the filesys
    thy_run_fs[-1].create(mod_thy_info[1:4])
    thy_run_path = thy_run_fs[-1].path(mod_thy_info[1:4])
    run_fs = autofile.fs.run(thy_run_path)

    # Set the displacement vectors and displace geometry
    disp_len = kickoff_size * phycon.ANG2BOHR
    if kickoff_backward:
        disp_len *= -1
    disp_xyzs = numpy.multiply(disp_xyzs, disp_len)

    print('geo test in kickoff_saddle:')
    print(automol.geom.string(geo), disp_xyzs)
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


# Old fake conf function that may not be useful
def fake_conf(mod_thy_info, filesystem, inf=()):
    """ generate data to be used for a fake well I think?
    """
    cnf_save_fs = filesystem[5]
    cnf_run_fs = filesystem[4]
    thy_save_fs = filesystem[3]
    run_fs = filesystem[-1]
    thy_save_path = thy_save_fs[-1].path(mod_thy_info[1:4])
    geo = thy_save_fs[-1].file.geometry.read(mod_thy_info[1:4])
    if inf:
        inf_obj, ene = inf
    else:
        ene = thy_save_fs[-1].file.energy.read(mod_thy_info[1:4])
        inf_obj = run_fs[0].file.info.read()
    tors_range_dct = {}
    cinf_obj = autofile.schema.info.conformer_trunk(0, tors_range_dct)
    cinf_obj.nsamp = 1
    cnf_save_fs = autofile.fs.conformer(thy_save_path)
    cnf_save_fs[0].create()
    cnf_run_fs[0].create()
    cnf_save_fs[0].file.info.write(cinf_obj)
    cnf_run_fs[0].file.info.write(cinf_obj)
    locs_lst = cnf_save_fs[-1].existing()
    if not locs_lst:
        cid = autofile.schema.generate_new_conformer_id()
        locs = [cid]
    else:
        locs = locs_lst[0]
    cnf_save_fs[-1].create(locs)
    cnf_run_fs[-1].create(locs)
    cnf_save_fs[-1].file.geometry_info.write(
        inf_obj, locs)
    cnf_run_fs[-1].file.geometry_info.write(
        inf_obj, locs)
    # method = inf_obj.method
    cnf_save_fs[-1].file.energy.write(ene, locs)
    cnf_run_fs[-1].file.energy.write(ene, locs)
    cnf_save_fs[-1].file.geometry.write(geo, locs)
    cnf_run_fs[-1].file.geometry.write(geo, locs)
# old way to check torsions of new geom after imag check
# # Check the torsons to see if ref matches (needed?)
# tors_names = automol.geom.zmatrix_torsion_coordinate_names(geo)
# locs_lst = cnf_save_fs[-1].existing()
# if locs_lst:
#     saved_geo = cnf_save_fs[-1].file.geometry.read(
#         locs_lst[0])
#     saved_tors = automol.geom.zmatrix_torsion_coordinate_names(
#         saved_geo)
#     if tors_names != saved_tors:
#         print("New reference geometry doesn't match original",
#               " reference geometry")
#         print('Removing original conformer save data')
#         cnf_run_fs.remove()
#         cnf_save_fs.remove()
