""" es_runners for initial geometry optimization
"""

import numpy
import automol
import elstruct
import autofile
from routines.es._routines import conformer
from routines.es import runner as es_runner
from lib import filesys
from lib import structure
from lib.phydat import phycon


def reference_geometry(
        spc_dct_i, thy_info, ini_thy_info,
        thy_run_fs, thy_save_fs,
        ini_thy_save_fs,
        cnf_run_fs, cnf_save_fs,
        run_fs,
        kickoff_size=0.1, kickoff_backward=False, overwrite=False):
    """ determine what to use as the reference geometry for all future runs
    If ini_thy_info refers to geometry dictionary then use that,
    otherwise values are from a hierarchy of:
    running level of theory, input level of theory, inchis.
    From the hierarchy an optimization is performed followed by a check for
    an imaginary frequency and then a conformer file system is set up.
    """
    ret = None

    if run_fs[0].file.info.exists([]):
        inf_obj = run_fs[0].file.info.read([])
        if inf_obj.status == autofile.schema.RunStatus.RUNNING:
            print('Reference geometry already running')
            return ret
    else:
        [prog, method, basis, _] = thy_info
        status = autofile.schema.RunStatus.RUNNING
        inf_obj = autofile.schema.info.run(
            job='', prog=prog, version='version', method=method, basis=basis,
            status=status)
        run_fs[0].file.info.write(inf_obj, [])

    # print('initializing geometry in reference_geometry')
    geo = None
    try:
        # Check to see if geometry should be obtained from dictionary
        spc_info = [spc_dct_i['ich'], spc_dct_i['chg'], spc_dct_i['mul']]
        if 'input_geom' in ini_thy_info:
            geom_obj = spc_dct_i['geo_obj']
            geo_init = geom_obj
            overwrite = True
            print('Found initial geometry from geometry dictionary')
        else:
            # Check to see if geo already exists at running_theory
            if thy_save_fs[-1].file.geometry.exists(thy_info[1:4]):
                thy_path = thy_save_fs[-1].path(thy_info[1:4])
                print('Getting reference geometry from {}'.format(thy_path))
                geo = thy_save_fs[-1].file.geometry.read(thy_info[1:4])
            if not geo:
                if ini_thy_save_fs:
                    geo_exists = ini_thy_save_fs[-1].file.geometry.exists(
                        ini_thy_info[1:4])
                    if geo_exists:
                        # If not, Compute geo at running_theory, using geo from
                        # initial_level as the starting point
                        # or from inchi is no initial level geometry
                        thy_path = ini_thy_save_fs[-1].path(
                            ini_thy_info[1:4])
                        geo_init = ini_thy_save_fs[-1].file.geometry.read(
                            ini_thy_info[1:4])
                    elif 'geo_obj' in spc_dct_i:
                        geo_init = spc_dct_i['geo_obj']
                        print('Getting geometry from geom dictionary')
                    else:
                        # print('Getting reference geometry from inchi',
                        #       spc_info[0])
                        geo_init = automol.inchi.geometry(spc_info[0])
                        print('Got reference geometry from inchi', spc_info[0])
                elif 'geo_obj' in spc_dct_i:
                    geo_init = spc_dct_i['geo_obj']
                    print('Getting geometry from geom dictionary')
                else:
                    geo_init = automol.inchi.geometry(spc_info[0])
                    print('Getting reference geometry from inchi')
        # Optimize from initial geometry to get reference geometry
        if not geo:
            _, opt_script_str, _, opt_kwargs = es_runner.par.run_qchem_par(
                *thy_info[0:2])
            params = {
                'spc_info': spc_info,
                'run_fs': run_fs,
                'thy_run_fs': thy_run_fs,
                'script_str': opt_script_str,
                'overwrite': overwrite,
                'thy_info': thy_info,
                'ini_geo': geo_init}
            geo, inf = run_initial_geometry_opt(**params, **opt_kwargs)
            thy_save_fs[-1].create(thy_info[1:4])
            thy_save_path = thy_save_fs[-1].path(thy_info[1:4])
            ncp = len(
                automol.graph.connected_components(
                    automol.geom.graph(geo)))
            if not automol.geom.is_atom(geo) and ncp < 2:
                geo, hess = remove_imag(
                    spc_dct_i, geo, thy_info, thy_run_fs,
                    run_fs, kickoff_size,
                    kickoff_backward,
                    overwrite=overwrite)

                tors_names = automol.geom.zmatrix_torsion_coordinate_names(geo)
                locs_lst = cnf_save_fs[-1].existing()
                if locs_lst:
                    saved_geo = cnf_save_fs[-1].file.geometry.read(
                        locs_lst[0])
                    saved_tors = automol.geom.zmatrix_torsion_coordinate_names(
                        saved_geo)
                    if tors_names != saved_tors:
                        print("New reference geometry doesn't match original",
                              " reference geometry")
                        print('Removing original conformer save data')
                        cnf_run_fs.remove()
                        cnf_save_fs.remove()

                print('Saving reference geometry')
                print(" - Save path: {}".format(thy_save_path))
                thy_save_fs[-1].file.hessian.write(hess, thy_info[1:4])

            thy_save_fs[-1].file.geometry.write(geo, thy_info[1:4])
            ncp = len(
                automol.graph.connected_components(
                    automol.geom.graph(geo)))
            if ncp < 2:
                zma = automol.geom.zmatrix(geo)
                thy_save_fs[-1].file.zmatrix.write(zma, thy_info[1:4])
                print('\nObtaining a single conformer using',
                      'the MonteCarlo conformer sampling routine')
                geo_path = thy_save_fs[0].path(thy_info[1:4])
                print('Sampling done using geom from {}'.format(geo_path))
                conformer.single_conformer(
                    zma, spc_info, thy_info,
                    thy_save_fs, cnf_run_fs, cnf_save_fs,
                    overwrite, saddle=False, dist_info=())
            else:
                print("Cannot create zmatrix for disconnected species")
                fake_conf(thy_info, filesys, inf)

        if geo:
            inf_obj.status = autofile.schema.RunStatus.SUCCESS
            run_fs[0].file.info.write(inf_obj, [])
        else:
            inf_obj.status = autofile.schema.RunStatus.FAILURE
            run_fs[0].file.info.write(inf_obj, [])

    except IOError:
        inf_obj.status = autofile.schema.RunStatus.FAILURE
        run_fs[0].file.info.write(inf_obj, [])

    return geo


def run_initial_geometry_opt(
        spc_info, thy_info, run_fs, thy_run_fs,
        script_str, overwrite, ini_geo, **kwargs):
    """ generate initial geometry via optimization from either reference
    geometries or from inchi
    """

    # set up the filesystem
    thy_run_fs[-1].create(thy_info[1:4])
    thy_run_path = thy_run_fs[-1].path(thy_info[1:4])
    # check if geometry has already been saved
    # if not call the electronic structure optimizer
    ncp1 = len(automol.graph.connected_components(automol.geom.graph(ini_geo)))
    if ncp1 < 2:
        geom = automol.geom.zmatrix(ini_geo)
    else:
        geom = ini_geo
    run_fs = autofile.fs.run(thy_run_path)
    es_runner.run_job(
        job=elstruct.Job.OPTIMIZATION,
        script_str=script_str,
        run_fs=run_fs,
        geom=geom,
        spc_info=spc_info,
        thy_info=thy_info,
        overwrite=overwrite,
        **kwargs,
    )
    ret = es_runner.read_job(job=elstruct.Job.OPTIMIZATION, run_fs=run_fs)
    geo = None
    inf = None
    if ret:
        print('Succesful reference geometry optimization')
        inf_obj, _, out_str = ret
        prog = inf_obj.prog
        geo = elstruct.reader.opt_geometry(prog, out_str)
        ncp2 = len(automol.graph.connected_components(automol.geom.graph(geo)))
        if ncp2 >= 2:
            method = inf_obj.method
            ene = elstruct.reader.energy(prog, method, out_str)
            inf = [inf_obj, ene]
    return geo, inf


def remove_imag(
        spc_dct_i, geo, thy_info, thy_run_fs, run_fs, kickoff_size=0.1,
        kickoff_backward=False,
        overwrite=False):
    """ if there is an imaginary frequency displace geometry along the imaginary
    mode and then reoptimize
    """

    print('The initial geometries will be checked for imaginary frequencies')
    spc_info = filesys.inf.get_spc_info(spc_dct_i)
    runinf = es_runner.par.run_qchem_par(
        *thy_info[0:2])
    [script_str, opt_script_str, kwargs, opt_kwargs] = runinf

    imag, geo, disp_xyzs, hess = run_check_imaginary(
        spc_info, geo, thy_info, thy_run_fs, script_str,
        overwrite, **kwargs)
    chk_idx = 0
    while imag and chk_idx < 5:
        chk_idx += 1
        print('Attemptin gkick off along mode, attempt {}...'.format(chk_idx))

        geo = run_kickoff_saddle(
            geo, disp_xyzs, spc_info, thy_info, run_fs, thy_run_fs,
            opt_script_str, kickoff_size, kickoff_backward,
            opt_cart=True, **opt_kwargs)
        print('Removing faulty geometry from filesystem. Rerunning Hessian...')

        thy_run_path = thy_run_fs[-1].path(thy_info[1:4])
        run_fs = autofile.fs.run(thy_run_path)
        run_fs[-1].remove([elstruct.Job.HESSIAN])
        imag, geo, disp_xyzs, hess = run_check_imaginary(
            spc_info, geo, thy_info, thy_run_fs, script_str,
            overwrite, **kwargs)
    return geo, hess


def run_check_imaginary(
        spc_info, geo, thy_info, thy_run_fs, script_str,
        overwrite=False, **kwargs):
    """ check if species has an imaginary frequency
    """

    thy_run_fs[-1].create(thy_info[1:4])
    thy_run_path = thy_run_fs[-1].path(thy_info[1:4])

    run_fs = autofile.fs.run(thy_run_path)
    imag = False
    disp_xyzs = []
    hess = ((), ())
    if automol.geom.is_atom(geo):
        hess = ((), ())
    else:
        es_runner.run_job(
            job=elstruct.Job.HESSIAN,
            spc_info=spc_info,
            thy_info=thy_info,
            geom=geo,
            run_fs=run_fs,
            script_str=script_str,
            overwrite=overwrite,
            **kwargs,
            )
        ret = es_runner.read_job(job=elstruct.Job.HESSIAN, run_fs=run_fs)
        if ret:
            inf_obj, _, out_str = ret
            prog = inf_obj.prog
            hess = elstruct.reader.hessian(prog, out_str)

            if hess:
                imag = False
                _, _, imag_freq, _ = structure.vib.projrot_freqs(
                    geo, hess, thy_run_path)
                if imag_freq:
                    imag = True

                # Mode for now set the imaginary frequency check to -100:
                # Should decrease once freq projector functions properly
                if imag:
                    imag = True
                    print('Imaginary mode found:')
                    norm_coos = elstruct.util.normal_coordinates(
                        geo, hess, project=True)
                    im_norm_coo = numpy.array(norm_coos)[:, 0]
                    disp_xyzs = numpy.reshape(im_norm_coo, (-1, 3))
    return imag, geo, disp_xyzs, hess


def run_kickoff_saddle(
        geo, disp_xyzs, spc_info, thy_info, run_fs, thy_run_fs,
        opt_script_str, kickoff_size=0.1, kickoff_backward=False,
        opt_cart=True, **kwargs):
    """ kickoff from saddle to find connected minima
    """
    thy_run_fs[-1].create(thy_info[1:4])
    thy_run_path = thy_run_fs[-1].path(thy_info[1:4])
    run_fs = autofile.fs.run(thy_run_path)
    disp_len = kickoff_size * phycon.ANG2BOHR
    if kickoff_backward:
        disp_len *= -1
    disp_xyzs = numpy.multiply(disp_xyzs, disp_len)
    geo = automol.geom.displaced(geo, disp_xyzs)
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
        thy_info=thy_info,
        overwrite=True,
        **kwargs,
    )
    ret = es_runner.read_job(job=elstruct.Job.OPTIMIZATION, run_fs=run_fs)
    if ret:
        inf_obj, _, out_str = ret
        prog = inf_obj.prog
        geo = elstruct.reader.opt_geometry(prog, out_str)
    return geo


def save_initial_geometry(
        thy_info, run_fs, thy_run_fs, thy_save_fs):
    """ save the geometry from the initial optimization as a reference geometry
    """
    thy_run_fs[-1].create(thy_info[1:4])
    thy_run_path = thy_run_fs[-1].path(thy_info[1:4])
    run_fs = autofile.fs.run(thy_run_path)

    thy_save_fs[-1].create(thy_info[1:4])
    thy_save_path = thy_save_fs[-1].path(thy_info[1:4])

    ret = es_runner.read_job(job=elstruct.Job.OPTIMIZATION, run_fs=run_fs)
    if ret:
        print('Saving reference geometry...')
        print(" - Save path: {}".format(thy_save_path))

        inf_obj, _, out_str = ret
        prog = inf_obj.prog
        geo = elstruct.reader.opt_geometry(prog, out_str)
        zma = automol.geom.zmatrix(geo)
        thy_save_fs[-1].file.geometry.write(geo, thy_info[1:4])
        thy_save_fs[-1].file.zmatrix.write(zma, thy_info[1:4])


def fake_conf(thy_info, filesystem, inf=()):
    """ generate data to be used for a fake well I think?
    """
    cnf_save_fs = filesystem[5]
    cnf_run_fs = filesystem[4]
    thy_save_fs = filesystem[3]
    run_fs = filesystem[-1]
    thy_save_path = thy_save_fs[-1].path(thy_info[1:4])
    geo = thy_save_fs[-1].file.geometry.read(thy_info[1:4])
    if inf:
        inf_obj, ene = inf
    else:
        ene = thy_save_fs[-1].file.energy.read(thy_info[1:4])
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


def fake_geo_gen(tsk, thy_info, filesystem):
    """ generate data to be used for a fake well I think?
    """
    if 'conf' in tsk:
        fake_conf(thy_info, filesystem)
    if 'scan' in tsk:
        pass
    if 'tau' in tsk:
        pass
