""" drivers for initial geometry optimization
"""
import os
import numpy
from qcelemental import constants as qcc
import automol
import elstruct
import autofile
import moldr
import scripts
import projrot_io

WAVEN2KCAL = qcc.conversion_factor('wavenumber', 'kcal/mol')
EH2KCAL = qcc.conversion_factor('hartree', 'kcal/mol')

def reference_geometry(
        spcdct, thy_level, thy_run_fs, thy_save_fs, cnf_run_fs, cnf_save_fs,
        run_fs, ini_thy_level=[], ini_thy_save_fs=None, kickoff_size=0.1,
        kickoff_backward=False, projrot_script_str='RPHt.exe',
        overwrite=False):
    """ determine what to use as the reference geometry for all future runs
    If ini_thy_info refers to geometry dictionary then use that,
    otherwise values are from a hierarchy of:
    running level of theory, input level of theory, inchis.
    From the hierarchy an optimization is performed followed by a check for
    an imaginary frequency and then a conformer file system is set up.
    """

    print('initializing geometry')
    geo = None
    # Check to see if geometry should be obtained from dictionary
    geom_obj = spcdct['geoobj']
    spc_info = [spcdct['ich'], spcdct['chg'], spcdct['mul']]
    if 'input_geom' in ini_thy_level: # geo is to be read in from dictionary of goemetries
        geo_init = geom_obj
        overwrite = True
        print('found initial geometry from geometry dictionary')
    else:
    # Check to see if geo already exists at running_theory
        if thy_save_fs.leaf.file.geometry.exists(thy_level[1:4]):
            thy_path = thy_save_fs.leaf.path(thy_level[1:4])
            print('getting reference geometry from {}'.format(thy_path))
            geo = thy_save_fs.leaf.file.geometry.read(thy_level[1:4])
        if not geo:
            if ini_thy_save_fs:
                if ini_thy_save_fs.leaf.file.geometry.exists(ini_thy_level[1:4]):
                # If not, Compute geo at running_theory, using geo from
                # initial_level as the starting point
                # or from inchi is no initial level geometry
                    thy_path = ini_thy_save_fs.leaf.path(ini_thy_level[1:4])
                    print('getting reference geometry from {}'.format(thy_path))
                    geo_init = ini_thy_save_fs.leaf.file.geometry.read(ini_thy_level[1:4])
                else:
                    geo_init = automol.inchi.geometry(spc_info[0])
                    print('getting reference geometry from inchi')
            else:
                geo_init = automol.inchi.geometry(spc_info[0])
                print('getting reference geometry from inchi')
        # Optimize from initial geometry to get reference geometry
    if not geo:
        _, opt_script_str, _, opt_kwargs = moldr.util.run_qchem_par(*thy_level[0:2])
        params = {
            'spc_info': spc_info,
            'run_fs': run_fs,
            'thy_run_fs': thy_run_fs,
           # 'thy_save_fs': thy_save_fs,
            'script_str': opt_script_str,
            'overwrite': overwrite,
            'thy_level': thy_level,
            'geo_init': geo_init}
        geo = run_initial_geometry_opt(**params, **opt_kwargs)

        geo, hess = remove_imag(
            spcdct, geo, thy_level, thy_run_fs,
            run_fs, kickoff_size,
            kickoff_backward,
            projrot_script_str,
            overwrite=overwrite)

        tors_names = automol.geom.zmatrix_torsion_coordinate_names(geo)
#        ntaudof = len(automol.graph.rotational_bond_keys(gra, with_h_rotors=False))
#        nsamp = moldr.util.nsamp_init(nsamp_par, ntaudof)

        locs_lst = cnf_save_fs.leaf.existing()
        if locs_lst:
            saved_geo = cnf_save_fs.leaf.file.geometry.read(locs_lst[0])
            saved_tors_names = automol.geom.zmatrix_torsion_coordinate_names(saved_geo)
            if tors_names != saved_tors_names:
                print("new reference geometry doesn't match original reference geometry")
                print('removing original conformer save data')
                cnf_run_fs.remove()
                cnf_save_fs.remove()

#        thy_run_fs.leaf.create(thy_level[1:4])
#        thy_run_path = thy_run_fs.leaf.path(thy_level[1:4])
#        run_fs = autofile.fs.run(thy_run_path)
        thy_save_fs.leaf.create(thy_level[1:4])
        thy_save_path = thy_save_fs.leaf.path(thy_level[1:4])
        print('Saving reference geometry')
        print(" - Save path: {}".format(thy_save_path))
#       ret = moldr.driver.read_job(job=elstruct.Job.OPTIMIZATION, run_fs=run_fs)
#        inf_obj, _, out_str = ret
#        prog = inf_obj.prog
#        geo = elstruct.reader.opt_geometry(prog, out_str)
        zma = automol.geom.zmatrix(geo)
        thy_save_fs.leaf.file.geometry.write(geo, thy_level[1:4])
        thy_save_fs.leaf.file.zmatrix.write(zma, thy_level[1:4])
        thy_save_fs.leaf.file.hessian.write(hess, thy_level[1:4])

        scripts.es.run_single_mc(
            spc_info, thy_level,
            thy_save_fs, cnf_run_fs, cnf_save_fs,
            overwrite)

    return geo


def run_initial_geometry_opt(
        spc_info, thy_level, run_fs, thy_run_fs,
        script_str, overwrite, geo_init, **kwargs):
    """ generate initial geometry via optimization from either reference
    geometries or from inchi
    """
    # set up the filesystem
    thy_run_fs.leaf.create(thy_level[1:4])
    thy_run_path = thy_run_fs.leaf.path(thy_level[1:4])

    # check if geometry has already been saved
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
        print('Succesful reference geometry optimization')
        inf_obj, _, out_str = ret
        prog = inf_obj.prog
        geo = elstruct.reader.opt_geometry(prog, out_str)
        zma = automol.geom.zmatrix(geo)
    return geo


def remove_imag(
        spcdct, geo, thy_level, thy_run_fs, run_fs, kickoff_size=0.1,
        kickoff_backward=False,
        projrot_script_str='RPHt.exe',
        overwrite=False):
    """ if there is an imaginary frequency displace the geometry along the imaginary
    mode and then reoptimize
    """

    print('the initial geometries will be checked for imaginary frequencies')
    spc_info = scripts.es.get_spc_info(spcdct)
    script_str, opt_script_str, kwargs, opt_kwargs = moldr.util.run_qchem_par(*thy_level[0:2])

    imag, geo, disp_xyzs, hess = run_check_imaginary(
        spc_info, geo, thy_level, thy_run_fs, script_str,
        projrot_script_str,
        overwrite, **kwargs)
    chk_idx = 0
    while imag and chk_idx < 5:
        chk_idx += 1
        print('imaginary frequency detected, attempting to kick off')

        moldr.geom.run_kickoff_saddle(
            geo, disp_xyzs, spc_info, thy_level, run_fs, thy_run_fs,
            opt_script_str, kickoff_size, kickoff_backward,
            opt_cart=True, **opt_kwargs)
        print('removing saddlepoint hessian')

        thy_run_path = thy_run_fs.leaf.path(thy_level[1:4])
        run_fs = autofile.fs.run(thy_run_path)
        run_fs.leaf.remove([elstruct.Job.HESSIAN])
        imag, geo, disp_xyzs, hess = run_check_imaginary(
            spc_info, geo, thy_level, thy_run_fs, script_str,
            projrot_script_str,
            overwrite, **kwargs)
    return geo, hess

def run_check_imaginary(
        spc_info, geo, thy_level, thy_run_fs, script_str,
        projrot_script_str='RPHt.exe',
        overwrite=False, **kwargs):
    """ check if species has an imaginary frequency
    """
    thy_run_fs.leaf.create(thy_level[1:4])
    thy_run_path = thy_run_fs.leaf.path(thy_level[1:4])

    run_fs = autofile.fs.run(thy_run_path)
    imag = False
    disp_xyzs = []
    hess = ()
    if automol.geom.is_atom(geo):
        hess = {}
    else:
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
                if automol.geom.is_linear(geo):
                    freqs = elstruct.util.harmonic_frequencies(geo, hess, project=False)
                else:
                    freqs = projrot_frequencies(
                        geo, hess, thy_level, thy_run_fs, projrot_script_str)

    # mode for now set the imaginary frequency check to -100:
    # Ultimately should decrease once frequency projector is functioning properly
                if min(freqs) < 0:
                    imag = True
                    print('Imaginary mode found:')
                    norm_coos = elstruct.util.normal_coordinates(
                        geo, hess, project=True)
                    im_norm_coo = numpy.array(norm_coos)[:, 0]
                    disp_xyzs = numpy.reshape(im_norm_coo, (-1, 3))
    return imag, geo, disp_xyzs, hess


def run_kickoff_saddle(
        geo, disp_xyzs, spc_info, thy_level, run_fs, thy_run_fs,
        opt_script_str, kickoff_size=0.1, kickoff_backward=False,
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
        script_str=opt_script_str,
        run_fs=run_fs,
        geom=geom,
        spc_info=spc_info,
        thy_level=thy_level,
        overwrite=True,
        **kwargs,
    )


def save_initial_geometry(
        thy_level, run_fs, thy_run_fs, thy_save_fs):
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


def projrot_frequencies(geo, hess, thy_level, thy_run_fs, projrot_script_str='RPHt.exe'):
    """ Get the projected frequencies from projrot code
    """
    # Write the string for the ProjRot input
    thy_run_fs.leaf.create(thy_level[1:4])
    thy_run_path = thy_run_fs.leaf.path(thy_level[1:4])

    coord_proj = 'cartesian'
    grad = ''
    rotors_str = ''
    projrot_inp_str = projrot_io._write.write_rpht_input(
        geo, grad, hess, rotors_str=rotors_str,
        coord_proj=coord_proj)

    proj_file_path = os.path.join(thy_run_path, 'RPHt_input_data.dat')
    with open(proj_file_path, 'w') as proj_file:
        proj_file.write(projrot_inp_str)

    proj_file_path = thy_run_path
    moldr.util.run_script(projrot_script_str, proj_file_path)
   
    rtproj_freqs, _ = projrot_io._read.read_rpht_output(
        proj_file_path+'/RTproj_freq.dat')
    rthrproj_freqs, _ = projrot_io._read.read_rpht_output(
        proj_file_path+'/hrproj_freq.dat')
    print('Projection test')
    print(rtproj_freqs)
    proj_freqs = rtproj_freqs
    return proj_freqs
