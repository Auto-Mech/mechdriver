""" drivers
"""
import os
import functools
import numpy
from qcelemental import constants as qcc
from qcelemental import periodictable as ptab
import projrot_io
import automol
import elstruct
import thermo
import autofile
import moldr
from moldr import runner
import mess_io

WAVEN2KCAL = qcc.conversion_factor('wavenumber', 'kcal/mol')
EH2KCAL = qcc.conversion_factor('hartree', 'kcal/mol')

def run_initial_geometry_opt(
        spc_info, theory_level, run_prefix, save_prefix,
        script_str, overwrite, geo_init, **kwargs):
    """ generate initial geometry via optimization from either reference
    geometries or from inchi
    """
    # set up the filesystem
    orb_restr = moldr.util.orbital_restriction(spc_info, theory_level)
    thy_level = theory_level[1:3]
    thy_level.append(orb_restr)

    thy_run_fs = autofile.fs.theory(run_prefix)
    thy_run_fs.leaf.create(thy_level)
    thy_run_path = thy_run_fs.leaf.path(thy_level)

    thy_save_fs = autofile.fs.theory(save_prefix)
    thy_save_fs.leaf.create(thy_level)
    thy_save_path = thy_save_fs.leaf.path(thy_level)

    # check if geometry has already been saved
    if thy_save_fs.leaf.file.geometry.exists(thy_level):
        geo = thy_save_fs.leaf.file.geometry.read(thy_level)
    # or if geometry has already been run
    elif thy_run_fs.leaf.file.geometry.exists(thy_level):
        geo = thy_run_fs.leaf.file.geometry.read(thy_level)
        zma = automol.geom.zmatrix(geo)
        thy_save_fs.leaf.file.geometry.write(geo, thy_level)
        thy_save_fs.leaf.file.zmatrix.write(zma, thy_level)
    else:
        # if not call the electronic structure optimizer
        zma = automol.geom.zmatrix(geo_init)
        run_job(
            job=elstruct.Job.OPTIMIZATION,
            geom=zma,
            spc_info=spc_info,
            theory_level=theory_level,
            prefix=thy_run_path,
            script_str=script_str,
            overwrite=overwrite,
            **kwargs,
        )
        ret = read_job(job=elstruct.Job.OPTIMIZATION, prefix=thy_run_path)
        geo = None
        if ret:
            print('Saving reference geometry')
            print(" - Save path: {}".format(thy_save_path))
            inf_obj, _, out_str = ret
            prog = inf_obj.prog
            geo = elstruct.reader.opt_geometry(prog, out_str)
            zma = automol.geom.zmatrix(geo)
            thy_save_fs.leaf.file.geometry.write(geo, thy_level)
            thy_save_fs.leaf.file.zmatrix.write(zma, thy_level)
    return geo


def run_check_imaginary(
        spc_info, theory_level, run_prefix, save_prefix, script_str,
        overwrite, **kwargs):
    """ check if species has an imaginary frequency
    """
    orb_restr = moldr.util.orbital_restriction(spc_info, theory_level)
    thy_level = theory_level[1:3]
    thy_level.append(orb_restr)

    thy_run_fs = autofile.fs.theory(run_prefix)
    thy_run_fs.leaf.create(thy_level)
    thy_run_path = thy_run_fs.leaf.path(thy_level)

    thy_save_fs = autofile.fs.theory(save_prefix)
    thy_save_fs.leaf.create(thy_level)
    thy_save_path = thy_save_fs.leaf.path(thy_level)
    print('thy_save_path in check imag:', thy_save_path)
    print('thy_level:', thy_level)
    print('save_prefix:', save_prefix)

    ret = read_job(job=elstruct.Job.OPTIMIZATION, prefix=thy_run_path)
    geo = None
    imag = False
    disp_xyzs = []
    if ret:
        print('Checking for saddle')
        inf_obj, _, out_str = ret
        prog = inf_obj.prog
        geo = elstruct.reader.opt_geometry(prog, out_str)

        if thy_save_fs.leaf.file.hessian.exists(thy_level):
            hess = thy_save_fs.leaf.file.hessian.read(thy_level)
        else:
            hess = None
        ret = read_job(job=elstruct.Job.HESSIAN, prefix=thy_run_path)
        if not hess:
            if not ret:
                run_job(
                    job=elstruct.Job.HESSIAN,
                    spc_info=spc_info,
                    theory_level=theory_level,
                    geom=geo,
                    prefix=thy_run_path,
                    script_str=script_str,
                    overwrite=overwrite,
                    **kwargs,
                    )
                ret = read_job(job=elstruct.Job.HESSIAN, prefix=thy_run_path)
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
        geo, disp_xyzs, spc_info, theory_level, run_path,
        script_str, kickoff_size=0.1, kickoff_backward=False,
        opt_cart=True, **kwargs):
    """ kickoff from saddle to find connected minima
    """
    print('kickoff from saddle')
    disp_len = kickoff_size * qcc.conversion_factor('angstrom', 'bohr')
    if kickoff_backward:
        disp_len *= -1
    disp_xyzs = numpy.multiply(disp_xyzs, disp_len)
    geo = automol.geom.displaced(geo, disp_xyzs)
    if opt_cart:
        geom = geo
    else:
        geom = automol.geom.zmatrix(geo)
    run_job(
        job=elstruct.Job.OPTIMIZATION,
        geom=geom,
        spc_info=spc_info,
        theory_level=theory_level,
        prefix=run_path,
        script_str=script_str,
        overwrite=True,
        **kwargs,
    )


def save_initial_geometry(
        spc_info, theory_level, run_prefix, save_prefix):
    """ save the geometry from the initial optimization as a reference geometry
    """
    orb_restr = moldr.util.orbital_restriction(spc_info, theory_level)
    thy_level = theory_level[1:3]
    thy_level.append(orb_restr)

    thy_run_fs = autofile.fs.theory(run_prefix)
    thy_run_fs.leaf.create(thy_level)
    thy_run_path = thy_run_fs.leaf.path(thy_level)

    thy_save_fs = autofile.fs.theory(save_prefix)
    thy_save_fs.leaf.create(thy_level)
    thy_save_path = thy_save_fs.leaf.path(thy_level)

    ret = read_job(job=elstruct.Job.OPTIMIZATION, prefix=thy_run_path)
    if ret:
        print('Saving reference geometry')
        print(" - Save path: {}".format(thy_save_path))

        inf_obj, _, out_str = ret
        prog = inf_obj.prog
        geo = elstruct.reader.opt_geometry(prog, out_str)
        zma = automol.geom.zmatrix(geo)
        thy_save_fs.leaf.file.geometry.write(geo, thy_level)
        thy_save_fs.leaf.file.zmatrix.write(zma, thy_level)


def conformer_sampling(
        spc_info, theory_level, run_prefix, save_prefix,
        script_str, overwrite, saddle=False, nsamp_par=(False, 3, 3, 1, 50, 50), tors_names='',
        **kwargs):
    """ Find the minimum energy conformer by optimizing from nsamp random
    initial torsional states
    """
    orb_restr = moldr.util.orbital_restriction(spc_info, theory_level)
    thy_level = theory_level[1:3]
    thy_level.append(orb_restr)

    thy_run_fs = autofile.fs.theory(run_prefix)
    thy_run_fs.leaf.create(thy_level)
    thy_run_path = thy_run_fs.leaf.path(thy_level)

    thy_save_fs = autofile.fs.theory(save_prefix)
    thy_save_fs.leaf.create(thy_level)
    thy_save_path = thy_save_fs.leaf.path(thy_level)

    run_path = thy_run_path
    save_path = thy_save_path
    geo = thy_save_fs.leaf.file.geometry.read(thy_level)
    tors_names = automol.geom.zmatrix_torsion_coordinate_names(geo)
    zma = automol.geom.zmatrix(geo)
    tors_ranges = automol.zmatrix.torsional_sampling_ranges(
        zma, tors_names)
    tors_range_dct = dict(zip(tors_names, tors_ranges))
    ich = spc_info[0]
    gra = automol.inchi.graph(ich)
    ntaudof = len(automol.graph.rotational_bond_keys(gra, with_h_rotors=False))
    nsamp = moldr.util.nsamp_init(nsamp_par, ntaudof)

    save_conformers(
        run_prefix=run_path,
        save_prefix=save_path,
    )
#    print('zma test in conformer sampling')
#    print(automol.zmatrix.string(zma))

    run_conformers(
        zma=zma,
        spc_info=spc_info,
        theory_level=theory_level,
        nsamp=nsamp,
        tors_range_dct=tors_range_dct,
        run_prefix=run_path,
        save_prefix=save_path,
        script_str=script_str,
        overwrite=overwrite,
        **kwargs,
    )

    save_conformers(
        run_prefix=run_path,
        save_prefix=save_path,
    )

    # save information about the minimum energy conformer in top directory
    cnf_save_fs = autofile.fs.conformer(save_path)
    min_cnf_locs = moldr.util.min_energy_conformer_locators(save_path)
#    print('thy_save_fs:', thy_save_fs)
#    print('thy_level:', thy_level)
    if min_cnf_locs:
        geo = cnf_save_fs.leaf.file.geometry.read(min_cnf_locs)
        zma = cnf_save_fs.leaf.file.zmatrix.read(min_cnf_locs)
        assert automol.zmatrix.almost_equal(zma, automol.geom.zmatrix(geo))
        thy_save_fs.leaf.file.geometry.write(geo, thy_level)
        thy_save_fs.leaf.file.zmatrix.write(zma, thy_level)


def ts_conformer_sampling(
        spc_info, geo, zma, tors_names, theory_level, run_prefix, save_prefix,
        script_str, overwrite, saddle=True, nsamp_par=(False, 3, 3, 1, 50, 50), 
        **kwargs):
    """ Find the minimum energy conformer by optimizing from nsamp random
    initial torsional states
    """
    orb_restr = moldr.util.orbital_restriction(spc_info, theory_level)
    thy_level = theory_level[1:3]
    thy_level.append(orb_restr)

    thy_run_fs = autofile.fs.theory(run_prefix)
    thy_run_fs.leaf.create(thy_level)
    thy_run_path = thy_run_fs.leaf.path(thy_level)

    thy_save_fs = autofile.fs.theory(save_prefix)
    thy_save_fs.leaf.create(thy_level)
    thy_save_path = thy_save_fs.leaf.path(thy_level)
    ts_run_fs = autofile.fs.ts(thy_run_path)
    ts_run_fs.trunk.create()
    ts_run_path = ts_run_fs.trunk.path()
    run_path = ts_run_path
#        run_path = run_prefix

    ts_save_fs = autofile.fs.ts(thy_save_path)
    ts_save_fs.trunk.create()
    ts_save_path = ts_save_fs.trunk.path()
    save_path = ts_save_path
#        save_path = ts_save_path
#    geo = ts_save_fs.trunk.file.geometry.read()

#    zma = automol.geom.zmatrix(geo)
    tors_ranges = automol.zmatrix.torsional_sampling_ranges(
        zma, tors_names)
    tors_range_dct = dict(zip(tors_names, tors_ranges))
    ich = spc_info[0]
    ntaudof = len(tors_names)
    nsamp = moldr.util.nsamp_init(nsamp_par, ntaudof)

    ts_save_conformers(
        run_prefix=run_path,
        save_prefix=save_path,
    )
#    print('zma test in conformer sampling')
#    print(automol.zmatrix.string(zma))

    run_conformers(
        zma=zma,
        spc_info=spc_info,
        theory_level=theory_level,
        nsamp=nsamp,
        tors_range_dct=tors_range_dct,
        run_prefix=run_path,
        save_prefix=save_path,
        script_str=script_str,
        overwrite=overwrite,
        **kwargs,
    )

    ts_save_conformers(
        run_prefix=run_path,
        save_prefix=save_path,
    )

    # save information about the minimum energy conformer in top directory
    cnf_save_fs = autofile.fs.conformer(save_path)
    min_cnf_locs = moldr.util.min_energy_conformer_locators(save_path)
#    print('thy_save_fs:', thy_save_fs)
#    print('thy_level:', thy_level)
    if min_cnf_locs:
        geo = cnf_save_fs.leaf.file.geometry.read(min_cnf_locs)
        zma = cnf_save_fs.leaf.file.zmatrix.read(min_cnf_locs)

        ts_save_fs.trunk.file.geometry.write(geo)
        ts_save_fs.trunk.file.zmatrix.write(zma)


def run_conformers(
        zma, spc_info, theory_level, nsamp, tors_range_dct,
        run_prefix, save_prefix, script_str, overwrite,
        **kwargs):
    """ run sampling algorithm to find conformers
    """
    if not tors_range_dct:
        print("No torsional coordinates. Setting nsamp to 1.")
        nsamp = 1

    print('Number of samples requested:', nsamp)

    cnf_run_fs = autofile.fs.conformer(run_prefix)
    cnf_save_fs = autofile.fs.conformer(save_prefix)
    cnf_save_fs.trunk.create()

    vma = automol.zmatrix.var_(zma)
    if cnf_save_fs.trunk.file.vmatrix.exists():
        existing_vma = cnf_save_fs.trunk.file.vmatrix.read()
        assert vma == existing_vma
    cnf_save_fs.trunk.file.vmatrix.write(vma)
    idx = 0
    nsamp0 = nsamp
    inf_obj = autofile.system.info.conformer_trunk(0, tors_range_dct)
    while True:
        if cnf_save_fs.trunk.file.info.exists():
            inf_obj_s = cnf_save_fs.trunk.file.info.read()
            nsampd = inf_obj_s.nsamp
        elif cnf_run_fs.trunk.file.info.exists():
            inf_obj_r = cnf_run_fs.trunk.file.info.read()
            nsampd = inf_obj_r.nsamp
        else:
            nsampd = 0

        nsamp = nsamp0 - nsampd
        if nsamp <= 0:
            print('Reached requested number of samples. '
                  'Conformer search complete.')
            break
        else:
            print("    New nsamp requested is {:d}.".format(nsamp))

            samp_zma, = automol.zmatrix.samples(zma, 1, tors_range_dct)
            cid = autofile.system.generate_new_conformer_id()
            locs = [cid]

            cnf_run_fs.leaf.create(locs)
            cnf_run_path = cnf_run_fs.leaf.path(locs)

            idx += 1
            print("Run {}/{}".format(idx, nsamp0))
            run_job(
                job=elstruct.Job.OPTIMIZATION,
                script_str=script_str,
                prefix=cnf_run_path,
                geom=samp_zma,
                spc_info=spc_info,
                theory_level=theory_level,
                overwrite=overwrite,
                **kwargs
            )

            nsampd += 1
            inf_obj.nsamp = nsampd
            cnf_save_fs.trunk.file.info.write(inf_obj)
            cnf_run_fs.trunk.file.info.write(inf_obj)


def save_conformers(run_prefix, save_prefix):
    """ save the conformers that have been found so far
    """
    cnf_run_fs = autofile.fs.conformer(run_prefix)
    cnf_save_fs = autofile.fs.conformer(save_prefix)

    locs_lst = cnf_save_fs.leaf.existing()
    seen_geos = [cnf_save_fs.leaf.file.geometry.read(locs)
                 for locs in locs_lst]
    seen_enes = [cnf_save_fs.leaf.file.energy.read(locs)
                 for locs in locs_lst]

    if not cnf_run_fs.trunk.exists():
        print("No conformers to save. Skipping...")
    else:
        for locs in cnf_run_fs.leaf.existing():
            run_path = cnf_run_fs.leaf.path(locs)
            print("Reading from conformer run at {}".format(run_path))

            ret = read_job(job=elstruct.Job.OPTIMIZATION, prefix=run_path)
            if ret:
                inf_obj, inp_str, out_str = ret
                prog = inf_obj.prog
                method = inf_obj.method
                ene = elstruct.reader.energy(prog, method, out_str)
                geo = elstruct.reader.opt_geometry(prog, out_str)
                zma = automol.geom.zmatrix(geo)
                gra = automol.geom.graph(geo)

                if len(automol.graph.connected_components(gra)) > 1:
                    print(" - Geometry is disconnected.. Skipping...")
                else:
                    unique = True

                    for idx, geoi in enumerate(seen_geos):
                        enei = seen_enes[idx]
                        etol = 1.e-6
                        if automol.geom.almost_equal_coulomb_spectrum(
                                geo, geoi, rtol=1e-2):
                            if abs(ene-enei) < etol:
                                unique = False

                    if not unique:
                        print(" - Geometry is not unique. Skipping...")
                    else:
                        save_path = cnf_save_fs.leaf.path(locs)
                        print(" - Geometry is unique. Saving...")
                        print(" - Save path: {}".format(save_path))

                        cnf_save_fs.leaf.create(locs)
                        cnf_save_fs.leaf.file.geometry_info.write(
                            inf_obj, locs)
                        cnf_save_fs.leaf.file.geometry_input.write(
                            inp_str, locs)
                        cnf_save_fs.leaf.file.energy.write(ene, locs)
                        cnf_save_fs.leaf.file.geometry.write(geo, locs)
                        cnf_save_fs.leaf.file.zmatrix.write(zma, locs)

                seen_geos.append(geo)
                seen_enes.append(ene)

        # update the conformer trajectory file
        locs_lst = cnf_save_fs.leaf.existing()
        if locs_lst:
            enes = [cnf_save_fs.leaf.file.energy.read(locs)
                    for locs in locs_lst]
            geos = [cnf_save_fs.leaf.file.geometry.read(locs)
                    for locs in locs_lst]

            traj = []
            for ene, geo in sorted(zip(enes, geos), key=lambda x: x[0]):
                comment = 'energy: {:>15.10f}'.format(ene)
                traj.append((comment, geo))

            traj_path = cnf_save_fs.trunk.file.trajectory.path()
            print("Updating conformer trajectory file at {}".format(traj_path))
            cnf_save_fs.trunk.file.trajectory.write(traj)


def ts_save_conformers(run_prefix, save_prefix):
    """ save the conformers that have been found so far
    """
    cnf_run_fs = autofile.fs.conformer(run_prefix)
    cnf_save_fs = autofile.fs.conformer(save_prefix)

    locs_lst = cnf_save_fs.leaf.existing()
    seen_geos = [cnf_save_fs.leaf.file.geometry.read(locs)
                 for locs in locs_lst]
    seen_enes = [cnf_save_fs.leaf.file.energy.read(locs)
                 for locs in locs_lst]

    if not cnf_run_fs.trunk.exists():
        print("No conformers to save. Skipping...")
    else:
        for locs in cnf_run_fs.leaf.existing():
            run_path = cnf_run_fs.leaf.path(locs)
            print("Reading from conformer run at {}".format(run_path))

            ret = read_job(job=elstruct.Job.OPTIMIZATION, prefix=run_path)
            if ret:
                inf_obj, inp_str, out_str = ret
                prog = inf_obj.prog
                method = inf_obj.method
                ene = elstruct.reader.energy(prog, method, out_str)
                geo = elstruct.reader.opt_geometry(prog, out_str)
                zma = elstruct.reader.opt_zmatrix(prog, out_str)

#                zma = automol.geom.zmatrix(geo)
                unique = True

                for idx, geoi in enumerate(seen_geos):
                    enei = seen_enes[idx]
                    etol = 1.e-6
                    if automol.geom.almost_equal_coulomb_spectrum(
                            geo, geoi, rtol=1e-2):
                        if abs(ene-enei) < etol:
                            unique = False
                if not unique:
                    print(" - Geometry is not unique. Skipping...")
                else:
                    save_path = cnf_save_fs.leaf.path(locs)
                    print(" - Geometry is unique. Saving...")
                    print(" - Save path: {}".format(save_path))

                    cnf_save_fs.leaf.create(locs)
                    cnf_save_fs.leaf.file.geometry_info.write(
                        inf_obj, locs)
                    cnf_save_fs.leaf.file.geometry_input.write(
                        inp_str, locs)
                    cnf_save_fs.leaf.file.energy.write(ene, locs)
                    cnf_save_fs.leaf.file.geometry.write(geo, locs)
                    cnf_save_fs.leaf.file.zmatrix.write(zma, locs)

                seen_geos.append(geo)
                seen_enes.append(ene)

        # update the conformer trajectory file
        locs_lst = cnf_save_fs.leaf.existing()
        if locs_lst:
            enes = [cnf_save_fs.leaf.file.energy.read(locs)
                    for locs in locs_lst]
            geos = [cnf_save_fs.leaf.file.geometry.read(locs)
                    for locs in locs_lst]

            traj = []
            for ene, geo in sorted(zip(enes, geos), key=lambda x: x[0]):
                comment = 'energy: {:>15.10f}'.format(ene)
                traj.append((comment, geo))

            traj_path = cnf_save_fs.trunk.file.trajectory.path()
            print("Updating conformer trajectory file at {}".format(traj_path))
            cnf_save_fs.trunk.file.trajectory.write(traj)


# TODO: rewrite this
def run_single_point_energy(
        geo, spc_info, theory_level, run_prefix, save_prefix,
        script_str, overwrite, **kwargs):
    """ Find the energy for the minimum energy conformer
    """
    orb_restr = moldr.util.orbital_restriction(spc_info, theory_level)
    thy_level = theory_level[1:3]
    thy_level.append(orb_restr)

    sp_run_fs = autofile.fs.single_point(run_prefix)
    sp_run_fs.leaf.create(thy_level)
    sp_run_path = sp_run_fs.leaf.path(thy_level)

    sp_save_fs = autofile.fs.single_point(save_prefix)
    sp_save_fs.leaf.create(thy_level)
    sp_save_path = sp_save_fs.leaf.path(thy_level)

    run_job(
        job='energy',
        script_str=script_str,
        prefix=sp_run_path,
        geom=geo,
        spc_info=spc_info,
        theory_level=theory_level,
        overwrite=overwrite,
        **kwargs,
    )

    ret = read_job(
        job='energy',
        prefix=sp_run_path,
    )

    if ret is not None:
        inf_obj, inp_str, out_str = ret

        print(" - Reading energy from output...")
        ene = elstruct.reader.energy(inf_obj.prog, inf_obj.method, out_str)

        print(" - Saving energy...")
        print(" - Save path: {}".format(sp_save_path))
        sp_save_fs.leaf.file.energy.write(ene, thy_level)
        sp_save_fs.leaf.file.input.write(inp_str, thy_level)
        sp_save_fs.leaf.file.info.write(inf_obj, thy_level)


def run_minimum_energy_gradient(
        spc_info, theory_level, run_prefix, save_prefix,
        script_str, overwrite, **kwargs):
    """ Find the gradient for the minimum energy conformer
    """
    cnf_run_fs = autofile.fs.conformer(run_prefix)
    cnf_save_fs = autofile.fs.conformer(save_prefix)

    min_cnf_locs = moldr.util.min_energy_conformer_locators(save_prefix)
    if min_cnf_locs:
        min_cnf_run_path = cnf_run_fs.leaf.path(min_cnf_locs)
        min_cnf_save_path = cnf_save_fs.leaf.path(min_cnf_locs)
        geo = cnf_save_fs.leaf.file.geometry.read(min_cnf_locs)
        print('Minimum energy conformer gradient')
        run_job(
            job='gradient',
            script_str=script_str,
            prefix=min_cnf_run_path,
            geom=geo,
            spc_info=spc_info,
            theory_level=theory_level,
            overwrite=overwrite,
            **kwargs,
        )

        ret = read_job(
            job='gradient',
            prefix=min_cnf_run_path,
        )

        if ret is not None:
            inf_obj, inp_str, out_str = ret

            print(" - Reading gradient from output...")
            grad = elstruct.reader.gradient(inf_obj.prog, out_str)

            print(" - Saving gradient...")
            print(" - Save path: {}".format(min_cnf_save_path))
            cnf_save_fs.leaf.file.gradient_info.write(inf_obj, min_cnf_locs)
            cnf_save_fs.leaf.file.gradient_input.write(inp_str, min_cnf_locs)
            cnf_save_fs.leaf.file.gradient.write(grad, min_cnf_locs)


def run_minimum_energy_hessian(
        spc_info, theory_level, run_prefix, save_prefix,
        script_str, overwrite, **kwargs):
    """ Find the hessian for the minimum energy conformer
    """
    cnf_run_fs = autofile.fs.conformer(run_prefix)
    cnf_save_fs = autofile.fs.conformer(save_prefix)

    min_cnf_locs = moldr.util.min_energy_conformer_locators(save_prefix)
    if min_cnf_locs:
        min_cnf_run_path = cnf_run_fs.leaf.path(min_cnf_locs)
        min_cnf_save_path = cnf_save_fs.leaf.path(min_cnf_locs)
        geo = cnf_save_fs.leaf.file.geometry.read(min_cnf_locs)
        print('Minimum energy conformer hessian')
        run_job(
            job='hessian',
            script_str=script_str,
            prefix=min_cnf_run_path,
            geom=geo,
            spc_info=spc_info,
            theory_level=theory_level,
            overwrite=overwrite,
            **kwargs
        )

        ret = read_job(
            job='hessian',
            prefix=min_cnf_run_path,
        )

        if ret is not None:
            inf_obj, inp_str, out_str = ret

            if automol.geom.is_atom(geo):
                freqs = ()
            else:
                print(" - Reading hessian from output...")
                hess = elstruct.reader.hessian(inf_obj.prog, out_str)
                freqs = elstruct.util.harmonic_frequencies(
                    geo, hess, project=False)
                print('Freqs test')
                print(freqs)
                freqs = elstruct.util.harmonic_frequencies(
                    geo, hess, project=True)

                print(" - Saving hessian...")
                print(" - Save path: {}".format(min_cnf_save_path))
                cnf_save_fs.leaf.file.hessian_info.write(inf_obj, min_cnf_locs)
                cnf_save_fs.leaf.file.hessian_input.write(inp_str, min_cnf_locs)
                cnf_save_fs.leaf.file.hessian.write(hess, min_cnf_locs)
                cnf_save_fs.leaf.file.harmonic_frequencies.write(
                    freqs, min_cnf_locs)


def run_minimum_energy_vpt2(
        spc_info, theory_level, run_prefix, save_prefix,
        script_str, overwrite, **kwargs):
    """  Run vpt2 for the minimum energy conformer
    """
    cnf_run_fs = autofile.fs.conformer(run_prefix)
    cnf_save_fs = autofile.fs.conformer(save_prefix)

    min_cnf_locs = moldr.util.min_energy_conformer_locators(save_prefix)
    if min_cnf_locs:
        min_cnf_run_path = cnf_run_fs.leaf.path(min_cnf_locs)
        min_cnf_save_path = cnf_save_fs.leaf.path(min_cnf_locs)
        geo = cnf_save_fs.leaf.file.geometry.read(min_cnf_locs)
        print('Minimum energy conformer vpt2')
        run_job(
            job='vpt2',
            script_str=script_str,
            prefix=min_cnf_run_path,
            geom=geo,
            spc_info=spc_info,
            theory_level=theory_level,
            overwrite=overwrite,
            **kwargs
        )

        ret = read_job(
            job='vpt2',
            prefix=min_cnf_run_path,
        )

        if ret is not None:
            inf_obj, inp_str, out_str = ret
            prog = inf_obj.prog

            vpt2_dict = elstruct.reader.vpt2(prog, out_str)
            anh_freqs = vpt2_dict['freqs']
            anh_zpe = vpt2_dict['zpe']
            xmat = vpt2_dict['xmat']
            vibrot_mat = vpt2_dict['vibrot']
            cent_dist = vpt2_dict['cent_dist']

            print(" - Reading anharmonic data from output...")
            # replace following with vpt2 information
#            print(" - Saving hessian...")
#            print(" - Save path: {}".format(min_cnf_save_path))
            cnf_save_fs.leaf.file.vpt2_info.write(inf_obj, min_cnf_locs)
            cnf_save_fs.leaf.file.vpt2_input.write(inp_str, min_cnf_locs)
            cnf_save_fs.leaf.file.anh_freqs.write(anh_freqs, min_cnf_locs)
#            cnf_save_fs.leaf.file.zpe.write(zpe, min_cnf_locs)
#            cnf_save_fs.leaf.file.harmonic_frequencies.write(freqs, min_cnf_locs)
            cnf_save_fs.leaf.file.anh_freqs.write(anh_freqs, min_cnf_locs)
            cnf_save_fs.leaf.file.xmat.write(xmat, min_cnf_locs)
#            cnf_save_fs.leaf.file.vibrot.write(vibrot, min_cnf_locs)
            cnf_save_fs.leaf.file.cent_dist.write(cent_dist, min_cnf_locs)


def run_conformer_gradients(
        spc_info, theory_level, run_prefix, save_prefix,
        script_str, overwrite, **kwargs):
    """ Determine the gradient for each of the conformers
    """
    cnf_run_fs = autofile.fs.conformer(run_prefix)
    cnf_save_fs = autofile.fs.conformer(save_prefix)

    cnf_locs_lst = cnf_save_fs.leaf.existing()
    for locs in cnf_locs_lst:
        cnf_run_path = cnf_run_fs.leaf.path(locs)
        cnf_save_path = cnf_save_fs.leaf.path(locs)
        geo = cnf_save_fs.leaf.file.geometry.read(locs)

        print('Running conformer gradient')
        run_job(
            job='gradient',
            script_str=script_str,
            prefix=cnf_run_path,
            geom=geo,
            spc_info=spc_info,
            theory_level=theory_level,
            overwrite=overwrite,
            **kwargs,
        )

        ret = read_job(
            job='gradient',
            prefix=cnf_run_path,
        )

        if ret is not None:
            inf_obj, inp_str, out_str = ret

            print(" - Reading gradient from output...")
            grad = elstruct.reader.gradient(inf_obj.prog, out_str)

            print(" - Saving gradient...")
            print(" - Save path: {}".format(cnf_save_path))
            cnf_save_fs.leaf.file.gradient_info.write(inf_obj, locs)
            cnf_save_fs.leaf.file.gradient_input.write(inp_str, locs)
            cnf_save_fs.leaf.file.gradient.write(grad, locs)


def run_conformer_hessians(
        spc_info, theory_level, run_prefix, save_prefix,
        script_str, overwrite, **kwargs):
    """ Determine the hessian for each of the conformers
    """
    cnf_run_fs = autofile.fs.conformer(run_prefix)
    cnf_save_fs = autofile.fs.conformer(save_prefix)

    cnf_locs_lst = cnf_save_fs.leaf.existing()
    for locs in cnf_locs_lst:
        cnf_run_path = cnf_run_fs.leaf.path(locs)
        cnf_save_path = cnf_save_fs.leaf.path(locs)
        geo = cnf_save_fs.leaf.file.geometry.read(locs)

        print('Running conformer hessian')
        run_job(
            job='hessian',
            script_str=script_str,
            prefix=cnf_run_path,
            geom=geo,
            spc_info=spc_info,
            theory_level=theory_level,
            overwrite=overwrite,
            **kwargs,
        )

        ret = read_job(
            job='hessian',
            prefix=cnf_run_path,
        )

        if ret is not None:
            inf_obj, inp_str, out_str = ret

            if automol.geom.is_atom(geo):
                freqs = ()
            else:
                print(" - Reading hessian from output...")
                hess = elstruct.reader.hessian(inf_obj.prog, out_str)
                freqs = elstruct.util.harmonic_frequencies(
                    geo, hess, project=False)
#                print('Conformer Freqs test')
#                print(freqs)

                print(" - Saving hessian...")
                print(" - Save path: {}".format(cnf_save_path))
                cnf_save_fs.leaf.file.hessian_info.write(inf_obj, locs)
                cnf_save_fs.leaf.file.hessian_input.write(inp_str, locs)
                cnf_save_fs.leaf.file.hessian.write(hess, locs)
                cnf_save_fs.leaf.file.harmonic_frequencies.write(freqs, locs)


# d. hindered rotor scans
def hindered_rotor_scans(
        spc_info, theory_level, run_prefix, save_prefix,
        script_str, overwrite, scan_increment=30., **opt_kwargs):
    """ Perform 1d scans over each of the torsional coordinates
    """
    cnf_run_fs = autofile.fs.conformer(run_prefix)
    cnf_save_fs = autofile.fs.conformer(save_prefix)

    min_cnf_locs = moldr.util.min_energy_conformer_locators(save_prefix)
    if min_cnf_locs:
        min_cnf_run_path = cnf_run_fs.leaf.path(min_cnf_locs)
        min_cnf_save_path = cnf_save_fs.leaf.path(min_cnf_locs)
        geo = cnf_save_fs.leaf.file.geometry.read(min_cnf_locs)
        zma = cnf_save_fs.leaf.file.zmatrix.read(min_cnf_locs)
        val_dct = automol.zmatrix.values(zma)
        tors_names = automol.geom.zmatrix_torsion_coordinate_names(geo)
        tors_linspaces = automol.zmatrix.torsional_scan_linspaces(
            zma, tors_names, scan_increment)
        tors_grids = [
            numpy.linspace(*linspace) + val_dct[name]
            for name, linspace in zip(tors_names, tors_linspaces)]

        for tors_name, tors_grid in zip(tors_names, tors_grids):
            run_scan(
                zma=zma,
                spc_info=spc_info,
                theory_level=theory_level,
                grid_dct={tors_name: tors_grid},
                run_prefix=min_cnf_run_path,
                save_prefix=min_cnf_save_path,
                script_str=script_str,
                overwrite=overwrite,
                **opt_kwargs,
            )

            print('min_cnf_save_path in hindered_rotor_scan')
            print(min_cnf_save_path)

            save_scan(
                run_prefix=min_cnf_run_path,
                save_prefix=min_cnf_save_path,
                coo_names=[tors_name],
            )


def tau_sampling(
        spc_info, theory_level, run_prefix, save_prefix,
        script_str, overwrite, nsamp_par, **opt_kwargs):
    """ Sample over torsions optimizing all other coordinates
    """
    orb_restr = moldr.util.orbital_restriction(spc_info, theory_level)
    thy_level = theory_level[1:3]
    thy_level.append(orb_restr)

    thy_run_fs = autofile.fs.theory(run_prefix)
    thy_run_fs.leaf.create(thy_level)
    thy_run_path = thy_run_fs.leaf.path(thy_level)

    thy_save_fs = autofile.fs.theory(save_prefix)
    thy_save_fs.leaf.create(thy_level)
    thy_save_path = thy_save_fs.leaf.path(thy_level)

    save_tau(
        run_prefix=thy_run_path,
        save_prefix=thy_save_path,
    )

    geo = thy_save_fs.leaf.file.geometry.read(thy_level)
    zma = automol.geom.zmatrix(geo)
    tors_names = automol.geom.zmatrix_torsion_coordinate_names(geo)
    tors_ranges = automol.zmatrix.torsional_sampling_ranges(
        zma, tors_names)
    tors_range_dct = dict(zip(tors_names, tors_ranges))
    ich = spc_info[0]
    gra = automol.inchi.graph(ich)
    ntaudof = len(automol.graph.rotational_bond_keys(gra, with_h_rotors=False))
    nsamp = moldr.util.nsamp_init(nsamp_par, ntaudof)

    run_tau(
        zma=zma,
        spc_info=spc_info,
        theory_level=theory_level,
        nsamp=nsamp,
        tors_range_dct=tors_range_dct,
        run_prefix=thy_run_path,
        save_prefix=thy_save_path,
        script_str=script_str,
        overwrite=overwrite,
        **opt_kwargs,
    )

    save_tau(
        run_prefix=thy_run_path,
        save_prefix=thy_save_path,
    )


def run_tau(
        zma, spc_info, theory_level, nsamp, tors_range_dct,
        run_prefix, save_prefix, script_str, overwrite, **kwargs):
    """ run sampling algorithm to find tau dependent geometries
    """
    if not tors_range_dct:
        print("No torsional coordinates. Setting nsamp to 1.")
        nsamp = 1

    tau_run_fs = autofile.fs.tau(run_prefix)
    tau_save_fs = autofile.fs.tau(save_prefix)
    tau_save_fs.trunk.create()

    vma = automol.zmatrix.var_(zma)
    if tau_save_fs.trunk.file.vmatrix.exists():
        existing_vma = tau_save_fs.trunk.file.vmatrix.read()
        assert vma == existing_vma
    tau_save_fs.trunk.file.vmatrix.write(vma)
    idx = 0
    nsamp0 = nsamp
    inf_obj = autofile.system.info.tau_trunk(0, tors_range_dct)
    while True:
        if tau_save_fs.trunk.file.info.exists():
            inf_obj_s = tau_save_fs.trunk.file.info.read()
            nsampd = inf_obj_s.nsamp
        elif tau_save_fs.trunk.file.info.exists():
            inf_obj_r = tau_save_fs.trunk.file.info.read()
            nsampd = inf_obj_r.nsamp
        else:
            nsampd = 0

        nsamp = nsamp0 - nsampd
        if nsamp <= 0:
            print('Reached requested number of samples. '
                  'Tau sampling complete.')
            break
        else:
            print("    New nsamp is {:d}.".format(nsamp))

            samp_zma, = automol.zmatrix.samples(zma, 1, tors_range_dct)
            tid = autofile.system.generate_new_tau_id()
            locs = [tid]

            tau_run_fs.leaf.create(locs)
            run_path = tau_run_fs.leaf.path(locs)

            idx += 1
            print("Run {}/{}".format(idx, nsamp0))
            run_job(
                job=elstruct.Job.OPTIMIZATION,
                script_str=script_str,
                prefix=run_path,
                geom=samp_zma,
                spc_info=spc_info,
                theory_level=theory_level,
                overwrite=overwrite,
                frozen_coordinates=tors_range_dct.keys(),
                **kwargs
            )

            nsampd += 1
            inf_obj.nsamp = nsampd
            tau_save_fs.trunk.file.info.write(inf_obj)
            tau_run_fs.trunk.file.info.write(inf_obj)


def save_tau(run_prefix, save_prefix):
    """ save the tau dependent geometries that have been found so far
    """
    tau_run_fs = autofile.fs.tau(run_prefix)
    tau_save_fs = autofile.fs.tau(save_prefix)

    saved_geos = [tau_save_fs.leaf.file.geometry.read(locs)
                  for locs in tau_save_fs.leaf.existing()]

    if not tau_run_fs.trunk.exists():
        print("No tau geometries to save. Skipping...")
    else:
        for locs in tau_run_fs.leaf.existing():
            run_path = tau_run_fs.leaf.path(locs)
            print("Reading from tau run at {}".format(run_path))

            ret = read_job(job=elstruct.Job.OPTIMIZATION, prefix=run_path)
            if ret:
                inf_obj, inp_str, out_str = ret
                prog = inf_obj.prog
                method = inf_obj.method
                ene = elstruct.reader.energy(prog, method, out_str)

                geo = elstruct.reader.opt_geometry(prog, out_str)

                save_path = tau_save_fs.leaf.path(locs)
                print(" - Saving...")
                print(" - Save path: {}".format(save_path))

                tau_save_fs.leaf.create(locs)
                tau_save_fs.leaf.file.geometry_info.write(inf_obj, locs)
                tau_save_fs.leaf.file.geometry_input.write(inp_str, locs)
                tau_save_fs.leaf.file.energy.write(ene, locs)
                tau_save_fs.leaf.file.geometry.write(geo, locs)

                saved_geos.append(geo)

        # update the tau trajectory file
        locs_lst = tau_save_fs.leaf.existing()
        if locs_lst:
            enes = [tau_save_fs.leaf.file.energy.read(locs)
                    for locs in locs_lst]
            geos = [tau_save_fs.leaf.file.geometry.read(locs)
                    for locs in locs_lst]

            traj = []
            for ene, geo in sorted(zip(enes, geos), key=lambda x: x[0]):
                comment = 'energy: {:>15.10f}'.format(ene)
                traj.append((comment, geo))

            traj_path = tau_save_fs.trunk.file.trajectory.path()
            print("Updating tau trajectory file at {}".format(traj_path))
            tau_save_fs.trunk.file.trajectory.write(traj)


def run_tau_gradients(
        spc_info, theory_level, run_prefix, save_prefix,
        script_str, overwrite, **kwargs):
    """ Determine gradients for all tau dependent optimized geometries
    """
    tau_run_fs = autofile.fs.tau(run_prefix)
    tau_save_fs = autofile.fs.tau(save_prefix)

    for locs in tau_save_fs.leaf.existing():
        geo = tau_save_fs.leaf.file.geometry.read(locs)
        print('Running tau gradient')
        tau_run_path = tau_run_fs.leaf.path(locs)
        tau_save_path = tau_save_fs.leaf.path(locs)

        run_job(
            job='gradient',
            script_str=script_str,
            prefix=tau_run_path,
            geom=geo,
            spc_info=spc_info,
            theory_level=theory_level,
            overwrite=overwrite,
            **kwargs,
        )

        ret = read_job(
            job='gradient',
            prefix=tau_run_path,
        )

        if ret is not None:
            inf_obj, inp_str, out_str = ret

            print(" - Reading gradient from output...")
            grad = elstruct.reader.gradient(inf_obj.prog, out_str)

            print(" - Saving gradient...")
            print(" - Save path: {}".format(tau_save_path))
            tau_save_fs.leaf.file.gradient_info.write(inf_obj, locs)
            tau_save_fs.leaf.file.gradient_input.write(inp_str, locs)
            tau_save_fs.leaf.file.gradient.write(grad, locs)


def run_tau_hessians(
        spc_info, theory_level, run_prefix, save_prefix,
        script_str, overwrite, **kwargs):
    """ Determine hessians for all tau dependent optimized geometries
    """
    tau_run_fs = autofile.fs.tau(run_prefix)
    tau_save_fs = autofile.fs.tau(save_prefix)

    for locs in tau_save_fs.leaf.existing():
        geo = tau_save_fs.leaf.file.geometry.read(locs)
        print('Running tau Hessian')
        tau_run_path = tau_run_fs.leaf.path(locs)
        tau_save_path = tau_run_fs.leaf.path(locs)
        run_job(
            job='hessian',
            script_str=script_str,
            prefix=tau_run_path,
            geom=geo,
            spc_info=spc_info,
            theory_level=theory_level,
            overwrite=overwrite,
            **kwargs,
        )

        ret = read_job(
            job='hessian',
            prefix=tau_run_path,
        )

        if ret is not None:
            inf_obj, inp_str, out_str = ret

            print(" - Reading hessian from output...")
            hess = elstruct.reader.hessian(inf_obj.prog, out_str)

            print(" - Saving hessian...")
            print(" - Save path: {}".format(tau_save_path))
            tau_save_fs.leaf.file.hessian_info.write(inf_obj, locs)
            tau_save_fs.leaf.file.hessian_input.write(inp_str, locs)
            tau_save_fs.leaf.file.hessian.write(hess, locs)


def tau_pf_write(
        name, save_prefix,
        run_grad=False, run_hess=False):
    """ Print out data fle for partition function evaluation 
    """
    cnf_save_fs = autofile.fs.conformer(save_prefix)
    min_cnf_locs = moldr.util.min_energy_conformer_locators(save_prefix)
    if min_cnf_locs:
        ene_ref = cnf_save_fs.leaf.file.energy.read(min_cnf_locs)
        print('ene_ref')
        print(ene_ref)

    tau_save_fs = autofile.fs.tau(save_prefix)
    evr = name+'\n'
    # cycle through saved tau geometries
    idx = 0
    for locs in tau_save_fs.leaf.existing():
        geo = tau_save_fs.leaf.file.geometry.read(locs)
        ene = tau_save_fs.leaf.file.energy.read(locs)
        ene = (ene - ene_ref)*qcc.conversion_factor('hartree', 'kcal/mol')
        ene_str = autofile.file.write.energy(ene)
        geo_str = autofile.file.write.geometry(geo)

        idx += 1
        idx_str = str(idx)

        evr += 'Sampling point'+idx_str+'\n'
        evr += 'Energy'+'\n'
        evr += ene_str+'\n'
        evr += 'Geometry'+'\n'
        evr += geo_str+'\n'
        if run_grad:
            grad = tau_save_fs.leaf.file.gradient.read(locs)
            grad_str = autofile.file.write.gradient(grad)
            evr += 'Gradient'+'\n'
            evr += grad_str
        if run_hess:
            hess = tau_save_fs.leaf.file.hessian.read(locs)
            hess_str = autofile.file.write.hessian(hess)
            evr += 'Hessian'+'\n'
            evr += hess_str+'\n'

    file_name = os.path.join(save_prefix, 'TAU', 'tau.out')
    with open(file_name, 'w') as tau_file:
        tau_file.write(evr)

    temp_list = [300., 500., 750., 1000., 1500.]
    for temp in temp_list:
        sumq = 0.
        sum2 = 0.
        idx = 0
        print('integral convergence for T = ', temp)
        for locs in tau_save_fs.leaf.existing():
            idx += 1
            ene = tau_save_fs.leaf.file.energy.read(locs)
            ene = (ene - ene_ref)*qcc.conversion_factor('hartree', 'kcal/mol')
            tmp = numpy.exp(-ene*349.7/(0.695*temp))
            sumq = sumq + tmp
            sum2 = sum2 + tmp**2
            sigma = numpy.sqrt(
                (abs(sum2/float(idx)-(sumq/float(idx))**2))/float(idx))
            print(sumq/float(idx), sigma, 100.*sigma*float(idx)/sumq, idx)


def run_scan(
        zma, spc_info, theory_level, grid_dct, run_prefix,
        save_prefix, script_str, overwrite, update_guess=True,
        reverse_sweep=True, **kwargs):
    """ run constrained optimization scan
    """
    if len(grid_dct) > 1:
        raise NotImplementedError

    scn_run_fs = autofile.fs.scan(run_prefix)
    scn_save_fs = autofile.fs.scan(save_prefix)

    vma = automol.zmatrix.var_(zma)
    if scn_save_fs.trunk.file.vmatrix.exists():
        existing_vma = scn_save_fs.trunk.file.vmatrix.read()
        assert vma == existing_vma

    # for now, running only one-dimensional hindered rotor scans
    ((coo_name, grid_vals),) = grid_dct.items()
    scn_save_fs.branch.create([[coo_name]])

    if scn_save_fs.branch.file.info.exists([[coo_name]]):
        inf_obj = scn_save_fs.branch.file.info.read([[coo_name]])
        existing_grid_dct = dict(inf_obj.grids)
        existing_grid_vals = existing_grid_dct[coo_name]
        assert (numpy.shape(grid_vals) == numpy.shape(existing_grid_vals) and
                numpy.allclose(grid_vals, existing_grid_vals))

    inf_obj = autofile.system.info.scan_branch({coo_name: grid_vals})
    scn_save_fs.branch.file.info.write(inf_obj, [[coo_name]])

    npoint = len(grid_vals)
    grid_idxs = tuple(range(npoint))

    for grid_idx in grid_idxs:
        scn_run_fs.leaf.create([[coo_name], [grid_idx]])

    prefixes = tuple(scn_run_fs.leaf.path([[coo_name], [grid_idx]])
                     for grid_idx in grid_idxs)

    _run_1d_scan(
        script_str=script_str,
        prefixes=prefixes,
        guess_zma=zma,
        coo_name=coo_name,
        grid_idxs=grid_idxs,
        grid_vals=grid_vals,
        spc_info=spc_info,
        theory_level=theory_level,
        overwrite=overwrite,
        update_guess=update_guess,
        **kwargs
    )

    if reverse_sweep:
        _run_1d_scan(
            script_str=script_str,
            prefixes=list(reversed(prefixes)),
            guess_zma=zma,
            coo_name=coo_name,
            grid_idxs=list(reversed(grid_idxs)),
            grid_vals=list(reversed(grid_vals)),
            spc_info=spc_info,
            theory_level=theory_level,
            overwrite=overwrite,
            update_guess=update_guess,
            **kwargs
        )
 

def run_multiref_scan(
        formula, high_mul, zma, spc_info, theory_level, grid_dct, 
        run_prefix, save_prefix, script_str, overwrite, update_guess=True,
        **kwargs):
    """ run constrained optimization scan
    """
    if len(grid_dct) > 1:
#    if len(grid_dct) > 1 or len(grid_dct2) > 1:
        raise NotImplementedError

    electron_count = automol.formula._formula.electron_count(formula)
    # this is only for 2e,2o case
    closed_orb = electron_count//2 - 1
    occ_orb = electron_count//2 + 1
    # end of 2e,2o case
    two_spin = spc_info[2]-1
    chg = spc_info[1]
    cas_options = [
        elstruct.option.specify(elstruct.Option.Scf.MAXITER_, 40),
        elstruct.option.specify(elstruct.Option.Casscf.OCC_, occ_orb),
        elstruct.option.specify(elstruct.Option.Casscf.CLOSED_, closed_orb),
        elstruct.option.specify(elstruct.Option.Casscf.WFN_, electron_count, 1, two_spin, chg)
        ]

    scn_run_fs = autofile.fs.scan(run_prefix)
    scn_save_fs = autofile.fs.scan(save_prefix)

    vma = automol.zmatrix.var_(zma)
    if scn_save_fs.trunk.file.vmatrix.exists():
        existing_vma = scn_save_fs.trunk.file.vmatrix.read()
        assert vma == existing_vma

    # for now, running only one-dimensional hindered rotor scans

    print('length test:', len(grid_dct))
    ((coo_name, grid_vals),) = grid_dct.items()
    print('coo_name test:', coo_name)
    scn_save_fs.branch.create([[coo_name]])
    if scn_save_fs.branch.file.info.exists([[coo_name]]):
        inf_obj = scn_save_fs.branch.file.info.read([[coo_name]])
        print('inf_obj test:', inf_obj)
        existing_grid_dct = dict(inf_obj.grids)
        existing_grid_vals = existing_grid_dct[coo_name]
        print('grid_vals test:', grid_vals, existing_grid_vals)
        assert (numpy.shape(grid_vals) == numpy.shape(existing_grid_vals) and
                numpy.allclose(grid_vals, existing_grid_vals))

    inf_obj = autofile.system.info.scan_branch({coo_name: grid_vals})
    scn_save_fs.branch.file.info.write(inf_obj, [[coo_name]])
    charge = spc_info[1]
    mul = spc_info[2]
    method = theory_level[1]
    basis = theory_level[2]
    prog = theory_level[0]
    orb_restr = moldr.util.orbital_restriction(spc_info, theory_level)
    #orb_restr = theory_level[3]
    theory_level[0] = 'molpro'
    prog = 'molpro'
    theory_level[1] = 'caspt2'
    _, script_str, _, kwargs = moldr.util.run_qchem_par(theory_level[0], theory_level[1])
    kwargs['casscf_options'] = cas_options
    kwargs['mol_options'] = ['nosym']
    print('orb_restr test:', orb_restr)

    guess_str1 = elstruct.writer.energy(
        geom=automol.zmatrix.set_values(zma, {coo_name: grid_vals[0]}),
        charge=charge,
        mult=high_mul,
        method='hf',
        basis=basis,
        prog=prog,
        orb_restricted=orb_restr,
        **kwargs)
    guess_str1 += '\n\n'
    guess_str1 = '\n'.join(guess_str1.splitlines()[2:])

    guess_str2 = elstruct.writer.energy(
        geom=automol.zmatrix.set_values(zma, {coo_name: grid_vals[0]}),
        charge=charge,
        mult=mul,
        method='casscf',
        basis=basis,
        prog=prog,
        orb_restricted=orb_restr,
        **kwargs)
    guess_str2 += '\n\n'
    guess_str2 = '\n'.join(guess_str2.splitlines()[2:])
    
    guess_str = guess_str1 + guess_str2
    guess_lines = guess_str.splitlines()
    kwargs['gen_lines'] = guess_lines

    npoint = len(grid_vals)
    grid_idxs = tuple(range(npoint))

    for grid_idx in grid_idxs:
        scn_run_fs.leaf.create([[coo_name], [grid_idx]])
    prefixes = tuple(scn_run_fs.leaf.path([[coo_name], [grid_idx]])
                     for grid_idx in grid_idxs)
    print('update_guess0 test:', update_guess)
    _run_1d_scan(
        script_str=script_str,
        prefixes=prefixes,
        guess_zma=zma,
        coo_name=coo_name,
        grid_idxs=grid_idxs,
        grid_vals=grid_vals,
        spc_info=spc_info,
        theory_level=theory_level,
        overwrite=overwrite,
        update_guess=update_guess,
        **kwargs
    )

#    ((coo_name, grid_vals2),) = grid_dct2.items()
#    scn_save_fs.branch.create([[coo_name]])
#    if scn_save_fs.branch.file.info.exists([[coo_name]]):
#        inf_obj = scn_save_fs.branch.file.info.read([[coo_name]])
#        existing_grid_dct = dict(inf_obj.grids)
#        existing_grid_vals2 = existing_grid_dct[coo_name]
#        print('grid_vals test 2:', grid_vals2, existing_grid_vals2)
#        assert (numpy.shape(grid_vals2) == numpy.shape(existing_grid_vals2) and
#                numpy.allclose(grid_vals2, existing_grid_vals2))
#
#    inf_obj = autofile.system.info.scan_branch({coo_name: grid_vals2})
#    scn_save_fs.branch.file.info.write(inf_obj, [[coo_name]])
#    charge = spc_info[1]
#    mult = spc_info[2]
#    method = theory_level[1]
#    basis = theory_level[2]
#    prog = theory_level[0]
#    orb_restr = theory_level[3]
#    prog = 'molpro'
#
#    if len(grid_dct2) > 1:
#        raise NotImplementedError
#
#    npoint2 = len(grid_vals2)
#    grid_idxs2 = tuple(range(npoint2))
#
#    for grid_idx in grid_idxs2:
#        scn_run_fs.leaf.create([[coo_name], [grid_idx]])
#
#    prefixes2 = tuple(scn_run_fs.leaf.path([[coo_name], [grid_idx]])
#                     for grid_idx in grid_idxs2)
#
#    _run_1d_scan(
#        script_str=script_str,
#        prefixes=prefixes2,
#        guess_zma=zma,
#        coo_name=coo_name,
#        grid_idxs=grid_idxs2,
#        grid_vals=grid_vals2,
#        spc_info=spc_info,
#        theory_level=theory_level,
#        overwrite=overwrite,
#        update_guess=update_guess,
#        **kwargs
#    )
  

def _run_1d_scan(
        script_str, prefixes, guess_zma, coo_name, grid_idxs, grid_vals,
        spc_info, theory_level, overwrite, errors=(),
        options_mat=(), retry_failed=True, update_guess=True, **kwargs):
    npoints = len(grid_idxs)
    assert len(grid_vals) == len(prefixes) == npoints
    for grid_idx, grid_val, prefix in zip(grid_idxs, grid_vals, prefixes):
        print("Point {}/{}".format(grid_idx+1, npoints))
        zma = automol.zmatrix.set_values(guess_zma, {coo_name: grid_val})
        print('zma test:', zma)

        run_job(
            job=elstruct.Job.OPTIMIZATION,
            script_str=script_str,
            prefix=prefix,
            geom=zma,
            spc_info=spc_info,
            theory_level=theory_level,
            overwrite=overwrite,
            frozen_coordinates=[coo_name],
            errors=errors,
            options_mat=options_mat,
            retry_failed=retry_failed,
            **kwargs
        )

        ret = read_job(job=elstruct.Job.OPTIMIZATION, prefix=prefix)
        # print('ret test:', ret)
        print('update_guess test:', update_guess)
        if update_guess and ret is not None:
            inf_obj, _, out_str = ret
            print('inf_obj test:', inf_obj)
            print('out_str test:', out_str)
            prog = inf_obj.prog
            method = inf_obj.method
            guess_zma = elstruct.reader.opt_zmatrix(prog, out_str)
            print('guess_zma test:', guess_zma)


def save_scan(run_prefix, save_prefix, coo_names):
    """ save the scans that have been run so far
    """
    print(run_prefix)
    print(save_prefix)
    print(coo_names)
    if len(coo_names) > 1:
        raise NotImplementedError

    scn_run_fs = autofile.fs.scan(run_prefix)
    scn_save_fs = autofile.fs.scan(save_prefix)

    if not scn_run_fs.branch.exists([coo_names]):
        print("No scan to save. Skipping...")
    else:
        locs_lst = []
        for locs in scn_run_fs.leaf.existing([coo_names]):
            run_path = scn_run_fs.leaf.path(locs)
            print("Reading from scan run at {}".format(run_path))

            ret = read_job(job=elstruct.Job.OPTIMIZATION, prefix=run_path)
            if ret:
                inf_obj, inp_str, out_str = ret
                prog = inf_obj.prog
                method = inf_obj.method
                ene = elstruct.reader.energy(prog, method, out_str)
                geo = elstruct.reader.opt_geometry(prog, out_str)
                zma = elstruct.reader.opt_zmatrix(prog, out_str)

                save_path = scn_save_fs.leaf.path(locs)
                print(" - Saving...")
                print(" - Save path: {}".format(save_path))

                scn_save_fs.leaf.create(locs)
                scn_save_fs.leaf.file.geometry_info.write(inf_obj, locs)
                scn_save_fs.leaf.file.geometry_input.write(inp_str, locs)
                scn_save_fs.leaf.file.energy.write(ene, locs)
                scn_save_fs.leaf.file.geometry.write(geo, locs)
                scn_save_fs.leaf.file.zmatrix.write(zma, locs)

                locs_lst.append(locs)

        if locs_lst:
            idxs_lst = [locs[-1] for locs in locs_lst]
            enes = [scn_save_fs.leaf.file.energy.read(locs)
                    for locs in locs_lst]
            geos = [scn_save_fs.leaf.file.geometry.read(locs)
                    for locs in locs_lst]

            traj = []
            for idxs, ene, geo in zip(idxs_lst, enes, geos):
                comment = 'energy: {:>15.10f}, grid idxs: {}'.format(ene, idxs)
                traj.append((comment, geo))

            traj_path = scn_save_fs.branch.file.trajectory.path([coo_names])
            print("Updating scan trajectory file at {}".format(traj_path))
            scn_save_fs.branch.file.trajectory.write(traj, [coo_names])


# centralized job runners and readers
JOB_ERROR_DCT = {
    elstruct.Job.ENERGY: elstruct.Error.SCF_NOCONV,
    elstruct.Job.GRADIENT: elstruct.Error.SCF_NOCONV,
    elstruct.Job.HESSIAN: elstruct.Error.SCF_NOCONV,
    elstruct.Job.OPTIMIZATION: elstruct.Error.OPT_NOCONV,
}

JOB_RUNNER_DCT = {
    elstruct.Job.ENERGY: functools.partial(
        runner.options_matrix_run, elstruct.writer.energy),
    elstruct.Job.GRADIENT: functools.partial(
        runner.options_matrix_run, elstruct.writer.gradient),
    elstruct.Job.HESSIAN: functools.partial(
        runner.options_matrix_run, elstruct.writer.hessian),
    elstruct.Job.OPTIMIZATION: runner.options_matrix_optimization,
}


def run_job(
        job, script_str, prefix,
        geom, spc_info, theory_level, 
        errors=(), options_mat=(), retry_failed=True, feedback=False,
        frozen_coordinates=(), freeze_dummy_atoms=True, overwrite=False,
        **kwargs):
    """ run an elstruct job by name
    """
    assert job in JOB_RUNNER_DCT
    assert job in JOB_ERROR_DCT

    run_fs = autofile.fs.run(prefix)
    run_fs.leaf.create([job])
    run_path = run_fs.leaf.path([job])

    if overwrite:
        do_run = True
        print(" - Running {} job at {}".format(job, run_path))
    else:
        if not run_fs.leaf.file.info.exists([job]):
            do_run = True
            print(" - Running {} job at {}".format(job, run_path))
        else:
            inf_obj = run_fs.leaf.file.info.read([job])
            if inf_obj.status == autofile.system.RunStatus.FAILURE:
                print(" - Found failed {} job at {}".format(job, run_path))
                if retry_failed:
                    print(" - Retrying...")
                    do_run = True
                else:
                    do_run = False
            else:
                do_run = False
                if inf_obj.status == autofile.system.RunStatus.SUCCESS:
                    print(" - Found completed {} job at {}"
                          .format(job, run_path))
                else:
                    print(" - Found running {} job at {}"
                          .format(job, run_path))
                    print(" - Skipping...")

    # I think the len(geom) test was supposed to stop the call for gradients and hessians for atoms
    # It didn't work properly for two reasons:
    # (i) if the geom is passed as a z-matrix, len(geom) is 2 not 1
    # (ii) it was checking for UPPERCASE job words and its or logic made it always true
    # (ii) if I fix those problems then if you turn off the optimization, then the conformer directories are not created
    # which creates evern more problems. For now I just turned these tests off, since gaussians seems to run fine for atom
    # gradients and hessians, etc.
#    print('do_run_0 test:', do_run)
#    print('geom test:',geom)
#    print(len(geom))
#    if len(geom) == 1:
#        print(len(geom))
#        print(job)
#        if job in ('hessian', 'optimization', 'gradient', 'vpt2'):
#            print(job)
#            do_run = False
#            print('do_run_1 test:', do_run)

#    print ('do_run test:', do_run)
    if do_run:
        # create the run directory
        status = autofile.system.RunStatus.RUNNING
        prog = theory_level[0]
        method = theory_level[1]
        basis = theory_level[2]
        inf_obj = autofile.system.info.run(
            job=job, prog=prog, method=method, basis=basis, status=status)
        inf_obj.utc_start_time = autofile.system.info.utc_time()
        run_fs.leaf.file.info.write(inf_obj, [job])

        runner = JOB_RUNNER_DCT[job]
        if job == elstruct.Job.OPTIMIZATION:
            runner = functools.partial(
                runner, feedback=feedback,
                frozen_coordinates=frozen_coordinates,
                freeze_dummy_atoms=freeze_dummy_atoms)

        orb_restr = moldr.util.orbital_restriction(spc_info, theory_level)
        inp_str, out_str = runner(
            script_str, run_path, geom=geom, chg=spc_info[1],
            mul=spc_info[2], method=theory_level[1], basis=theory_level[2],
            orb_restricted=orb_restr, prog=theory_level[0],
            errors=errors, options_mat=options_mat, **kwargs
        )

        inf_obj.utc_end_time = autofile.system.info.utc_time()
        prog = inf_obj.prog

        if is_successful_output(out_str, job, prog):
            run_fs.leaf.file.output.write(out_str, [job])
            print(" - Run succeeded.")
            status = autofile.system.RunStatus.SUCCESS
        else:
            print(" - Run failed.")
            status = autofile.system.RunStatus.FAILURE
        inf_obj.status = status
        run_fs.leaf.file.info.write(inf_obj, [job])
        run_fs.leaf.file.input.write(inp_str, [job])
        print('finished run_job')


def read_job(job, prefix):
    """ read from an elstruct job by name
    """
    ret = None

    run_fs = autofile.fs.run(prefix)
    run_path = run_fs.leaf.path([job])

    print(" - Reading from {} job at {}".format(job, run_path))
    if not run_fs.leaf.file.output.exists([job]):
        print(" - No output file found. Skipping...")
    else:
        assert run_fs.leaf.file.info.exists([job])
        assert run_fs.leaf.file.input.exists([job])
        inf_obj = run_fs.leaf.file.info.read([job])
        inp_str = run_fs.leaf.file.input.read([job])
        out_str = run_fs.leaf.file.output.read([job])
        prog = inf_obj.prog

        if is_successful_output(out_str, job, prog):
            print(" - Found successful output. Reading...")
            ret = (inf_obj, inp_str, out_str)
        else:
            print(" - Output has an error message. Skipping...")

    return ret


def is_successful_output(out_str, job, prog):
    """ is this a successful output string?
    """
    assert job in JOB_ERROR_DCT
    error = JOB_ERROR_DCT[job]

    ret = False
    if elstruct.reader.has_normal_exit_message(prog, out_str):
        if not elstruct.reader.has_error_message(prog, error, out_str):
            ret = True

    return ret


