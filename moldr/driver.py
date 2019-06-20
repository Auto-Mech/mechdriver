""" drivers
"""
# import itertools
import functools
import os
import warnings
import numpy
from qcelemental import constants as qcc
import automol
import elstruct
import autofile
from autofile import SFS
from moldr import optsmat


def run_conformers(ich, charge, mult, method, basis, orb_restricted,
                   nsamp, run_prefix, save_prefix, script_str, prog,
                   **kwargs):
    """ run sampling algorithm to find conformers
    """
    geo = automol.inchi.geometry(ich)
    zma = automol.geom.zmatrix(geo)
    tors_names = automol.geom.zmatrix_torsion_coordinate_names(geo)
    tors_range_vals = automol.zmatrix.torsional_sampling_ranges(
        zma, tors_names)
    tors_ranges = dict(zip(tors_names, tors_range_vals))

    if not tors_ranges:
        print("No torsional coordinates. Setting nsamp to 1.")
        nsamp = 1

    root_specs = (ich, charge, mult, method, basis, orb_restricted)

    # check for a previously saved run
    vma = automol.zmatrix.var_(zma)
    if SFS.conf_trunk.dir.exists(save_prefix, root_specs):
        _vma = SFS.conf_trunk.file.vmatrix.read(save_prefix, root_specs)
        assert vma == _vma
        inf_obj = SFS.conf_trunk.file.info.read(save_prefix, root_specs)
        nsamp = max(nsamp - inf_obj.nsamp, 0)
        print("Found previous saved run. Adjusting nsamp.")
        print("    New nsamp is {:d}.".format(nsamp))
    else:
        SFS.conf_trunk.dir.create(save_prefix, root_specs)
        inf_obj = autofile.system.info.conformer_trunk(
            nsamp=0, tors_ranges=tors_ranges)
    SFS.conf_trunk.file.vmatrix.write(vma, save_prefix, root_specs)
    SFS.conf_trunk.file.info.write(inf_obj, save_prefix, root_specs)

    # update the number of samples

    inp_zmas = automol.zmatrix.samples(zma, nsamp, tors_ranges)

    cids = tuple(autofile.system.generate_new_conformer_id()
                 for _ in range(nsamp))

    print()
    print("Running optimizations in run directories.")
    job = 'optimization'
    for idx, (cid, inp_zma) in enumerate(zip(cids, inp_zmas)):
        specs = root_specs + (cid, job)

        if not SFS.conf_run.dir.exists(run_prefix, specs):
            run_path = SFS.conf_run.dir.path(run_prefix, specs)
            print("Starting run {}/{} at {}".format(idx+1, nsamp, run_path))

            status = "running"
            run_inf_obj = autofile.system.info.run(
                job=job, prog=prog, method=method, basis=basis, status=status)
            run_inf_obj.utc_start_time = autofile.system.info.utc_time()

            SFS.conf_run.dir.create(run_prefix, specs)
            SFS.conf_run.file.info.write(run_inf_obj, run_prefix, specs)

#            inp_str, out_str = feedback_optimization(
#                script_str, run_path,
#                geom=inp_zma, charge=charge, mult=mult, method=method,
#                basis=basis, prog=prog,
#                **kwargs)
            options_mat = [
                [{'scf_options': (elstruct.Option.Scf.Guess.CORE,)},
                 {'scf_options': (elstruct.Option.Scf.Guess.HUCKEL,)},
                 {'scf_options': (
                     elstruct.option.specify(elstruct.Option.Scf.DIIS_, True),
                     elstruct.Option.Scf.Guess.HUCKEL,)}],
                [{'job_options': (elstruct.Option.Opt.Coord.ZMATRIX,)},
                 {'job_options': (
                     elstruct.option.specify(elstruct.Option.Opt.MAXITER_, 10),
                     elstruct.Option.Opt.Coord.ZMATRIX,)}],
            ]
            errors = [
                elstruct.Error.SCF_NOCONV,
                elstruct.Error.OPT_NOCONV,
            ]
#            inp_str, out_str = robust_feedback_opt(
            inp_str, out_str = robust_feedback_opt(
                script_str, run_path,
                geom=inp_zma, charge=charge, mult=mult, method=method,
                basis=basis, prog=prog,
                errors=errors, options_mat=options_mat, **kwargs)

            run_inf_obj.utc_end_time = autofile.system.info.utc_time()

            if elstruct.reader.has_normal_exit_message(prog, out_str):
                SFS.conf_run.file.output.write(out_str, run_prefix, specs)
                status = "succeeded"
            else:
                status = "failed"
            run_inf_obj.status = status
            SFS.conf_run.file.info.write(run_inf_obj, run_prefix, specs)
            SFS.conf_run.file.input.write(inp_str, run_prefix, specs)

            print("Run {}/{} {} at {}".format(idx+1, nsamp, status, run_path))


def save_conformers(ich, charge, mult, method, basis, orb_restricted,
                    run_prefix, save_prefix):
    """ save the conformers that have been found so far
    """
#    if not SFS.conf.dir.exists(run_prefix, root_specs):
#        print ('conformer directory does not exist')
#    else:
#        SFS.conf_trunk.dir.create(run_prefix, root_specs)
    root_specs = (ich, charge, mult, method, basis, orb_restricted)
    print('save_conf_test', run_prefix, root_specs)
    if SFS.conf_trunk.dir.exists(run_prefix, root_specs):
        run_conf_specs_lst = SFS.conf.dir.existing(run_prefix, root_specs)
        saved_conf_specs_lst = SFS.conf.dir.existing(save_prefix, root_specs)

        conf_specs_lst = []
        ene_lst = []
        geo_lst = []
        inp_str_lst = []
        inf_obj_lst = []

        print()
        print("Reading optimizations from run directories.")
        print(root_specs)
        run_specs = ('optimization',)
        for conf_specs in run_conf_specs_lst:
            specs = root_specs + conf_specs + run_specs

            run_path = SFS.conf_run.dir.path(run_prefix, specs)
            print("Reading from run at {}".format(run_path))

            if SFS.conf_run.file.output.exists(run_prefix, specs):
                inf_obj = SFS.conf_run.file.info.read(run_prefix, specs)
                inp_str = SFS.conf_run.file.input.read(run_prefix, specs)
                out_str = SFS.conf_run.file.output.read(run_prefix, specs)
                prog = inf_obj.prog
                if elstruct.reader.has_normal_exit_message(prog, out_str):
                    ene = elstruct.reader.energy(prog, method, out_str)
                    geo = elstruct.reader.opt_geometry(prog, out_str)

                    print(automol.geom.coulomb_spectrum(geo))
                    # save the information to a list
                    conf_specs_lst.append(conf_specs)
                    inf_obj_lst.append(inf_obj)
                    inp_str_lst.append(inp_str)
                    ene_lst.append(ene)
                    geo_lst.append(geo)

        seen_geo_lst = []
        for conf_specs in saved_conf_specs_lst:
            specs = root_specs + conf_specs
            geo = SFS.conf.file.geometry.read(save_prefix, specs)
            seen_geo_lst.append(geo)

        print("Writing unique conformer information to save directories.")
        idxs = automol.geom.argunique_coulomb_spectrum(
            geo_lst, seen_geos=seen_geo_lst, rtol=1e-3)
        for idx in idxs:
            conf_specs = conf_specs_lst[idx]
            inf_obj = inf_obj_lst[idx]
            inp_str = inp_str_lst[idx]
            ene = ene_lst[idx]
            geo = geo_lst[idx]

            specs = root_specs + conf_specs
            save_path = SFS.conf.dir.path(save_prefix, specs)
            print("Saving values from run at {}".format(save_path))

            SFS.conf.dir.create(save_prefix, specs)
            SFS.conf.file.geometry_info.write(inf_obj, save_prefix, specs)
            SFS.conf.file.geometry_input.write(inp_str, save_prefix, specs)
            SFS.conf.file.energy.write(ene, save_prefix, specs)
            SFS.conf.file.geometry.write(geo, save_prefix, specs)

        # update the number of samples
        nsamp_new = len(conf_specs_lst)
        trunk_inf_obj = SFS.conf_trunk.file.info.read(save_prefix, root_specs)
        trunk_inf_obj.nsamp += nsamp_new
        SFS.conf_trunk.file.info.write(trunk_inf_obj, save_prefix, root_specs)


def run_scan(ich, charge, mult, method, basis, orb_restricted, cid,
             run_prefix, save_prefix, script_str, prog, scan_incr=30.,
             # ncoords,
             **kwargs):
    """ run a scan
    """
    root_specs = (ich, charge, mult, method, basis, orb_restricted, cid)
    print('just before exists assert')
    print(cid)
    print(save_prefix)
    print(root_specs)
    print(SFS.conf.file.geometry.path(save_prefix, root_specs))
#    assert SFS.conf.file.geometry.exists(save_prefix, root_specs)
    if not SFS.conf.file.geometry.exists(save_prefix, root_specs):
        print('file does not exist')
    else:
        geo = SFS.conf.file.geometry.read(save_prefix, root_specs)
        zma = automol.geom.zmatrix(geo)

        vma = automol.zmatrix.var_(zma)
        if SFS.scan_trunk.dir.exists(save_prefix, root_specs):
            _vma = SFS.scan_trunk.file.vmatrix.read(save_prefix, root_specs)
            assert vma == _vma
        SFS.scan_trunk.dir.create(save_prefix, root_specs)
        SFS.scan_trunk.file.vmatrix.write(vma, save_prefix, root_specs)

        print(root_specs)
        print(run_prefix, save_prefix, script_str, prog, kwargs)
        print(SFS.scan_trunk.dir.path(save_prefix, root_specs))

        tors_names = automol.geom.zmatrix_torsion_coordinate_names(geo)
        increment = scan_incr*qcc.conversion_factor('degree', 'radian')
        tors_linspace_vals = automol.zmatrix.torsional_scan_grids(
            zma, tors_names, increment)
        tors_linspaces = dict(zip(tors_names, tors_linspace_vals))

        job = 'optimization'
        for tors_name, linspace in tors_linspaces.items():
            branch_specs = root_specs + ([tors_name],)
            inf_obj = autofile.system.info.scan_branch({tors_name: linspace})

            SFS.scan_branch.dir.create(save_prefix, branch_specs)
            SFS.scan_branch.file.info.write(inf_obj, save_prefix, branch_specs)

            last_zma = zma

            grid = numpy.linspace(*linspace)
            npoint = len(grid)
            for grid_idx, grid_val in enumerate(grid):
                specs = branch_specs + ([grid_idx], job)

                if not SFS.scan_run.dir.exists(run_prefix, specs):
                    run_path = SFS.scan_run.dir.path(run_prefix, specs)
                    print("Starting run {}/{} at {}"
                          .format(grid_idx+1, npoint, run_path))
                    inp_zma = automol.zmatrix.set_values(
                        last_zma, {tors_name: grid_val})

                    SFS.scan_run.dir.create(run_prefix, specs)

                    status = "running"
                    run_inf_obj = autofile.system.info.run(
                        job=job, prog=prog, method=method, basis=basis,
                        status=status)
                    run_inf_obj.utc_start_time = (
                        autofile.system.info.utc_time())

                    SFS.scan_run.dir.create(run_prefix, specs)
                    SFS.scan_run.file.info.write(run_inf_obj, run_prefix, specs)

    #                print (inp_zma)
    #                inp_str, out_str = feedback_optimization(
    #                    script_str, run_path,
    #                    geom=inp_zma, charge=charge, mult=mult, method=method,
    #                    basis=basis, prog=prog,
    #                    frozen_coordinates=[tors_name],
    #                    **kwargs)
    #                inp_str,out_str = robust_feedback_opt(
                    options_mat = [
                        [{'scf_options': (elstruct.Option.Scf.Guess.CORE,)},
                         {'scf_options': (elstruct.Option.Scf.Guess.HUCKEL,)},
                         {'scf_options': (
                             elstruct.option.specify(elstruct.Option.Scf.DIIS_, True),
                             elstruct.Option.Scf.Guess.HUCKEL,)}],
                        [{'job_options': (elstruct.Option.Opt.Coord.ZMATRIX,)},
                         {'job_options': (
                             elstruct.option.specify(elstruct.Option.Opt.MAXITER_, 10),
                             elstruct.Option.Opt.Coord.ZMATRIX,)}],
                    ]
                    errors = [
                        elstruct.Error.SCF_NOCONV,
                        elstruct.Error.OPT_NOCONV,
                    ]
                    inp_str, out_str = robust_feedback_opt(
                        script_str, run_path,
                        geom=inp_zma, charge=charge, mult=mult, method=method,
                        basis=basis, prog=prog,
                        frozen_coordinates=[tors_name],
                        errors=errors, options_mat=options_mat, **kwargs)

                    run_inf_obj.utc_end_time = autofile.system.info.utc_time()


                    if elstruct.reader.has_normal_exit_message(prog, out_str):
                        SFS.scan_run.file.output.write(out_str, run_prefix, specs)
                        status = "succeeded"
                        last_zma = elstruct.reader.opt_zmatrix(prog, out_str)
                    else:
                        status = "failed"

                    run_inf_obj.status = status
                    SFS.scan_run.file.info.write(run_inf_obj, run_prefix, specs)
                    SFS.scan_run.file.input.write(inp_str, run_prefix, specs)
                    print("Run {}/{} {} at {}"
                          .format(grid_idx+1, npoint, status, run_path))


#def save_scan(ich, charge, mult, method, basis, orb_restricted, cid,
#              run_prefix, save_prefix):
#    """ save geometries and energies from a scan
#    """
#    root_specs = (ich, charge, mult, method, basis, orb_restricted, cid)
#
#    run_conf_specs_lst = SFS.conf.dir.existing(run_prefix, root_specs)
#    saved_conf_specs_lst = SFS.conf.dir.existing(save_prefix, root_specs)
#    saved_scan_lst = SFS.conf.
#    conf_specs_lst = []
#    ene_lst = []
#    geo_lst = []
#    inp_str_lst = []
#    inf_obj_lst = []
#
#    print()
#    print("Reading torsional scans from run directories.")
#    print(root_specs)
#
#    assert SFS.conf.file.geometry.exists(save_prefix, root_specs)
#    geo = SFS.conf.file.geometry.read(save_prefix, root_specs)
#    zma = automol.geom.zmatrix(geo)
#
#    vma = automol.zmatrix.var_(zma)
#    if SFS.scan_trunk.dir.exists(save_prefix, root_specs):
#        _vma = SFS.scan_trunk.file.vmatrix.read(save_prefix, root_specs)
#        assert vma == _vma
#    else:
#        SFS.scan_trunk.dir.create(save_prefix, root_specs)
#    SFS.scan_trunk.file.vmatrix.write(vma, save_prefix, root_specs)
#
#    print(root_specs)
#    print(run_prefix, save_prefix, script_str, prog, kwargs)
#    print(SFS.scan_trunk.dir.path(save_prefix, root_specs))
#
#    tors_names = automol.geom.zmatrix_torsion_coordinate_names(geo)
#    tors_linspace_vals = automol.zmatrix.torsional_scan_grids(zma, tors_names)
#    tors_linspaces = dict(zip(tors_names, tors_linspace_vals))
#
#    job = 'optimization'
#    for tors_name, linspace in tors_linspaces.items():
#        branch_specs = root_specs + ([tors_name],)
#        inf_obj = autofile.system.info.scan_branch({tors_name: linspace})
#
#        SFS.scan_branch.dir.create(save_prefix, branch_specs)
#        SFS.scan_branch.file.info.write(inf_obj, save_prefix, branch_specs)
#
#        last_zma = zma
#
#        grid = numpy.linspace(*linspace)
#        npoint = len(grid)
#        for grid_idx, grid_val in enumerate(grid):
#            specs = branch_specs + ([grid_idx], job)
#
#            if not SFS.scan_run.dir.exists(run_prefix, specs):
#                run_path = SFS.scan_run.dir.path(run_prefix, specs)
#                print("Starting run {}/{} at {}"
#                      .format(grid_idx+1, npoint, run_path))
#                inp_zma = automol.zmatrix.set_values(
#                    last_zma, {tors_name: grid_val})
#
#                SFS.scan_run.dir.create(run_prefix, specs)
#
#                run_inf_obj = autofile.system.info.run(
#                    job=job, prog=prog, method=method, basis=basis)
#                run_inf_obj.utc_start_time = autofile.system.info.utc_time()
#
#                SFS.scan_run.dir.create(run_prefix, specs)
#                SFS.scan_run.file.info.write(run_inf_obj, run_prefix, specs)
#
#                inp_str, out_str = feedback_optimization(
#                    script_str, run_path,
#                    prog, method, basis, inp_zma, mult, charge,
#                    frozen_coordinates=[tors_name],
#                    **kwargs)
#
#                run_inf_obj.utc_end_time = autofile.system.info.utc_time()
#
#                SFS.scan_run.file.info.write(run_inf_obj, run_prefix, specs)
#                SFS.scan_run.file.input.write(inp_str, run_prefix, specs)
#
#                status = "failed"
#                if elstruct.reader.has_normal_exit_message(prog, out_str):
#                    SFS.scan_run.file.output.write(out_str, run_prefix, specs)
#                    status = "succeeded"
#
#                    last_zma = elstruct.reader.opt_zmatrix(prog, out_str)
#
#                print("Run {}/{} {} at {}"
#                      .format(grid_idx+1, npoint, status, run_path))
#


def run_tau(ich, charge, mult, method, basis, orb_restricted,
            nsamp, run_prefix, save_prefix, script_str, prog,
            **kwargs):
    """ run sampling algorithm to find taus
    """
    geo = automol.inchi.geometry(ich)
    zma = automol.geom.zmatrix(geo)
    tors_names = automol.geom.zmatrix_torsion_coordinate_names(geo)
    tors_range_vals = automol.zmatrix.torsional_sampling_ranges(zma, tors_names)
    tors_ranges = dict(zip(tors_names, tors_range_vals))

    if not tors_ranges:
        print("No torsional coordinates. Setting nsamp to 1.")
        nsamp = 1

    root_specs = (ich, charge, mult, method, basis, orb_restricted)

    # check for a previously saved run
    vma = automol.zmatrix.var_(zma)
    if SFS.tau_trunk.dir.exists(save_prefix, root_specs):
        _vma = SFS.tau_trunk.file.vmatrix.read(save_prefix, root_specs)
        assert vma == _vma
        inf_obj = SFS.tau_trunk.file.info.read(save_prefix, root_specs)
        nsamp = max(nsamp - inf_obj.nsamp, 0)
        print("Found previous saved run. Adjusting nsamp.")
        print("    New nsamp is {:d}.".format(nsamp))
    else:
        SFS.tau_trunk.dir.create(save_prefix, root_specs)
        inf_obj = autofile.system.info.tau_trunk(
            nsamp=0, tors_ranges=tors_ranges)
    SFS.tau_trunk.file.vmatrix.write(vma, save_prefix, root_specs)
    SFS.tau_trunk.file.info.write(inf_obj, save_prefix, root_specs)

    # update the number of samples

    inp_zmas = automol.zmatrix.samples(zma, nsamp, tors_ranges)

    cids = tuple(autofile.system.generate_new_conformer_id()
                 for _ in range(nsamp))

    print()
    print("Running optimizations in run directories.")
    job = 'optimization'
    for idx, (cid, inp_zma) in enumerate(zip(cids, inp_zmas)):
        specs = root_specs + (cid, job)

        if not SFS.tau_run.dir.exists(run_prefix, specs):
            run_path = SFS.tau_run.dir.path(run_prefix, specs)
            print("Starting run {}/{} at {}".format(idx+1, nsamp, run_path))

            status = "running"
            run_inf_obj = autofile.system.info.run(
                job=job, prog=prog, method=method, basis=basis, status=status)
            run_inf_obj.utc_start_time = autofile.system.info.utc_time()

            SFS.tau_run.dir.create(run_prefix, specs)
            SFS.tau_run.file.info.write(run_inf_obj, run_prefix, specs)

#            inp_str, out_str = feedback_optimization(
#                script_str, run_path,
#                geom=inp_zma, charge=charge, mult=mult, method=method,
#                basis=basis, prog=prog, frozen_coordinates=tors_names ,
#                **kwargs)
#            inp_str,out_str = robust_feedback_opt(
            options_mat = [
                [{'scf_options': (elstruct.Option.Scf.Guess.CORE,)},
                 {'scf_options': (elstruct.Option.Scf.Guess.HUCKEL,)},
                 {'scf_options': (
                     elstruct.option.specify(elstruct.Option.Scf.DIIS_, True),
                     elstruct.Option.Scf.Guess.HUCKEL,)}],
                [{'job_options': (elstruct.Option.Opt.Coord.ZMATRIX,)},
                 {'job_options': (
                     elstruct.option.specify(elstruct.Option.Opt.MAXITER_, 10),
                     elstruct.Option.Opt.Coord.ZMATRIX,)}],
            ]
            errors = [
                elstruct.Error.SCF_NOCONV,
                elstruct.Error.OPT_NOCONV,
            ]
            inp_str, out_str = robust_feedback_opt(
                script_str, run_path,
                geom=inp_zma, charge=charge, mult=mult, method=method,
                basis=basis, prog=prog, frozen_coordinates=tors_names, 
                errors=errors, options_mat=options_mat, **kwargs)

            run_inf_obj.utc_end_time = autofile.system.info.utc_time()

            if elstruct.reader.has_normal_exit_message(prog, out_str):
                SFS.tau_run.file.output.write(out_str, run_prefix, specs)
                status = "succeeded"
            else:
                status = "failed"
            run_inf_obj.status = status
            SFS.tau_run.file.info.write(run_inf_obj, run_prefix, specs)
            SFS.tau_run.file.input.write(inp_str, run_prefix, specs)
            print("Run {}/{} {} at {}".format(idx+1, nsamp, status, run_path))


def save_tau(ich, charge, mult, method, basis, orb_restricted,
             run_prefix, save_prefix):
    """ save the taus that have been found so far
    """
    root_specs = (ich, charge, mult, method, basis, orb_restricted)
    run_tau_specs_lst = SFS.tau.dir.existing(run_prefix, root_specs)
    saved_tau_specs_lst = SFS.tau.dir.existing(save_prefix, root_specs)

    tau_specs_lst = []
    ene_lst = []
    geo_lst = []
    inp_str_lst = []
    inf_obj_lst = []

    print()
    print("Reading optimizations from run directories.")
    print(root_specs)
    print (run_tau_specs_lst)
    run_specs = ('optimization',)
    for tau_specs in run_tau_specs_lst:
        specs = root_specs + tau_specs + run_specs

        run_path = SFS.tau_run.dir.path(run_prefix, specs)
        print("Reading from run at {}".format(run_path))

        if SFS.tau_run.file.output.exists(run_prefix, specs):
            inf_obj = SFS.tau_run.file.info.read(run_prefix, specs)
            inp_str = SFS.tau_run.file.input.read(run_prefix, specs)
            out_str = SFS.tau_run.file.output.read(run_prefix, specs)
            prog = inf_obj.prog
            if elstruct.reader.has_normal_exit_message(prog, out_str):
                ene = elstruct.reader.energy(prog, method, out_str)
                geo = elstruct.reader.opt_geometry(prog, out_str)

            print(automol.geom.coulomb_spectrum(geo))
            save_specs = root_specs + tau_specs
            save_path = SFS.tau.dir.path(save_prefix, save_specs)
            print("Saving values from run at {}".format(save_path))

            SFS.tau.dir.create(save_prefix, save_specs)
            SFS.tau.file.geometry_info.write(inf_obj, save_prefix, save_specs)
            SFS.tau.file.geometry_input.write(inp_str, save_prefix, save_specs)
            SFS.tau.file.energy.write(ene, save_prefix, save_specs)
            SFS.tau.file.geometry.write(geo, save_prefix, save_specs)

    # update the number of samples
    nsamp_new = len(tau_specs_lst)
    trunk_inf_obj = SFS.tau_trunk.file.info.read(save_prefix, root_specs)
    trunk_inf_obj.nsamp += nsamp_new
    SFS.tau_trunk.file.info.write(trunk_inf_obj, save_prefix, root_specs)


def run_tau_hessian(ich, charge, mult, method, basis, orb_restricted, 
                    run_prefix, save_prefix, script_str, prog, vignore=1e10,
                    **kwargs):
    """ save the taus that have been found so far
    """
    root_specs = (ich, charge, mult, method, basis, orb_restricted)
    saved_tau_specs_lst = SFS.tau.dir.existing(save_prefix, root_specs)

    tau_specs_lst = []
    ene_lst = []
    geo_lst = []
    hessian_lst = []
    inp_str_lst = []
    inf_obj_lst = []

    print()
    print("Reading optimizations from run directories as prelude to hessians.")
    print(root_specs)
    run_specs = ('hessian',)

    for tau_specs in saved_tau_specs_lst:
        specs = root_specs + tau_specs
        geo = SFS.tau.file.geometry.read(save_prefix, specs)
        ene = SFS.tau.file.energy.read(save_prefix, specs)
        if ene < vignore:
            geo_lst.append(geo)
            ene_lst.append(ene)
            tau_specs_lst.append(tau_specs)
    for idx, (tau_specs, geo) in enumerate(zip(tau_specs_lst, geo_lst)):
        specs = root_specs + tau_specs + run_specs 
        run_path = SFS.tau_run.dir.path(run_prefix, specs)
        if not SFS.tau_run.dir.exists(run_prefix, specs):
            do_run = True
        else:
            run_inf_obj = SFS.tau_run.file.info.read(run_prefix, specs) 
            do_run = run_inf_obj.status == "failed"
            if do_run:
                print ('Found failed hessian at {}'.format(run_path))
                SFS.tau_run.dir.remove (run_prefix,specs)
            else:
                print ('Skipping hessian at {}'.format(run_path))

        if do_run:
            run_path = SFS.tau_run.dir.path(run_prefix, specs)
            print("Starting run {} at {}".format(idx+1, run_path))

            status = "running"
            run_inf_obj = autofile.system.info.run(
                job='hessian', prog=prog, method=method, basis=basis, status=status)
            run_inf_obj.utc_start_time = autofile.system.info.utc_time()

            SFS.tau_run.dir.create(run_prefix, specs)
            SFS.tau_run.file.info.write(run_inf_obj, run_prefix, specs)

            options_mat = [
                [{'scf_options': (elstruct.Option.Scf.Guess.CORE,)},
                 {'scf_options': (elstruct.Option.Scf.Guess.HUCKEL,)},
                 {'scf_options': (
                     elstruct.option.specify(elstruct.Option.Scf.DIIS_, True),
                     elstruct.Option.Scf.Guess.HUCKEL,)}],
            ]
            errors = [
                elstruct.Error.SCF_NOCONV,
            ]
            inp_str, out_str = robust_run(
                elstruct.writer.hessian,
                script_str, run_path,
                geom=geo, charge=charge, mult=mult, method=method,
                basis=basis, prog=prog,
                errors=errors, options_mat=options_mat, **kwargs)

            run_inf_obj.utc_end_time = autofile.system.info.utc_time()

            if elstruct.reader.has_normal_exit_message(prog, out_str):
                SFS.tau_run.file.output.write(out_str, run_prefix, specs)
                status = "succeeded"
            else:
                status = "failed"

            run_inf_obj.status = status
            SFS.tau_run.file.info.write(run_inf_obj, run_prefix, specs)
            SFS.tau_run.file.input.write(inp_str, run_prefix, specs)
            print("Run {} {} at {}".format(idx+1, status, run_path))


def robust_run(input_writer, script_str, run_dir,
               geom, charge, mult, method, basis, prog,
               errors=(), options_mat=(),
               **kwargs):
    """ try several sets of options to generate an output file
    :returns: the input string, the output string, and the run directory
    :rtype: (str, str, str)
    """
    assert len(errors) == len(options_mat)

    try_idx = 0
    kwargs_dct = dict(kwargs)
    while not optsmat.is_exhausted(options_mat):
        if try_idx > 0:
            kwargs_dct = optsmat.updated_kwargs(kwargs, options_mat)
        try_dir_name = 'try{:d}'.format(try_idx)
        try_dir_path = os.path.join(run_dir, try_dir_name)
        assert not os.path.exists(try_dir_path)
        os.mkdir(try_dir_path)

        # filter out the warnings from the trial runs
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            input_str, output_str = elstruct.run.direct(
                input_writer, script_str, try_dir_path,
                geom=geom, charge=charge, mult=mult, method=method,
                basis=basis, prog=prog, **kwargs_dct)

        error_vals = [
            elstruct.reader.has_error_message(prog, error, output_str)
            for error in errors]

        if not any(error_vals):
            break

        try_idx += 1
        row_idx = error_vals.index(True)
        options_mat = optsmat.advance(row_idx, options_mat)

    if (any(error_vals) or not
            elstruct.reader.has_normal_exit_message(prog, output_str)):
        warnings.resetwarnings()
        warnings.warn("elstruct robust run failed; last run was in {}"
                      .format(run_dir))

    return input_str, output_str


def feedback_optimization(script_str, run_dir,
                          geom, charge, mult, method, basis, prog,
                          ntries=3, **kwargs):
    """ retry an optimization from the last (unoptimized) structure
    """
    assert automol.geom.is_valid(geom) or automol.zmatrix.is_valid(geom)
    is_zmat = automol.zmatrix.is_valid(geom)
    read_geom_ = (elstruct.reader.opt_geometry_(prog) if not is_zmat else
                  elstruct.reader.opt_zmatrix_(prog))
    has_noconv_error_ = functools.partial(
        elstruct.reader.has_error_message, prog, elstruct.Error.OPT_NOCONV)

    for try_idx in range(ntries):
        try_dir_name = 'try{:d}'.format(try_idx)
        try_dir_path = os.path.join(run_dir, try_dir_name)
        assert not os.path.exists(try_dir_path)
        os.mkdir(try_dir_path)

        # filter out the warnings from the trial runs
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            input_str, output_str = elstruct.run.direct(
                elstruct.writer.optimization, script_str, try_dir_path,
                geom=geom, charge=charge, mult=mult, method=method,
                basis=basis, prog=prog, **kwargs)

        if has_noconv_error_(output_str):
            geom = read_geom_(output_str)
        else:
            break

    if has_noconv_error_(output_str):
        warnings.resetwarnings()
        warnings.warn("elstruct feedback optimization failed; "
                      "last try was in {}".format(run_dir))

    return input_str, output_str


def robust_feedback_opt(script_str, run_dir,
                        geom, charge, mult, method, basis, prog,
                        errors=(), options_mat=(), ntries=3, **kwargs):
    """ try several sets of options to generate an output file
        retry an optimization from the last (unoptimized) structure
    """
    assert automol.geom.is_valid(geom) or automol.zmatrix.is_valid(geom)
    is_zmat = automol.zmatrix.is_valid(geom)
    read_geom_ = (elstruct.reader.opt_geometry_(prog) if not is_zmat else
                  elstruct.reader.opt_zmatrix_(prog))
    has_no_opt_conv_error_ = functools.partial(
        elstruct.reader.has_error_message, prog, elstruct.Error.OPT_NOCONV)

    for try_idx in range(ntries):
        try_dir_name = 'try{:d}'.format(try_idx)
        try_dir_path = os.path.join(run_dir, try_dir_name)
        assert not os.path.exists(try_dir_path)
        os.mkdir(try_dir_path)

        # filter out the warnings from the trial runs
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            input_str, output_str = robust_run(
                elstruct.writer.optimization, script_str, try_dir_path,
                geom=geom, charge=charge, mult=mult, method=method,
                basis=basis, prog=prog, errors=errors, options_mat=options_mat,
                **kwargs)

        if has_no_opt_conv_error_(output_str):
            geom = read_geom_(output_str)
        else:
            break

    if has_no_opt_conv_error_(output_str):
        warnings.resetwarnings()
        warnings.warn("elstruct feedback optimization failed; "
                      "last try was in {}".format(run_dir))

    return input_str, output_str
