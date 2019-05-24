""" molecule driver routines
"""
import elstruct
import automol
import autodir
import moldr.run
import moldr.read


def conformers(nsamp, script_str, run_prefix, save_prefix,
               # elstruct robust run arguments
               prog, method, basis, geo, mult, charge,
               errors=(), options_mat=(),
               **kwargs):
    """ optimize from randomly sampled torions to find unique conformers
    """
    nsamp_tot = nsamp
    seen_geos = ()

    if autodir.conf.has_base_information_file(save_prefix):
        inf_obj = autodir.conf.read_base_information_file(save_prefix)
        nsamp = max(nsamp_tot - inf_obj.nsamp, 0)
        nsamp_tot = max(nsamp_tot, inf_obj.nsamp)
        print("Found previous saved run. Adjusting `nsamp`.")
        print("nsamp done={:d}".format(inf_obj.nsamp))
        print("nsamp left={:d}".format(nsamp))

        rids = autodir.conf.identifiers(save_prefix)
        seen_geos = [
            autodir.conf.read_geometry_file(save_prefix, rid) for rid in rids]

    if nsamp == 0:
        print("Updating trajectory file at {}"
              .format(autodir.conf.trajectory_file_path(save_prefix)))
        autodir.conf.update_trajectory_file(save_prefix)
        return

    # generate the sample z-matrices
    zma = automol.geom.zmatrix(geo)
    print (zma)
    tors_names = automol.geom.zmatrix_torsion_coordinate_names(geo)
    print (tors_names)
    inp_zmas = automol.zmatrix.tors.samples(zma, nsamp, tors_names)

    # generate random ids for each
    rids = tuple(autodir.id_.identifier() for _ in range(nsamp))

    # run the optimizations
    run_prefixes = [
        autodir.conf.run_directory_path(run_prefix, rid) for rid in rids]
    for rid in rids:
        autodir.conf.create_run_directory(run_prefix, rid)

    moldr.run.geometries(
        geoms=inp_zmas,
        run_prefixes=run_prefixes,
        run_name=autodir.conf.OPT_RUN_NAME,
        script_str=script_str,
        input_writer=elstruct.writer.optimization,
        prog=prog,
        method=method,
        basis=basis,
        mult=mult,
        charge=charge,
        errors=errors,
        options_mat=options_mat,
        **kwargs,
    )

    # read out the optimized z-matrices and energies
    inp_strs, inf_objs, vals_lst = moldr.read.values(
        run_prefixes=run_prefixes,
        run_name=autodir.conf.OPT_RUN_NAME,
        output_readers=(
            elstruct.reader.energy_(prog, method),
            elstruct.reader.opt_zmatrix_(prog),
        ),
    )
    enes, zmas = zip(*vals_lst)
    geos = tuple(map(automol.zmatrix.geometry, zmas))

    # filter out the non-unique geometries
    filter_idxs = automol.geom.argunique_coulomb_spectrum(
        geos, seen_geos=seen_geos, rtol=7e-2)
    nuniq = len(filter_idxs)

    # save them
    autodir.conf.create_base(save_prefix)
    base_inf_obj = autodir.conf.base_information(nsamp=nsamp_tot)
    autodir.conf.write_base_information_file(save_prefix, base_inf_obj)

    if not filter_idxs:
        print("No new conformers found.")

    for idx, filter_idx in enumerate(filter_idxs):
        rid = rids[filter_idx]
        dir_path = autodir.conf.directory_path(save_prefix, rid)
        print("Saving new conformer {}/{} at {}"
              .format(idx+1, nuniq, dir_path))

        autodir.conf.create(save_prefix, rid)
        inp_str = inp_strs[filter_idx]
        inf_obj = inf_objs[filter_idx]
        ene = enes[filter_idx]
        geo = geos[filter_idx]

        autodir.conf.write_input_file(save_prefix, rid, inp_str)
        autodir.conf.write_information_file(save_prefix, rid, inf_obj)
        autodir.conf.write_energy_file(save_prefix, rid, ene)
        autodir.conf.write_geometry_file(save_prefix, rid, geo)

    print("Updating trajectory file at {}"
          .format(autodir.conf.trajectory_file_path(save_prefix)))
    autodir.conf.update_trajectory_file(save_prefix)


def add_conformer_gradients(script_str, run_prefix, save_prefix,
                            # elstruct robust run arguments
                            prog, method, basis, mult, charge,
                            errors=(), options_mat=(),
                            **kwargs):
    """ determine gradients for conformers
    """
    # filter out the ones that already have gradients
    rids = [rid for rid in autodir.conf.identifiers(save_prefix)
            if not autodir.conf.has_gradient_file(save_prefix, rid)]

    if rids:
        geos = [autodir.conf.read_geometry_file(save_prefix, rid)
                for rid in rids]

        # run the gradient calculations
        run_prefixes = [
            autodir.conf.run_directory_path(run_prefix, rid) for rid in rids]

        moldr.run.geometries(
            geoms=geos,
            run_prefixes=run_prefixes,
            run_name=autodir.conf.GRAD_RUN_NAME,
            script_str=script_str,
            input_writer=elstruct.writer.gradient,
            prog=prog,
            method=method,
            basis=basis,
            mult=mult,
            charge=charge,
            errors=errors,
            options_mat=options_mat,
            **kwargs,
        )

        # save the gradient values
        inp_strs, inf_objs, vals_lst = moldr.read.values(
            run_prefixes=run_prefixes,
            run_name=autodir.conf.GRAD_RUN_NAME,
            output_readers=(
                elstruct.reader.gradient_(prog),
            ),
        )

        grads, = zip(*vals_lst)
        ngrad = len(grads)

        for idx, (rid, inp_str, inf_obj, grad) in enumerate(
                zip(rids, inp_strs, inf_objs, grads)):
            dir_path = autodir.conf.directory_path(save_prefix, rid)
            print("Saving gradient {}/{} at {}"
                  .format(idx+1, ngrad, dir_path))
            autodir.conf.write_gradient_information_file(
                save_prefix, rid, inf_obj)
            autodir.conf.write_gradient_input_file(save_prefix, rid, inp_str)
            autodir.conf.write_gradient_file(save_prefix, rid, grad)


def add_conformer_hessians(script_str, run_prefix, save_prefix,
                           # elstruct robust run arguments
                           prog, method, basis, mult, charge,
                           errors=(), options_mat=(),
                           **kwargs):
    """ determine hessians for conformers
    """
    # filter out the ones that already have hessians
    rids = [rid for rid in autodir.conf.identifiers(save_prefix)
            if not autodir.conf.has_hessian_file(save_prefix, rid)]

    if rids:
        geos = [autodir.conf.read_geometry_file(save_prefix, rid)
                for rid in rids]

        # run the hessian calculations
        run_prefixes = [
            autodir.conf.run_directory_path(run_prefix, rid) for rid in rids]

        moldr.run.geometries(
            geoms=geos,
            run_prefixes=run_prefixes,
            run_name=autodir.conf.HESS_RUN_NAME,
            script_str=script_str,
            input_writer=elstruct.writer.hessian,
            prog=prog,
            method=method,
            basis=basis,
            mult=mult,
            charge=charge,
            errors=errors,
            options_mat=options_mat,
            **kwargs,
        )

        # save the hessian values
        inp_strs, inf_objs, vals_lst = moldr.read.values(
            run_prefixes=run_prefixes,
            run_name=autodir.conf.HESS_RUN_NAME,
            output_readers=(
                elstruct.reader.hessian_(prog),
            ),
        )

        hess_lst, = zip(*vals_lst)
        nhess = len(hess_lst)

        for idx, (rid, inp_str, inf_obj, hess) in enumerate(
                zip(rids, inp_strs, inf_objs, hess_lst)):
            dir_path = autodir.conf.directory_path(save_prefix, rid)
            print("Saving hessian {}/{} at {}"
                  .format(idx+1, nhess, dir_path))
            autodir.conf.write_hessian_information_file(
                save_prefix, rid, inf_obj)
            autodir.conf.write_hessian_input_file(save_prefix, rid, inp_str)
            autodir.conf.write_hessian_file(save_prefix, rid, hess)
