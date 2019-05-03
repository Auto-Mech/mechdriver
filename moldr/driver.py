""" drivers
"""
import automol
import elstruct
import elcarro
import autofile
from autofile import fs


def run_conformers(ich, mult, method, basis, orb_restricted,
                   run_prefix, save_prefix, nsamp, script_str, prog,
                   **kwargs):
    """ run sampling algorithm to find conformers
    """
    geo = automol.inchi.geometry(ich)
    zma = automol.geom.zmatrix(geo)
    tors_names = automol.geom.zmatrix_torsion_coordinate_names(geo)
    tors_range_vals = automol.zmatrix.tors.sampling_ranges(zma, tors_names)
    tors_ranges = dict(zip(tors_names, tors_range_vals))

    if not tors_ranges:
        print("No torsional coordinates. Setting nsamp to 1.")
        nsamp = 1

    root_specs = (ich, mult, method, basis, orb_restricted)

    # check for a previously saved run
    vma = automol.zmatrix.var_(zma)
    if fs.conformer_trunk.dir.exists(save_prefix, root_specs):
        _vma = fs.conformer_trunk.file.vmatrix.read(save_prefix, root_specs)
        assert vma == _vma
        inf_obj = fs.conformer_trunk.file.info.read(save_prefix, root_specs)
        nsamp = max(nsamp - inf_obj.nsamp, 0)
        print("Found previous saved run. Adjusting nsamp.")
        print("    New nsamp is {:d}.".format(nsamp))
    else:
        fs.conformer_trunk.dir.create(save_prefix, root_specs)
        inf_obj = autofile.system.info.conformer_trunk(
            nsamp=0, tors_ranges=tors_ranges)
    fs.conformer_trunk.file.vmatrix.write(vma, save_prefix, root_specs)
    fs.conformer_trunk.file.info.write(inf_obj, save_prefix, root_specs)

    # update the number of samples

    inp_zmas = automol.zmatrix.tors.samples(zma, nsamp, tors_ranges)

    cids = tuple(autofile.system.generate_new_conformer_id()
                 for _ in range(nsamp))

    print()
    print("Running optimizations in run directories.")
    job = 'optimization'
    for idx, (cid, inp_zma) in enumerate(zip(cids, inp_zmas)):
        specs = root_specs + (cid, job)

        if not fs.conformer_run.dir.exists(run_prefix, specs):
            run_path = fs.conformer_run.dir.path(run_prefix, specs)
            print("Starting run {}/{} at {}".format(idx+1, nsamp, run_path))

            run_inf_obj = autofile.system.info.run(
                job=job, prog=prog, method=method, basis=basis)
            run_inf_obj.utc_start_time = autofile.system.info.utc_time()

            fs.conformer_run.dir.create(run_prefix, specs)
            fs.conformer_run.file.info.write(run_inf_obj, run_prefix, specs)

            inp_str, out_str = elcarro.feedback_optimization(
                script_str, run_path,
                prog, method, basis, inp_zma, mult, charge=0,
                **kwargs)

            run_inf_obj.utc_end_time = autofile.system.info.utc_time()

            fs.conformer_run.file.info.write(run_inf_obj, run_prefix, specs)
            fs.conformer_run.file.input.write(inp_str, run_prefix, specs)

            status = "succeeded"
            if elstruct.reader.has_normal_exit_message(prog, out_str):
                fs.conformer_run.file.output.write(out_str, run_prefix, specs)
                status = "failed"

            print("Run {}/{} {} at {}".format(idx+1, nsamp, status, run_path))


def save_conformers(ich, mult, method, basis, orb_restricted,
                    run_prefix, save_prefix):
    """ save the conformers that have been found so far
    """
    root_specs = (ich, mult, method, basis, orb_restricted)
    run_conf_specs_lst = fs.conformer.dir.existing(run_prefix, root_specs)
    saved_conf_specs_lst = fs.conformer.dir.existing(save_prefix, root_specs)

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

        run_path = fs.conformer_run.dir.path(run_prefix, specs)
        print("Reading from run at {}".format(run_path))

        if fs.conformer_run.file.output.exists(run_prefix, specs):
            inf_obj = fs.conformer_run.file.info.read(run_prefix, specs)
            inp_str = fs.conformer_run.file.input.read(run_prefix, specs)
            out_str = fs.conformer_run.file.output.read(run_prefix, specs)
            prog = inf_obj.prog
            if elstruct.reader.has_normal_exit_message(prog, out_str):
                ene = elstruct.reader.energy(prog, method, out_str)
                geo = elstruct.reader.opt_geometry(prog, out_str)

                # save the information to a list
                conf_specs_lst.append(conf_specs)
                inf_obj_lst.append(inf_obj)
                inp_str_lst.append(inp_str)
                ene_lst.append(ene)
                geo_lst.append(geo)

    seen_geo_lst = []
    for conf_specs in saved_conf_specs_lst:
        specs = root_specs + conf_specs
        geo = fs.conformer.file.geometry.read(save_prefix, specs)
        seen_geo_lst.append(geo)

    print("Writing unique confomrer information to save directories.")
    idxs = automol.geom.argunique_coulomb_spectrum(
        geo_lst, seen_geos=seen_geo_lst, rtol=7e-2)
    for idx in idxs:
        conf_specs = conf_specs_lst[idx]
        inf_obj = inf_obj_lst[idx]
        inp_str = inp_str_lst[idx]
        ene = ene_lst[idx]
        geo = geo_lst[idx]

        specs = root_specs + conf_specs
        save_path = fs.conformer.dir.path(save_prefix, specs)
        print("Saving values from run at {}".format(save_path))

        fs.conformer.dir.create(save_prefix, specs)
        fs.conformer.file.geometry_info.write(inf_obj, save_prefix, specs)
        fs.conformer.file.geometry_input.write(inp_str, save_prefix, specs)
        fs.conformer.file.energy.write(ene, save_prefix, specs)
        fs.conformer.file.geometry.write(geo, save_prefix, specs)

    # update the number of samples
    nsamp_new = len(conf_specs_lst)
    trunk_inf_obj = fs.conformer_trunk.file.info.read(save_prefix, root_specs)
    trunk_inf_obj.nsamp += nsamp_new
    fs.conformer_trunk.file.info.write(trunk_inf_obj, save_prefix, root_specs)
