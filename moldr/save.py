""" save functions
"""
import os
import autodir
import elstruct
# import autoinf
import automol


def conformers(run_prefix, save_prefix):
    """ save conformers to a filesystem
    """
    # read in the run information to get the nsamp
    assert autodir.run.has_base_information_file(run_prefix)
    inf_obj = autodir.run.read_base_information_file(run_prefix)
    assert inf_obj.complete

    prog = inf_obj.prog
    method = inf_obj.method
    nsamp_new = inf_obj.function_info.nsamp
    nsamp_old = 0

    # read in the save information to get the original nsamp
    if autodir.conf.has_base_information_file(save_prefix):
        inf_obj = autodir.conf.read_base_information_file(save_prefix)
        nsamp_old = inf_obj.nsamp

    nsamp_tot = nsamp_old + nsamp_new

    if nsamp_new > 0:
        rids, inp_strs, out_strs, (enes, geos) = read_run(
            prefix=run_prefix,
            output_readers=(
                elstruct.reader.energy_(prog, method),
                elstruct.reader.opt_geometry_(prog),
            ),
        )

        # grab only the unique ones
        idxs = automol.geom.argunique_coulomb_spectrum(geos, rtol=7e-2)

        rids = list(map(rids.__getitem__, idxs))
        inp_strs = list(map(inp_strs.__getitem__, idxs))
        out_strs = list(map(out_strs.__getitem__, idxs))
        enes = list(map(enes.__getitem__, idxs))
        geos = list(map(geos.__getitem__, idxs))

        autodir.conf.create_base(save_prefix)
        base_inf_obj = autodir.conf.base_information(nsamp=nsamp_tot)
        autodir.conf.write_base_information_file(save_prefix, base_inf_obj)

        for rid, ene, geo in zip(rids, enes, geos):
            autodir.conf.create(save_prefix, rid)
            autodir.conf.write_energy_file(save_prefix, rid, ene)
            autodir.conf.write_geometry_file(save_prefix, rid, geo)


def read_run(prefix, output_readers):
    """ read several values from each output in a run directory
    """
    base_pth = autodir.run.base_path(prefix)
    assert os.path.exists(base_pth)
    rids = autodir.id_.directory_identifiers_at(base_pth)

    inp_strs = []
    out_strs = []
    reads_lst = []
    for rid in rids:
        inp_str = autodir.run.read_input_file(prefix, rid)
        out_str = autodir.run.read_output_file(prefix, rid)
        out_vals = tuple(output_reader(out_str)
                         for output_reader in output_readers)
        inp_strs.append(inp_str)
        out_strs.append(out_str)
        reads_lst.append(out_vals)

    inp_strs = tuple(inp_strs)
    out_strs = tuple(out_strs)
    vals_lst = tuple(zip(*reads_lst))

    return rids, inp_strs, out_strs, vals_lst
