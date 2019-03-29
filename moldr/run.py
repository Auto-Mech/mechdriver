""" run functions
"""
import automol
import autodir
import autoinf
import elstruct


def torsion_sampling_information(nsamp):
    """ torsion sampling information object
    """
    inf_obj = autoinf.Info(nsamp=nsamp)
    assert autoinf.matches_function_signature(inf_obj,
                                              torsion_sampling_information)
    return inf_obj


def torsion_sampling(nsamp,
                     # elstruct robust run arguments
                     script_str, input_writer, prefix,
                     prog, method, basis, geo, mult, charge,
                     errors, options_mat,
                     **kwargs):
    """ sample torsions
    """
    autodir.run.create_base(prefix)

    # for now, require cartesian inputs
    assert automol.geom.is_valid(geo)
    ich = automol.geom.inchi(geo)

    # write/update the information file
    autodir.run.create_base(prefix)
    complete = False
    base_inf_obj = autodir.run.base_information(
        function=torsion_sampling.__name__,
        function_info=torsion_sampling_information(nsamp=nsamp),
        job=input_writer.__name__,
        prog=prog,
        method=method,
        basis=basis,
        inchi=ich,
        complete=complete,
    )
    autodir.run.write_base_information_file(prefix, base_inf_obj)
    autodir.run.write_run_script(prefix, script_str)

    # generate the sample z-matrices
    zma = automol.geom.zmatrix(geo)
    tors_names = automol.geom.zmatrix_torsion_coordinate_names(geo)
    zmas = automol.zmatrix.tors.samples(zma, nsamp, tors_names)

    # run the input files for each geometry
    nzmas = len(zmas)
    rids = tuple(autodir.id_.identifier() for _ in range(nzmas))

    for zma, rid in zip(zmas, rids):
        run_path = autodir.run.directory_path(prefix, rid)
        autodir.run.create(prefix, rid)
        inp_str, out_str = elstruct.run.robust(
            script_str, run_path, input_writer,
            prog, method, basis, zma, mult, charge,
            errors, options_mat,
            **kwargs)

        autodir.run.write_input_file(prefix, rid, inp_str)
        autodir.run.write_output_file(prefix, rid, out_str)

    # now that we're done, set the `complete` flag in the info file to True
    base_inf_obj.complete = True
    autodir.run.write_base_information_file(prefix, base_inf_obj)
