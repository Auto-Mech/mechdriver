""" runners
"""
import os
import autodir
import elstruct


def geometries(geoms, run_prefixes, run_name, job_name, runner_,
               # runner arguments
               script_str, prog, method, basis, mult, charge,
               **kwargs):
    """ run jobs for multiple geometries of a single species
    """
    nruns = len(run_prefixes)
    assert len(geoms) == nruns
    for idx, (geom, run_prefix) in enumerate(zip(geoms, run_prefixes)):
        run_path = autodir.run.directory_path(run_prefix, run_name)
        run_path = os.path.abspath(run_path)
        autodir.run.create(run_prefix, run_name)
        print("Starting run {}/{} at {}".format(idx+1, nruns, run_path))

        inf_obj = autodir.run.information(
            job=job_name,
            prog=prog,
            method=method,
            basis=basis,
        )
        autodir.run.add_start_time_to_information(inf_obj)
        autodir.run.write_information_file(run_prefix, run_name, inf_obj)

        inp_str, out_str = runner_(
            script_str, run_path,
            prog, method, basis, geom, mult, charge,
            **kwargs)

        if elstruct.reader.has_normal_exit_message(prog, out_str):
            autodir.run.add_end_time_to_information(inf_obj)
            autodir.run.write_information_file(run_prefix, run_name, inf_obj)
            autodir.run.write_input_file(run_prefix, run_name, inp_str)
            autodir.run.write_output_file(run_prefix, run_name, out_str)
