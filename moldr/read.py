""" readers
"""
import os
import autodir


def values(run_prefixes, run_name, output_readers):
    """ read values from each output file
    """
    nruns = len(run_prefixes)
    inp_strs = []
    inf_objs = []
    vals_lst = []
    for idx, run_prefix in enumerate(run_prefixes):
        run_path = autodir.run.directory_path(run_prefix, run_name)
        run_path = os.path.abspath(run_path)
        print("Reading from run {}/{} at {}".format(idx+1, nruns, run_path))

        inp_str = autodir.run.read_input_file(run_prefix, run_name)
        inf_obj = autodir.run.read_information_file(run_prefix, run_name)
        out_str = autodir.run.read_output_file(run_prefix, run_name)
        vals = tuple(output_reader(out_str)
                     for output_reader in output_readers)
        inp_strs.append(inp_str)
        inf_objs.append(inf_obj)
        vals_lst.append(vals)

    inp_strs = tuple(inp_strs)
    inf_objs = tuple(inf_objs)
    vals_lst = tuple(vals_lst)
    return inp_strs, inf_objs, vals_lst
