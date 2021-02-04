""" conformer new functions
"""

def _calc_nsampd(cnf_save_fs, cnf_run_fs):
    """ Determine the number of samples completed
    """

    cnf_save_fs[0].create()
    inf_obj = autofile.schema.info_objects.conformer_trunk(0)
    if cnf_save_fs[0].file.info.exists():
        inf_obj_s = cnf_save_fs[0].file.info.read()
        nsampd = inf_obj_s.nsamp
    elif cnf_run_fs[0].file.info.exists():
        inf_obj_r = cnf_run_fs[0].file.info.read()
        nsampd = inf_obj_r.nsamp
    else:
        nsampd = 0

    return nsampd
