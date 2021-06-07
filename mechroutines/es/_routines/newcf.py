""" conformer new functions
"""

# Samples functions
def calc_namp(tors_names, nsamp_par, cnf_save_fs, cnf_run_fs):
    """ Determine the number of samples to od
    """
    
    tors_ranges = tuple((0, 2*numpy.pi) for tors in tors_names)
    tors_range_dct = dict(zip(tors_names, tors_ranges))
    nsamp = util.nsamp_init(nsamp_par, len(tors_names))
    ioprinter.debug_message('tors_names', tors_names)
    ioprinter.debug_message('tors_range_dct', tors_range_dct)
    if not tors_range_dct:
        ioprinter.info_message(
            " - No torsional coordinates. Setting nsamp to 1.")
        nsamp = 1

    return nsamp


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
