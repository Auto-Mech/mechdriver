""" utilites
"""

import numpy
import automol
from mechlib.amech_io.printer import info_message


def calc_nsamp(tors_names, nsamp_par, zma, zrxn=None):
    """ Determine the number of samples to od
    """
    ntaudof = None
    if not any(tors_names):
        info_message(
            " - No torsional coordinates. Setting nsamp to 1.")
        nsamp = 1
        tors_range_dct = None
    else:
        if zrxn is None:
            gra = automol.zmat.graph(zma)
        else:
            gra = automol.reac.ts_graph(zrxn)

        ntaudof = len(automol.graph.rotational_bond_keys(
            gra, with_chx_rotors=False))
        nsamp = nsamp_init(nsamp_par, ntaudof)

        tors_ranges = tuple((0, 2*numpy.pi) for tors in tors_names)
        tors_range_dct = dict(zip(tors_names, tors_ranges))

    return nsamp, tors_range_dct


def calc_nsampd(cnf_save_fs, cnf_run_fs, rid=None):
    """ Determine the number of samples completed
    """

    if rid is None:
        cnf_save_fs[0].create()
        if cnf_save_fs[0].file.info.exists():
            inf_obj_s = cnf_save_fs[0].file.info.read()
            nsampd = inf_obj_s.nsamp
        elif cnf_run_fs[0].file.info.exists():
            inf_obj_r = cnf_run_fs[0].file.info.read()
            nsampd = inf_obj_r.nsamp
        else:
            nsampd = 0
    else:
        cnf_save_fs[1].create([rid])
        if cnf_save_fs[1].file.info.exists([rid]):
            inf_obj_s = cnf_save_fs[1].file.info.read([rid])
            nsampd = inf_obj_s.nsamp
        elif cnf_run_fs[1].file.info.exists([rid]):
            inf_obj_r = cnf_run_fs[1].file.info.read([rid])
            nsampd = inf_obj_r.nsamp
        else:
            nsampd = 0

    return nsampd


def nsamp_init(nsamp_par, ntaudof):
    """ determine nsamp for given species"""
    if nsamp_par[0]:
        nsamp = min(nsamp_par[1] + nsamp_par[2] * nsamp_par[3]**ntaudof,
                    nsamp_par[4])
        # print('Setting nsamp using formula: min(A+B*C**n')
    else:
        nsamp = nsamp_par[5]
    return nsamp
