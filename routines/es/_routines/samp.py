"""
  Sampling routines
"""

import automol


def torsional_sampling(constraints=('tors', 'none')):
    """ Run a torsional sampling procedure
    """


def _set_coordinates():
    """ Set up the coordinates to be frozen
    """

def _generate_samp_zma(zma, tors_range_dct, nsamp):
    """ Generate a sample ZMA
    """
    if nsampd > 0:
        samp_zma, = automol.zmatrix.samples(zma, 1, tors_range_dct)
    else:
        samp_zma = zma

    return samp_zma


def _set_tors_range(tors_name_grps):
    """ Set the range for all of the torsions
    """
    tors_ranges = automol.zmatrix.torsional_sampling_ranges(tors_name_grps)
    tors_range_dct = dict(zip(
        tuple(grp[0] for grp in tors_name_grps), tors_ranges))

    return tors_range_dct


def _set_nsamp(nsamp_par, ich, tors_range_dct, saddle):
    """ Determine the number of samples to run
    """

    if tors_range_dct:

        # Generate nsamp using num_tors and parameters
        if not saddle:
            gra = automol.inchi.graph(ich)
            ntaudof = len(
                automol.graph.rotational_bond_keys(gra, with_h_rotors=False))
        else:
            ntaudof = len((tors_range_dct.keys()))

        # Set nsamp
        if nsamp_par[0]:
            nsamp = min(nsamp_par[1] + nsamp_par[2] * nsamp_par[3]**ntaudof,
                        nsamp_par[4])
        else:
            nsamp = nsamp_par[5]

    else:
        print("No torsional coordinates. Setting nsamp to 1.")
        nsamp = 1

    return nsamp


def _find_nsampd(save_fs, run_fs):
    """ Determine the number of coordinates that have been sampled
        using info files from the filesystem
    """

    if save_fs[0].file.info.exists():
        inf_obj_s = save_fs[0].file.info.read()
        nsampd = inf_obj_s.nsamp
    elif run_fs[0].file.info.exists():
        inf_obj_r = run_fs[0].file.info.read()
        nsampd = inf_obj_r.nsamp
    else:
        nsampd = 0

    return nsampd
