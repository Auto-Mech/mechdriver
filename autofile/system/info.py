""" Info objects
"""
import numbers
import numpy
import autofile.info
from autofile.system._util import utc_time as _utc_time


def conformer_trunk(nsamp, tors_ranges):
    """ conformer trunk information

    :param nsamp: the number of samples
    :type nsamp: int
    :param tors_ranges: sampling ranges [(start, end)] for each torsional
        coordinate, by z-matrix coordinate name
    :type tors_ranges: dict[str: (float, float)]
    """
    tors_range_dct = dict(tors_ranges)
    for key, rng in tors_range_dct.items():
        tors_range_dct[key] = (rng[0]*180./numpy.pi, rng[1]*180./numpy.pi)

    assert all(isinstance(key, str) and len(rng) == 2
               and all(isinstance(x, numbers.Real) for x in rng)
               for key, rng in tors_range_dct.items())

    tors_ranges = autofile.info.Info(**tors_range_dct)
    assert isinstance(nsamp, numbers.Integral)
    inf_obj = autofile.info.Info(nsamp=nsamp, tors_ranges=tors_ranges)
    assert autofile.info.matches_function_signature(inf_obj, conformer_trunk)
    return inf_obj


def tau_trunk(nsamp, tors_ranges):
    """ tau trunk information

    :param nsamp: the number of samples
    :type nsamp: int
    :param tors_ranges: sampling ranges [(start, end)] for each torsional
        coordinate, by z-matrix coordinate name
    :type tors_ranges: dict[str: (float, float)]
    """
    tors_range_dct = dict(tors_ranges)
    for key, rng in tors_range_dct.items():
        tors_range_dct[key] = (rng[0]*180./numpy.pi, rng[1]*180./numpy.pi)

    assert all(isinstance(key, str) and len(rng) == 2
               and all(isinstance(x, numbers.Real) for x in rng)
               for key, rng in tors_range_dct.items())

    tors_ranges = autofile.info.Info(**tors_range_dct)
    assert isinstance(nsamp, numbers.Integral)
    inf_obj = autofile.info.Info(nsamp=nsamp, tors_ranges=tors_ranges)
    assert autofile.info.matches_function_signature(inf_obj, tau_trunk)
    return inf_obj


def scan_branch(grids):
    """ scan trunk information

    :param grids: sampling grids, [val1, val2, ...], for each coordinate,
        by coordinate name
    :type grids: dict[str: list[float]]
    """
    grid_dct = dict(grids)
    # note:renormalization of angle ranges needs to be updated for 2D grids.
    for key, rng in grid_dct.items():
        if 'R' not in key:
            grid_dct[key] = rng*180./numpy.pi

    assert all(isinstance(key, str) and numpy.ndim(vals) == 1
               and all(isinstance(x, numbers.Real) for x in vals)
               for key, vals in grid_dct.items())

    grids = autofile.info.Info(**grid_dct)
    inf_obj = autofile.info.Info(grids=grids)
    assert autofile.info.matches_function_signature(inf_obj, scan_branch)
    return inf_obj


def vpt2_trunk(fermi):
    """ vpt2 trunk information

    :param fermi: description of fermi resonance treatment
    :type fermi: str
    """

    assert isinstance(fermi, str)
    inf_obj = autofile.info.Info(fermi=fermi)
    assert autofile.info.matches_function_signature(inf_obj, vpt2_trunk)
    return inf_obj


def lennard_jones(potential, nsamp,
                  method, basis, program, version):
    """ energy transfer trunk """
    inf_obj = autofile.info.Info(potential=potential, nsamp=nsamp,
                                 method=method, basis=basis,
                                 program=program, version=version)
    assert autofile.info.matches_function_signature(
        inf_obj, lennard_jones)
    return inf_obj


class RunStatus():
    """ run statuses """
    RUNNING = "running"
    SUCCESS = "succeeded"
    FAILURE = "failed"


def run(job, prog, version, method, basis, status, utc_start_time=None,
        utc_end_time=None):
    """ run information
    """
    inf_obj = autofile.info.Info(
        job=job,
        prog=prog,
        version=version,
        method=method,
        basis=basis,
        status=status,
        utc_start_time=utc_start_time,
        utc_end_time=utc_end_time,
    )
    assert autofile.info.matches_function_signature(inf_obj, run)
    return inf_obj


def utc_time():
    """ current run time
    """
    return _utc_time()
