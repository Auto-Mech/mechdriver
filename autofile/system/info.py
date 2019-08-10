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


class RunStatus():
    """ run statuses """
    RUNNING = "running"
    SUCCESS = "succeeded"
    FAILURE = "failed"


def run(job, prog, method, basis, status, utc_start_time=None,
#def run(job, theory_level, status, utc_start_time=None,
        utc_end_time=None):
    """ run information
    """
#    prog = theory_level[0]
#    method = theory_level[1]
#    basis = theory_level[2]
    inf_obj = autofile.info.Info(
        job=job,
        prog=prog,
        method=method,
        basis=basis,
        status=status,
        utc_start_time=utc_start_time,
        utc_end_time=utc_end_time,
    )
    print('inf_obj in run:', inf_obj)
    print('run in run:', run)
    assert autofile.info.matches_function_signature(inf_obj, run)
    return inf_obj


def utc_time():
    """ current run time
    """
    return _utc_time()
