""" Info objects
"""
import numbers
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


def scan_branch(tors_linspaces):
    """ scan trunk information

    :param tors_linspaces: sampling linspaces [(start, end, num)] for each
        torsional coordinate, by z-matrix coordinate name
    :type tors_linspaces: dict[str: (float, float, int)]
    """
    tors_linspace_dct = dict(tors_linspaces)

    assert all(isinstance(key, str) and len(lsp) == 3
               and all(isinstance(x, numbers.Real) for x in lsp[:2])
               and isinstance(lsp[2], numbers.Integral)
               for key, lsp in tors_linspace_dct.items())

    tors_linspaces = autofile.info.Info(**tors_linspace_dct)
    inf_obj = autofile.info.Info(tors_linspaces=tors_linspaces)
    assert autofile.info.matches_function_signature(inf_obj, scan_branch)
    return inf_obj


class RunStatus():
    """ run statuses """
    RUNNING = "running"
    SUCCESS = "succeeded"
    FAILURE = "failed"


def run(job, prog, method, basis, status, utc_start_time=None,
        utc_end_time=None):
    """ run information
    """
    inf_obj = autofile.info.Info(
        job=job,
        prog=prog,
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
