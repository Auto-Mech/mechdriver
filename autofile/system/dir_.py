""" DataSeriesDirs
"""
try:
    from inspect import getfullargspec as function_argspec
except ImportError:
    from inspect import getargspec as function_argspec
import automol
from autofile.system import map_
from autofile.system import file_
from autofile.system import model


SPEC_FILE_PREFIX = 'dir'


def species_trunk(root_dsdir=None):
    """ species trunk DataSeriesDir
    """
    _map = _pack_arguments(map_.species_trunk)
    nspecs = _count_arguments(map_.species_trunk)
    return model.DataSeriesDir(map_=_map, nspecs=nspecs, depth=1,
                               root_dsdir=root_dsdir)


def species_leaf(root_dsdir=None):
    """ species leaf DataSeriesDir
    """
    spec_dfile = file_.data_series_specifier(
        file_prefix=SPEC_FILE_PREFIX,
        map_dct_={
            'inchi': lambda specs: specs[0],
            'charge': lambda specs: specs[1],
            'multiplicity': lambda specs: specs[2],
            'smiles': lambda specs: automol.inchi.smiles(specs[0])},
        spec_keys=['inchi', 'charge', 'multiplicity'])

    _map = _pack_arguments(map_.species_leaf)
    nspecs = _count_arguments(map_.species_leaf)
    return model.DataSeriesDir(map_=_map, nspecs=nspecs, depth=5,
                               spec_dfile=spec_dfile,
                               root_dsdir=root_dsdir)


def reaction_trunk(root_dsdir=None):
    """ reaction trunk DataSeriesDir
    """
    _map = _pack_arguments(map_.reaction_trunk)
    nspecs = _count_arguments(map_.reaction_trunk)
    return model.DataSeriesDir(map_=_map, nspecs=nspecs, depth=1,
                               root_dsdir=root_dsdir)


def reaction_leaf(root_dsdir=None):
    """ reaction leaf DataSeriesDir
    """
    spec_dfile = file_.data_series_specifier(
        file_prefix=SPEC_FILE_PREFIX,
        map_dct_={
            'inchis': lambda specs: specs[0],
            'charges': lambda specs: specs[1],
            'multiplicities': lambda specs: specs[2],
            'ts_multiplicity': lambda specs: specs[3],
            'smiles': lambda specs: [
                list(map(automol.inchi.smiles, specs[0][0])),
                list(map(automol.inchi.smiles, specs[0][1]))],
        },
        spec_keys=['inchis', 'charges', 'multiplicities', 'ts_multiplicity'])

    _map = _pack_arguments(map_.reaction_leaf)
    nspecs = _count_arguments(map_.reaction_leaf)
    return model.DataSeriesDir(map_=_map, nspecs=nspecs, depth=11,
                               spec_dfile=spec_dfile,
                               root_dsdir=root_dsdir)


def theory_leaf(root_dsdir=None):
    """ theory leaf DataSeriesDir
    """
    spec_dfile = file_.data_series_specifier(
        file_prefix=SPEC_FILE_PREFIX,
        map_dct_={
            'method': lambda specs: specs[0],
            'basis': lambda specs: specs[1],
            'orb_restricted': lambda specs: specs[2]},
        spec_keys=['method', 'basis', 'orb_restricted'])

    _map = _pack_arguments(map_.theory_leaf)
    nspecs = _count_arguments(map_.theory_leaf)
    return model.DataSeriesDir(map_=_map, nspecs=nspecs, depth=1,
                               spec_dfile=spec_dfile,
                               root_dsdir=root_dsdir)


def run_trunk(root_dsdir=None):
    """ run trunk DataSeriesDir
    """
    _map = _pack_arguments(map_.run_trunk)
    nspecs = _count_arguments(map_.run_trunk)
    return model.DataSeriesDir(map_=_map, nspecs=nspecs, depth=1,
                               root_dsdir=root_dsdir)


def run_leaf(root_dsdir=None):
    """ run leaf DataSeriesDir
    """
    spec_dfile = file_.data_series_specifier(
        file_prefix=SPEC_FILE_PREFIX,
        map_dct_={'job': lambda specs: specs[0]},
        spec_keys=['job'])

    _map = _pack_arguments(map_.run_leaf)
    nspecs = _count_arguments(map_.run_leaf)
    return model.DataSeriesDir(map_=_map, nspecs=nspecs, depth=1,
                               spec_dfile=spec_dfile,
                               root_dsdir=root_dsdir, removable=True)


def subrun_leaf(root_dsdir=None):
    """ subrun leaf DataSeriesDir
    """
    spec_dfile = file_.data_series_specifier(
        file_prefix=SPEC_FILE_PREFIX,
        map_dct_={'macro_idx': lambda specs: specs[0],
                  'micro_idx': lambda specs: specs[1]},
        spec_keys=['macro_idx', 'micro_idx'])

    _map = _pack_arguments(map_.subrun_leaf)
    nspecs = _count_arguments(map_.subrun_leaf)
    return model.DataSeriesDir(map_=_map, nspecs=nspecs, depth=1,
                               spec_dfile=spec_dfile,
                               root_dsdir=root_dsdir)


def conformer_trunk(root_dsdir=None):
    """ conformer trunk DataSeriesDir
    """
    _map = _pack_arguments(map_.conformer_trunk)
    nspecs = _count_arguments(map_.conformer_trunk)
    return model.DataSeriesDir(map_=_map, nspecs=nspecs, depth=1,
                               root_dsdir=root_dsdir)


def conformer_leaf(root_dsdir=None):
    """ conformer leaf DataSeriesDir
    """
    spec_dfile = file_.data_series_specifier(
        file_prefix=SPEC_FILE_PREFIX,
        map_dct_={'conformer_id': lambda specs: specs[0]},
        spec_keys=['conformer_id'])

    _map = _pack_arguments(map_.conformer_leaf)
    nspecs = _count_arguments(map_.conformer_leaf)
    return model.DataSeriesDir(map_=_map, nspecs=nspecs, depth=1,
                               spec_dfile=spec_dfile,
                               root_dsdir=root_dsdir)


def single_point_trunk(root_dsdir=None):
    """ single point trunk DataSeriesDir
    """
    _map = _pack_arguments(map_.single_point_trunk)
    nspecs = _count_arguments(map_.single_point_trunk)
    return model.DataSeriesDir(map_=_map, nspecs=nspecs, depth=1,
                               root_dsdir=root_dsdir)


def scan_trunk(root_dsdir=None):
    """ scan trunk DataSeriesDir
    """
    _map = _pack_arguments(map_.scan_trunk)
    nspecs = _count_arguments(map_.scan_trunk)
    return model.DataSeriesDir(map_=_map, nspecs=nspecs, depth=1,
                               root_dsdir=root_dsdir)


def scan_branch(root_dsdir=None):
    """ scan branch DataSeriesDir
    """
    spec_dfile = file_.data_series_specifier(
        file_prefix=SPEC_FILE_PREFIX,
        map_dct_={'tors_names': lambda specs: specs[0]},
        spec_keys=['tors_names'])

    _map = _pack_arguments(map_.scan_branch)
    nspecs = _count_arguments(map_.scan_branch)
    return model.DataSeriesDir(map_=_map, nspecs=nspecs, depth=1,
                               spec_dfile=spec_dfile,
                               root_dsdir=root_dsdir)


def scan_leaf(root_dsdir=None):
    """ scan leaf DataSeriesDir
    """
    spec_dfile = file_.data_series_specifier(
        file_prefix=SPEC_FILE_PREFIX,
        map_dct_={'grid_idxs': lambda specs: specs[0]},
        spec_keys=['grid_idxs'])

    _map = _pack_arguments(map_.scan_leaf)
    nspecs = _count_arguments(map_.scan_leaf)
    return model.DataSeriesDir(map_=_map, nspecs=nspecs, depth=1,
                               spec_dfile=spec_dfile,
                               root_dsdir=root_dsdir)


def single_point_leaf(root_dsdir=None):
    """ single-point leaf DataSeriesDir
    """
    return theory_leaf(root_dsdir=root_dsdir)


def tau_trunk(root_dsdir=None):
    """ tau trunk DataSeriesDir
    """
    _map = _pack_arguments(map_.tau_trunk)
    nspecs = _count_arguments(map_.tau_trunk)
    return model.DataSeriesDir(map_=_map, nspecs=nspecs, depth=1,
                               root_dsdir=root_dsdir)


def tau_leaf(root_dsdir=None):
    """ tau leaf DataSeriesDir
    """
    spec_dfile = file_.data_series_specifier(
        file_prefix=SPEC_FILE_PREFIX,
        map_dct_={'conformer_id': lambda specs: specs[0]},
        spec_keys=['conformer_id'])

    _map = _pack_arguments(map_.tau_leaf)
    nspecs = _count_arguments(map_.tau_leaf)
    return model.DataSeriesDir(map_=_map, nspecs=nspecs, depth=1,
                               spec_dfile=spec_dfile,
                               root_dsdir=root_dsdir)


# helpers
def _pack_arguments(function):
    """ generate an equivalent function that takes all of its arguments packed
    into a sequence
    """
    def _function(args=()):
        return function(*args)
    return _function


def _count_arguments(function):
    """ conut the number of arguments that a function takes in
    """
    argspec = function_argspec(function)
    return len(argspec.args)
