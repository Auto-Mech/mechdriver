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
    nlocs = _count_arguments(map_.species_trunk)
    return model.DataSeriesDir(map_=_map, nlocs=nlocs, depth=1,
                               root_dsdir=root_dsdir)

def ts_trunk(root_dsdir=None):
    """ ts trunk DataSeriesDir
    """
    _map = _pack_arguments(map_.ts_trunk)
    nlocs = _count_arguments(map_.ts_trunk)
    return model.DataSeriesDir(map_=_map, nlocs=nlocs, depth=1,
                               root_dsdir=root_dsdir)


def species_leaf(root_dsdir=None):
    """ species leaf DataSeriesDir
    """
    loc_dfile = file_.locator(
        file_prefix=SPEC_FILE_PREFIX,
        map_dct_={
            'inchi': lambda locs: locs[0],
            'charge': lambda locs: locs[1],
            'multiplicity': lambda locs: locs[2],
            'smiles': lambda locs: automol.inchi.smiles(locs[0])},
        loc_keys=['inchi', 'charge', 'multiplicity'])

    _map = _pack_arguments(map_.species_leaf)
    nlocs = _count_arguments(map_.species_leaf)
    return model.DataSeriesDir(map_=_map, nlocs=nlocs, depth=5,
                               loc_dfile=loc_dfile,
                               root_dsdir=root_dsdir)


def reaction_trunk(root_dsdir=None):
    """ reaction trunk DataSeriesDir
    """
    _map = _pack_arguments(map_.reaction_trunk)
    nlocs = _count_arguments(map_.reaction_trunk)
    return model.DataSeriesDir(map_=_map, nlocs=nlocs, depth=1,
                               root_dsdir=root_dsdir)


def reaction_leaf(root_dsdir=None):
    """ reaction leaf DataSeriesDir
    """
    loc_dfile = file_.locator(
        file_prefix=SPEC_FILE_PREFIX,
        map_dct_={
            'inchis': lambda locs: locs[0],
            'charges': lambda locs: locs[1],
            'multiplicities': lambda locs: locs[2],
            'ts_multiplicity': lambda locs: locs[3],
            'smiles': lambda locs: [
                list(map(automol.inchi.smiles, locs[0][0])),
                list(map(automol.inchi.smiles, locs[0][1]))],
        },
        loc_keys=['inchis', 'charges', 'multiplicities', 'ts_multiplicity'])

    _map = _pack_arguments(map_.reaction_leaf)
    nlocs = _count_arguments(map_.reaction_leaf)
    return model.DataSeriesDir(map_=_map, nlocs=nlocs, depth=11,
                               loc_dfile=loc_dfile,
                               root_dsdir=root_dsdir)


def theory_leaf(root_dsdir=None):
    """ theory leaf DataSeriesDir
    """
    loc_dfile = file_.locator(
        file_prefix=SPEC_FILE_PREFIX,
        map_dct_={
            'method': lambda locs: locs[0],
            'basis': lambda locs: locs[1],
            'orb_restricted': lambda locs: locs[2]},
        loc_keys=['method', 'basis', 'orb_restricted'])

    _map = _pack_arguments(map_.theory_leaf)
    nlocs = _count_arguments(map_.theory_leaf)
    return model.DataSeriesDir(map_=_map, nlocs=nlocs, depth=1,
                               loc_dfile=loc_dfile,
                               root_dsdir=root_dsdir)


def run_trunk(root_dsdir=None):
    """ run trunk DataSeriesDir
    """
    _map = _pack_arguments(map_.run_trunk)
    nlocs = _count_arguments(map_.run_trunk)
    return model.DataSeriesDir(map_=_map, nlocs=nlocs, depth=1,
                               root_dsdir=root_dsdir)


def run_leaf(root_dsdir=None):
    """ run leaf DataSeriesDir
    """
    loc_dfile = file_.locator(
        file_prefix=SPEC_FILE_PREFIX,
        map_dct_={'job': lambda locs: locs[0]},
        loc_keys=['job'])

    _map = _pack_arguments(map_.run_leaf)
    nlocs = _count_arguments(map_.run_leaf)
    return model.DataSeriesDir(map_=_map, nlocs=nlocs, depth=1,
                               loc_dfile=loc_dfile,
                               root_dsdir=root_dsdir, removable=True)


def subrun_leaf(root_dsdir=None):
    """ subrun leaf DataSeriesDir
    """
    loc_dfile = file_.locator(
        file_prefix=SPEC_FILE_PREFIX,
        map_dct_={'macro_idx': lambda locs: locs[0],
                  'micro_idx': lambda locs: locs[1]},
        loc_keys=['macro_idx', 'micro_idx'])

    _map = _pack_arguments(map_.subrun_leaf)
    nlocs = _count_arguments(map_.subrun_leaf)
    return model.DataSeriesDir(map_=_map, nlocs=nlocs, depth=1,
                               loc_dfile=loc_dfile,
                               root_dsdir=root_dsdir)


def conformer_trunk(root_dsdir=None):
    """ conformer trunk DataSeriesDir
    """
    _map = _pack_arguments(map_.conformer_trunk)
    nlocs = _count_arguments(map_.conformer_trunk)
    return model.DataSeriesDir(map_=_map, nlocs=nlocs, depth=1,
                               root_dsdir=root_dsdir)


def conformer_leaf(root_dsdir=None):
    """ conformer leaf DataSeriesDir
    """
    loc_dfile = file_.locator(
        file_prefix=SPEC_FILE_PREFIX,
        map_dct_={'conformer_id': lambda locs: locs[0]},
        loc_keys=['conformer_id'])

    _map = _pack_arguments(map_.conformer_leaf)
    nlocs = _count_arguments(map_.conformer_leaf)
    return model.DataSeriesDir(map_=_map, nlocs=nlocs, depth=1,
                               loc_dfile=loc_dfile,
                               root_dsdir=root_dsdir)


def single_point_trunk(root_dsdir=None):
    """ single point trunk DataSeriesDir
    """
    _map = _pack_arguments(map_.single_point_trunk)
    nlocs = _count_arguments(map_.single_point_trunk)
    return model.DataSeriesDir(map_=_map, nlocs=nlocs, depth=1,
                               root_dsdir=root_dsdir)


def scan_trunk(root_dsdir=None):
    """ scan trunk DataSeriesDir
    """
    _map = _pack_arguments(map_.scan_trunk)
    nlocs = _count_arguments(map_.scan_trunk)
    return model.DataSeriesDir(map_=_map, nlocs=nlocs, depth=1,
                               root_dsdir=root_dsdir)


def scan_branch(root_dsdir=None):
    """ scan branch DataSeriesDir
    """
    loc_dfile = file_.locator(
        file_prefix=SPEC_FILE_PREFIX,
        map_dct_={'tors_names': lambda locs: locs[0]},
        loc_keys=['tors_names'])

    _map = _pack_arguments(map_.scan_branch)
    nlocs = _count_arguments(map_.scan_branch)
    return model.DataSeriesDir(map_=_map, nlocs=nlocs, depth=1,
                               loc_dfile=loc_dfile,
                               root_dsdir=root_dsdir)


def scan_leaf(root_dsdir=None):
    """ scan leaf DataSeriesDir
    """
    loc_dfile = file_.locator(
        file_prefix=SPEC_FILE_PREFIX,
        map_dct_={'grid_idxs': lambda locs: locs[0]},
        loc_keys=['grid_idxs'])

    _map = _pack_arguments(map_.scan_leaf)
    nlocs = _count_arguments(map_.scan_leaf)
    return model.DataSeriesDir(map_=_map, nlocs=nlocs, depth=1,
                               loc_dfile=loc_dfile,
                               root_dsdir=root_dsdir)


def single_point_leaf(root_dsdir=None):
    """ single-point leaf DataSeriesDir
    """
    return theory_leaf(root_dsdir=root_dsdir)


def tau_trunk(root_dsdir=None):
    """ tau trunk DataSeriesDir
    """
    _map = _pack_arguments(map_.tau_trunk)
    nlocs = _count_arguments(map_.tau_trunk)
    return model.DataSeriesDir(map_=_map, nlocs=nlocs, depth=1,
                               root_dsdir=root_dsdir)


def tau_leaf(root_dsdir=None):
    """ tau leaf DataSeriesDir
    """
    loc_dfile = file_.locator(
        file_prefix=SPEC_FILE_PREFIX,
        map_dct_={'conformer_id': lambda locs: locs[0]},
        loc_keys=['conformer_id'])

    _map = _pack_arguments(map_.tau_leaf)
    nlocs = _count_arguments(map_.tau_leaf)
    return model.DataSeriesDir(map_=_map, nlocs=nlocs, depth=1,
                               loc_dfile=loc_dfile,
                               root_dsdir=root_dsdir)


def build_trunk(root_dsdir=None):
    """ build trunk DataSeriesDir
    """
    loc_dfile = file_.locator(
        file_prefix=SPEC_FILE_PREFIX,
        map_dct_={'head': lambda locs: locs[0]},
        loc_keys=['head'])

    _map = _pack_arguments(map_.build_trunk)
    nlocs = _count_arguments(map_.build_trunk)
    return model.DataSeriesDir(map_=_map, nlocs=nlocs, depth=1,
                               loc_dfile=loc_dfile,
                               root_dsdir=root_dsdir)


def build_leaf(root_dsdir=None):
    """ build leaf DataSeriesDir
    """
    loc_dfile = file_.locator(
        file_prefix=SPEC_FILE_PREFIX,
        map_dct_={'num': lambda locs: locs[0]},
        loc_keys=['num'])

    _map = _pack_arguments(map_.build_leaf)
    nlocs = _count_arguments(map_.build_leaf)
    return model.DataSeriesDir(map_=_map, nlocs=nlocs, depth=1,
                               loc_dfile=loc_dfile,
                               root_dsdir=root_dsdir, removable=True)


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
