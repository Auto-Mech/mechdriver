""" DataSeriess
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


def species_trunk(prefix, root_ds=None):
    """ species trunk DataSeries
    """
    _map = _pack_arguments(map_.species_trunk)
    nlocs = _count_arguments(map_.species_trunk)
    return model.DataSeries(prefix, map_=_map, nlocs=nlocs, depth=1,
                            root_ds=root_ds)


def species_leaf(prefix, root_ds=None):
    """ species leaf DataSeries
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
    return model.DataSeries(prefix, map_=_map, nlocs=nlocs, depth=5,
                            loc_dfile=loc_dfile, root_ds=root_ds)


def theory_leaf(prefix, root_ds=None):
    """ theory leaf DataSeries
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
    return model.DataSeries(prefix, map_=_map, nlocs=nlocs, depth=1,
                            loc_dfile=loc_dfile, root_ds=root_ds)


def conformer_trunk(prefix, root_ds=None):
    """ conformer trunk DataSeries
    """
    _map = _pack_arguments(map_.conformer_trunk)
    nlocs = _count_arguments(map_.conformer_trunk)
    return model.DataSeries(prefix, map_=_map, nlocs=nlocs, depth=1,
                            root_ds=root_ds)


def conformer_leaf(prefix, root_ds=None):
    """ conformer leaf DataSeries
    """
    loc_dfile = file_.locator(
        file_prefix=SPEC_FILE_PREFIX,
        map_dct_={'conformer_id': lambda locs: locs[0]},
        loc_keys=['conformer_id'])

    _map = _pack_arguments(map_.conformer_leaf)
    nlocs = _count_arguments(map_.conformer_leaf)
    return model.DataSeries(prefix, map_=_map, nlocs=nlocs, depth=1,
                            loc_dfile=loc_dfile, root_ds=root_ds)


def single_point_trunk(prefix, root_ds=None):
    """ single point trunk DataSeries
    """
    _map = _pack_arguments(map_.single_point_trunk)
    nlocs = _count_arguments(map_.single_point_trunk)
    return model.DataSeries(prefix, map_=_map, nlocs=nlocs, depth=1,
                            root_ds=root_ds)


def single_point_leaf(prefix, root_ds=None):
    """ single-point leaf DataSeries
    """
    return theory_leaf(prefix, root_ds=root_ds)


def high_spin_trunk(prefix, root_ds=None):
    """ high spin, single point trunk DataSeries
    """
    _map = _pack_arguments(map_.high_spin_trunk)
    nlocs = _count_arguments(map_.high_spin_trunk)
    return model.DataSeries(prefix, map_=_map, nlocs=nlocs, depth=1,
                            root_ds=root_ds)


def high_spin_leaf(prefix, root_ds=None):
    """ high-spin, single-point leaf DataSeries
    """
    return theory_leaf(prefix, root_ds=root_ds)


def scan_trunk(prefix, root_ds=None):
    """ scan trunk DataSeries
    """
    _map = _pack_arguments(map_.scan_trunk)
    nlocs = _count_arguments(map_.scan_trunk)
    return model.DataSeries(prefix, map_=_map, nlocs=nlocs, depth=1,
                            root_ds=root_ds)


def scan_branch(prefix, root_ds=None):
    """ scan branch DataSeries
    """
    loc_dfile = file_.locator(
        file_prefix=SPEC_FILE_PREFIX,
        map_dct_={'coo_names': lambda locs: locs[0]},
        loc_keys=['coo_names'])

    _map = _pack_arguments(map_.scan_branch)
    nlocs = _count_arguments(map_.scan_branch)
    return model.DataSeries(prefix, map_=_map, nlocs=nlocs, depth=1,
                            loc_dfile=loc_dfile, root_ds=root_ds)


def scan_leaf(prefix, root_ds=None):
    """ scan leaf DataSeries
    """
    loc_dfile = file_.locator(
        file_prefix=SPEC_FILE_PREFIX,
        # this should actually say coo_vals, but leave for compatibility
        map_dct_={'grid_idxs': lambda locs: locs[0]},
        loc_keys=['grid_idxs'])

    _map = _pack_arguments(map_.scan_leaf)
    nlocs = _count_arguments(map_.scan_leaf)
    return model.DataSeries(prefix, map_=_map, nlocs=nlocs, depth=1,
                            loc_dfile=loc_dfile, root_ds=root_ds)


def cscan_trunk(prefix, root_ds=None):
    """ constrained scan trunk DataSeries
    """
    _map = _pack_arguments(map_.cscan_trunk)
    nlocs = _count_arguments(map_.cscan_trunk)
    return model.DataSeries(prefix, map_=_map, nlocs=nlocs, depth=1,
                            root_ds=root_ds)


def cscan_branch1(prefix, root_ds=None):
    """ constrained scan branch 1 DataSeries
    """
    loc_dfile = file_.locator(
        file_prefix=SPEC_FILE_PREFIX,
        map_dct_={'coo_names': lambda locs: locs[0]},
        loc_keys=['coo_names'])

    _map = _pack_arguments(map_.cscan_branch1)
    nlocs = _count_arguments(map_.cscan_branch1)
    return model.DataSeries(prefix, map_=_map, nlocs=nlocs, depth=1,
                            loc_dfile=loc_dfile, root_ds=root_ds)


def cscan_branch2(prefix, root_ds=None):
    """ constrained scan branch 2 DataSeries
    """
    loc_dfile = file_.locator(
        file_prefix=SPEC_FILE_PREFIX,
        map_dct_={'coo_vals': lambda locs: locs[0]},
        loc_keys=['coo_vals'])

    _map = _pack_arguments(map_.cscan_branch2)
    nlocs = _count_arguments(map_.cscan_branch2)
    return model.DataSeries(prefix, map_=_map, nlocs=nlocs, depth=1,
                            loc_dfile=loc_dfile, root_ds=root_ds)


def cscan_leaf(prefix, root_ds=None):
    """ constrained scan branch 2 DataSeries
    """

    def _round_values(val_dct):
        names = list(val_dct.keys())
        vals = [float(round(val, 2)) for val in val_dct.values()]
        val_dct = dict(zip(names, vals))
        return val_dct

    loc_dfile = file_.locator(
        file_prefix=SPEC_FILE_PREFIX,
        map_dct_={'cons_coo_vals': lambda locs: _round_values(locs[0])},
        loc_keys=['cons_coo_vals'])

    _map = _pack_arguments(map_.cscan_leaf)
    nlocs = _count_arguments(map_.cscan_leaf)
    return model.DataSeries(prefix, map_=_map, nlocs=nlocs, depth=1,
                            loc_dfile=loc_dfile, root_ds=root_ds)


def tau_trunk(prefix, root_ds=None):
    """ tau trunk DataSeries
    """
    _map = _pack_arguments(map_.tau_trunk)
    nlocs = _count_arguments(map_.tau_trunk)
    return model.DataSeries(prefix, map_=_map, nlocs=nlocs, depth=1,
                            root_ds=root_ds)


def energy_transfer_trunk(prefix, root_ds=None):
    """ energy transfer trunk DataSeries
    """
    _map = _pack_arguments(map_.energy_transfer_trunk)
    nlocs = _count_arguments(map_.energy_transfer_trunk)
    tmp = model.DataSeries(prefix, map_=_map, nlocs=nlocs, depth=1,
                            root_ds=root_ds)
    return model.DataSeries(prefix, map_=_map, nlocs=nlocs, depth=1,
                            root_ds=root_ds)


def energy_transfer_branch(prefix, root_ds=None):
    """ energy transfer branch DataSeries
    """
    return species_leaf(prefix, root_ds=root_ds)


def energy_transfer_leaf(prefix, root_ds=None):
    """ energy transfer leaf DataSeries
    """
    return theory_leaf(prefix, root_ds=root_ds)


def tau_leaf(prefix, root_ds=None):
    """ tau leaf DataSeries
    """
    loc_dfile = file_.locator(
        file_prefix=SPEC_FILE_PREFIX,
        map_dct_={'conformer_id': lambda locs: locs[0]},
        loc_keys=['conformer_id'])

    _map = _pack_arguments(map_.tau_leaf)
    nlocs = _count_arguments(map_.tau_leaf)
    return model.DataSeries(prefix, map_=_map, nlocs=nlocs, depth=1,
                            loc_dfile=loc_dfile, root_ds=root_ds)


def reaction_trunk(prefix, root_ds=None):
    """ reaction trunk DataSeries
    """
    _map = _pack_arguments(map_.reaction_trunk)
    nlocs = _count_arguments(map_.reaction_trunk)
    return model.DataSeries(prefix, map_=_map, nlocs=nlocs, depth=1,
                            root_ds=root_ds)


def reaction_leaf(prefix, root_ds=None):
    """ reaction leaf DataSeries
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
    return model.DataSeries(prefix, map_=_map, nlocs=nlocs, depth=11,
                            loc_dfile=loc_dfile, root_ds=root_ds)


def ts_trunk(prefix, root_ds=None):
    """ ts trunk DataSeries
    """
    _map = _pack_arguments(map_.ts_trunk)
    nlocs = _count_arguments(map_.ts_trunk)
    return model.DataSeries(prefix, map_=_map, nlocs=nlocs, depth=1,
                            root_ds=root_ds)


def direction_leaf(prefix, root_ds=None):
    """ direction leaf  DataSeries
    """
    _map = _pack_arguments(map_.direction_leaf)
    nlocs = _count_arguments(map_.direction_leaf)
    return model.DataSeries(prefix, map_=_map, nlocs=nlocs, depth=1,
                            root_ds=root_ds)


def run_trunk(prefix, root_ds=None):
    """ run trunk DataSeries
    """
    _map = _pack_arguments(map_.run_trunk)
    nlocs = _count_arguments(map_.run_trunk)
    return model.DataSeries(prefix, map_=_map, nlocs=nlocs, depth=1,
                            root_ds=root_ds)


def run_leaf(prefix, root_ds=None):
    """ run leaf DataSeries
    """
    loc_dfile = file_.locator(
        file_prefix=SPEC_FILE_PREFIX,
        map_dct_={'job': lambda locs: locs[0]},
        loc_keys=['job'])

    _map = _pack_arguments(map_.run_leaf)
    nlocs = _count_arguments(map_.run_leaf)
    return model.DataSeries(prefix, map_=_map, nlocs=nlocs, depth=1,
                            loc_dfile=loc_dfile, root_ds=root_ds,
                            removable=True)


def subrun_leaf(prefix, root_ds=None):
    """ subrun leaf DataSeries
    """
    loc_dfile = file_.locator(
        file_prefix=SPEC_FILE_PREFIX,
        map_dct_={'macro_idx': lambda locs: locs[0],
                  'micro_idx': lambda locs: locs[1]},
        loc_keys=['macro_idx', 'micro_idx'])

    _map = _pack_arguments(map_.subrun_leaf)
    nlocs = _count_arguments(map_.subrun_leaf)
    return model.DataSeries(prefix, map_=_map, nlocs=nlocs, depth=1,
                            loc_dfile=loc_dfile, root_ds=root_ds)


def build_trunk(prefix, root_ds=None):
    """ build trunk DataSeries
    """
    loc_dfile = file_.locator(
        file_prefix=SPEC_FILE_PREFIX,
        map_dct_={'head': lambda locs: locs[0]},
        loc_keys=['head'])

    _map = _pack_arguments(map_.build_trunk)
    nlocs = _count_arguments(map_.build_trunk)
    return model.DataSeries(prefix, map_=_map, nlocs=nlocs, depth=1,
                            loc_dfile=loc_dfile, root_ds=root_ds)


def build_leaf(prefix, root_ds=None):
    """ build leaf DataSeries
    """
    loc_dfile = file_.locator(
        file_prefix=SPEC_FILE_PREFIX,
        map_dct_={'num': lambda locs: locs[0]},
        loc_keys=['num'])

    _map = _pack_arguments(map_.build_leaf)
    nlocs = _count_arguments(map_.build_leaf)
    return model.DataSeries(prefix, map_=_map, nlocs=nlocs, depth=1,
                            loc_dfile=loc_dfile, root_ds=root_ds,
                            removable=True)


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
