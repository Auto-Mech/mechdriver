""" DataSeriess
"""

from inspect import getfullargspec as function_argspec
import automol
from autofile import model
from autofile.schema import loc_maps
from autofile.schema import data_files


SPEC_FILE_PREFIX = 'dir'


# DataSeries for species-specific layers
def species_trunk(prefix, root_ds=None):
    """ species trunk DataSeries

    :param prefix: path to trunk
    :type prefix: string
    :return: dataseries filesystem object for trunk
    :type: DataSeries
    """
    _map = _pack_arguments(loc_maps.species_trunk)
    nlocs = _count_arguments(loc_maps.species_trunk)
    return model.DataSeries(prefix, map_=_map, nlocs=nlocs, depth=1,
                            root_ds=root_ds)


def species_leaf(prefix, root_ds=None):
    """ species leaf DataSeries

    :param prefix: path to leaf
    :type prefix: string
    :return: dataseries filesystem object for leaf
    :type: DataSeries
    """
    loc_dfile = data_files.locator(
        file_prefix=SPEC_FILE_PREFIX,
        map_dct_={
            'inchi': lambda locs: locs[0],
            'charge': lambda locs: locs[1],
            'multiplicity': lambda locs: locs[2],
            'smiles': lambda locs: automol.chi.smiles(locs[0])},
        loc_keys=['inchi', 'charge', 'multiplicity'])

    _map = _pack_arguments(loc_maps.species_leaf)
    nlocs = _count_arguments(loc_maps.species_leaf)
    return model.DataSeries(prefix, map_=_map, nlocs=nlocs, depth=5,
                            loc_dfile=loc_dfile, root_ds=root_ds)


# DataSeries for reaction-specific layers
def reaction_trunk(prefix, root_ds=None):
    """ reaction trunk DataSeries

    :param prefix: path to trunk
    :type prefix: string
    :return: dataseries filesystem object for trunk
    :type: DataSeries
    """
    _map = _pack_arguments(loc_maps.reaction_trunk)
    nlocs = _count_arguments(loc_maps.reaction_trunk)
    return model.DataSeries(prefix, map_=_map, nlocs=nlocs, depth=1,
                            root_ds=root_ds)


def reaction_leaf(prefix, root_ds=None):
    """ reaction leaf DataSeries

    :param prefix: path to leaf
    :type prefix: string
    :return: dataseries filesystem object for leaf
    :type: DataSeries
    """
    loc_dfile = data_files.locator(
        file_prefix=SPEC_FILE_PREFIX,
        map_dct_={
            'inchis': lambda locs: locs[0],
            'charges': lambda locs: locs[1],
            'multiplicities': lambda locs: locs[2],
            'ts_multiplicity': lambda locs: locs[3],
            'smiles': lambda locs: [
                list(map(automol.chi.smiles, locs[0][0])),
                list(map(automol.chi.smiles, locs[0][1]))],
        },
        loc_keys=['inchis', 'charges', 'multiplicities', 'ts_multiplicity'])

    _map = _pack_arguments(loc_maps.reaction_leaf)
    nlocs = _count_arguments(loc_maps.reaction_leaf)
    return model.DataSeries(prefix, map_=_map, nlocs=nlocs, depth=11,
                            loc_dfile=loc_dfile, root_ds=root_ds)


def transition_state_trunk(prefix, root_ds=None):
    """ transition state trunk DataSeries

    :param prefix: path to trunk
    :type prefix: string
    :return: dataseries filesystem object for trunk
    :type: DataSeries
    """
    _map = _pack_arguments(loc_maps.transition_state_trunk)
    nlocs = _count_arguments(loc_maps.transition_state_trunk)
    return model.DataSeries(prefix, map_=_map, nlocs=nlocs, depth=1,
                            root_ds=root_ds)


def transition_state_leaf(prefix, root_ds=None):
    """ transition state leaf DataSeries

    :param prefix: path to leaf
    :type prefix: string
    :return: dataseries filesystem object for leaf
    :type: DataSeries
    """
    loc_dfile = data_files.locator(
        file_prefix=SPEC_FILE_PREFIX,
        map_dct_={'num': lambda locs: locs[0]},
        loc_keys=['num'])

    _map = _pack_arguments(loc_maps.transition_state_leaf)
    nlocs = _count_arguments(loc_maps.transition_state_leaf)
    return model.DataSeries(prefix, map_=_map, nlocs=nlocs, depth=1,
                            loc_dfile=loc_dfile, root_ds=root_ds)


# DataSeries for layers used by both species and reaction file systems
def theory_leaf(prefix, root_ds=None):
    """ theory leaf DataSeries

    :param prefix: path to leaf
    :type prefix: string
    :return: dataseries filesystem object for leaf
    :type: DataSeries
    """

    loc_dfile = data_files.locator(
        file_prefix=SPEC_FILE_PREFIX,
        map_dct_={
            'method': lambda locs: locs[0],
            'basis': lambda locs: locs[1],
            'orb_type': lambda locs: locs[2]},
        loc_keys=['method', 'basis', 'orb_type'])

    _map = _pack_arguments(loc_maps.theory_leaf)
    nlocs = _count_arguments(loc_maps.theory_leaf)
    return model.DataSeries(prefix, map_=_map, nlocs=nlocs, depth=1,
                            loc_dfile=loc_dfile, root_ds=root_ds)


def conformer_trunk(prefix, root_ds=None):
    """ conformer trunk DataSeries

    :param prefix: path to trunk
    :type prefix: string
    :return: dataseries filesystem object for trunk
    :type: DataSeries
    """
    _map = _pack_arguments(loc_maps.conformer_trunk)
    nlocs = _count_arguments(loc_maps.conformer_trunk)
    return model.DataSeries(prefix, map_=_map, nlocs=nlocs, depth=1,
                            root_ds=root_ds)


def conformer_branch(prefix, root_ds=None):
    """ conformer branch DataSeries

    :param prefix: path to branch
    :type prefix: string
    :return: dataseries filesystem object for branch
    :type: DataSeries
    """
    loc_dfile = data_files.locator(
        file_prefix=SPEC_FILE_PREFIX,
        map_dct_={'ring_id': lambda locs: locs[0]},
        loc_keys=['ring_id'])

    _map = _pack_arguments(loc_maps.conformer_branch)
    nlocs = _count_arguments(loc_maps.conformer_branch)
    return model.DataSeries(prefix, map_=_map, nlocs=nlocs, depth=1,
                            loc_dfile=loc_dfile, root_ds=root_ds)


def conformer_leaf(prefix, root_ds=None):
    """ conformer leaf DataSeries

    :param prefix: path to leaf
    :type prefix: string
    :return: dataseries filesystem object for leaf
    :type: DataSeries
    """
    loc_dfile = data_files.locator(
        file_prefix=SPEC_FILE_PREFIX,
        map_dct_={
            'conformer_id': lambda locs: locs[0]
        },
        loc_keys=['conformer_id'])

    _map = _pack_arguments(loc_maps.conformer_leaf)
    nlocs = _count_arguments(loc_maps.conformer_leaf)
    return model.DataSeries(prefix, map_=_map, nlocs=nlocs, depth=1,
                            loc_dfile=loc_dfile, root_ds=root_ds)


def single_point_trunk(prefix, root_ds=None):
    """ single point trunk DataSeries

    :param prefix: path to trunk
    :type prefix: string
    :return: dataseries filesystem object for trunk
    :type: DataSeries
    """
    _map = _pack_arguments(loc_maps.single_point_trunk)
    nlocs = _count_arguments(loc_maps.single_point_trunk)
    return model.DataSeries(prefix, map_=_map, nlocs=nlocs, depth=1,
                            root_ds=root_ds)


def single_point_leaf(prefix, root_ds=None):
    """ single-point leaf DataSeries

    :param prefix: path to leaf
    :type prefix: string
    :return: dataseries filesystem object for leaf
    :type: DataSeries
    """
    return theory_leaf(prefix, root_ds=root_ds)


def high_spin_trunk(prefix, root_ds=None):
    """ high spin, single point trunk DataSeries

    :param prefix: path to trunk
    :type prefix: string
    :return: dataseries filesystem object for trunk
    :type: DataSeries
    """
    _map = _pack_arguments(loc_maps.high_spin_trunk)
    nlocs = _count_arguments(loc_maps.high_spin_trunk)
    return model.DataSeries(prefix, map_=_map, nlocs=nlocs, depth=1,
                            root_ds=root_ds)


def high_spin_leaf(prefix, root_ds=None):
    """ high-spin, single-point leaf DataSeries

    :param prefix: path to leaf
    :type prefix: string
    :return: dataseries filesystem object for leaf
    :type: DataSeries
    """
    return theory_leaf(prefix, root_ds=root_ds)


def symmetry_trunk(prefix, root_ds=None):
    """ symmetric-conformer trunk DataSeries

    :param prefix: path to trunk
    :type prefix: string
    :return: dataseries filesystem object for trunk
    :type: DataSeries
    """
    _map = _pack_arguments(loc_maps.symmetry_trunk)
    nlocs = _count_arguments(loc_maps.symmetry_trunk)
    return model.DataSeries(prefix, map_=_map, nlocs=nlocs, depth=1,
                            root_ds=root_ds)


def zmatrix_trunk(prefix, root_ds=None):
    """ zmatrix trunk DataSeries

    :param prefix: path to trunk
    :type prefix: string
    :return: dataseries filesystem object for trunk
    :type: DataSeries
    """
    _map = _pack_arguments(loc_maps.zmatrix_trunk)
    nlocs = _count_arguments(loc_maps.zmatrix_trunk)
    return model.DataSeries(prefix, map_=_map, nlocs=nlocs, depth=1,
                            root_ds=root_ds)


def zmatrix_leaf(prefix, root_ds=None):
    """ zmatrix leaf DataSeries

    :param prefix: path to leaf
    :type prefix: string
    :return: dataseries filesystem object for leaf
    :type: DataSeries
    """
    loc_dfile = data_files.locator(
        file_prefix=SPEC_FILE_PREFIX,
        map_dct_={'num': lambda locs: locs[0]},
        loc_keys=['num'])

    _map = _pack_arguments(loc_maps.zmatrix_leaf)
    nlocs = _count_arguments(loc_maps.zmatrix_leaf)
    return model.DataSeries(prefix, map_=_map, nlocs=nlocs, depth=1,
                            loc_dfile=loc_dfile, root_ds=root_ds)


def scan_trunk(prefix, root_ds=None):
    """ scan trunk DataSeries

    :param prefix: path to trunk
    :type prefix: string
    :return: dataseries filesystem object for trunk
    :type: DataSeries
    """
    _map = _pack_arguments(loc_maps.scan_trunk)
    nlocs = _count_arguments(loc_maps.scan_trunk)
    return model.DataSeries(prefix, map_=_map, nlocs=nlocs, depth=1,
                            root_ds=root_ds)


def scan_branch(prefix, root_ds=None):
    """ scan branch DataSeries

    :param prefix: path to branch
    :type prefix: string
    :return: dataseries filesystem object for branch
    :type: DataSeries
    """
    loc_dfile = data_files.locator(
        file_prefix=SPEC_FILE_PREFIX,
        map_dct_={'coo_names': lambda locs: locs[0]},
        loc_keys=['coo_names'])

    _map = _pack_arguments(loc_maps.scan_branch)
    nlocs = _count_arguments(loc_maps.scan_branch)
    return model.DataSeries(prefix, map_=_map, nlocs=nlocs, depth=1,
                            loc_dfile=loc_dfile, root_ds=root_ds)


def scan_leaf(prefix, root_ds=None):
    """ scan leaf DataSeries

    :param prefix: path to leaf
    :type prefix: string
    :return: dataseries filesystem object for leaf
    :type: DataSeries
    """
    loc_dfile = data_files.locator(
        file_prefix=SPEC_FILE_PREFIX,
        # this should actually say coo_vals, but leave for compatibility
        map_dct_={'grid_idxs': lambda locs: locs[0]},
        loc_keys=['grid_idxs'])

    _map = _pack_arguments(loc_maps.scan_leaf)
    nlocs = _count_arguments(loc_maps.scan_leaf)
    return model.DataSeries(prefix, map_=_map, nlocs=nlocs, depth=1,
                            loc_dfile=loc_dfile, root_ds=root_ds)


def cscan_trunk(prefix, root_ds=None):
    """ constrained scan trunk DataSeries

    :param prefix: path to trunk
    :type prefix: string
    :return: dataseries filesystem object for trunk
    :type: DataSeries
    """
    _map = _pack_arguments(loc_maps.cscan_trunk)
    nlocs = _count_arguments(loc_maps.cscan_trunk)
    return model.DataSeries(prefix, map_=_map, nlocs=nlocs, depth=1,
                            root_ds=root_ds)


def cscan_branch1(prefix, root_ds=None):
    """ constrained scan branch 2 DataSeries

    :param prefix: path to branch
    :type prefix: string
    :return: dataseries filesystem object for branch
    :type: DataSeries
    """

    def _round_values(val_dct):
        names = list(val_dct.keys())
        vals = [float(round(val, 2)) for val in val_dct.values()]
        val_dct = dict(zip(names, vals))
        return val_dct

    loc_dfile = data_files.locator(
        file_prefix=SPEC_FILE_PREFIX,
        map_dct_={'cons_coo_vals': lambda locs: _round_values(locs[0])},
        loc_keys=['cons_coo_vals'])

    _map = _pack_arguments(loc_maps.cscan_branch1)
    nlocs = _count_arguments(loc_maps.cscan_branch1)
    return model.DataSeries(prefix, map_=_map, nlocs=nlocs, depth=1,
                            loc_dfile=loc_dfile, root_ds=root_ds)


def cscan_branch2(prefix, root_ds=None):
    """ constrained scan branch 1 DataSeries

    :param prefix: path to branch
    :type prefix: string
    :return: dataseries filesystem object for branch
    :type: DataSeries
    """
    loc_dfile = data_files.locator(
        file_prefix=SPEC_FILE_PREFIX,
        map_dct_={'coo_names': lambda locs: locs[0]},
        loc_keys=['coo_names'])

    _map = _pack_arguments(loc_maps.cscan_branch2)
    nlocs = _count_arguments(loc_maps.cscan_branch2)
    return model.DataSeries(prefix, map_=_map, nlocs=nlocs, depth=1,
                            loc_dfile=loc_dfile, root_ds=root_ds)


def cscan_leaf(prefix, root_ds=None):
    """ constrained scan leaf DataSeries

    :param prefix: path to leaf
    :type prefix: string
    :return: dataseries filesystem object for leaf
    :type: DataSeries
    """
    loc_dfile = data_files.locator(
        file_prefix=SPEC_FILE_PREFIX,
        map_dct_={'coo_vals': lambda locs: locs[0]},
        loc_keys=['coo_vals'])

    _map = _pack_arguments(loc_maps.cscan_leaf)
    nlocs = _count_arguments(loc_maps.cscan_leaf)
    return model.DataSeries(prefix, map_=_map, nlocs=nlocs, depth=1,
                            loc_dfile=loc_dfile, root_ds=root_ds)


def tau_trunk(prefix, root_ds=None):
    """ tau trunk DataSeries

    :param prefix: path to trunk
    :type prefix: string
    :return: dataseries filesystem object for trunk
    :type: DataSeries
    """
    _map = _pack_arguments(loc_maps.tau_trunk)
    nlocs = _count_arguments(loc_maps.tau_trunk)
    return model.DataSeries(prefix, map_=_map, nlocs=nlocs, depth=1,
                            root_ds=root_ds)


def tau_leaf(prefix, root_ds=None):
    """ tau leaf DataSeries

    :param prefix: path to leaf
    :type prefix: string
    :return: dataseries filesystem object for leaf
    :type: DataSeries
    """
    loc_dfile = data_files.locator(
        file_prefix=SPEC_FILE_PREFIX,
        map_dct_={'conformer_id': lambda locs: locs[0]},
        loc_keys=['conformer_id'])

    _map = _pack_arguments(loc_maps.tau_leaf)
    nlocs = _count_arguments(loc_maps.tau_leaf)
    return model.DataSeries(prefix, map_=_map, nlocs=nlocs, depth=1,
                            loc_dfile=loc_dfile, root_ds=root_ds)


def energy_transfer_trunk(prefix, root_ds=None):
    """ energy transfer trunk DataSeries

    :param prefix: path to trunk
    :type prefix: string
    :return: dataseries filesystem object for trunk
    :type: DataSeries
    """
    _map = _pack_arguments(loc_maps.energy_transfer_trunk)
    nlocs = _count_arguments(loc_maps.energy_transfer_trunk)
    return model.DataSeries(prefix, map_=_map, nlocs=nlocs, depth=1,
                            root_ds=root_ds)


def energy_transfer_branch(prefix, root_ds=None):
    """ energy transfer branch DataSeries

    :param prefix: path to branch
    :type prefix: string
    :return: dataseries filesystem object for branch
    :type: DataSeries
    """
    return species_leaf(prefix, root_ds=root_ds)


def energy_transfer_leaf(prefix, root_ds=None):
    """ energy transfer leaf DataSeries

    :param prefix: path to leaf
    :type prefix: string
    :return: dataseries filesystem object for leaf
    :type: DataSeries
    """
    return theory_leaf(prefix, root_ds=root_ds)


def vrctst_trunk(prefix, root_ds=None):
    """ vrctst trunk DataSeries

    :param prefix: path to trunk
    :type prefix: string
    :return: dataseries filesystem object for trunk
    :type: DataSeries
    """
    _map = _pack_arguments(loc_maps.vrctst_trunk)
    nlocs = _count_arguments(loc_maps.vrctst_trunk)
    return model.DataSeries(prefix, map_=_map, nlocs=nlocs, depth=1,
                            root_ds=root_ds)


def vrctst_leaf(prefix, root_ds=None):
    """ vrctst leaf DataSeries

    :param prefix: path to leaf
    :type prefix: string
    :return: dataseries filesystem object for leaf
    :type: DataSeries
    """
    loc_dfile = data_files.locator(
        file_prefix=SPEC_FILE_PREFIX,
        map_dct_={'num': lambda locs: locs[0]},
        loc_keys=['num'])

    _map = _pack_arguments(loc_maps.vrctst_leaf)
    nlocs = _count_arguments(loc_maps.vrctst_leaf)
    return model.DataSeries(prefix, map_=_map, nlocs=nlocs, depth=1,
                            loc_dfile=loc_dfile, root_ds=root_ds)


# DataSeries specific to the run file system
def run_trunk(prefix, root_ds=None):
    """ run trunk DataSeries

    :param prefix: path to trunk
    :type prefix: string
    :return: dataseries filesystem object for trunk
    :type: DataSeries
    """
    _map = _pack_arguments(loc_maps.run_trunk)
    nlocs = _count_arguments(loc_maps.run_trunk)
    return model.DataSeries(prefix, map_=_map, nlocs=nlocs, depth=1,
                            root_ds=root_ds)


def run_leaf(prefix, root_ds=None):
    """ run leaf DataSeries

    :param prefix: path to leaf
    :type prefix: string
    :return: dataseries filesystem object for leaf
    :type: DataSeries
    """
    loc_dfile = data_files.locator(
        file_prefix=SPEC_FILE_PREFIX,
        map_dct_={'job': lambda locs: locs[0]},
        loc_keys=['job'])

    _map = _pack_arguments(loc_maps.run_leaf)
    nlocs = _count_arguments(loc_maps.run_leaf)
    return model.DataSeries(prefix, map_=_map, nlocs=nlocs, depth=1,
                            loc_dfile=loc_dfile, root_ds=root_ds,
                            removable=True)


def subrun_leaf(prefix, root_ds=None):
    """ subrun leaf DataSeries

    :param prefix: path to leaf
    :type prefix: string
    :return: dataseries filesystem object for leaf
    :type: DataSeries
    """
    loc_dfile = data_files.locator(
        file_prefix=SPEC_FILE_PREFIX,
        map_dct_={'macro_idx': lambda locs: locs[0],
                  'micro_idx': lambda locs: locs[1]},
        loc_keys=['macro_idx', 'micro_idx'])

    _map = _pack_arguments(loc_maps.subrun_leaf)
    nlocs = _count_arguments(loc_maps.subrun_leaf)
    return model.DataSeries(prefix, map_=_map, nlocs=nlocs, depth=1,
                            loc_dfile=loc_dfile, root_ds=root_ds)


def build_trunk(prefix, root_ds=None):
    """ build trunk DataSeries

    :param prefix: path to trunk
    :type prefix: string
    :return: dataseries filesystem object for trunk
    :type: DataSeries
    """
    loc_dfile = data_files.locator(
        file_prefix=SPEC_FILE_PREFIX,
        map_dct_={'head': lambda locs: locs[0]},
        loc_keys=['head'])

    _map = _pack_arguments(loc_maps.build_trunk)
    nlocs = _count_arguments(loc_maps.build_trunk)
    return model.DataSeries(prefix, map_=_map, nlocs=nlocs, depth=1,
                            loc_dfile=loc_dfile, root_ds=root_ds)


def build_branch(prefix, root_ds=None):
    """ build trunk DataSeries

    :param prefix: path to trunk
    :type prefix: string
    :return: dataseries filesystem object for trunk
    :type: DataSeries
    """
    loc_dfile = data_files.locator(
        file_prefix=SPEC_FILE_PREFIX,
        map_dct_={'fml_str': lambda locs: locs[0]},
        loc_keys=['fml_str'])

    _map = _pack_arguments(loc_maps.build_branch)
    nlocs = _count_arguments(loc_maps.build_branch)
    return model.DataSeries(prefix, map_=_map, nlocs=nlocs, depth=1,
                            loc_dfile=loc_dfile, root_ds=root_ds)


def build_leaf(prefix, root_ds=None):
    """ build leaf DataSeries

    :param prefix: path to leaf
    :type prefix: string
    :return: dataseries filesystem object for leaf
    :type: DataSeries
    """
    loc_dfile = data_files.locator(
        file_prefix=SPEC_FILE_PREFIX,
        map_dct_={'num': lambda locs: locs[0]},
        loc_keys=['num'])

    _map = _pack_arguments(loc_maps.build_leaf)
    nlocs = _count_arguments(loc_maps.build_leaf)
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

    _function.__name__ = function.__name__
    return _function


def _count_arguments(function):
    """ conut the number of arguments that a function takes in
    """
    argspec = function_argspec(function)
    return len(argspec.args)
