""" implements the option matrix data structure for robust runs
opts = sequence of strings fed into an elstruct *_options argument
opt = an individual string from `opts`
opts_dct = a dictionary of `optn`s by argument keyword
opts_row = a list of `optn_dct`s for cycling to fix a particular error
opts_mat = a list of `optn_row`s for each error
note: options (`opts_dct`) are a subset of keyword arguments (`kwargs_dct`)
"""
import itertools
try:
    from collections.abc import Sequence as _Sequence
except ImportError:
    from collections import Sequence as _Sequence


def is_exhausted(opts_mat):
    """ is this options matrix exhausted?
    the matrix is exhausted if any of its rows are
    """
    return any(not opts_row for opts_row in opts_mat)


def advance(row_idx, opts_mat):
    """ advance one row from the options matrix
    """
    assert not is_exhausted(opts_mat)
    assert row_idx < len(opts_mat)

    opts_mat = list(map(list, opts_mat))
    opts_mat[row_idx].pop(0)
    opts_mat = tuple(map(tuple, opts_mat))

    return opts_mat


def updated_kwargs(kwargs_dct, opts_mat):
    """ update `kwargs_dct` with the current set of options
    merges the first entries from each options row into `kwargs_dct`
    """
    assert not is_exhausted(opts_mat)

    opts_col = [opts_row[0] for opts_row in opts_mat]

    # prohibit conflicting options across rows
    all_keys = list(itertools.chain(*opts_col))
    assert len(all_keys) == len(set(all_keys))

    for opts_dct in opts_col:
        kwargs_dct = _update_kwargs(kwargs_dct, opts_dct)

    return kwargs_dct


def _update_kwargs(kwargs_dct, opts_dct):
    """ update a kwargs dictionary with a dictionary of options
    where an option has already been set in kwargs, the values from `opts_dct`
    are appended to the existing ones; otherwise this is a regular dictionary
    update
    """
    kwargs_dct = dict(kwargs_dct).copy()

    for key, opts in opts_dct.items():
        # sanity check: make sure the options are valid
        #assert key.endswith('_options')
        #assert _is_nonstring_sequence(opts)

        if key in kwargs_dct:
            #assert _is_nonstring_sequence(kwargs_dct[key])
            kwargs_dct[key] = tuple(itertools.chain(kwargs_dct[key], opts))
        else:
            kwargs_dct[key] = opts

    return kwargs_dct


def _is_nonstring_sequence(obj):
    return (isinstance(obj, _Sequence)
            and not isinstance(obj, (str, bytes, bytearray)))
