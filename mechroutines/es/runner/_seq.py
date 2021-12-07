""" Handles sequences of electronic structure jobs to deal with error handling.

    Implements the option matrix data structure for robust runs
opts = sequence of strings fed into an elstruct *_options argument
opt = an individual string from `opts`
opts_dct = a dictionary of `optn`s by argument keyword
opts_row = a list of `optn_dct`s for cycling to fix a particular error
opts_mat = a list of `optn_row`s for each error
note: options (`opts_dct`) are a subset of keyword arguments (`kwargs_dct`)

    The options matrix is used to update `kwargs` dictionaries that are passed
    to elstruct package functions to write electronic structure input files.
"""

import itertools
import warnings
from collections.abc import Sequence as _Sequence
import automol
import elstruct
import autofile


# FUNCTIONS FOR HANDLING THE SEQUENCE OF OPTIONS
def options_matrix_optimization(script_str, prefix,
                                geo, chg, mul, method, basis, prog,
                                errors=(), options_mat=(), feedback=False,
                                frozen_coordinates=(),
                                freeze_dummy_atoms=True,
                                **kwargs):
    """ try several sets of options to generate an output file

        :param script_str: BASH submission script for electronic structure job
        :type script_str: str
        :param prefix:
        :type prefix:
        :param geo: molecular geometry or Z-Matrix
        :type geo:
        :param chg: electric charge
        :type chg: int
        :param mul: spin-multiplicity
        :type mul: int
        :param method: name of the electronic structure method
        :type method: str
        :param basis: name of the basis set
        :type basis: str
        :param prog: name of the electronic structure program
        :type prog: str
        :param errors: list of error message types to search output for
        :type errors: tuple(str)
        :param options_mat: varopis options to run job with
        :type options_mat: tuple(dict[str: str])
        :param feedback: update geom with job from previous sequence
        :type feedback: bool
        :param frozen_coordinates: Z-matrix coordinate names to freeze in opts
        :type frozen_coordinates: tuple(str)
        :param freeze_dummy_atoms: freeze any coords defined by dummy atoms
        :type freeze_dummy_atoms: bool
        :param kwargs:
        :type:
        :returns: the input string and the output string
        :rtype: (str, str)
    """
    assert len(errors) == len(options_mat)

    subrun_fs = autofile.fs.subrun(prefix)
    max_macro_idx, _ = max(subrun_fs[-1].existing(), default=(-1, -1))
    macro_idx = max_macro_idx + 1
    if macro_idx == 26:
        macro_idx = 0
    micro_idx = 0

    if freeze_dummy_atoms and automol.zmat.is_valid(geo):
        frozen_coordinates = (tuple(frozen_coordinates) +
                              automol.zmat.dummy_coordinate_names(geo))

    kwargs_ = dict(kwargs)

    # Initialize loop geo
    step_geo = geo
    while True:
        subrun_fs[-1].create([macro_idx, micro_idx])
        path = subrun_fs[-1].path([macro_idx, micro_idx])
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            inp_str, out_str = elstruct.run.direct(
                elstruct.writer.optimization, script_str, path,
                geo=step_geo, charge=chg, mult=mul, method=method,
                basis=basis, prog=prog, frozen_coordinates=frozen_coordinates,
                **kwargs_)

        # List any errors found in the output
        errs_found = [err for err in errors
                      if elstruct.reader.has_error_message(prog, err, out_str)]

        # Hacky nonsense, need new system for errors
        # Failure: Break if MCSCF Failure found that cnnout be fixed
        if elstruct.reader.has_error_message(
             prog, elstruct.Error.MCSCF_NOCONV, out_str):
            print("elstruct robust run failed; "
                  "unfixable MCSCF convergence issues")
            break

        if elstruct.reader.has_error_message(
             prog, elstruct.Error.LIN_DEP_BASIS, out_str):
            if automol.zmat.is_valid(step_geo):
                step_geo = automol.zmat.geometry(step_geo)
                frozen_coordinates = ()
                print('fail for linear dependence, trying a geom')
                continue
            else:
                print('linear dependence issue with geom, breaking')
                break

        # Break if successful
        if not any(errs_found):
            # success
            break

        if not is_exhausted(options_mat):
            # try again
            micro_idx += 1
            error_row_idx = errors.index(errs_found[0])
            print('options mat errs test', errors, errs_found, error_row_idx)
            kwargs_ = updated_kwargs(kwargs, options_mat)
            options_mat = advance(error_row_idx, options_mat)
            if feedback:
                # Try and get ZMA, then geo
                # if neither present use geo from prev. step (for weird errs)
                geo = elstruct.reader.opt_geometry(prog, out_str)
                if automol.zmat.is_valid(step_geo):
                    dum_key_dct = automol.zmat.dummy_key_dictionary(step_geo)
                    geo_wdum = automol.geom.insert_dummies(geo, dum_key_dct)
                    geo = automol.zmat.from_geometry(step_geo, geo_wdum)
                if geo is not None:
                    step_geo = geo
        else:
            # failure
            warnings.resetwarnings()
            warnings.warn("elstruct robust run failed; "
                          f"last run was in {path}")
            break

    return inp_str, out_str


def options_matrix_run(input_writer, script_str, prefix,
                       geo, chg, mul, method, basis, prog,
                       errors=(), options_mat=(),
                       **kwargs):
    """ try several sets of options to generate an output file

        :param script_str: Shell submission script for electronic structure job
        :type script_str: str
        :param prefix:
        :type prefix:
        :param geo: molecular geometry or Z-Matrix
        :type geo:
        :param chg: electric charge
        :type chg: int
        :param mul: spin-multiplicity
        :type mul: int
        :param method: name of the electronic structure method
        :type method: str
        :param basis: name of the basis set
        :type basis: str
        :param prog: name of the electronic structure program
        :type prog: str
        :param errors: list of error message types to search output for
        :type errors: tuple(str)
        :param options_mat: varopis options to run job with
        :type options_mat: tuple(dict[str: str])

    :returns: the input string and the output string
    :rtype: (str, str)
    """

    assert len(errors) == len(options_mat)

    subrun_fs = autofile.fs.subrun(prefix)
    max_macro_idx, _ = max(subrun_fs[-1].existing(), default=(-1, -1))
    macro_idx = max_macro_idx + 1
    if macro_idx == 26:
        macro_idx = 0
    micro_idx = 0

    kwargs_ = dict(kwargs)
    while True:
        subrun_fs[-1].create([macro_idx, micro_idx])
        path = subrun_fs[-1].path([macro_idx, micro_idx])
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            inp_str, out_str = elstruct.run.direct(
                input_writer, script_str, path,
                geo=geo, charge=chg, mult=mul, method=method,
                basis=basis, prog=prog, **kwargs_)

        # List any errors found in the output
        errs_found = [err for err in errors
                      if elstruct.reader.has_error_message(prog, err, out_str)]

        # Break if successful
        if not any(errs_found):
            # success
            break

        if not is_exhausted(options_mat):
            # try again
            micro_idx += 1
            error_row_idx = errors.index(errs_found[0])
            print('errs test', errors, errs_found, error_row_idx)
            kwargs_ = updated_kwargs(kwargs, options_mat)
            options_mat = advance(error_row_idx, options_mat)
        else:
            # failure
            warnings.resetwarnings()
            warnings.warn("elstruct robust run failed; "
                          f"last run was in {path}")
            break

    return inp_str, out_str


# OPTIONS MATRIX IMPLEMENTATION
def is_exhausted(opts_mat):
    """ Assess if the options matrix has no remaining option
        rows remaining inside of the matrix by assessing if
        any row is empty.

        :param opts_mat: options matrix
        :type opts_mat: tuple(tuple(dict))
        :rtype: bool
    """
    return any(not opts_row for opts_row in opts_mat)


def advance(row_idx, opts_mat):
    """ Advance by one the list of options that are defined
        in the given row of an options matrix. Does so removing
        the first element of options list of that row.

        :param row_idx: idx of row containing options list to advance
        :type row_idx: int
        :param opts_mat: options matrix
        :type opts_mat: tuple(tuple(dict))
        :rtype: tuple(tuple(dict))
    """

    assert not is_exhausted(opts_mat)
    assert row_idx < len(opts_mat)

    opts_mat = list(map(list, opts_mat))
    opts_mat[row_idx].pop(0)
    opts_mat = tuple(map(tuple, opts_mat))

    return opts_mat


def updated_kwargs(kwargs_dct, opts_mat):
    """ Update a `kwargs_dct` with the current options in the options matrix
        with the current set of options. Does so by merging the first
        entries from each options row into `kwargs_dct`

        :param kwargs_dct: dictionary of elstruct arguments
        :type kwargs_dct: dict[str:]
        :param opts_mat: options matrix
        :type opts_mat: tuple(tuple(dict))
        :rtype: dict[str:]
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
    """ Update a kwargs dictionary with a dictionary of options,
        where an option has already been set in kwargs. The values from
        `opts_dct` are appended to the existing ones; otherwise
        this is a regular dictionary.

        :param opts_dct: dictionary of options within matrix row
        :type opts_mat: tuple(tuple(dict))
    """

    kwargs_dct = dict(kwargs_dct).copy()

    for key, opts in opts_dct.items():
        # sanity check: make sure the options are valid
        # assert key.endswith('_options')
        # assert _is_nonstring_sequence(opts)

        if key in kwargs_dct:
            # assert _is_nonstring_sequence(kwargs_dct[key])
            kwargs_dct[key] = tuple(itertools.chain(kwargs_dct[key], opts))
        else:
            kwargs_dct[key] = opts

    return kwargs_dct


def _is_nonstring_sequence(obj):
    """ Assess if the sequence is composed of non-string objects.

        :param obj: sequence object to assess
        :rtype: bool
    """
    return (isinstance(obj, _Sequence)
            and not isinstance(obj, (str, bytes, bytearray)))
