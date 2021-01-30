""" Handles sequences of electronic structure jobs to deal with error handling

    Uses implements the option matrix data structure for robust runs
opts = sequence of strings fed into an elstruct *_options argument
opt = an individual string from `opts`
opts_dct = a dictionary of `optn`s by argument keyword
opts_row = a list of `optn_dct`s for cycling to fix a particular error
opts_mat = a list of `optn_row`s for each error
note: options (`opts_dct`) are a subset of keyword arguments (`kwargs_dct`)
"""

import itertools
import warnings
try:
    from collections.abc import Sequence as _Sequence
except ImportError:
    from collections import Sequence as _Sequence
import automol
import elstruct
import autofile
from autoparse import pattern as app
from autoparse import find as apf


# FUNCTIONS FOR HANDLING THE SEQUENCE OF OPTIONS
def options_matrix_optimization(script_str, prefix,
                                geo, chg, mul, method, basis, prog,
                                errors=(), options_mat=(), feedback=False,
                                frozen_coordinates=(),
                                freeze_dummy_atoms=True,
                                **kwargs):
    """ try several sets of options to generate an output file

    :returns: the input string and the output string
    :rtype: (str, str)
    """
    assert len(errors) == len(options_mat)

    subrun_fs = autofile.fs.subrun(prefix)
    max_macro_idx, _ = max(subrun_fs[-1].existing(), default=(-1, -1))
    macro_idx = max_macro_idx + 1
    micro_idx = 0

    if freeze_dummy_atoms and automol.zmat.is_valid(geo):
        frozen_coordinates = (tuple(frozen_coordinates) +
                              automol.zmat.dummy_coordinate_names(geo))

    kwargs_ = dict(kwargs)
    while True:
        subrun_fs[-1].create([macro_idx, micro_idx])
        path = subrun_fs[-1].path([macro_idx, micro_idx])

        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            inp_str, out_str = elstruct.run.direct(
                elstruct.writer.optimization, script_str, path,
                geo=geo, charge=chg, mult=mul, method=method,
                basis=basis, prog=prog, frozen_coordinates=frozen_coordinates,
                **kwargs_)

        error_vals = [elstruct.reader.has_error_message(prog, error, out_str)
                      for error in errors]

        # Kill the while loop if we Molpro error signaling a hopeless point
        # When an MCSCF WF calculation fails to converge at some step in opt
        # it is not clear how to save the optimization, so we give up on opt
        fail_pattern = app.one_of_these([
            app.escape('The problem occurs in Multi'),
            app.escape('The problem occurs in cipro')
        ])
        if apf.has_match(fail_pattern, out_str, case=False):
            break

        if not any(error_vals):
            # success
            break
        if not is_exhausted(options_mat):
            # try again
            micro_idx += 1
            error_row_idx = error_vals.index(True)
            kwargs_ = updated_kwargs(kwargs, options_mat)
            options_mat = advance(error_row_idx, options_mat)
            if feedback:
                geo = (elstruct.reader.opt_zmatrix(prog, out_str)
                       if automol.zmat.is_valid(geo) else
                       elstruct.reader.opt_geometry(prog, out_str))
        else:
            # failure
            warnings.resetwarnings()
            warnings.warn("elstruct robust run failed; "
                          "last run was in, {}".format(path))
            break

    return inp_str, out_str


def options_matrix_run(input_writer, script_str, prefix,
                       geo, chg, mul, method, basis, prog,
                       errors=(), options_mat=(),
                       **kwargs):
    """ try several sets of options to generate an output file

    :returns: the input string and the output string
    :rtype: (str, str)
    """
    assert len(errors) == len(options_mat)

    subrun_fs = autofile.fs.subrun(prefix)
    max_macro_idx, _ = max(subrun_fs[-1].existing(), default=(-1, -1))
    macro_idx = max_macro_idx + 1
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

        error_vals = [elstruct.reader.has_error_message(prog, error, out_str)
                      for error in errors]

        if not any(error_vals):
            # success
            break
        if not is_exhausted(options_mat):
            # try again
            micro_idx += 1
            error_row_idx = error_vals.index(True)
            kwargs_ = updated_kwargs(kwargs, options_mat)
            options_mat = advance(error_row_idx, options_mat)
        else:
            # failure
            warnings.resetwarnings()
            warnings.warn("elstruct robust run failed; "
                          "last run was in, {}".format(path))
            break

    return inp_str, out_str


def molpro_opts_mat(spc_info, geo):
    """ prepare the errors and options mat to perform successive
        single-point energy calculations in Molpro when the RHF fails to
        converge. This currently only works for doublets.
    """

    # Get the nelectrons, spins, and orbitals for the wf card
    formula = automol.geom.formula(geo)
    elec_count = automol.formula.electron_count(formula)
    two_spin = spc_info[2] - 1
    num_act_elc = two_spin
    num_act_orb = num_act_elc
    closed_orb = (elec_count - num_act_elc) // 2
    occ_orb = closed_orb + num_act_orb

    # Build the strings UHF and CASSCF wf card and set the errors and options
    uhf_str = (
        "{{uhf,maxit=300;wf,{0},1,{1};orbprint,3}}"
    ).format(elec_count, two_spin)
    cas_str = (
        "{{casscf,maxit=40;"
        "closed,{0};occ,{1};wf,{2},1,{3};canonical;orbprint,3}}"
    ).format(closed_orb, occ_orb, elec_count, two_spin)

    errors = [elstruct.Error.SCF_NOCONV]
    options_mat = [
        [{'gen_lines': {2: [uhf_str]}},
         {'gen_lines': {2: [cas_str]}},
         {'gen_lines': {2: [cas_str]}}
         ]
    ]

    return errors, options_mat


# OPTIONS MATRIX IMPLEMENTATION
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
        # assert key.endswith('_options')
        # assert _is_nonstring_sequence(opts)

        if key in kwargs_dct:
            # assert _is_nonstring_sequence(kwargs_dct[key])
            kwargs_dct[key] = tuple(itertools.chain(kwargs_dct[key], opts))
        else:
            kwargs_dct[key] = opts

    return kwargs_dct


def _is_nonstring_sequence(obj):
    return (isinstance(obj, _Sequence)
            and not isinstance(obj, (str, bytes, bytearray)))