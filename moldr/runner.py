""" elstruct runners (formerly elcarro)
"""
import warnings
import automol
import elstruct
import autofile
import moldr.optsmat


def options_matrix_optimization(script_str, prefix,
                                geom, chg, mul, method, basis, prog,
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
    max_macro_idx, _ = max(subrun_fs.leaf.existing(), default=(-1, -1))
    macro_idx = max_macro_idx + 1
    micro_idx = 0
    read_geom_ = (elstruct.reader.opt_zmatrix_(prog)
                  if automol.zmatrix.is_valid(geom) else
                  elstruct.reader.opt_geometry_(prog))

    if freeze_dummy_atoms and automol.zmatrix.is_valid(geom):
        frozen_coordinates = (tuple(frozen_coordinates) +
                              automol.zmatrix.dummy_coordinate_names(geom))

    kwargs_ = dict(kwargs)
    while True:
        subrun_fs.leaf.create([macro_idx, micro_idx])
        path = subrun_fs.leaf.path([macro_idx, micro_idx])

        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            inp_str, out_str = elstruct.run.direct(
                elstruct.writer.optimization, script_str, path,
                geom=geom, charge=chg, mult=mul, method=method,
                basis=basis, prog=prog, frozen_coordinates=frozen_coordinates,
                **kwargs_)

        error_vals = [elstruct.reader.has_error_message(prog, error, out_str)
                      for error in errors]

        if not any(error_vals):
            # success
            break
        elif not moldr.optsmat.is_exhausted(options_mat):
            # try again
            micro_idx += 1
            error_row_idx = error_vals.index(True)
            kwargs_ = moldr.optsmat.updated_kwargs(kwargs, options_mat)
            options_mat = moldr.optsmat.advance(error_row_idx, options_mat)
            if feedback:
                geom = read_geom_(out_str)
        else:
            # failure
            warnings.resetwarnings()
            warnings.warn("elstruct robust run failed; "
                          "last run was in, {}".format(path))
            break

    return inp_str, out_str


def options_matrix_run(input_writer, script_str, prefix,
                       geom, chg, mul, method, basis, prog,
                       errors=(), options_mat=(),
                       **kwargs):
    """ try several sets of options to generate an output file

    :returns: the input string and the output string
    :rtype: (str, str)
    """
    assert len(errors) == len(options_mat)

    subrun_fs = autofile.fs.subrun(prefix)
    max_macro_idx, _ = max(subrun_fs.leaf.existing(), default=(-1, -1))
    macro_idx = max_macro_idx + 1
    micro_idx = 0

    kwargs_ = dict(kwargs)
    while True:
        subrun_fs.leaf.create([macro_idx, micro_idx])
        path = subrun_fs.leaf.path([macro_idx, micro_idx])

        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            inp_str, out_str = elstruct.run.direct(
                input_writer, script_str, path,
                geom=geom, charge=chg, mult=mul, method=method,
                basis=basis, prog=prog, **kwargs_)

        error_vals = [elstruct.reader.has_error_message(prog, error, out_str)
                      for error in errors]

        if not any(error_vals):
            # success
            break
        elif not moldr.optsmat.is_exhausted(options_mat):
            # try again
            micro_idx += 1
            error_row_idx = error_vals.index(True)
            kwargs_ = moldr.optsmat.updated_kwargs(kwargs, options_mat)
            options_mat = moldr.optsmat.advance(error_row_idx, options_mat)
        else:
            # failure
            warnings.resetwarnings()
            warnings.warn("elstruct robust run failed; "
                          "last run was in, {}".format(path))
            break

    return inp_str, out_str
