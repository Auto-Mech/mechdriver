""" elstruct runners (formerly elcarro)
"""
import warnings
import numpy
from qcelemental import constants as qcc
import automol
import elstruct
import autofile
import moldr.optsmat


def options_matrix_optimization(script_str, prefix,
                                geom, charge, mult, method, basis, prog,
                                errors=(), options_mat=(), feedback=False,
                                frozen_coordinates=(),
                                freeze_dummy_atoms=True,
                                saddle=False,
                                kickoff_saddle=False,
                                **kwargs):
    """ try several sets of options to generate an output file

    :returns: the input string and the output string
    :rtype: (str, str)
    """
    assert len(errors) == len(options_mat)

    success = False

    subrun_ds = autofile.system.series.subrun_leaf()
    max_macro_idx, _ = max(subrun_ds.dir.existing(prefix), default=(-1, -1))
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
        subrun_ds.dir.create(prefix, (macro_idx, micro_idx))
        path = subrun_ds.dir.path(prefix, (macro_idx, micro_idx))

        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            inp_str, out_str = elstruct.run.direct(
                elstruct.writer.optimization, script_str, path,
                geom=geom, charge=charge, mult=mult, method=method,
                basis=basis, prog=prog, frozen_coordinates=frozen_coordinates,
                saddle=saddle, **kwargs_)

        error_vals = [elstruct.reader.has_error_message(prog, error, out_str)
                      for error in errors]

        if not any(error_vals):
            # success
            success = True
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

    if success and not saddle and kickoff_saddle:
        geo = elstruct.reader.opt_geometry(prog, out_str)
        _, hess_out_str = options_matrix_run(
            elstruct.writer.hessian, script_str, prefix, geo, charge, mult,
            method, basis, prog)
        hess = elstruct.reader.hessian(prog, hess_out_str)
        freqs = elstruct.util.harmonic_frequencies(geo, hess, project=True)
        norm_coos = elstruct.util.normal_coordinates(geo, hess, project=True)

        # if there's an imaginary frequency, try again after displacing along
        # the mode
        if freqs[0] < -1e-1:
            im_norm_coo = numpy.array(norm_coos)[:, 0]
            disp_len = 0.1 * qcc.conversion_factor('angstrom', 'bohr')
            disp_xyzs = numpy.reshape(im_norm_coo, (-1, 3)) * disp_len
            geo = automol.geom.displace(geo, disp_xyzs)
            if automol.geom.is_valid(geom):
                geom = geo
            else:
                assert automol.zmatrix.is_valid(geom)
                vma = automol.zmatrix.var_(geom)
                geom = automol.vmatrix.zmatrix_from_geometry(vma, geo)

            options_matrix_optimization(
                script_str, prefix,
                geom, charge, mult, method, basis, prog,
                errors=errors, options_mat=options_mat, feedback=feedback,
                frozen_coordinates=frozen_coordinates,
                freeze_dummy_atoms=freeze_dummy_atoms,
                saddle=saddle,
                kickoff_saddle=False,
                **kwargs)

    return inp_str, out_str


def options_matrix_run(input_writer, script_str, prefix,
                       geom, charge, mult, method, basis, prog,
                       errors=(), options_mat=(),
                       **kwargs):
    """ try several sets of options to generate an output file

    :returns: the input string and the output string
    :rtype: (str, str)
    """
    assert len(errors) == len(options_mat)

    subrun_ds = autofile.system.series.subrun_leaf()
    max_macro_idx, _ = max(subrun_ds.dir.existing(prefix), default=(-1, -1))
    macro_idx = max_macro_idx + 1
    micro_idx = 0

    kwargs_ = dict(kwargs)
    while True:
        subrun_ds.dir.create(prefix, (macro_idx, micro_idx))
        path = subrun_ds.dir.path(prefix, (macro_idx, micro_idx))

        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            inp_str, out_str = elstruct.run.direct(
                input_writer, script_str, path,
                geom=geom, charge=charge, mult=mult, method=method,
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
