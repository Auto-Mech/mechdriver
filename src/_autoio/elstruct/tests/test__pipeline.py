""" test elstruct writer/run/reader pipelines
"""
import warnings
import tempfile
import numpy
import automol
import elstruct


SCRIPT_DCT = {
    'cfour2': None,
    'gaussian09': None,
    'gaussian03': None,
    'gaussian16': None,
    'molpro2015': None,
    'molpro2021': None,
    'mrcc2018': None,
    'nwchem6': None,
    'orca4': None,
    'psi4': "#!/usr/bin/env bash\n"
            "psi4 -i run.inp -o run.out >> stdout.log &> stderr.log",
}

# Certain methods, though nominally supported, crash and cannot be tested
BAD_METHODS = (
    ('psi4', 'ccsd'),
    ('psi4', 'ccsd(t)'),
)


def test__energy():
    """ test the energy pipeline
    """
    basis = '6-31g'
    geo = (('O', (0.0, 0.0, -0.110)),
           ('H', (0.0, -1.635, 0.876)),
           ('H', (-0.0, 1.635, 0.876)))
    mult_vals = [1, 2]
    charge_vals = [0, 1]

    for prog in elstruct.writer.programs():
        # print(prog)
        for method in elstruct.program_methods(prog):
            # print(method)
            for mult, charge in zip(mult_vals, charge_vals):
                for orb_type in (
                        elstruct.program_method_orbital_types(
                            prog, method, singlet=(mult == 1))):

                    vals = _test_pipeline(
                        script_str=SCRIPT_DCT[prog],
                        writer=elstruct.writer.energy,
                        readers=(
                            elstruct.reader.energy,
                        ),
                        args=(prog, geo, charge, mult, method, basis),
                        kwargs={'orb_type': orb_type, 'memory': 8},
                        error=elstruct.Error.SCF_NOCONV,
                        error_kwargs={'scf_options': [
                            elstruct.option.specify(
                                elstruct.Option.Scf.MAXITER_, 2)
                        ]},
                    )
                    # Print the value for Psi4 since it was run and read
                    if prog == elstruct.par.Program.PSI4:
                        print('ene\n', vals)


def test__gradient():
    """ test the gradient pipeline
    """
    basis = 'sto-3g'
    geo = (('O', (0.0, 0.0, -0.110)),
           ('H', (0.0, -1.635, 0.876)),
           ('H', (-0.0, 1.635, 0.876)))
    mult_vals = [1, 2]
    charge_vals = [0, 1]

    for prog in elstruct.writer.gradient_programs():
        methods = list(elstruct.program_nondft_methods(prog))
        dft_methods = list(elstruct.program_dft_methods(prog))
        if dft_methods:
            methods.append(numpy.random.choice(dft_methods))

        for method in methods:
            if (prog, method) in BAD_METHODS:
                continue
            for mult, charge in zip(mult_vals, charge_vals):
                for orb_type in (
                        elstruct.program_method_orbital_types(
                            prog, method, singlet=(mult == 1))):

                    vals = _test_pipeline(
                        script_str=SCRIPT_DCT[prog],
                        writer=elstruct.writer.gradient,
                        readers=(
                            elstruct.reader.energy,
                            elstruct.reader.gradient,
                        ),
                        args=(prog, geo, charge, mult, method, basis),
                        kwargs={'orb_type': orb_type, 'memory': 8},
                    )
                    # Print the value for Psi4 since it was run and read
                    if prog == elstruct.par.Program.PSI4:
                        print('grad\n', vals)


def test__hessian():
    """ test the hessian pipeline
    """
    basis = 'sto-3g'
    geo = (('O', (0.0, 0.0, -0.110)),
           ('H', (0.0, -1.635, 0.876)),
           ('H', (-0.0, 1.635, 0.876)))
    mult_vals = [1, 2]
    charge_vals = [0, 1]

    for prog in elstruct.writer.hessian_programs():
        methods = list(elstruct.program_nondft_methods(prog))
        dft_methods = list(elstruct.program_dft_methods(prog))
        if dft_methods:
            methods.append(numpy.random.choice(dft_methods))

        for method in methods:
            # Psi4 segfaults for ccsd
            if (prog, method) in BAD_METHODS:
                continue
            for mult, charge in zip(mult_vals, charge_vals):
                for orb_type in (
                        elstruct.program_method_orbital_types(
                            prog, method, singlet=(mult == 1))):

                    vals = _test_pipeline(
                        script_str=SCRIPT_DCT[prog],
                        writer=elstruct.writer.hessian,
                        readers=(
                            elstruct.reader.energy,
                            elstruct.reader.hessian,
                        ),
                        args=(prog, geo, charge, mult, method, basis),
                        kwargs={'orb_type': orb_type, 'memory': 8},
                    )
                    # Print the value for Psi4 since it was run and read
                    if prog == elstruct.par.Program.PSI4:
                        print('hess\n', vals)


def test__optimization():
    """ test elstruct optimization writes and reads
    """
    method = 'hf'
    basis = 'sto-3g'
    geo = (('C', (None, None, None), (None, None, None), (None, None, None)),
           ('O', (0, None, None), ('R1', None, None), (2.6, None, None)),
           ('H', (0, 1, None), ('R2', 'A2', None), (2.0, 1.9, None)),
           ('H', (0, 1, 2), ('R3', 'A3', 'D3'), (2.0, 1.9, 2.1)),
           ('H', (0, 1, 2), ('R4', 'A4', 'D4'), (2.0, 1.9, 4.1)),
           ('H', (1, 0, 2), ('R5', 'A5', 'D5'), (1.8, 1.8, 5.2)))
    mult = 1
    charge = 0
    orb_type = 'R'
    frozen_coordinates = ('R5', 'A5', 'D3')
    ref_frozen_values = (1.8, 1.8, 2.1)
    for prog in elstruct.writer.optimization_programs():
        script_str = SCRIPT_DCT[prog]

        # MRCC2018 does not support constrained optimizations
        if prog != 'mrcc2018':
            opt_kwargs = {'orb_type': orb_type,
                          'memory': 8,
                          'frozen_coordinates': frozen_coordinates}
        else:
            opt_kwargs = {'orb_type': orb_type, 'memory': 8}

        vals = _test_pipeline(
            script_str=script_str,
            writer=elstruct.writer.optimization,
            readers=(
                elstruct.reader.energy,
                elstruct.reader.opt_geometry,
                elstruct.reader.opt_zmatrix,
                elstruct.reader.opt_zmatrices,
            ),
            args=(prog, geo, charge, mult, method, basis),
            kwargs=opt_kwargs,
            error=elstruct.Error.OPT_NOCONV,
            error_kwargs={'job_options': [
                elstruct.option.specify(
                    elstruct.Option.Opt.MAXITER_, 2)
            ]},
        )

        if script_str is not None:
            # check that the frozen coordinates didn't change
            zma = vals[-1]
            val_dct = automol.zmat.value_dictionary(zma)
            frozen_values = tuple(
                map(val_dct.__getitem__, frozen_coordinates))
            assert numpy.allclose(
                frozen_values, ref_frozen_values, rtol=1e-4)

        # Print the value for Psi4 since it was run and read
        if prog == elstruct.par.Program.PSI4:
            print('geom\n', vals)


def _test_pipeline(script_str, writer, readers,
                   args, kwargs, error=None, error_kwargs=None):
    """ pipe
    """

    read_vals = []
    prog, method = args[0], args[4]

    # for programs with no run test, ensure input file generated
    _ = writer(*args, **kwargs)
    if script_str is not None:
        script_str = SCRIPT_DCT[prog]
        run_dir = tempfile.mkdtemp()
        run_dir = 'tmp'
        inp_str, out_str = elstruct.run.direct(
            writer, script_str, run_dir, *args, **kwargs)
        assert elstruct.reader.has_normal_exit_message(prog, out_str)

        for i, reader in enumerate(readers):
            if i == 0:
                val = reader(prog, method, out_str)
            else:
                val = reader(prog, out_str)
            read_vals.append(val)

        if error is not None:
            run_dir = tempfile.mkdtemp()
            assert not elstruct.reader.has_error_message(prog, error,
                                                         out_str)

            err_kwargs = kwargs.copy()
            err_kwargs.update(error_kwargs)

            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                _, err_out_str = elstruct.run.direct(
                    writer, script_str, run_dir, *args, **err_kwargs)

            assert elstruct.reader.has_error_message(prog, error,
                                                     err_out_str)
    return read_vals


if __name__ == '__main__':
    test__energy()
    test__gradient()
    test__hessian()
    test__optimization()
