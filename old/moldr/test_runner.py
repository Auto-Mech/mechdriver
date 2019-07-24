""" test moldr.runner
"""
import tempfile
import pytest
import elstruct
import moldr
import moldr.optsmat

SCRIPT_DCT = {
    'psi4': "#!/usr/bin/env bash\n"
            "psi4 -i run.inp -o run.out >> stdout.log &> stderr.log",
    'g09': None,
    'molpro': None,
}

# SCRIPT_DCT = {
#     'psi4': "#!/usr/bin/env bash\n"
#             "psi4 -i run.inp -o run.out >> stdout.log &> stderr.log",
#     'g09': "#!/usr/bin/env bash\n"
#            "g09 run.inp run.out >> stdout.log &> stderr.log",
#     'molpro': "#!/usr/bin/env bash\n"
#               "molpro run.inp -o run.out >> stdout.log &> stderr.log",
# }


def test__feedback_optimization():
    """ test runner.feedback_optimization
    """
    method = 'hf'
    basis = 'sto-3g'
    geom = ((('C', (None, None, None), (None, None, None)),
             ('O', (0, None, None), ('R1', None, None)),
             ('H', (0, 1, None), ('R2', 'A2', None)),
             ('H', (0, 1, 2), ('R3', 'A3', 'D3')),
             ('H', (0, 1, 2), ('R4', 'A4', 'D4')),
             ('H', (1, 0, 2), ('R5', 'A5', 'D5'))),
            {'R1': 2.6, 'R2': 2.0, 'A2': 1.9,
             'R3': 2.0, 'A3': 1.9, 'D3': 2.1,
             'R4': 2.0, 'A4': 1.9, 'D4': 4.1,
             'R5': 1.8, 'A5': 1.8, 'D5': 5.2})
    mult = 1
    charge = 0
    orb_restricted = True
    frozen_coordinates = ('R5', 'A5', 'D3')

    for prog in elstruct.writer.optimization_programs():
        script_str = SCRIPT_DCT[prog]
        if script_str is not None:
            run_dir = tempfile.mkdtemp()
            print(prog, run_dir)
            moldr.runner.feedback_optimization(
                script_str, run_dir,
                geom, charge, mult, method, basis, prog,
                orb_restricted=orb_restricted,
                frozen_coordinates=frozen_coordinates,
                job_options=[
                    elstruct.option.specify(elstruct.Option.Opt.MAXITER_, 4)],
            )


def test__options_matrix_run():
    """ test runner.options_matrix_run
    """
    input_writer = elstruct.writer.optimization
    method = 'hf'
    basis = 'sto-3g'
    geom = ((('O', (None, None, None), (None, None, None)),
             ('H', (0, None, None), ('R1', None, None)),
             ('H', (0, 1, None), ('R2', 'A2', None))),
            {'R1': 1.83114, 'R2': 1.83115, 'A2': 1.81475845})
    mult = 1
    charge = 0
    kwargs = {
        'scf_options': [
            elstruct.option.specify(elstruct.Option.Scf.MAXITER_, 2)],
        'job_options': [
            elstruct.option.specify(elstruct.Option.Opt.MAXITER_, 2)]}

    errors = [elstruct.Error.SCF_NOCONV, elstruct.Error.OPT_NOCONV]
    options_mat = [
        [{'scf_options': [elstruct.Option.Scf.Guess.CORE]},
         {'scf_options': [
             elstruct.option.specify(elstruct.Option.Scf.MAXITER_, 100)]}],
        [{'job_options': [
            elstruct.option.specify(elstruct.Option.Opt.MAXITER_, 100)]}]]
    bad_options_mat = [[{'scf_options': [elstruct.Option.Scf.Guess.CORE]}],
                       [{'job_options': []}]]

    for prog in elstruct.writer.optimization_programs():
        print()
        print(prog, method)
        if prog in SCRIPT_DCT:
            script_str = SCRIPT_DCT[prog]
            if script_str is not None:
                run_dir = tempfile.mkdtemp()
                print(run_dir)
                _, out_str = moldr.runner.options_matrix_run(
                    input_writer, script_str, run_dir,
                    geom, charge, mult, method, basis, prog,
                    errors=errors, options_mat=options_mat, **kwargs
                )

                assert elstruct.reader.has_normal_exit_message(prog, out_str)

                ene = elstruct.reader.energy(prog, method, out_str)
                zma = elstruct.reader.opt_zmatrix(prog, out_str)
                print(ene)
                print(zma)

                # make sure it throws a warning when it fails
                with pytest.warns(UserWarning):
                    moldr.runner.options_matrix_run(
                        input_writer, script_str, run_dir,
                        geom, charge, mult, method, basis, prog,
                        errors=errors, options_mat=bad_options_mat, **kwargs
                    )


if __name__ == '__main__':
    # test__options_matrix_run()
    test__feedback_optimization()
