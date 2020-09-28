""" Set runtime parameters
"""

import automol
import elstruct


def qchem_params(prog, method):  # , saddle=False):
    """ dictionary of parameters for different electronic structure codes
    """
    if prog == 'gaussian09':
        sp_script_str = ("#!/usr/bin/env bash\n"
                         "g09 run.inp run.out >> stdout.log &> stderr.log")
        opt_script_str = sp_script_str
        kwargs = {
            'memory': 20,
            'machine_options': ['%NProcShared=10'],
            'gen_lines': {1: ['# int=ultrafine']},
        }
        opt_kwargs = {
            'memory': 20,
            'machine_options': ['%NProcShared=10'],
            'gen_lines': {1: ['# int=ultrafine']},
            # 'gen_lines': {1: ['# int=superfine']},
            'feedback': True,
            # 'job_options': ['verytight'],
            # 'job_options': ['verytight'],
            'errors': [
                elstruct.Error.OPT_NOCONV
            ],
            'options_mat': [
                [{},
                 {},
                 {},
                 {'job_options': ['calcfc']},
                 {'job_options': ['calcfc']},
                 {'job_options': ['calcall']}]
            ],
        }

    if prog == 'psi4':
        sp_script_str = (
            "#!/usr/bin/env bash\n"
            "psi4 -i run.inp -o run.out >> stdout.log &> stderr.log"
        )
        opt_script_str = sp_script_str
        kwargs = {}
        opt_kwargs = {}

    if prog == 'molpro2015':
        sp_script_str = (
            "#!/usr/bin/env bash\n"
            "molpro -n 6 run.inp -o run.out >> stdout.log &> stderr.log"
        )
        if method == 'caspt2':
            opt_script_str = (
                "#!/usr/bin/env bash\n"
                "molpro -n 8 run.inp -o run.out "
                "--nouse-logfile --no-xml-output >> "
                "stdout.log &> stderr.log"
                # "molpro -n 8 run.inp -o run.out >> stdout.log &> stderr.log"
            )
        else:
            opt_script_str = (
                "#!/usr/bin/env bash\n"
                "molpro --mppx -n 12 run.inp -o run.out "
                "--nouse-logfile --no-xml-output >> "
                # "molpro --mppx -n 12 run.inp -o run.out >> "
                "stdout.log &> stderr.log"
            )
        if method in ('caspt2', 'caspt2c'):
            kwargs = {
                'memory': 10,
                'corr_options': ['shift=0.2'],
                'mol_options': ['nosym'],
            }
            opt_kwargs = {
                'memory': 5,
                'corr_options': ['shift=0.2'],
                'mol_options': ['nosym'],
                'feedback': True,
                'errors': [
                    elstruct.Error.OPT_NOCONV
                ],
                'options_mat': [
                    [{'job_options': ['numhess=0']},
                     {'job_options': ['numhess=10']},
                     {'job_options': ['numhess=1']}]
                ],
            }
        elif method == 'caspt2i':
            kwargs = {
                'memory': 10,
                'corr_options': ['shift=0.2', 'ipea=0.25'],
                'mol_options': ['nosym'],
            }
            opt_kwargs = {
                'memory': 5,
                'corr_options': ['shift=0.2', 'ipea=0.25'],
                'mol_options': ['nosym'],
                'feedback': True,
                'errors': [
                    elstruct.Error.OPT_NOCONV
                ],
                'options_mat': [
                    [{'job_options': ['numhess=0']},
                     {'job_options': ['numhess=10']},
                     {'job_options': ['numhess=1']}]
                ],
            }
        else:
            kwargs = {
                'memory': 20,
                'corr_options': ['maxit=100'],
                'mol_options': ['nosym'],
            }
            opt_kwargs = {
                'memory': 5,
                'mol_options': ['nosym'],
                'corr_options': ['maxit=100'],
                'feedback': True,
                'errors': [
                ],
                'options_mat': [
                    [{'job_options': ['numhess=0']},
                     {'job_options': ['numhess=10']},
                     {'job_options': ['numhess=1']}]
                ],
            }

    return sp_script_str, opt_script_str, kwargs, opt_kwargs
