""" electronic structure parameter list
"""
import elstruct


def run_qchem_par(prog):
    """ dictionary of parameters for different electronic structure codes
    """

    if prog == 'g09':
        script_str = ("#!/usr/bin/env bash\n"
                      "g09 run.inp run.out >> stdout.log &> stderr.log")
        opt_script_str = script_str
        kwargs = {
            'memory': 20,
            'machine_options': ['%NProcShared=15'],
            'gen_lines': ['# int=ultrafine'],
        }
        opt_kwargs = {
            'memory': 20,
            'machine_options': ['%NProcShared=15'],
            'gen_lines': ['# int=ultrafine'],
            'feedback': True,
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
        script_str = ("#!/usr/bin/env bash\n"
                      "psi4 -i run.inp -o run.out >> stdout.log &> stderr.log")
        opt_script_str = script_str
        kwargs = {}
        opt_kwargs = {}

    if prog == 'molpro':
        script_str = ("#!/usr/bin/env bash\n"
                      "molpro -n 8 run.inp -o run.out >> stdout.log &> stderr.log")
        opt_script_str = ("#!/usr/bin/env bash\n"
                          "molpro --mppx -n 12 run.inp -o run.out >> stdout.log &> stderr.log")
        kwargs = {
            'memory': 10,
        }
        opt_kwargs = {
            'memory': 10,
        }

    if prog == 'qchem':
        script_str = ("#!/usr/bin/env bash\n"
                      "molpro -i run.inp -o run.out >> stdout.log &> stderr.log")
        opt_script_str = script_str
        kwargs = {
            'memory': 10,
        }
        opt_kwargs = {}

    if prog == 'cfour':
        script_str = ("#!/usr/bin/env bash\n"
                      "molpro -i run.inp -o run.out >> stdout.log &> stderr.log")
        opt_script_str = script_str
        kwargs = {
            'memory': 10,
        }
        opt_kwargs = {}

    if prog == 'orca':
        script_str = ("#!/usr/bin/env bash\n"
                      "molpro -i run.inp -o run.out >> stdout.log &> stderr.log")
        opt_script_str = script_str
        kwargs = {
            'memory': 10,
        }
        opt_kwargs = {}

    return  script_str, opt_script_str, kwargs, opt_kwargs
