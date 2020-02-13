""" Set runtime parameters
"""

import automol
import elstruct


def run_qchem_par(prog, method):  # , saddle=False):
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
        # if spin:
        #     kwargs['scf_options'] = [
        #         elstruct.option.specify(
        #             elstruct.Option.Scf.Guess.MIX)]
        #     opt_kwargs['scf_options'] = [
        #         elstruct.option.specify(
        #             elstruct.Option.Scf.Guess.MIX)]

        # doesn't quite work for gen_lines case right now
        # if saddle:
        #    opt_kwargs = {
        #        'memory': 20,
        #        'machine_options': ['%NProcShared=10'],
        #        'gen_lines': {1: ['# int=ultrafine']},
        #        'feedback': True,
        #        'errors': [
        #            elstruct.Error.OPT_NOCONV
        #        ],
        #        'options_mat': [
        #            [{'gen_lines': [{1: ['# int=superfine']}]},
        #            {'gen_lines': [{1: ['# int=superfine']}]},
        #            {'gen_lines': [{1: ['# int=superfine']}]},
        #            {'gen_lines': [{1: ['# int=superfine']}]}]
        #        ],
        #    }

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
            "molpro -n 4 run.inp -o run.out >> stdout.log &> stderr.log"
        )
        if method == 'caspt2':
            opt_script_str = (
                "#!/usr/bin/env bash\n"
                "molpro -n 8 run.inp -o run.out",
                "--nouse-logfile --no-xml-output >> "
                "stdout.log &> stderr.log"
                # "molpro -n 8 run.inp -o run.out >> stdout.log &> stderr.log"
            )
        else:
            opt_script_str = (
                "#!/usr/bin/env bash\n"
                "molpro --mppx -n 12 run.inp -o run.out",
                "--nouse-logfile --no-xml-output >> "
                # "molpro --mppx -n 12 run.inp -o run.out >> "
                "stdout.log &> stderr.log"
            )
        if method in ('caspt2', 'caspt2c'):
            kwargs = {
                'memory': 15,
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
                'memory': 16,
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

    # if prog == 'qchem':
    #     sp_script_str = (
    #         "#!/usr/bin/env bash\n"
    #         "molpro -i run.inp -o run.out >> stdout.log &> stderr.log")
    #     opt_script_str = sp_script_str
    #     kwargs = {
    #         'memory': 20,
    #     }
    #     opt_kwargs = {}

    # if prog == 'cfour':
    #     sp_script_str = (
    #         "#!/usr/bin/env bash\n"
    #         "molpro -i run.inp -o run.out >> stdout.log &> stderr.log")
    #     opt_script_str = sp_script_str
    #     kwargs = {
    #         'memory': 20,
    #     }
    #     opt_kwargs = {}

    # if prog == 'orca':
    #     sp_script_str = (
    #         "#!/usr/bin/env bash\n"
    #         "molpro -i run.inp -o run.out >> stdout.log &> stderr.log")
    #     opt_script_str = sp_script_str
    #     kwargs = {
    #         'memory': 20,
    #     }
    #     opt_kwargs = {}

    return sp_script_str, opt_script_str, kwargs, opt_kwargs


def set_molpro_options_mat(spc_info, geo):
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
