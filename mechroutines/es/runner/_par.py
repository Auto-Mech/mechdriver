""" Set runtime options and job submission strings for
    electronic structure package calculations.

    The options will be determined by the user input and
    will correspond to the default otherwise.

    Constructs `kwargs` dictionaries that are passed
    to elstruct package functions to write electronic structure input files.

    Constructs Shell submission sctings
"""

import elstruct
import automol
from autorun import SCRIPT_DCT


def qchem_params(method_dct, job=None):
    """ Build the kwargs dictionary and Shell submission script string to
        be used to write and run the electronic structure job.

        :param method_dct:
        :type method_dct: dict[str: obj]
        :param job: elstronic structure calculation
        :type job: str
        :rtype: (dict[str:tuple(str)], str)
    """

    prog = method_dct.get('program', None)

    # Build the defaul values
    ret = INI_PARAM_BUILD_DCT[prog](method_dct, job=job)

    # Alter with the input method_dct (a massive pain...)
    # opt_kwargs.update(method_dct)
    # kwargs.update(method_dct)

    return ret


def _gaussian(method_dct, job=None):
    """ Build kwargs dictionary and Shell submission script for Gaussian jobs.

        :param method_dct:
        :type method_dct: dict[str: obj]
        :param job: elstronic structure calculation
        :type job: str
        :rtype: (dict[str:tuple(str)], str)
    """

    # Set the options
    nprocs = method_dct.get('nprocs', 10)
    memory = method_dct.get('memory', 20)
    nprocs = nprocs if nprocs is not None else 10
    memory = memory if memory is not None else 20

    method = method_dct.get('method')

    # Build the submission script string
    script_str = SCRIPT_DCT['gaussian09']

    # Build the options dictionary
    machine_options = ['%NProcShared={}'.format(nprocs)]

    gen_lines = method_dct.get('gen_lines', {})
    if not gen_lines:
        if elstruct.par.Method.is_dft(method):
            gen_lines = {1: ['# int=ultrafine']}

    kwargs = {
        'memory': memory,
        'machine_options': machine_options,
        'gen_lines': gen_lines,
    }

    if job == elstruct.Job.OPTIMIZATION:
        kwargs.update({
            'feedback': True,
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
        })
    elif job in (elstruct.Job.IRCF, elstruct.Job.IRCR):
        irc_pts = method_dct.get('nirc', 10)
        irc_step = 3
        irc_dir = 'forward' if job == elstruct.Job.IRCF else 'reverse'
        kwargs.update({
            'job_options': [
                'calcfc',
                'stepsize={}'.format(irc_step),
                'maxpoints={}'.format(irc_pts),
                irc_dir
            ]
        })

    return script_str, kwargs


def _molpro(method_dct, job=None):
    """ Build kwargs dictionary and Shell submission script for Molpro jobs

        :param method_dct:
        :type method_dct: dict[str: obj]
        :param job: elstronic structure calculation
        :type job: str
        :rtype: (dict[str:tuple(str)], str)
    """

    # Pull stuff from the method_dct
    method = method_dct.get('method')
    if method in ('caspt2', 'caspt2c', 'caspt2i'):
        nprocs = method_dct.get('nprocs', 4)
        memory = method_dct.get('memory', 10)
        nprocs = nprocs if nprocs is not None else 4
        memory = memory if memory is not None else 10
    else:
        nprocs = method_dct.get('nprocs', 4)
        memory = method_dct.get('memory', 20)
        nprocs = nprocs if nprocs is not None else 4
        memory = memory if memory is not None else 10

    # Build the script string
    if method in ('caspt2c', 'caspt2i'):
        script_str = SCRIPT_DCT['molpro2015_mppx'].format(nprocs)
    else:
        script_str = SCRIPT_DCT['molpro2015'].format(nprocs)

    # Build the kwargs
    kwargs = {
        'memory': memory,
        # 'mol_options': ['no_symmetry'],
        'mol_options': ['nosym'],
    }

    corr_options = ['maxit=100']
    if method in ('caspt2', 'caspt2c', 'caspt2i'):
        corr_options.append('shift=0.2')
        if method == 'caspt2i':
            corr_options.append('ipea=0.25')
    kwargs.update({'corr_options': corr_options})

    if job == elstruct.Job.OPTIMIZATION:
        kwargs.update({
            'feedback': True,
            'errors': [
                elstruct.Error.OPT_NOCONV
            ],
            'options_mat': [
                [{'job_options': ['numhess=0']},
                 {'job_options': ['numhess=10']},
                 {'job_options': ['numhess=1']}]
            ],
        })

    return script_str, kwargs


def _psi4(method_dct, job=None):
    """ Build kwargs dictionary and Shell submission script for Psi4 jobs.

        :param method_dct:
        :type method_dct: dict[str: obj]
        :param job: elstronic structure calculation
        :type job: str
        :rtype: (dict[str:tuple(str)], str)
    """

    # Job unneeded for now
    _, _ = method_dct, job
    method = method_dct.get('method')
    memory = method_dct.get('memory', 20)

    # Build the submission script string
    script_str = SCRIPT_DCT['psi4']

    # Build the options dictionary
    kwargs = {
        'memory': memory,
        # 'mol_options': ['no_symmetry']
        'mol_options': ['symmetry c1']
    }
    if job == elstruct.Job.OPTIMIZATION:
        kwargs.update({
            'feedback': False,
            'errors': [
                elstruct.Error.SYMM_NOFIND
            ],
            'options_mat': [
                [{'job_options': ['set intrafrag_step_limit 0.08']},
                 {'job_options': ['set intrafrag_step_limit 0.03']}]
            ],
        })

    if elstruct.par.Method.is_dft(method):
        kwargs.update({
            'scf_options': [
                # 'set dft_spherical_points 590',
                # 'set dft_radial_points 99',
                'set dft_basis_tolerance 1.0E-11',
                # 'set dft_pruning_scheme robust'
                'set opt_coordinates delocalized',  # for scan
                'set ensure_bt_convergence true'
            ]})

    return script_str, kwargs


INI_PARAM_BUILD_DCT = {
    elstruct.par.Program.GAUSSIAN09: _gaussian,
    elstruct.par.Program.GAUSSIAN16: _gaussian,
    elstruct.par.Program.MOLPRO2015: _molpro,
    elstruct.par.Program.PSI4: _psi4,
}


# Unique constructors for programs
def molpro_opts_mat(spc_info, geo):
    """ prepare the errors and options mat to perform successive
        single-point energy calculations in Molpro when the RHF fails to
        converge. This currently only works for doublets.

        :param spc_info:
        :type spc_info:
        :param geo: input molecular geometry object
        :type geo: automol.geom object
        :rtype: (tuple(str), tuple(dict[str: str]))
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
