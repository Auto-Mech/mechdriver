""" Set runtime options and job submission strings for
    electronic structure package calculations.
"""

import elstruct
from autorun import SCRIPT_DCT


def qchem_params(method_dct, job=None):
    """ Build kwargs and script string

        Function assumes a single-point energy job is none
        is given.

        :param method_dct: 
        :type method_dct: dict[str:
    """

    prog = method_dct.get('program', None)

    # Build the defaul values
    ret = INI_PARAM_BUILD_DCT[prog](method_dct, job=job)

    # Alter with the input method_dct (a massive pain...)
    # opt_kwargs.update(method_dct)
    # kwargs.update(method_dct)

    return ret


def _gaussian(method_dct, job=None):
    """ Build params for Gaussian
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
    """ Build params for Gaussian
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
    """ Build params for Gaussian
    """

    # Job unneeded for now
    _, _ = method_dct, job
    memory = method_dct.get('memory', 20)

    # Build the submission script string
    script_str = SCRIPT_DCT['psi4']

    # Build the options dictionary
    kwargs = {
        'memory': memory,
        # 'mol_options': ['no_symmetry']
        'mol_options': ['symmetry c1']
    }

    return script_str, kwargs


INI_PARAM_BUILD_DCT = {
    elstruct.par.Program.GAUSSIAN09: _gaussian,
    elstruct.par.Program.GAUSSIAN16: _gaussian,
    elstruct.par.Program.MOLPRO2015: _molpro,
    elstruct.par.Program.PSI4: _psi4,
}
