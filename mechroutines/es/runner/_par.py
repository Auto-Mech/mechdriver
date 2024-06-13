""" Set runtime options and job submission strings for
    electronic structure package calculations.

    The options will be determined by the user input and
    will correspond to the default otherwise.

    Constructs `kwargs` dictionaries that are passed
    to elstruct package functions to write electronic structure input files.

    Constructs BASH submission sctings
"""

import elstruct
import automol
from autorun import SCRIPT_DCT


def qchem_params(method_dct, job=None, geo=None, spc_info=None):
    """ Build the kwargs dictionary and BASH submission script string to
        be used to write and run the electronic structure job.

        :param method_dct:
        :type method_dct: dict[str: obj]
        :param job: elstronic structure calculation
        :type job: str
        :rtype: (dict[str:tuple(str)], str)
    """

    prog = method_dct.get('program', None)

    # Build the defaul values
    ret = INI_PARAM_BUILD_DCT[prog](
        method_dct, prog,
        job=job, geo=geo, spc_info=spc_info)

    # Alter with the input method_dct (a massive pain...)
    # opt_kwargs.update(method_dct)
    # kwargs.update(method_dct)

    return ret


def _gaussian(method_dct, prog, job=None, geo=None, spc_info=None):
    """ Build kwargs dictionary and BASH submission script for Gaussian jobs.

        :param method_dct:
        :type method_dct: dict[str: obj]
        :param job: elstronic structure calculation
        :type job: str
        :rtype: (dict[str:tuple(str)], str)
    """

    _, _ = geo, spc_info

    # Set the options
    nprocs = method_dct.get('nprocs', 9)
    memory = method_dct.get('mem', 20)
    nprocs = nprocs if nprocs is not None else 9
    memory = memory if memory is not None else 20

    method = method_dct.get('method')

    # Build the submission script string
    script_str = SCRIPT_DCT[prog]

    # Build the options dictionary
    machine_options = [f'%NProcShared={nprocs}']

    gen_lines = method_dct.get('gen_lines', {})
    if not gen_lines:
        if elstruct.par.Method.is_dft(method):
            # gen_lines = {1: ['# int=superfine']}
            gen_lines = {1: ['# int='+ method_dct.get('grid')]}

    if job == 'tightfreq':
        gen_lines = {1: ['# int=superfine']}
        job = elstruct.Job.HESSIAN

    kwargs = {
        'memory': memory,
        'machine_options': machine_options,
        'gen_lines': gen_lines,
    }
    if job == 'tightopt':
        kwargs.update({
            'gen_lines': {1: ['# int=superfine']},
            'job_options': ['Tight'],
            'feedback': True,
            'errors': [
                elstruct.Error.OPT_NOCONV
            ],
            'options_mat': [
                [{'job_options': ['Tight']},
                 {'job_options': ['Tight']},
                 {'job_options': ['Tight']},
                 {'job_options': ['Tight', 'calcfc']},
                 {'job_options': ['Tight', 'calcfc']},
                 {'job_options': ['Tight', 'calcall']}],
            ],
        })

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
                f'stepsize={irc_step}',
                f'maxpoints={irc_pts}',
                irc_dir
            ]
        })

    return script_str, kwargs


def _molpro(method_dct, prog, job=None, geo=None, spc_info=None):
    """ Build kwargs dictionary and BASH submission script for Molpro jobs

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
        memory = method_dct.get('mem', 20)
        econv = method_dct.get('econv', 1.0e-6)
        gconv = method_dct.get('gconv', 3.0e-4)
        nprocs = nprocs if nprocs is not None else 4
        memory = memory if memory is not None else 10
        econv = econv if econv is not None else 1.0e-6
        gconv = gconv if gconv is not None else 3.0e-4
    else:
        nprocs = method_dct.get('nprocs', 4)
        memory = method_dct.get('mem', 20)
        econv = method_dct.get('econv', 1.0e-6)
        gconv = method_dct.get('gconv', 3.0e-4)
        nprocs = nprocs if nprocs is not None else 4
        memory = memory if memory is not None else 20
        econv = econv if econv is not None else 1.0e-6
        gconv = gconv if gconv is not None else 3.0e-4

    scf_econv_line = f'energy={econv:.1E}'.replace('E', 'd')
    corr_econv_line = f'energy={econv:.1E}'.replace('E', 'd')
    gconv_line = f'gradient={gconv:.1E}'.replace('E', 'd')

    # Build the script string
    prog = prog+'_mppx' if method_dct['mppx'] else prog
    script_str = SCRIPT_DCT[prog].format(nprocs)

    # Set glob thresholds line
    thrsh_line = f'gthresh,orbital={econv:.1E}'.replace('E', 'd')
    if method_dct.get('tight_integral'):
        thrsh_line += ',oneint=1.0d-16,twoint=1.0d-16,compress=1.0d-13'

    # Build the kwargs
    kwargs = {
        'memory': memory,
        'mol_options': ['nosym'],
        'scf_options': [scf_econv_line, 'maxit=150']
        # 'scf_options': [scf_econv_line, 'so-sci', 'maxit=150']
        # 'scf_options': [scf_econv_line, 'maxit=5']
    }

    corr_options = [corr_econv_line, 'maxit=100']
    if method in ('caspt2', 'caspt2c', 'caspt2i'):
        corr_options.append('shift=0.2')
        if method == 'caspt2i':
            corr_options.append('ipea=0.25')
    kwargs.update({'corr_options': corr_options})

    if job == elstruct.Job.OPTIMIZATION:
        kwargs.update({
            'gen_lines': {1: [thrsh_line]},
            'job_options': [gconv_line],
            'feedback': True,
            'errors': [
                elstruct.Error.OPT_NOCONV
            ],
            'options_mat': [
                [{'job_options': [gconv_line, 'numhess=0']},
                 {'job_options': [gconv_line, 'numhess=10']},
                 {'job_options': [gconv_line, 'numhess=1']}]
            ],
        })
    else:
        # Get the nelectrons, spins, and orbitals for the wf card
        formula = automol.geom.formula(geo)
        elec_count = automol.form.electron_count(formula)
        two_spin = spc_info[2] - 1
        num_act_elc = two_spin
        num_act_orb = num_act_elc
        closed_orb = (elec_count - num_act_elc) // 2
        occ_orb = closed_orb + num_act_orb

        # Build strings UHF and CASSCF wf card and set the errors and options
        uhf_str = (
            f"{{uhf,maxit=150;wf,{elec_count},1,{two_spin};orbprint,3}}"
        )
        cas_str = (
            "{casscf,maxit=40;"
            f"closed,{closed_orb};occ,{occ_orb};"
            f"wf,{elec_count},1,{two_spin};canonical;orbprint,3}}"
        )

        kwargs.update({
            'gen_lines': {1: [thrsh_line]},
            'feedback': True,
            'errors': [
                elstruct.Error.SCF_NOCONV
            ],
            'options_mat': [
                [{'gen_lines': {1: [thrsh_line], 2: [uhf_str]}},
                 {'gen_lines': {1: [thrsh_line], 2: [cas_str]}},
                 {'gen_lines': {1: [thrsh_line], 2: [cas_str]}}]
             ],
        })

    return script_str, kwargs


def _psi4(method_dct, prog, job=None, geo=None, spc_info=None):
    """ Build kwargs dictionary and BASH submission script for Psi4 jobs

        :param method_dct:
        :type method_dct: dict[str: obj]
        :param job: elstronic structure calculation
        :type job: str
        :rtype: (dict[str:tuple(str)], str)
    """

    _, _ = geo, spc_info

    # Job unneeded for now
    method = method_dct.get('method')
    nprocs = method_dct.get('nprocs', 8)
    memory = method_dct.get('mem', 10)
    nprocs = nprocs if nprocs is not None else 8
    memory = memory if memory is not None else 10

    # Build the submission script string
    script_str = SCRIPT_DCT[prog]

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
                # 'set opt_coordinates delocalized',  # for scan
                # 'set ensure_bt_convergence true'
            ]})

    return script_str, kwargs


def _qchem(method_dct, prog, job=None, geo=None, spc_info=None):
    """ Build kwargs dictionary and BASH submission script for Gaussian jobs.

        :param method_dct:
        :type method_dct: dict[str: obj]
        :param job: elstronic structure calculation
        :type job: str
        :rtype: (dict[str:tuple(str)], str)
    """

    _, _ = geo, spc_info

    # Set the options
    nprocs = method_dct.get('nprocs', 8)
    memory = method_dct.get('mem', 20)
    nprocs = nprocs if nprocs is not None else 8
    memory = memory if memory is not None else 20

    method = method_dct.get('method')

    # Build the submission script string
    script_str = SCRIPT_DCT[prog].format(nprocs)

    kwargs = {
        'memory': memory,
    }

    return script_str, kwargs


def _orca(method_dct, prog, job=None, geo=None, spc_info=None):
    """ Build kwargs dictionary and BASH submission script for Orca jobs.

        :param method_dct:
        :type method_dct: dict[str: obj]
        :param job: elstronic structure calculation
        :type job: str
        :rtype: (dict[str:tuple(str)], str)
    """

    _, _ = geo, spc_info

    # Set the options
    nprocs = method_dct.get('nprocs', 8)
    memory = method_dct.get('mem', 20)
    nprocs = nprocs if nprocs is not None else 8
    memory = memory if memory is not None else 20

    method = method_dct.get('method')

    # Build the submission script string
    script_str = SCRIPT_DCT[prog].format(nprocs)

    kwargs = {
        'memory': memory,
    }

    return script_str, kwargs


INI_PARAM_BUILD_DCT = {
    elstruct.Program.GAUSSIAN09: _gaussian,
    elstruct.Program.GAUSSIAN16: _gaussian,
    elstruct.Program.MOLPRO2021: _molpro,
    elstruct.Program.MOLPRO2015: _molpro,
    elstruct.Program.PSI4: _psi4,
    elstruct.Program.QCHEM5: _qchem,
    elstruct.Program.ORCA4: _orca,
}
