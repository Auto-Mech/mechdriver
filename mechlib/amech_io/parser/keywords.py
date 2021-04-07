"""
  Libraries to check for allowed and supported keywords
"""

from phydat import phycon


# Run Keywords
RUN_INP_DEFAULT_DCT = {
    'inp_mech': (('chemkin'), 'chemkin'),
    'inp_spc': (('csv',), 'csv'),
    'out_mech': (('chemkin'), 'chemkin'),
    'out_spc': (('csv',), 'csv'),
    'print_mech': ((True, False), False),
    'print_debug': ((True, False), False),
    'run_prefix': (None, None),
    'save_prefix': (None, None)
}


# Electronic Structure Tasks
ES_TSK_SUPPORTED_DCT = {
    'spc': [
        'init_geom',
        'conf_pucker', 'conf_samp', 'conf_energy', 'conf_grad', 'conf_hess',
        'conf_vpt2', 'conf_prop', 'conf_opt',
        'hr_scan', 'hr_energy', 'hr_grad', 'hr_hess', 'hr_vpt2', 'hr_reopt',
        'tau_samp', 'tau_energy', 'tau_grad', 'tau_hess'],
    'ts': [
        'find_ts',
        'find_sadpt', 'find_molrad_vtst', 'find_radrad_vtst', 'find_vrctst',
        'conf_pucker', 'conf_samp', 'conf_energy', 'conf_grad', 'conf_hess',
        'conf_vpt2', 'conf_prop',
        'hr_scan', 'hr_energy', 'hr_grad', 'hr_hess', 'hr_vpt2', 'hr_reopt',
        'rpath_scan', 'rpath_energy', 'rpath_grad', 'rpath_hess'],
    'vdw': [
        'find',
        'conf_pucker', 'conf_samp', 'conf_energy', 'conf_grad', 'conf_hess']
}
ES_TSK_KEYWORDS_SUPPORTED_DCT = {
    'init_geom': ['runlvl', 'inplvl', 'retryfail', 'overwrite'],
    'find_ts': ['runlvl', 'inplvl', 'rxndirn',
                'var_splvl1', 'var_splvl2', 'var_scnlvl',
                'nobarrier', 'retryfail', 'overwrite'],
    'find_sadpt': ['runlvl', 'inplvl', 'rxndirn',
                   'nobarrier', 'retryfail', 'overwrite'],
    'find_molrad_vtst': ['runlvl', 'inplvl', 'rxndirn',
                         'var_splvl1', 'var_splvl2', 'var_scnlvl',
                         'nobarrier', 'retryfail', 'overwrite'],
    'find_radrad_vtst': ['runlvl', 'inplvl', 'rxndirn',
                         'var_splvl1', 'var_splvl2', 'var_scnlvl',
                         'nobarrier', 'pot_thresh',
                         'retryfail', 'overwrite'],
    'find_vrctst': ['runlvl', 'inplvl', 'rxndirn',
                    'var_splvl1', 'var_splvl2', 'var_scnlvl',
                    'nobarrier', 'retryfail', 'overwrite'],
    'conf_pucker': ['runlvl', 'inplvl', 'cnf_range', 'retryfail', 'overwrite', 'resave'],
    'conf_samp': ['runlvl', 'inplvl', 'cnf_range', 'retryfail', 'overwrite', 'resave'],
    'conf_energy': ['runlvl', 'inplvl', 'cnf_range', 'retryfail', 'overwrite'],
    'conf_grad': ['runlvl', 'inplvl', 'cnf_range', 'retryfail', 'overwrite'],
    'conf_hess': ['runlvl', 'inplvl', 'cnf_range', 'retryfail', 'overwrite'],
    'conf_vpt2': ['runlvl', 'inplvl', 'cnf_range', 'retryfail', 'overwrite'],
    'conf_prop': ['runlvl', 'inplvl', 'cnf_range', 'retryfail', 'overwrite'],
    'conf_opt': ['runlvl', 'inplvl', 'cnf_range', 'retryfail', 'overwrite'],
    'hr_scan': ['runlvl', 'inplvl', 'tors_model', 'resamp_min',
                'retryfail', 'overwrite'],
    'hr_grad': ['runlvl', 'inplvl', 'tors_model',
                'retryfail', 'overwrite'],
    'hr_hess': ['runlvl', 'inplvl', 'tors_model',
                'retryfail', 'overwrite'],
    'hr_energy': ['runlvl', 'inplvl', 'tors_model',
                  'retryfail', 'overwrite'],
    'hr_vpt2': ['runlvl', 'inplvl', 'tors_model',
                'retryfail', 'overwrite'],
    'hr_reopt': ['runlvl', 'inplvl', 'tors_model',
                 'retryfail', 'overwrite', 'hrthresh'],
    'tau_samp': ['runlvl', 'inplvl', 'retryfail', 'overwrite'],
    'tau_energy': ['runlvl', 'inplvl', 'retryfail', 'overwrite'],
    'tau_grad': ['runlvl', 'inplvl', 'retryfail', 'overwrite'],
    'tau_hess': ['runlvl', 'inplvl', 'hessmax', 'retryfail', 'overwrite'],
    'rpath_scan': ['runlvl', 'inplvl', 'rxncoord', 'retryfail', 'overwrite'],
    'rpath_energy': ['runlvl', 'inplvl', 'rxncoord', 'retryfail', 'overwrite'],
    'rpath_grad': ['runlvl', 'inplvl', 'rxncoord', 'retryfail', 'overwrite'],
    'rpath_hess': ['runlvl', 'inplvl', 'rxncoord', 'retryfail', 'overwrite'],
}
TRANS_TSK_SUPPORTED_DCT = {
    'spc': ['onedmin']
}
TRANS_TSK_KEYWORDS_SUPPORTED_DCT = {
    'onedmin': ['runlvl', 'inplvl', 'bath', 'pot',
                'njobs', 'nsamp',
                'smin', 'smax', 'conf',
                'retryfail', 'overwrite']
}
PRNT_TSK_SUPPORTED_DCT = {
    'spc': [
        'freqs', 'energy', 'geo', 'zmatrix', 'enthalpy', 'coeffs'],
    'ts': [
        'freqs'],
    'vdw': [
        'freqs']}
PRNT_TSK_KEYWORDS_SUPPORTED_DCT = {
    'freqs': ['geolvl', 'proplvl', 'nconfs', 'econfs', 'scale'],
    'energy': ['geolvl', 'proplvl', 'nconfs', 'econfs'],
    'geo': ['geolvl', 'proplvl', 'nconfs', 'econfs'],
    'zmatrix': ['geolvl', 'proplvl', 'nconfs', 'econfs'],
    'enthalpy': ['geolvl', 'proplvl', 'nconfs', 'econfs'],
    'coeffs': [],
    }


# TASK DICTIONARIES
# es tsk: (object type, (allowed values), default)  # use functions for weird
ES_TSK_DCT = {
    # all tasks
    'runlvl': (str, (), None),
    'inplvl': (str, (), None),
    'var_splvl1': (str, (), None),
    'var_splvl2': (str, (), None),
    'var_scnlvl': (str, (), None),
    'retryfail': (bool, (True, False), True),
    'overwrite': (bool, (True, False), False),
    # conformer tasks
    'cnf_range': (str, _cnf_string, 'min'),
    # tau tasks
    'hessmax': (int, (), 1000),
    # hr tasks
    'tors_model': (str, ('1dhr', '1dhrf', '1dhrfa', 'mdhr', 'mdhrv'), '1dhr'),
    'resamp_min': (bool, (True, False), False),
    'hrthresh': (float, (), -0.5),
    'potthresh': (float, (), 0.3),
    # reaction
    'rxncoord': (str, ('irc', 'auto'), 'auto'),
    'nobarrier': (str, ('pst', 'rpvtst', 'vrctst'), None)
}
TRANS_TSK_DCT = {
    'pot': (str, ('sphere',), 'lj_12_6'),
    'njobs': (int, (), 1),
    'nsamp': (int, (), 1),
    'smin': (float, (), 2.0),
    'smax': (float, (), 6.0),
    'conf': (str, ('sphere',), 'sphere'),
    # prob redundant with es
    'retryfail': (bool, (True, False), True),
    'overwrite': (bool, (True, False), False)
}
PRNT_TSK_DCT = {
    'geolvl': None,
    'proplvl': None,
    'nconfs': 'min',
    'econfs': 'min',
    'scale': None
    }


# Species keywords
SPC_REQUIRED_KEYWORDS = [
    'mult',
    'charge',
]
SPC_SUPPORTED_KEYWORDS = [
    'charge',
    'mult',
    'geom',
    'mc_nsamp',
    'tau_nsamp',
    'hind_inc',
    'tors_names',
    'elec_levels',
    'sym_factor',
    'inchi',
    'smiles',
    'geom',
    'kickoff',
    'ts_search',
    'active',
    'zma_idx',
    # etrans
    'smin',
    'smax',
    'etrans_nsamp',
    'kt_pst',
    'temp_pst',
    'n_pst',
    'lj',
    'edown'
]
SPC_DEFAULT_DCT = {
    'kickoff': (0.1, False),
    'hind_inc': 30.0*phycon.DEG2RAD,
    'mc_nsamp': (True, 12, 1, 3, 100, 25),
    'tau_nsamp': (True, 12, 1, 3, 100, 25),
    'lj': None,
    'edown': None
}
TS_DEFAULT_DCT = {**SPC_DEFAULT_DCT, **{
    'pst_params': (1.0, 6),
    'rxndirn': 'forw',
    'kt_pst': 4.0e-10,
    'temp_pst': 300.0,
    'n_pst': 6.0
}}


# Theory Keywords
THY_REQUIRED_KEYWORDS = [
    'program',
    'method',
    'basis',
    'orb_res'
]
THY_SUPPORTED_KEYWORDS = [
    'program',
    'method',
    'basis',
    'orb_res',
    'ncycles',
    'mem',
    'nprocs',
    'econv',
    'gconv'
]

THY = {
    'program': (str),
    'method': (str),
    'basis': (str),
    'orb_res': (str, ('RR', 'UU', 'RU')),
    'ncycles': (int),
    'mem': (float),
    'nprocs': (int),
    'econv': (float),
    'gconv': (float)
}

# Model keywords
MOD_KIN_DEFAULT = {
    'pressures': (),
    'rate_temps': (),
    'thermo_temps': (),
    'rate_fit': {
        'fit_method': 'arrhenius',
        'pdep_temps': (500, 100),
        'pdep_tol': 20.0,
        'pdep_pval': 1.0,
        'pdep_low': None,
        'pdep_high': None,
        'arr_dbl_tol': 15.0,
        'troe_param_fit_list': ('ts1', 'ts2', 'ts3', 'alpha'),
    },
    'thermo_fit': {
        'ref_scheme': 'basic',
        'ref_enes': 'ANL0'
    }
}

MOD_PF_DEFAULT = {
    'ene': (('sp', 'composite'), 'sp'),
    'rot': (('rigid', 'vpt2'), 'rigid'),
    'vib': (('harm', 'vpt2', 'tau'), 'harm'),
    'tors': (('rigid', '1dhr', '1dhrf', '1dhrfa', 'mdhr', 'mdhrv', 'tau'), 'rigid'),
    'sym': (('none', 'sampling', '1dhr'), 'none'),
    'ts_nobar': (('pst', 'rpvtst', 'vrctst'), 'pst'),
    'ts_sadpt': (('fixed', 'pst', 'rpvtst', 'vrctst'), 'fixed'),
    'wells': (('fake', 'find', 'none'), 'fake'),
    'rwells': (('fake', 'find', 'none'), 'fake'),
    'pwells': (('fake', 'find', 'none'), 'fake'),
    'tunnel': (('none', 'eckart', 'sct'), 'eckart'),
    'etrans': (('none', 'estimate', 'read'), 'estimate')
}


# MISC
VRC_DCT = {
    'fortran_compiler': 'gfortran',
    'spc_name': 'mol',
    'memory': 4.0,
    'r1dists_lr': [8., 6., 5., 4.5, 4.],
    'r1dists_sr': [4., 3.8, 3.6, 3.4, 3.2, 3., 2.8, 2.6, 2.4, 2.2],
    'r2dists_sr': [4., 3.8, 3.6, 3.4, 3.2, 3., 2.8, 2.6, 2.4, 2.2],
    'd1dists': [0.01, 0.5, 1.],
    'd2dists': [0.01, 0.5, 1.],
    'conditions': {},
    'nsamp_max': 2000,
    'nsamp_min': 50,
    'flux_err': 10,
    'pes_size': 2,
    'exe_path': '/blues/gpfs/home/sjklipp/bin/molpro'
}


# Dictionary Builders and Checkers
def build_default():
    """ Way of building the default dcts for various things, this is
        for ES tasks
    """
    supp_keywrds = ES_TSK_KEYWORDS_SUPPORTED_DCT[tsk]  
    default_dct = dict(
        zip(keywrds, (NEW_ES_TSK_DCT[key][2] for key in keywrds)))

    return default_dct


def check_dictionary(inp_dct, chk_dct, section, dyn_vals=()):
    """ Check if the dictionary to see if it has the allowed vals
    """

    # if inp_dct is not None:  # check if nonempty to see if section undefined

    # Assess if user-defined keywords
    # (1) include requird keywords and (2) only define supported keywords
    inp_keys = set(inp_dct.keys())
    chk_keys = set(chk_dct.keys())
    unsupported_keys = inp_keys - chk_keys
    undefined_required_keys = chk_keys - inp_keys

    if unsupported_keys:
        print('User defined unsupported keywords in {}'.format(section))
        for key in unsupported_keys:
            print(key)
    if undefined_required_keys:
        print('User has not defined required keywords in {}'.format(section))
        for key in undefined_required_keys:
            print(key)

    # Assess if the keywords have the appropriate value
    for key, val in inp_dct.items():
        allowed_typ, allowed_vals, _ = chk_dct[key]

        if not isinstance(type(val), allowed_typ):
            print('val must be type {}'.format(allowed_typ))
        if allowed_vals:
            if val not in allowed_vals:
                print('val is {}, must be {}'.format(val, allowed_vals))


def check_thy_lvls(key_dct, method_dct):
    """ For specific tasks, we need to a second level of cheeck to ensure
        that values of keywords that correspond to blocks defined in either
        the thy or model dat files.

        :param key_dct:
        :type key_dct: dict[]
        :param method_dct: thy or mod dct 
            
    """
    method_keys = ['runlvl', 'inplvl', 'var_splvl1', 'var_splvl2', 'var_scnlvl']
    assert set(key_dct[key] for key in method_keys) <= set(method_dct.keys())


def check_lst(inp_lst, sup_lst):
    """ Check
    """
    if set(inp_lst) >= sup_lst:
        print('Unsupported keys')
