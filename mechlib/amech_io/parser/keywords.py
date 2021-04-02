"""
  Libraries to check for allowed and supported keywords
"""

from phydat import phycon


# Run Keywords
RUN_INP_REQUIRED_KEYWORDS = [
    'mech',
    'spc'
]
RUN_INP_SUPPORTED_KEYWORDS = [
    'mech',
    'spc',
    'run_prefix',
    'save_prefix',
    'print_mech',
    'print_debug'
]
RUN_INP_KEY_DCT = {
    'mech': ['chemkin'],
    'spc': ['csv'],
    'print_mech': [True, False],
    'print_debug': [True, False]
}
RUN_INP_DEFAULT_DCT = {
    'inp_mech': 'chemkin',
    'inp_spc': 'csv',
    'out_mech': 'chemkin',
    'out_spc': 'csv',
    'print_mech': False,
    'print_debug': False,
    'run_prefix': None,
    'save_prefix': None
}

RUN_SUPPORTED_KEYWORDS = [
    'es',
    'thermochem',
    'kinetics',
    'write_messpf',
    'run_messpf',
    'run_nasa',
    'write_messrate',
    'run_messrate',
    'run_fits',
    'transport',
    'print'
]

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

MODEL_PF_SUPPORTED_DCT = {
    'ene': ['sp', 'composite'],
    'rot': ['rigid', 'vpt2'],
    'vib': ['harm', 'vpt2', 'tau'],
    'tors': ['rigid', '1dhr', '1dhrf', '1dhrfa', 'mdhr', 'mdhrv', 'tau'],
    'sym': ['none', 'sampling', '1dhr'],
    'ts_barrierless': ['pst', 'rpvtst', 'vrctst'],
    'ts_sadpt': ['fixed', 'pst', 'rpvtst', 'vrctst'],
    'wells': ['fake', 'find', 'none'],
    'rwells': ['fake', 'find', 'none'],
    'pwells': ['fake', 'find', 'none'],
    'tunnel': ['none', 'eckart', 'sct'],
    'etrans': ['none', 'estimate', 'read']
}
MODEL_PF_DEFAULT_DCT = {
    'ene': 'sp',
    'rot': 'rigid',
    'vib': 'harm',
    'tors': 'rigid',
    'sym': 'none',
    'ts_nobar': 'pst',
    'ts_sadpt': 'fixed',
    'wells': 'fake',
    'tunnel': 'none',
    'etrans': 'estimate'
}

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



# Electronic Structure Tasks
ES_TSK_OBJ_SUPPORTED_LST = [
    'spc',
    'ts',
    'vdw'
]
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

# es tsk: (object type, (allowed values), default)  # use functions for weird
NEW_ES_TSK_DCT = {
    # all tasks
    'runlvl': (str, (), None),
    'inplvl': (str, (), None),
    'var_splvl1': (str, (), None),
    'var_splvl2': (str, (), None),
    'var_scnlvl': (str, (), None),
    'retryfail': (bool, (True, False), True),
    'overwrite': (bool, (True, False), False),
    # conformer tasks
    'cnf_range': (str, (), 'min'),
    # tau tasks
    'hessmax': (int, (), 1000),
    # hr tasks
    'tors_model': (str, ('1dhr', '1dhrf', '1dhrfa', 'mdhr'), '1dhr'),
    'resamp_min': (bool, (True, False), False),
    'hrthresh': (float, (), -0.5),
    # reaction
    # 'rxndirn': ['forw', 'back', 'exo'],
    'rxncoord': (str, ('irc', 'auto'), 'auto')
}


ES_TSK_KEYWORDS_VAL_SUPPORTED_DCT = {
    'tors_model': ['1dhr', '1dhrf', '1dhrfa', 'mdhr'],
    'cnf_range': ['min'],
    'nobarrier': ['pst', 'rpvtst', 'vrctst'],
    'retryfail': [True, False],
    'overwrite': [True, False],
    'resave': [True, False],
    'rxndirn': ['forw', 'back', 'exo'],
    'rxncoord': ['irc', 'auto'],
    'resamp_min': [True, False],
}
ES_TSK_KEYWORDS_DEFAULT_DCT = {
    'runlvl': None,
    'inplvl': None,
    'var_splvl1': None,
    'var_splvl2': None,
    'var_scnlvl': None,
    'tors_model': '1dhr',
    'cnf_range': 'min',
    'nobarrier': 'pst',
    'retryfail': True,
    'overwrite': False,
    'resave': False,
    'rxndirn': 'forw',
    'rxncoord': 'auto',
    'hessmax': 1000,
    'hrthresh': -0.5,
    'pot_thresh': 0.3
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

# ETrans
TRANS_TSK_SUPPORTED_DCT = {
    'spc': ['onedmin']
}
TRANS_TSK_KEYWORDS_SUPPORTED_DCT = {
    'onedmin': ['runlvl', 'inplvl', 'bath', 'pot',
                'njobs', 'nsamp',
                'smin', 'smax', 'conf',
                'retryfail', 'overwrite']
}
TRANS_TSK_KEYWORDS_VAL_SUPPORTED_DCT = {
    'retryfail': [True, False],
    'overwrite': [True, False],
}
TRANS_TSK_KEYWORDS_VAL_DEFAULT_DCT = {
    'pot': 'lj_12_6',
    'njobs': 1,
    'nsamp': 1,
    'smin': 2.0,
    'smax': 6.0,
    'conf': 'sphere',
    'retryfail': True,
    'overwrite': False,
}
# Print
PRNT_TSK_OBJ_SUPPORTED_LST = [
    'spc',
    'ts',
    'vdw'
]
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
PRNT_TSK_KEYWORDS_DEFAULT_DCT = {
    'geolvl': None,
    'proplvl': None,
    'nconfs': 'min',
    'econfs': 'min',
    'scale': None
    }
