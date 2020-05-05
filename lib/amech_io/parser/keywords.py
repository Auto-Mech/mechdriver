"""
  Libraries to check for allowed and supported keywords
"""

# Run Keywords
RUN_INP_REQUIRED_KEYWORDS = [
    'mech',
    'run_prefix',
    'save_prefix',
]
RUN_SUPPORTED_KEYWORDS = [
    'es',
    'thermochem',
    'kinetics',
    'write_messpf',
    'run_messpf',
    'run_nasa',
    'write_messrate',
    'run_messrate',
    'run_fits'
]

# ES Task keywords


# Model keywords
MODEL_PF_SUPPORTED_DCT = {
    'ene': ['sp', 'composite'],
    'rot': ['rigid', 'vpt2'],
    'vib': ['harm', 'vpt2', 'tau'],
    'tors': ['rigid', '1dhr', '1dhrf', 'mdhr', 'tau'],
    'sym': ['none', 'sampling', '1dhr'],
    'ts_barrierless': ['pst', 'rpvtst', 'vrctst'],
    'ts_sadpt': ['fixed', 'rpvtst'],
    'wells': ['fake', 'find'],
    'tunnel': ['none', 'eckart', 'sct']
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
    'tunnel': 'none'
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
        'conf_samp', 'conf_energy', 'conf_grad', 'conf_hess', 'conf_vpt2',
        'hr_scan', 'hr_energy', 'hr_grad', 'hr_hess',
        'tau_samp', 'tau_energy', 'tau_grad', 'tau_hess'],
    'ts': [
        'find_ts',
        'conf_samp', 'conf_energy', 'conf_grad', 'conf_hess', 'conf_vpt2',
        'hr_scan', 'hr_energy', 'hr_grad', 'hr_hess',
        'tau_samp', 'tau_energy', 'tau_grad', 'tau_hess',
        'irc_scan', 'irc_energy', 'irc_grad', 'irc_hess',
        'drp_samp', 'drp_energy', 'drp_grad', 'drp_hess'],
    'vdw': [
        'find',
        'conf_samp', 'conf_energy', 'conf_grad', 'conf_hess']
}
ES_TSK_KEYWORDS_SUPPORTED_DCT = {
    'init_geom': ['runlvl', 'inplvl', 'overwrite'],
    'find_ts': ['runlvl', 'inplvl', 'var_splvl1', 'var_splvl2', 'var_scnlvl',
                'nobarrier', 'overwrite'],
    'conf_samp': ['runlvl', 'inplvl', 'overwrite'],
    'conf_energy': ['runlvl', 'inplvl', 'zpve_min', 'overwrite'],
    'conf_grad': ['runlvl', 'inplvl', 'zpve_min', 'overwrite'],
    'conf_hess': ['runlvl', 'inplvl', 'zpve_min', 'overwrite'],
    'conf_vpt2': ['runlvl', 'inplvl', 'zpve_min', 'overwrite'],
    'hr_scan': ['runlvl', 'inplvl', 'frz_all_tors', 'ndim_tors', 'overwrite'],
    'hr_grad': ['runlvl', 'inplvl', 'frz_all_tors', 'ndim_tors', 'overwrite'],
    'hr_hess': ['runlvl', 'inplvl', 'frz_all_tors', 'ndim_tors', 'overwrite'],
    'hr_energy': ['runlvl', 'inplvl', 'frz_all_tors', 'ndim_tors',
                  'overwrite'],
    'tau_samp': ['runlvl', 'inplvl', 'overwrite'],
    'tau_energy': ['runlvl', 'inplvl', 'overwrite'],
    'tau_grad': ['runlvl', 'inplvl', 'overwrite'],
    'tau_hess': ['runlvl', 'inplvl', 'overwrite'],
    'irc_scan': ['runlvl', 'inplvl', 'overwrite'],
    'irc_energy': ['runlvl', 'inplvl', 'overwrite'],
    'irc_grad': ['runlvl', 'inplvl', 'overwrite'],
    'irc_hess': ['runlvl', 'inplvl', 'overwrite'],
    'drp_scan': ['runlvl', 'inplvl', 'overwrite'],
    'drp_energy': ['runlvl', 'inplvl', 'overwrite'],
    'drp_grad': ['runlvl', 'inplvl', 'overwrite'],
    'drp_hess': ['runlvl', 'inplvl', 'overwrite'],
}
ES_TSK_KEYWORDS_VAL_SUPPORTED_DCT = {
    'frz_all_tors': [True, False],
    'ndim_tors': ['1dhr', 'mdhr'],
    'zpve_min': [True, False],
    'nobarrier': ['pst', 'rpvtst', 'vrctst']
}
ES_TSK_KEYWORDS_DEFAULT_DCT = {
    'runlvl': None,
    'inplvl': None,
    'var_splvl1': None,
    'var_splvl2': None,
    'var_scnlvl': None,
    'frz_all_tors': False,
    'ndim_tors': '1dhr',
    'zpve_min': False,
    'nobarrier': 'pst',
    'overwrite': False
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
    'mc_tau',
    'hind_inc',
    'tors_names',
    'elec_levels',
    'sym_factor',
    'inchi',
    'smiles',
    'geom',
    'kickoff',
    'ts_search'
]
SPC_DEFAULT_DCT = {
    'kickoff': [0.1, False],
    'pst_params': [1.0, 6],
    'hind_inc': 30.0,
    'mc_nsamp': [True, 10, 1, 3, 50, 10]
}

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
