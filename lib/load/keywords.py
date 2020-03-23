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
    'sort',
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
    'tors': ['rigid', '1dhr', 'mdhr', 'tau'],
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
    'tunneling': 'none'
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
ES_TSK_OPTIONS_SUPPORTED_DCT = {
    'init_geom': ['overwrite'],
    'find_ts': ['overwrite'],
    'conf_samp': ['overwrite'],
    'conf_energy': ['zpve_min', 'overwrite'],
    'conf_grad': ['zpve_min', 'overwrite'],
    'conf_hess': ['zpve_min', 'overwrite'],
    'conf_vpt2': ['zpve_min', 'overwrite'],
    'hr_scan': ['frz_all_tors', 'mdhr', 'overwrite'],
    'hr_grad': ['frz_all_tors', 'mdhr', 'overwrite'],
    'hr_hess': ['frz_all_tors', 'mdhr', 'overwrite'],
    'hr_energy': ['frz_all_tors', 'mdhr', 'overwrite'],
    'tau_samp': ['overwrite'],
    'tau_energy': ['overwrite'],
    'tau_grad': ['overwrite'],
    'tau_hess': ['overwrite'],
    'irc_scan': ['overwrite'],
    'irc_energy': ['overwrite'],
    'irc_grad': ['overwrite'],
    'irc_hess': ['overwrite'],
    'drp_scan': ['overwrite'],
    'drp_energy': ['overwrite'],
    'drp_grad': ['overwrite'],
    'drp_hess': ['overwrite'],
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
    'elec_levels',
    'sym_factor',
    'inchi',
    'smiles',
    'geom',
    'kickoff'
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
