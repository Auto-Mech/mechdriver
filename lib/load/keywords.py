"""
  Libraries to check for allowed and supported keywords
"""

# Run Keywords
RUN_INP_REQUIRED_KEYWORDS = [
    'ids',
    'spc',
    'mech',
    'run_prefix',
    'save_prefix',
    'sort_rxns',
    'check_stereo',
    'rad_rad_sort'
]
RUN_SUPPORTED_KEYWORDS = [
    'sort',
    'es',
    'thermo',
    'rates',
    'fits',
    'poly'
]
OPTIONS_SUPPORTED_KEYWORDS = [
    'check_stereo',
    'overwrite',
    'ref_mols',
    'es_job_set',
    'hind_inc',
    'mc_nsamp'
]

# Model keywords
MODEL_SUPPORTED_DCT = {
    'ene': ['sp', 'composite'],
    'vib': ['harm'],
    'tors': ['rigid', '1dhr', 'mdhr', 'tau'],
    'sym': ['sampling', '1dhr'],
    'ts_barrierless': ['pst', 'vtst', 'vrctst'],
    'ts_sadpt': ['fixed', 'variational'],
    'wells': ['fake', 'find'],
    'tunneling': ['none', 'eckart', 'sct']
}
ES_TSK_SUPPORTED_LST = [
    'find_ts',
    'find_geom',
    'conf_samp',
    'conf_hess',
    'sym_samp',
    'hr_scan',
    'conf_vpt2',
    'conf_energy',
    'run_irc'
]

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
    'geom'
]

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
