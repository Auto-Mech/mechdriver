"""
  Libraries to check for allowed and supported keywords
"""

import sys
from phydat import phycon


# Run Keywords
RUN_INP_DCT = {
    'inp_mech': (('chemkin'), 'chemkin'),
    'inp_spc': (('csv',), 'csv'),
    'out_mech': (('chemkin'), 'chemkin'),
    'out_spc': (('csv',), 'csv'),
    'print_mech': ((True, False), False),
    'print_debug': ((True, False), False),
    'run_prefix': (None, None),
    'save_prefix': (None, None)
}

# HANDLE TASK KEYS

# Commonly useful task keyword lists
BASE = ('runlvl', 'inplvl', 'retryfail', 'overwrite')
MREF = ('var_splvl1', 'var_splvl2', 'var_scnlvl')
TRANS = ('bath', 'pot', 'njobs', 'nsamp', 'smin', 'smax', 'conf')
PRNT = ('geolvl', 'proplvl', 'nconfs', 'econfs')

# Determines what objects and keywords are allowed for tasks for ES,Trans,Print
# Need way to set required tsks
# Tasks: (allowed obj, allowed_keywords)
TSK_KEY_DCT = {
    # Electronic Structure Driver Tasks
    'init_geom': (('spc',), BASE),
    'find_ts': (('spc', 'ts'), BASE + MREF + ('nobarrier',)),
    'conf_pucker': (('spc', 'ts'), BASE + ('cnf_range',)),
    'conf_samp': (('spc', 'ts'), BASE + ('cnf_range', 'resave',)),
    'conf_energy': (('spc', 'ts'), BASE + ('cnf_range',)),
    'conf_grad': (('spc', 'ts'), BASE + ('cnf_range',)),
    'conf_hess': (('spc', 'ts'), BASE + ('cnf_range',)),
    'conf_vpt2': (('spc', 'ts'), BASE + ('cnf_range',)),
    'conf_prop': (('spc', 'ts'), BASE + ('cnf_range',)),
    'conf_opt': (('spc', 'ts'), BASE + ('cnf_range',)),
    'hr_scan': (('spc', 'ts'), BASE + ('tors_model', 'resamp_min',)),
    'hr_grad': (('spc', 'ts'), BASE + ('tors_model',)),
    'hr_hess': (('spc', 'ts'), BASE + ('tors_model',)),
    'hr_energy': (('spc', 'ts'), BASE + ('tors_model',)),
    'hr_vpt2': (('spc', 'ts'), BASE + ('tors_model',)),
    'hr_reopt': (('spc', 'ts'), BASE + ('tors_model',)),
    'tau_samp': (('spc', 'ts'), BASE),
    'tau_energy': (('spc', 'ts'), BASE),
    'tau_grad': (('spc', 'ts'), BASE),
    'tau_hess': (('spc', 'ts'), BASE + ('hessmax',)),
    'rpath_scan': (('ts',), BASE + ('rxncoord',)),
    'rpath_energy': (('ts',), BASE + ('rxncoord',)),
    'rpath_grad': (('ts',), BASE + ('rxncoord',)),
    'rpath_hess': (('ts',), BASE + ('rxncoord',)),
    # Transport Driver Tasks
    'onedmin': (('spc',), (BASE + TRANS)),
    # Process Driver Tasks
    'freqs': (('spc', 'ts', 'vdw'), PRNT + ('scale',)),
    'energy': (('spc',), PRNT),
    'geo': (('spc',), PRNT),
    'zmatrix': (('spc',), PRNT),
    'enthalpy': (('spc',), PRNT),
    'coeffs': (('spc',), ()),
}

# es tsk: (object type, (allowed values), default)  # use functions for weird
# maybe the required checks use if None given?
TSK_VAL_DCT = {
    # Common
    'runlvl': (str, (), None),
    'inplvl': (str, (), None),
    'var_splvl1': (str, (), None),
    'var_splvl2': (str, (), None),
    'var_scnlvl': (str, (), None),
    'resave': (bool, (True, False), True),
    'retryfail': (bool, (True, False), True),
    'overwrite': (bool, (True, False), False),
    # ES
    'cnf_range': (str, (), 'min'),   # change to econfs, nconfs
    'hessmax': (int, (), 1000),
    'tors_model': (str, ('1dhr', '1dhrf', '1dhrfa', 'mdhr', 'mdhrv'), '1dhr'),
    'resamp_min': (bool, (True, False), False),
    'hrthresh': (float, (), -0.5),
    'potthresh': (float, (), 0.3),
    'rxncoord': (str, ('irc', 'auto'), 'auto'),
    'nobarrier': (str, ('pst', 'rpvtst', 'vrctst'), None),
    # TRans
    'pot': (str, ('sphere',), 'lj_12_6'),
    'njobs': (int, (), 1),
    'nsamp': (int, (), 1),
    'smin': (float, (), 2.0),
    'smax': (float, (), 6.0),
    'conf': (str, ('sphere',), 'sphere'),
    # PRoc
    'geolvl': (str, (), None),
    'proplvl': (str, (), None),
    'nconfs': (str, (), 'min'),
    'econfs': (str, (), 'min'),
    'scale': (str, (), None),
}
# Have nconfs and econfs keywords and combine them to figure out which to use?

OTHER_TSK_KEY_DCT = {
    'write_mess': ('kin_model', 'spc_model', 'overwrite'),
    'run_mess': ('nprocs', 'inpname'),
    'run_fits': ('kin_model',),
}
OTHER_TSK_VAL_DCT = {
    'kin_model': (str, (), None),
    'spc_model': (str, (), None),
    'overwrite': (str, (), False),
    'nprocs': (int, (), 10),
    'inpname': (str, (), None)
}

# Species keywords
SPC_REQUIRED_KEYWORDS = [
    'mult',
    'charge',
]
SPC_SUPPORTED_KEYWORDS = [
    'hind_inc',
    'geom',
    'ts_search',
    'active',
    'zma_idx',
]
SPC_ALT_DCT = {
    'smin': (None, float, (), None),
    'smax': (None, float, (), None),
    'etrans_nsamp': (None, int, (), None),
    'lj': (None, list, (), None),
    'edown': (None, list, (), None)
}

SPC_DEFAULT_DCT = {
    # requied
    'mult': (True, int, (), None),
    'charge': (True, int, (), None),
    'inchi': (),
    'smiles': (),
    # other
    'tors_names': (),
    'elec_levels': (),
    'sym_factor': (),

    # ones I know need to be in
    'kickoff': (True, (0.1, False)),
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
    'n_pst': 6.0,
    'active': (False, str, (), None),
}}


# Theory Keywords
# rquired, Type, allowed, default,
# maybe set defaults using the qchem params script?
THY = {
    'program': (True, str, (), None),
    'method': (True, str, (), None),
    'basis': (True, str, (), None),
    'orb_res': (True, str, ('RR', 'UU', 'RU'), None),
    'ncycles': (False, int, (), None),
    'mem': (False, float, (), None),
    'nprocs': (False, int, (), None),
    'econv': (False, float, (), None),
    'gconv': (False, float, (), None)
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
    },
    'glob_etransfer': {
        'lj': 1.0,
        'alpha': 1.0
    }
}

MOD_PF_DEFAULT = {
    'ene': (str, ('sp', 'composite'), 'sp'),
    'rot': (str, ('rigid', 'vpt2'), 'rigid'),
    'vib': (str, ('harm', 'vpt2', 'tau'), 'harm'),
    'tors': (str, ('rigid', '1dhr', '1dhrf', '1dhrfa', 'mdhr', 'mdhrv', 'tau'),
             'rigid'),
    'sym': (str, ('none', 'sampling', '1dhr'), 'none'),
    'ts_nobar': (str, ('pst', 'rpvtst', 'vrctst'), 'pst'),
    'ts_sadpt': (str, ('fixed', 'pst', 'rpvtst', 'vrctst'), 'fixed'),
    'wells': (str, ('fake', 'find', 'none'), 'fake'),
    'rwells': (str, ('fake', 'find', 'none'), 'fake'),
    'pwells': (str, ('fake', 'find', 'none'), 'fake'),
    'tunnel': (str, ('none', 'eckart', 'sct'), 'eckart'),
    'etrans': (str, ('none', 'estimate', 'read'), 'estimate')
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


# Dictionary Builders
def build_tsk_default(tsk, tsk_key_dct, tsk_val_dct):
    """ Way of building the default dcts for various things, this is
        for ES tasks
    """

    # Set all of the keywords that are allowed for a task
    supp_keywrds = tsk_key_dct[tsk][1]

    # Now build a dct where all the keywords are defaulted to internal value
    default_dct = dict(
        zip(supp_keywrds, (tsk_val_dct[key][2] for key in supp_keywrds)))

    return default_dct


def build_spc_default(sadpt=False):
    """ Use the spc dict and build the default values from it.
        Have some way fo getting the TSs as well?
    """

    defaults = {}
    if not sadpt:
        for key, (_, _, val) in SPC_DEFAULT_DCT.items():
            if val is not None:
                defaults[key] = val
    else:
        for key, (_, _, val) in TS_DEFAULT_DCT.items():
            if val is not None:
                defaults[key] = val

    return defaults


# Dictionary Checkers
def check_dictionary(inp_dct, key_dct, val_dct, section):
    """ Check if the dictionary to see if it has the allowed vals
    """

    # if inp_dct is not None:  # check if nonempty to see if section undefined

    # Assess if user-defined keywords
    # (1) include requird keywords and (2) only define supported keywords
    inp_keys = set(inp_dct.keys())
    chk_keys = set(val_dct.keys())
    unsupported_keys = inp_keys - chk_keys
    undefined_required_keys = chk_keys - inp_keys

    # print('inp\n', inp_keys)
    # print('chk\n', chk_keys)
    # print('unsupport\n', unsupported_keys)
    # print('unrequired\n', undefined_required_keys)

    if unsupported_keys:
        print('User defined unsupported keywords in {}'.format(section))
        for key in unsupported_keys:
            print(key)
        sys.exit()
    # not correct, need new way to do this
    # if undefined_required_keys:
    #     print('User has not defined required keywords in {}'.format(section))
    #     for key in undefined_required_keys:
    #         print(key)
    #     sys.exit()

    # Assess if the keywords have the appropriate value
    for key, val in inp_dct.items():

        allowed_typ, allowed_vals, _ = val_dct[key]

        # fails if None hit, need some way of aviding this
        # maybe the required checks use if None given?
        if not isinstance(val, allowed_typ):
            print('val {} must be type {}'.format(val, allowed_typ))
            sys.exit()
        if allowed_vals:
            if val not in allowed_vals:
                print('val is {}, must be {}'.format(val, allowed_vals))
                sys.exit()


def check_thy_lvls(key_dct, method_dct, section=''):
    """ For specific tasks, we need to a second level of cheeck to ensure
        that values of keywords that correspond to blocks defined in either
        the thy or model dat files.

        :param key_dct:
        :type key_dct: dict[]
        :param method_dct: thy or mod dct
    """

    thy_defined_methods = set(method_dct.keys())

    for key in ('runlvl', 'inplvl', 'var_splvl1', 'var_splvl2', 'var_scnlvl'):
        val = key_dct.get(key)
        if val is not None:
            if val not in thy_defined_methods:
                print('User has not defined val in thy.dat')
                print(key, val)
                sys.exit()


def check_lst(inp_lst, sup_lst):
    """ Check
    """
    if set(inp_lst) >= sup_lst:
        print('Unsupported keys')
