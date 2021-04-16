""" Libraries of supported keywords and their corresponding values

    Also includes functionalities for constructing dictionaries
    of default values as well as assessing the validity of user input.
"""

import sys
import automol
from phydat import symm, eleclvl


# Run Keywords
RUN_INP_VAL_DCT = {
    'inp_mech': (str, ('chemkin'), 'chemkin'),
    'inp_spc': (str, ('csv',), 'csv'),
    'out_mech': (str, ('chemkin'), 'chemkin'),
    'out_spc': (str, ('csv',), 'csv'),
    'print_mech': (bool, (True, False), False),
    'print_debug': (bool, (True, False), False),
    'run_prefix': (str, None, None),
    'save_prefix': (str, None, None)
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
    # KTP/Therm
    'write_mess': ((), ('kin_model', 'spc_model', 'overwrite')),
    'run_mess': ((), ('nprocs', 'inpname')),
    'run_fits': ((), ('kin_model',)),
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
    # Trans
    'pot': (str, ('sphere',), 'lj_12_6'),
    'njobs': (int, (), 1),
    'nsamp': (int, (), 1),
    'smin': (float, (), 2.0),
    'smax': (float, (), 6.0),
    'conf': (str, ('sphere',), 'sphere'),
    # Proc
    'geolvl': (str, (), None),
    'proplvl': (str, (), None),
    'nconfs': (str, (), 'min'),
    'econfs': (str, (), 'min'),
    'scale': (str, (), None),
    # KTP/Therm
    'kin_model': (str, (), None),
    'spc_model': (str, (), None),
    'nprocs': (int, (), 10),
    'inpname': (str, (), None)
}
# Have nconfs and econfs keywords and combine them to figure out which to use?

SPC_VAL_DCT = {
    'mult': (int, (), None),
    'charge': (int, (), None),
    'inchi': (str, (), None),
    'smiles': (str, (), None),
    'tors_names': (str, (), None),
    'elec_levels': (str, (), None),
    'sym_factor': (str, (), None),
    'kickoff': (tuple, (), (True, (0.1, False))),
    'hind_inc': (float, (), 30.0),
    'mc_nsamp': (tuple, (), (True, 12, 1, 3, 100, 25)),
    'tau_nsamp': (tuple, (), (True, 12, 1, 3, 100, 25)),
    'smin': (float, (), None),
    'smax': (float, (), None),
    'etrans_nsamp': (int, (), None),
    'lj': (tuple, (), None),
    'edown': (tuple, list, (), None),
    'active': (tuple, (), None),
    'zma_idx': (int, (), 0)
}
TS_VAL_DCT = {
    'pst_params': (tuple, (), (1.0, 6)),
    'rxndirn': (str, (), 'forw'),
    'kt_pst': (float, (), 4.0e-10),
    'temp_pst': (float, (), 300.0),
    'n_pst': (float, (), 6.0),
    'active': (str, (), None),
    'ts_seatch': (str, (), 'sadpt'),
    'ts_idx': (int, (), 0)
}
TS_VAL_DCT.update(SPC_VAL_DCT)

# Theory Keywords
# rquired, Type, allowed, default,
# maybe set defaults using the qchem params script?
THY_VAL_DCT = {
    'program': (str, (), None),
    'method': (str, (), None),
    'basis': (str, (), None),
    'orb_res': (str, ('RR', 'UU', 'RU'), None),
    'ncycles': (int, (), None),
    'mem': (float, (), None),
    'nprocs': (int, (), None),
    'econv': (float, (), None),
    'gconv': (float, (), None)
}

# Model keywords
MODKIN_VAL_DEFAULT = {
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

MODPF_VAL_DCT = {
    'ene': (str, ('sp', 'composite'), 'sp'),
    'rot': (str, ('rigid', 'vpt2'), 'rigid'),
    'vib': (str, ('harm', 'vpt2', 'tau'), 'harm'),
    'tors': (str, ('rigid', '1dhr', '1dhrf', '1dhrfa', 'mdhr', 'mdhrv', 'tau'),
             'rigid'),
    'sym': (str, ('none', 'sampling', '1dhr'), 'none'),
    'ts_nobar': (str, ('pst', 'rpvtst', 'vrctst'), 'pst'),
    'ts_sadpt': (str, ('fixed', 'pst', 'rpvtst', 'vrctst'), 'fixed'),
    # 'wells': (str, ('fake', 'find', 'none'), 'fake'),
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


# Functions needed to build custom values
def elvl_symf(dct, ich, mul):
    """ set elec levels and sym factor
    """

    if 'elec_levels' not in dct:
        dct['elec_levels'] = eleclvl.DCT.get(
            (ich, mul), (0.0, mul))
    if 'sym_factor' not in dct:
        dct['sym_factor'] = symm.DCT.get(
            (ich, mul), 1.0)


# def active(dct):
#     """ maybe just read them (like the geometry dct)
#     """
#     # Add speciaized calls not in the default dct
#     if 'active' not in dct:
#         dct['active_space'] = None
#     else:
#         aspace = dct.get('active')
#         assert len(aspace) == 4, (
#             'active must be length 4: {}'.format(aspace)
#         )
#         wfn_file = aspace[3]
#         wfn_inp = os.path.join(os.path.join(job_path, 'inp/'+wfn_file))
#         if os.path.exists(wfn_inp):
#             wfn_str = ioformat.ptt.read_inp_str(job_path, wfn_inp)
#             print('Found file: {}. Reading file...'.format(wfn_file))
#         else:
#             wfn_str = None
#             print('No file: {}. Reading file...'.format(wfn_file))
#         dct['active_space'] = (
#             aspace[0], aspace[1], aspace[2], wfn_str)
#     return None


# Dictionary Builders
def defaults_from_val_dct(dct):
    """ Building the default dictionary for run inp dictionary

        prob works for spc as well
    """
    supp_keywrds = tuple(dct.keys())
    default_dct = dict(
        zip(supp_keywrds, (dct[key][2] for key in supp_keywrds)))

    return default_dct


def defaults_from_key_val_dcts(key, key_dct, val_dct):
    """ Way of building the default dcts for various things, this is
        for tasks blocks

        only works for tsks since it's info is in two dcts
    """

    # Set all of the keywords that are allowed for a task
    keywrds = key_dct[key][1]

    # Now build a dct where all the keywords are defaulted to internal value
    default_dct = dict(zip(keywrds, (val_dct[kwrd][2] for kwrd in keywrds)))

    return default_dct


# Dictionary Checkers
def check_val_dictionary1(inp_dct, val_dct, section):
    """ Check if the dictionary to see if it has the allowed vals
    """

    # if inp_dct is not None:  # check if nonempty to see if section undefined

    # Assess if user-defined keywords
    # (1) include requird keywords and (2) only define supported keywords
    inp_keys = set(inp_dct.keys())
    chk_keys = set(val_dct.keys())
    unsupported_keys = inp_keys - chk_keys
    # undefined_required_keys = chk_keys - inp_keys

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


def check_val_dictionary2(inp_dct, val_dct, section):
    """ check dct function 2
    """

    # Assess if the keywords have the appropriate value
    print('inpdct\n', inp_dct)
    print('valdct\n', val_dct)
    for key, val in inp_dct.items():

        allowed_typ, allowed_vals, _ = val_dct[key]

        # fails if None hit, need some way of aviding this
        # maybe the required checks use if None given?
        if not isinstance(val, allowed_typ):
            print('bad {}'.format(section))
            print('val {} must be type {}'.format(val, allowed_typ))
            sys.exit()
        if allowed_vals:
            if val not in allowed_vals:
                print('bad {}'.format(section))
                print('val is {}, must be {}'.format(val, allowed_vals))
                sys.exit()


# special checkers
def _new_check_dct(tsk_lsts, tsk_key_dct, tsk_val_dct, thy_dct):
    """ Loop over all of the tasks, add default keywords and parameters
        and assesses if all the input is valid
    """

    for tsk_lst in tsk_lsts:

        # Unpack the task
        [obj, tsk, keyword_dct] = tsk_lst

        # Build the dictionary of default values for task
        default_dct = defaults_from_val_dct(tsk, tsk_key_dct, tsk_val_dct)
        if obj not in tsk_key_dct[tsk][0]:  # correct
            print('tsk {}, not allowed for {}'.format(tsk, obj))
            print('')
            sys.exit()

        # Update the current task dct with the default
        new_key_dct = automol.util.dict_.right_update(default_dct, keyword_dct)

        # Check if the keyword values are allowed
        # need 2nd for anything that takes a string from the thy.dat file
        if check_val_dictionary1(new_key_dct, tsk_val_dct, 'ES_TSKS'):
            print('\n\nCHECK FAILED, QUITTING...')
            sys.exit()
        if check_val_dictionary2(new_key_dct, tsk_val_dct, 'ES_TSKS'):
            print('\n\nCHECK FAILED, QUITTING...')
            sys.exit()
        if check_thy_lvls(new_key_dct, thy_dct):
            print('\n\nCHECK FAILED, QUITTING...')
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
                print('User has not defined val in {}'.format(section))
                print(key, val)
                sys.exit()


def check_model_combinations(pf_dct):
    """ Check if a model combination is not implemented for PF routines
    """
    if pf_dct['vib'] == 'vpt2' and pf_dct['tors'] == '1dhr':
        print('*ERROR: VPT2 and 1DHR combination is not yet implemented')
        sys.exit()
    elif pf_dct['vib'] == 'vpt2' and pf_dct['tors'] == 'tau':
        print('*ERROR: VPT2 and TAU combination is not yet implemented')
        sys.exit()
