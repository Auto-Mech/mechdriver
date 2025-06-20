""" Parses the `run.dat` input file for MechDriver that specifices all
    of calculations to run for a given session of the code.

    Specifcally, looks for and parses several subsections:
        (1) `input` block: various input
        (2) `pes' block: idxs denoting what PESs in mech file to run
        (3) `spc` block: idxs denoting what species in .csv file to run
        (4) `els tasks` block: set of tasks for ESDriver
        (5) `therm tasks` block: set of tasks for ThermDriver
        (6) `ktp tasks` block: set of tasks for kTPDriver
        (7) `trans tasks` block: set of tasks for TransDriver
        (8) `proc tasks` block: set of tasks for ProcDriver

    Function parses the strings and converts them into formatted dictionaries
    that are passed to the sub-drivers of the code:
        ESDriver, ThermoDriver, kTPDriver, TransDriver, ProcDriver

    These dictionaries are built in three stages:
        (1) filled with user-specified options
        (2) default values not defined by the user are added, and
        (3) assessed that all keywordws and values are supported by the code.
"""

import sys
import ioformat
from mechlib.amech_io.printer import error_message
from mechlib.amech_io.parser._keywrd import defaults_from_val_dct
from mechlib.amech_io.parser._keywrd import defaults_from_key_val_dcts
from mechlib.amech_io.parser._keywrd import check_dct1
from mechlib.amech_io.parser._keywrd import check_thy_lvls
from mechlib.amech_io.parser._keywrd import right_update


# DICTIONARIES OF DEFAULTS #
# Run Keywords
RUN_INP_REQ = [
    'inp_mech', 'out_mech', 'inp_spc', 'out_spc', 'run_prefix', 'save_prefix']
RUN_INP_VAL_DCT = {
    'inp_mech': ((str,), ('chemkin'), 'chemkin'),
    'inp_spc': ((str,), ('csv',), 'csv'),
    'out_mech': ((str,), ('chemkin'), 'chemkin'),
    'out_spc': ((str,), ('csv',), 'csv'),
    'print_mech': ((bool,), (True, False), False),
    'print_debug': ((bool,), (True, False), False),
    'run_prefix': ((str,), (), None),
    'save_prefix': ((str,), (), None),
    'canonical': ((bool,), (True, False), False)
}

# HANDLE TASK KEYS

# Commonly useful task keyword lists
BASE = ('runlvl', 'inplvl', 'retryfail', 'overwrite')
MREF = ('var_splvl1', 'var_splvl2', 'var_scnlvl')
TRANS = ('njobs', 'nsamp', 'conf', 'cnf_range', 'sort')
PRNT = ('geolvl', 'proplvl', 'cnf_range', 'sort', 'nprocs', 'spc_model', 'kin_model')

# Supported object types for task (useful if task requestes 'all')
SUPP_OBJS = ('spc', 'ts')

# Determines what objects and keywords are allowed for tasks for ES,Trans,Print
# Need way to set required tsks
# Tasks: (allowed obj, allowed_keywords)
TSK_KEY_DCT = {
    # Electronic Structure Driver Tasks
    'init_geom': (('spc',), BASE),
    'find_ts': (('spc', 'ts'), BASE + MREF + ('nobarrier', 'varecof_nprocs')),
    'conf_pucker': (('spc', 'ts'), BASE + ('cnf_range', 'sort', 
                                           'algorithm','thresholds','eps','checks','rand_tors')),
    'conf_samp': (('spc', 'ts'), BASE + ('cnf_range', 'sort', 'resave',)),
    'conf_energy': (('spc', 'ts'), BASE + ('cnf_range', 'sort',)),
    'conf_grad': (('spc', 'ts'), BASE + ('cnf_range', 'sort',)),
    'conf_hess': (('spc', 'ts'), BASE + ('cnf_range', 'sort',)),
    'conf_vpt2': (('spc', 'ts'), BASE + ('cnf_range', 'sort',)),
    'conf_prop': (('spc', 'ts'), BASE + ('cnf_range', 'sort',)),
    'conf_opt': (('spc', 'ts'), BASE + ('cnf_range', 'sort', 'resave',)),
    'hr_scan': (('spc', 'ts'), BASE + ('tors_model', 'resamp_min',
                                       'cnf_range', 'sort',)),
    'hr_grad': (('spc', 'ts'), BASE + ('tors_model', 'cnf_range', 'sort',)),
    'hr_hess': (('spc', 'ts'), BASE + ('tors_model', 'cnf_range', 'sort',)),
    'hr_energy': (('spc', 'ts'), BASE + ('tors_model', 'cnf_range', 'sort',)),
    'hr_vpt2': (('spc', 'ts'), BASE + ('tors_model', 'cnf_range', 'sort',)),
    'hr_reopt': (('spc', 'ts'), BASE + ('tors_model', 'hrthresh',
                                        'cnf_range', 'sort',)),
    'tau_samp': (('spc', 'ts'), BASE + ('resave',)),
    'tau_energy': (('spc', 'ts'), BASE),
    'tau_grad': (('spc', 'ts'), BASE),
    'tau_hess': (('spc', 'ts'), BASE + ('hessmax',)),
    'rpath_scan': (('ts',), BASE + ('rxncoord',)),
    'rpath_energy': (('ts',), BASE + ('rxncoord',)),
    'rpath_grad': (('ts',), BASE + ('rxncoord',)),
    'rpath_hess': (('ts',), BASE + ('rxncoord',)),
    # Transport Driver Tasks
    'onedmin': (('spc',), (BASE + TRANS)),
    'write_transport': (('spc',), (BASE + TRANS)),
    # Process Driver Tasks
    'date': (('spc', 'ts', 'vdw'), PRNT),
    'freqs': (('spc', 'ts', 'vdw'), PRNT + ('scale',)),
    'energy': (('spc', 'ts'), PRNT),
    'weight': (('spc', 'ts'), PRNT),
    'gibbs': (('spc', 'ts'), PRNT),
    'geo': (('spc', 'ts'), PRNT),
    'molden': (('spc', 'ts'), PRNT),
    'sidata': (('spc', 'ts'), PRNT),
    'hess_json': (('spc', 'ts'), PRNT),
    'zmatrix': (('spc', 'ts'), PRNT),
    'torsions': (('spc', 'ts'), PRNT),
    'enthalpy': (('spc', 'ts'), PRNT),
    'pf': (('spc', 'ts'), PRNT),
    'messpf_inp': (('spc', 'ts'), PRNT),
    'coeffs': (('spc', 'ts'), ()),
    # KTP/Therm
    'write_mess': ((), ('kin_model', 'spc_model', 'overwrite',
                        'well_extension', 'well_lumping', 'mess_version',
                        'float_precision',
                        'cnf_range', 'sort', 'nprocs')),
    'run_mess': ((), ('kin_model', 'spc_model', 'nprocs',
                      'well_lumping', 'mess_version',
                      'cnf_range', 'sort')),
    'run_fits': ((), ('kin_model', 'spc_model',
                      'well_lumping', 'mess_version',
                      'combine',
                      'cnf_range', 'sort', 'nprocs',)),
}

# tsk: (object types, (allowed values), default)  # use functions for weird
# maybe the required checks use if None given?
TSK_VAL_DCT = {
    # Common
    'runlvl': ((str,), (), None),
    'inplvl': ((str,), (), None),
    'var_splvl1': ((str,), (), None),
    'var_splvl2': ((str,), (), None),
    'var_scnlvl': ((str,), (), None),
    'resave': ((bool,), (True, False), False),
    'retryfail': ((bool,), (True, False), True),
    'overwrite': ((bool,), (True, False), False),
    # ES
    'cnf_range': ((str,), (), 'min'),   # change to econfs, nconfs
    'sort': ((str,), (), None),
    'hessmax': ((int,), (), 1000),
    'tors_model': ((str,),
                   ('1dhr', '1dhrf', '1dhrfa', 'mdhr', 'mdhrv'), '1dhr'),
    'resamp_min': ((bool,), (True, False), False),
    'hrthresh': ((float,), (), -0.2),
    'potthresh': ((float,), (), 0.3),
    'rxncoord': ((str,), ('irc', 'auto'), 'auto'),
    'nobarrier': ((str,), ('pst', 'rpvtst', 'vrctst'), None),
    're_id': ((bool,), (True, False), False),
    'varecof_nprocs': ((int,), (), 10),
    #adl added arguments for ring puckering
    'algorithm': ((str,), 
                  ('crest','pucker','torsions','robust','torsions2','etkdg'), 'pucker'),
    'thresholds': ((str,), 
                  ('default','relaxed'), 'default'),
    'eps': ((float,), (), 0.2),
    'checks': ((int,), (1,2), 1),
    'rand_tors': ((int,), (), 2),
    # Trans
    'njobs': ((int,), (), 1),
    'nsamp': ((int,), (), 1),
    'conf': ((str,), ('sphere', 'min'), 'sphere'),
    # Proc
    'geolvl': ((str,), (), None),
    'proplvl': ((str,), (), None),
    'nconfs': ((str,), (), 'min'),
    'econfs': ((str,), (), 'min'),
    'scale': ((str,), (), None),
    # KTP/Therm
    'kin_model': ((str,), (), None),
    'spc_model': ((str,), (), None),
    'nprocs': ((int,), (), 9),
    'mess_version': ((str,), ('v1', 'v2'), 'v1'),
    'well_extension': ((bool,), (), False),
    'well_lumping': ((bool,), (), False),
    'combine': ((str,), ('stereo',), None),
    'linked_pes': ((tuple,), (), None),
    'float_precision': ((str,), ('double', 'quadruple'), 'double'),
}
# Have nconfs and econfs keywords and combine them to figure out which to use?


# INPUT PARSERS #
# Input Section
def input_dictionary(run_str):
    """ Parses the `input` block and builds a
        dictionary of keywords and their corresponding values.

        :param run_str: input string of the run.dat block
        :type run_str: str
        :rtype: dict[str: obj]
    """

    # Read the input block
    inp_block = ioformat.ptt.end_block(run_str, 'input', footer='input')
    inp_dct = ioformat.ptt.keyword_dct_from_block(inp_block)

    # Add defaults to the dictionary
    inp_dct = right_update(
        defaults_from_val_dct(RUN_INP_VAL_DCT), inp_dct)

    # Check the dictionary
    check_dct1(inp_dct, RUN_INP_VAL_DCT, RUN_INP_REQ, 'Run-Input')

    return inp_dct


# Chemistry objects
def chem_idxs(run_str):
    """  Parses the `pes` block of the run.dat file and
         builds a dictionary of the PESs and corresponding channels the
         user wishes to run.

         Parses the `spc` block of the run.dat file and
         builds a dictionary of the species the
         user wishes to run.

         May break if idx is given on two lines of string.

        :param run_str: string of the run.dat input file
        :type run_str: str
        :returns: ({pes_idx: list of channel_idxs}, {1: list of species idxs})
        :rtype: dict[str: tuple]
    """

    # PES idxs to run
    pes_block = ioformat.ptt.end_block(run_str, 'pes', footer='pes')

    if pes_block is not None:
        _pes_idxs = {}
        for line in pes_block.strip().splitlines():
            [pes_nums, chn_nums] = line.split(':')
            _pes_nums = ioformat.ptt.idx_lst_from_line(pes_nums)
            _chn_nums = ioformat.ptt.idx_lst_from_line(chn_nums)
            for idx in _pes_nums:
                _pes_idxs.update({idx-1: tuple(val-1 for val in _chn_nums)})
    else:
        _pes_idxs = None

    # SPC idxs to run
    spc_block = ioformat.ptt.end_block(run_str, 'spc', footer='spc')

    if spc_block is not None:
        _idxs = ()
        for line in spc_block.splitlines():
            _idxs += ioformat.ptt.idx_lst_from_line(line)
        _spc_idxs = {1: tuple(val-1 for val in _idxs)}
    else:
        _spc_idxs = None

    # Kill code if no idxs given
    if _pes_idxs is None and _spc_idxs is None:
        error_message('No pes or spc section given in run.dat file. Quitting')
        sys.exit()

    return _pes_idxs, _spc_idxs


# Driver Task Lists
def extract_task(tsk, tsk_lst):
    """ Searches for a task in the task lst and if found:
        the corresponding keywords and values will be returned

        Function only works if task is present in the list one time.

        :param tsk: task to extract information for
        :type tsk: str
        :param tsk_lst: list of tasks to run for some driver
        :type tsk_lst: tuple(tuple(str/dict))
        :rtype: tuple(str/dict)
    """

    tsk_inf = None
    for _tsk_inf in tsk_lst:
        if any(x == tsk for x in _tsk_inf):  # just looks in all pars
            tsk_inf = _tsk_inf
            break

    return tsk_inf


def tasks(run_str, thy_dct):
    """ runstr
    """

    # Read blocks and build user determined task lists`
    es_block = ioformat.ptt.end_block(run_str, 'els', footer='els')
    trans_block = ioformat.ptt.end_block(run_str, 'trans', footer='trans')
    therm_block = ioformat.ptt.end_block(run_str, 'thermo', footer='thermo')
    ktp_block = ioformat.ptt.end_block(run_str, 'ktp', footer='ktp')
    proc_block = ioformat.ptt.end_block(run_str, 'proc', footer='proc')

    # print('els\n', es_block)
    # print('therm\n', therm_block)
    # print('trans\n', trans_block)
    # print('proc\n', proc_block)

    es_tsks = _tsk_lst(es_block, 3)
    therm_tsks = _tsk_lst(therm_block, 2)
    ktp_tsks = _tsk_lst(ktp_block, 2)
    trans_tsks = _tsk_lst(trans_block, 3)
    proc_tsks = _tsk_lst(proc_block, 3)

    # Add defaults to each task as needed
    es_tsks = _tsk_defaults(es_tsks)
    therm_tsks = _tsk_defaults(therm_tsks)
    ktp_tsks = _tsk_defaults(ktp_tsks)
    trans_tsks = _tsk_defaults(trans_tsks)
    proc_tsks = _tsk_defaults(proc_tsks)

    # Assess each dictionary for correctness
    _check_tsks(es_tsks, thy_dct)
    _check_tsks(therm_tsks, thy_dct)
    _check_tsks(ktp_tsks, thy_dct)
    _check_tsks(trans_tsks, thy_dct)
    _check_tsks(proc_tsks, thy_dct)

    tsk_dct = {
        'es': es_tsks,
        'thermo': therm_tsks,
        'ktp': ktp_tsks,
        'trans': trans_tsks,
        'proc': proc_tsks
    }

    return tsk_dct


def _tsk_lst(tsk_str, num):
    """ Set the sequence of electronic structure tasks for a given
        species or PESs
    """

    # Build the task lists from the string
    if tsk_str is not None:
        tsks = []
        tsk_str = ioformat.remove_whitespace_from_string(tsk_str)
        for line in tsk_str.splitlines():
            _tsk = _split_line(line, num)
            tsks.append(_tsk)
        mod_tsks = tsks
        # mod_tsks = _expand_tsks(tsks) if num == 3 else tsks
    else:
        mod_tsks = None

    return mod_tsks


def _expand_tsks(tsks_lst):
    """ Loops over the driver task list and checks if each task is a
        macro-task that should be expanded into sub-tasks.

            Right now, it splits all obj tasks into spc and ts

        :param tsk_lst: list of tasks to run for some driver
        :type tsk_lst: tuple(tuple(str/dict))
        :rtype: tuple(str/dict)
    """

    mod_tsks_lst = []
    for tsk_lst in tsks_lst:
        [obj, tsk, dct] = tsk_lst
        objs = ['spc', 'ts'] if obj == 'all' else [obj]
        for obj in objs:
            mod_tsks_lst.append([obj, tsk, dct])

    return mod_tsks_lst


def _tsk_defaults(tsk_lst):
    """ Fill out the keyword dictionaries for various task lists with
        default values
    """

    if tsk_lst is not None:
        mod_tsk_lst = []
        for _tsk_lst in tsk_lst:
            keyword_dct = _tsk_lst[-1]
            tsk = _tsk_lst[:-1][-1]
            default_dct = defaults_from_key_val_dcts(
                tsk, TSK_KEY_DCT, TSK_VAL_DCT)
            new_key_dct = right_update(
                default_dct, keyword_dct)

            mod_lst = _tsk_lst[:-1] + [new_key_dct]
            mod_tsk_lst.append(mod_lst)
    else:
        mod_tsk_lst = None

    return mod_tsk_lst


def _check_tsks(tsk_lsts, thy_dct):
    """ Loop over all of the tasks, add default keywords and parameters
        and assesses if all the input is valid
    """

    if tsk_lsts is not None:

        for tsk_lst in tsk_lsts:

            # Unpack the task
            _tsk = tsk_lst[:-1]
            if len(_tsk) == 2:
                # Case(1): spc task keywords (ESDriver)
                obj, tsk = _tsk[0], _tsk[1]
                # Have to make lst to handle case where obj == 'all'
                obj_lst = SUPP_OBJS if obj == 'all' else (obj,)
                allowed_keys = TSK_KEY_DCT[tsk][0]
            else:
                # Case(2): task keywords (ThermoDriver, kTPDriver)
                tsk = _tsk[0]
                obj_lst = list(tsk_lst[1].keys())
                allowed_keys = TSK_KEY_DCT[tsk][1]
                
            for _obj in obj_lst:
                if _obj not in allowed_keys:
                    error_message(f'obj {_obj}, not allowed for {tsk}')
                    print(f'Allowed keys: {" ".join(allowed_keys)}')
                    sys.exit()

            key_dct = tsk_lst[-1]

            # Check if keyword values are allowed
            check_dct1(key_dct, TSK_VAL_DCT, (), 'Task')

            # Check keywords with thylvls as values use lvls defined in thy dct
            check_thy_lvls(key_dct, thy_dct)


def _split_line(line, num):
    """ Split a line
    """
    line = line.split()
    if num == 3:
        tsk, key_lst = line[:2], line[2:]
    elif num == 2:
        tsk, key_lst = line[:1], line[1:]
    key_dct = ioformat.ptt.keyword_dct_from_block('\n'.join(key_lst))

    return tsk + [key_dct]  # could convert to empty dct instead of None


# Check a bunch of stuff
def check_inputs(tsk_dct, pes_dct, pes_mod_dct, spc_mod_dct):
    """ Check if inputs placed that is required
    """

    # Check if a mechanism has been provided where required
    if tsk_dct['ktp'] or tsk_dct['thermo']:
        if pes_mod_dct is None:
            error_message(
                'kTPDriver or Thermo Requested. \n'
                ' However no kin model provided in models.dat\n'
                ' Exiting MechDriver...')
            sys.exit()
        if spc_mod_dct is None:
            error_message(
                'kTPDriver or Thermo Requested. \n'
                '  However no spc model provided in models.dat\n'
                '  Exiting MechDriver...')
            sys.exit()

    if tsk_dct['ktp']:
        if pes_dct is None:
            error_message(
                'kTPDriver Requested. \n'
                '  However no reaction channels provided in mechanism.dat\n'
                '  Exiting MechDriver...')
            sys.exit()
