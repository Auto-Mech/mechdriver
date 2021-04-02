""" Functions to handle reading the files of the various inputs

    Perhaps return the strings of all the sections that can be parsed out of the file
    # Read the blocks, perhaps return a dct of the strings?

    Check if the strings are none
"""

import sys
import automol
import autoparse.find as apf
import ioformat
# from mechlib.amech_io.parser.tsks import
from mechlib.amech_io.parser.keywords import THY_REQUIRED_KEYWORDS
from mechlib.amech_io.parser.keywords import THY_SUPPORTED_KEYWORDS
from mechlib.amech_io import printer as ioprinter


RUN_INP = 'inp/run.dat'
CSV_INP = 'inp/species.csv'
DAT_INP = 'inp/species.dat'
MECH_INP = 'inp/mechanism.dat'
SORT_INP = 'inp/sort.dat'
MODEL_INP = 'inp/models.dat'
THY_INP = 'inp/theory.dat'


def read_run(job_path):
    """ Parse the run.dat file
    """
   
    # Read the input string
    run_str = ioformat.pathtools.read_file(
        job_path, RUN_INP, remove_comments='#', remove_whitespace=True)
    ioprinter.reading('run.dat...', newline=1)  # Add a Found <file> to msg

    # Read the main blocks
    inp_block = ioformat.ptt.end_block(run_str, 'input', footer='input')

    pes_block = ioformat.ptt.end_block(run_str, 'pes', footer='pes')
    spc_block = ioformat.ptt.end_block(run_str, 'spc', footer='spc')

    es_tsks_block = ioformat.ptt.end_block(run_str, 'es', footer='es')
    trans_tsks_block = ioformat.ptt.end_block(run_str, 'transport', footer='transport')
    therm_tsks_block = ioformat.ptt.end_block(run_str, 'thermo', footer='thermo')
    ktp_tsks_block = ioformat.ptt.end_block(run_str, 'ktp', footer='ktp')
    proc_tsks_block = ioformat.ptt.end_block(run_str, 'proc', footer='proc')

    # Parse information in the blocks
    inp_key_dct = _keyword_dct(inp_block, {})
    print('inp_dct\n', inp_key_dct) 

    pes_idx_dct = _pes_idxs(pes_block)
    spc_idx_dct = _spc_idxs(spc_block)
    print('pes_dct\n', pes_idx_dct) 
    print('spc_dct\n', spc_idx_dct) 

    # es_tsk_lst = tsks.es_tsk_lst(es_tsk_str, thy_dct)
    # therm_tsk_lst = tsks.trans_tsk_lst(trans_tsk_str)
    # ktp_tsk_lst = tsks.trans_tsk_lst(trans_tsk_str)
    # trans_tsk_lst = tsks.trans_tsk_lst(trans_tsk_str)
    # proc_tsk_lst = tsks.proc_tsk_lst(proc_tsk_str)

    return None


def read_thy(job_path):
    """ Parse the theory.dat file
    """

    thy_str = ioformat.pathtools.read_file(
        job_path, THY_INP, remove_comments='#', remove_whitespace=True)
    ioprinter.reading('theory.dat...', newline=1)

    thy_blocks = ioformat.ptt.named_end_blocks(thy_str, 'level', footer='level')
    thy_dct = ioformat.ptt.keyword_dcts_from_lst(thy_blocks)
    print(thy_dct)

    # check_dictionary2(thy_dct, THY_REQUIRED_KEYWORDS THY_SUPPORTED_KEYWORDS)

    return thy_dct


def read_model(job_path):
    """ Parse the models.dat file
    """

    mod_str = ioformat.pathtools.read_file(
        job_path, MODEL_INP, remove_comments='#', remove_whitespace=True)
    ioprinter.reading('model.dat...', newline=1)

    kin_blocks = ioformat.ptt.named_end_blocks(mod_str, 'kin')
    spc_blocks = ioformat.ptt.named_end_blocks(mod_str, 'spc')

    kin_mod_dct = ioformat.ptt.keyword_dcts_from_lst(kin_blocks)
    spc_mod_dct = ioformat.ptt.keyword_dcts_from_lst(spc_blocks)

    print('kin_mod\n', kin_mod_dct)
    print('spc_mod\n', spc_mod_dct)
    # mix the global dicts
    # def combine_glob_and_spc_dct(glob_dct, spc_dct):


def read_spc(job_path):
    """ a
    """
    spc_str = ioformat.pathtools.read_file(job_path, CSV_INP)
    ioprinter.reading('species.csv...', newline=1)
    dat_str = ioformat.pathtools.read_file(job_path, DAT_INP)
    ioprinter.reading('species.dat...', newline=1)

    return spc_str


def read_mech(job_path):
    """Build the PES dct
    """

    # Read the string
    mech_str = ioformat.pathtools.read_file(
        job_path, MECH_INP, remove_comments='!', remove_whitespace=True)
    mech_info = util.read_mechanism_file(
        mech_str, mech_type, spc_dct, sort_rxns=sort_rxns)

    return mech_info


# Keyword dict build and check 
def _keyword_dct(inp_str, def_dct):
    """ merge a dictionary from the input
        :param inp_dct: dictionary of input keywords
        :param def_dct: dictionary of default values
                        (only has ones where defaults, won't have vals for 
                         keys that user must input like run_prefix)
    """
    inp_dct = ioformat.ptt.keyword_dct_from_string(inp_str)
    return automol.util.dict_.right_update(def_dct, inp_dct)


def _keyword_lst(inp_str, key_lst):
    """ build lst and check against supported lsts
        :param inp_dct: dictionary of input keywords
        :param def_dct: dictionary of default values
                        (only has ones where defaults, won't have vals for 
                         keys that user must input like run_prefix)
    """
    return ioformat.ptt.keyword_lst(inp_str)


def _pes_idxs(string):
    """  Build a dictionary of the PESs.
            {pes_idx: [chn_idxs]}
            breaks if pes_idx is given on two lines
    """
   
    run_pes = {}
    for line in string.strip().splitlines():
        [pes_nums, chn_nums] = line.split(':')
        pes_idxs = ioformat.ptt.parse_idxs(pes_nums)
        chn_idxs = ioformat.ptt.parse_idxs(chn_nums)
        for idx in pes_idxs:
            run_pes.update({idx: chn_idxs})

    return run_pes


def _spc_idxs(string):
    """  Build a dictionary of the PESs.
            {pes_idx: [chn_idxs]}
    """
   
    spc_idxs = ()
    for line in string.splitlines():
        spc_idxs += ioformat.ptt.parse_idxs(line)

    return {1: spc_idxs}


def _merge_glob(amech_dct):
    """ [change to pull out the glob dct and just use it to overwrite the spc dct
        or have parser return an additional one?
    """
    new_amech_dct = {}
    glob = amech_dct.get('global', {})
    for spc in (x for x in amech_dct if x != 'global'):
        new_amech_dct[spc] = right_update(glob, amech_dct[spc])

    return new_amech_dct


# Check the keyword dictionaries to see if proper things defined
def check_dictionary(inp_dct, chk_dct, section):
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
            if not val in allowed_vals:
                print('val is {}, must be {}'.format(val, allowed_vals))


def check_dictionary2(dct, req_keys, supp_keys):
    """ Make sure the theory dictionary keywords are all correct
    """
    req_keys_def = all(key in dct.keys() for key in req_keys)
    if not req_keys_def:
        print('*ERROR: Required keywords missing from thy section')
        print('level with issue: ', name)
        print('Required keys:')
        for key in THY_REQUIRED_KEYWORDS:
            print(key)
        sys.exit()
    sup_keys_def = all(key in supp_keys for key in dct.keys())
    if not sup_keys_def:
        print('*ERROR: unsupported required keywords missing from thy section')
        print('level with issue: ', name)
        print('Supported keys:')
        for key in THY_SUPPORTED_KEYWORDS:
            print(key)
        sys.exit()


def check_lst(inp_lst, sup_lst):
    """ Check 
    """
    if set(inp_lst) >= sup_lst:
        print('Unsupported keys')
