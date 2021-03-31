""" Functions to handle reading the files of the various inputs

    Perhaps return the strings of all the sections that can be parsed out of the file
    # Read the blocks, perhaps return a dct of the strings?

    Check if the strings are none
"""

import sys
import autoparse.find as apf
from ioformat import ptt
from mechlib.amech_io.parser.keywords import THY_REQUIRED_KEYWORDS
from mechlib.amech_io.parser.keywords import THY_SUPPORTED_KEYWORDS


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
    run_str = ioformat.ptt.read_inp_str(job_path, RUN_INP, remove_comments='#')
    ioprinter.reading('run.dat...', newline=1)  # Add a Found <file> to msg

    # Read the main blocks
    inp_block = _end_block(run_str, 'input')
    obj_block = _end_block(run_str, 'objs')
    job_block = _end_block(run_str, 'jobs')  # could read and set as neccessary
    es_tsks_block = _end_block(run_str, 'es_tsks')
    if es_tsks_block is None:
        print('no block')
    # trans_tsks_block = _end_block(run_str, 'trans_tsks')
    # therm_tsks_block = _end_block(run_str, 'therm_tsks')
    # ktp_tsks_block = _end_block(run_str, 'ktp_tsks')
    # print_tsks_block = _end_block(run_str, 'print_tsks')

    # Parse inp block
    key_dct = keyword_dct(inp_str, def_dct)

    # Parse objs blokcs
    run_dct['pes'] = get_pes_idxs(_paren_block(obj_str, 'pes'))
    run_dct['spc'] = get_spc_idxs(_paren_block(obj_str, 'spc'))
    # check the run_obj_dct 

    # Parse the job blocks
    key_lst = keyword_lst(inp_str, def_dct)

    # Parse tsks blokcs
    es_tsk_lst = tsks.es_tsk_lst(es_tsk_str, thy_dct)
    trans_tsk_lst = tsks.trans_tsk_lst(trans_tsk_str)
    prnt_tsk_lst = tsks.prnt_tsk_lst(prnt_tsk_str)

    # Check if needed strings exist

    if prnt_tsks_str is None:
        print('*ERROR: No "print_tsks" section defined in run.dat')
        # sys.exit()

    return None


def read_thy(job_path):
    """ Parse the theory.dat file
    """

    thy_str = ptt.read_inp_str(job_path, THY_INP, remove_comments='#')
    ioprinter.reading('theory.dat...', newline=1)

    thy_blocks = _named_end_blocks(string, 'level')
    if thy_blocks is not None:
        thy_methods = {}
        for section in thy_sections:
            name = section[0]
            keyword_dct = ptt.build_keyword_dct(thy_str)
            check_dictionary2(
                keyword_dct,
                THY_REQUIRED_KEYWORDS
                THY_SUPPORTED_KEYWORDS)
            thy_methods[name] = keyword_dct

    return thy_methods


def read_model(job_path):
    """ Parse the models.dat file
    """

    mod_str = ptt.read_inp_str(job_path, MOD_INP, remove_comments='#')
    ioprinter.reading('model.dat...', newline=1)

    pes_blocks = _named_end_blocks(mod_str, 'pes_model')
    spc_blocks = _named_end_blocks(mod_str, 'spc_model')


def read_spc(job_path):
    """ a
    """
    spc_str = ioformat.ptt.read_inp_str(job_path, CSV_INP)
    ioprinter.reading('species.csv...', newline=1)
    dat_str = ioformat.ptt.read_inp_str(job_path, DAT_INP)
    ioprinter.reading('species.dat...', newline=1)

    return spc_str


def _read_mech(job_path):
    """Build the PES dct
    """

    # Read the string
    mech_str = ptt.read_inp_str(job_path, MECH_INP, remove_comments='!')
    mech_info = util.read_mechanism_file(
        mech_str, mech_type, spc_dct, sort_rxns=sort_rxns)

    return mech_info


# Keyword dict build and check functions
def keyword_dct(inp_str, def_dct):
    """ merge a dictionary from the input
        :param inp_dct: dictionary of input keywords
        :param def_dct: dictionary of default values
                        (only has ones where defaults, won't have vals for 
                         keys that user must input like run_prefix)
    """
    inp_dct = ioformat.ptt.build_keyword_dct(inp_str)
    return automol.util.dict_.right_update(def_dct, inp_dct)


def keyword_lst(inp_str, key_lst):
    """ build lst and check against supported lsts
        :param inp_dct: dictionary of input keywords
        :param def_dct: dictionary of default values
                        (only has ones where defaults, won't have vals for 
                         keys that user must input like run_prefix)
    """
    return ioformat.ptt.build_keyword_lst(inp_str)


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
            print('val must be type {}'.format(allowed_typ)
        if allowed_vals:
            if not val in allowed_vals:
                print('val is {}, must be {}'.format(val, allowed_vals)


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


# Patterns of use
# maybe loop over different block types for keywords
def _end_block(string, header):
    """ A pattern for a certain block
        rtype: str
    """
    return ioformat.remove_whitespace(
        apf.first_capture(ioformat.ptt.end_section(header), string))


def _named_end_blocks(string, header):
    """ A pattern for a certian block
        rtype: dict[str: str]
    """ 
    caps = apf.all_captures(ptt.end_section_wname2(header), string)
    if caps is not None:
        caps = dict(zip((cap[0] for cap in caps), (cap[1] for cap in caps)))
    return caps


def _paren_block(string, header):
    """ A patter for a certain block
    """
    return apf.first_capture(ptt.paren_section(header), string)
