""" library of parser functions for the run.dat file
"""

import sys
import autoparse.find as apf
from lib.amech_io.parser import ptt
from lib.amech_io.parser import tsks
from lib.amech_io.parser.keywords import RUN_INP_REQUIRED_KEYWORDS
from lib.amech_io.parser.keywords import RUN_SUPPORTED_KEYWORDS
from lib.amech_io.cleaner import remove_whitespace


RUN_INP = 'inp/run.dat'


# PARSE THE INPUT SECTION OF THE FILE #
def build_run_inp_dct(job_path):
    """ Build a dictionary for all the theory keywords
    """

    # Read the input section
    run_str = ptt.read_inp_str(job_path, RUN_INP)
    keyword_dct = ptt.build_keyword_dct(inp_block(run_str))

    # Check if section specified fully and supported
    check_run_keyword_dct(keyword_dct)

    return keyword_dct


def inp_block(inp_str):
    """ Read the string that has the global model information
    """
    return remove_whitespace(
        apf.first_capture(ptt.end_section('input'), inp_str))


def check_run_keyword_dct(dct):
    """ Make sure the theory dictionary keywords are all correct
    """
    if not dct:
        print('*ERROR: No "input" section specified in run.dat')
        sys.exit()

    kdef = all(key in dct.keys() for key in RUN_INP_REQUIRED_KEYWORDS)
    if not kdef:
        print('*ERROR: Not all required keywords specified',
              'in run.dat "input" section')
        print('*Required keys:')
        for key in RUN_INP_REQUIRED_KEYWORDS:
            print(key)
        sys.exit()


# PARSE THE OBJ SECTION OF THE FILE #
def objects_dct(job_path):
    """ Get the sections for the run block
    """

    # Read the obj section
    run_str = ptt.read_inp_str(job_path, RUN_INP)
    obj_str = object_block(run_str)

    # Read the sections of the obj section
    pes_block_str = apf.first_capture(ptt.paren_section('pes'), obj_str)
    spc_block_str = apf.first_capture(ptt.paren_section('spc'), obj_str)

    # Check if the obj section has been specified
    check_obj_spec(obj_str, pes_block_str, spc_block_str)

    # Build the run dictionary
    run_dct = {}
    if pes_block_str is not None:
        run_dct['pes'] = get_pes_idxs(remove_whitespace(pes_block_str))
    else:
        run_dct['pes'] = []
    if spc_block_str is not None:
        run_dct['spc'] = get_spc_idxs(remove_whitespace(spc_block_str))
    else:
        run_dct['spc'] = []

    return run_dct


def object_block(inp_str):
    """ Read the string that has the global model information
    """
    return remove_whitespace(
        apf.first_capture(ptt.end_section('obj'), inp_str))


def get_pes_idxs(pes_str):
    """ Determine the indices corresponding to what PES and channels to run
    """
    run_pes = {}
    for line in pes_str.splitlines():
        try:
            [pes, chns] = line.strip().split(';')
            [pes_idx_str, pes_mod_str] = pes.strip().split()
            [chns_idx_str, chns_mod_str] = chns.strip().split()
            pes_lst = ptt.parse_idx_inp(pes_idx_str)
            chn_lst = ptt.parse_idx_inp(chns_idx_str)
            pes_mod = pes_mod_str.strip()
            chns_mod = chns_mod_str.strip()
            for pes in pes_lst:
                for chn in chn_lst:
                    run_pes[(pes, chn)] = (pes_mod, chns_mod)
        except:
            print('*ERROR: pes line formatted incorrectly, cannot be parsed')
            print('Line:   ', line)
            sys.exit()

    return run_pes


def get_spc_idxs(pes_str):
    """ Determine the indices corresponding to what species to run
    """
    run_spc = {}
    for line in pes_str.splitlines():
        try:
            [spc_idxs, proc1, proc2] = line.strip().split()
            spc_lst = ptt.parse_idx_inp(spc_idxs)
            proc1 = proc1.strip()
            proc2 = proc2.strip()
            for spc in spc_lst:
                run_spc[spc] = (proc1, proc2)
        except:
            print('*ERROR: spc line formatted incorrectly, cannot be parsed')
            print('Line:   ', line)
            sys.exit()

    return run_spc


def check_obj_spec(obj_str, pes_block_str, spc_block_str):
    """ Check the obj string
    """
    if obj_str is None:
        print('*ERROR: No "obj" section specified in run.dat')
        sys.exit()
    else:
        if pes_block_str is None and spc_block_str is None:
            print('*ERROR: No "pes" or "spc" requested in the "obj" section ',
                  'specified in run.dat')
            sys.exit()


# PARSE THE JOBS SECTION OF THE FILE #
def build_run_jobs_lst(job_path):
    """ Build a dictionary for all the theory keywords
    """

    # Read the jobs section
    job_str = ptt.read_inp_str(job_path, RUN_INP)
    keyword_lst = ptt.build_keyword_lst(jobs_block(job_str))

    # Check the jobs sectuib
    check_run_jobs_section(job_str, keyword_lst)

    return keyword_lst


def jobs_block(inp_str):
    """ Read the string that has the global model information
    """
    return remove_whitespace(
        apf.first_capture(ptt.end_section('jobs'), inp_str))


def set_thermodriver(run_jobs_lst):
    """ Set vars for the thermodriver run using the run_jobs_lst
    """
    write_messpf = False
    run_messpf = False
    run_nasa = False
    if 'thermochem' in run_jobs_lst:
        write_messpf = True
        run_messpf = True
        run_nasa = True
    else:
        if 'write_messpf' in run_jobs_lst:
            write_messpf = True
        if 'run_messpf' in run_jobs_lst:
            run_messpf = True
        if 'run_nasa' in run_jobs_lst:
            run_nasa = True

    return write_messpf, run_messpf, run_nasa


def set_ktpdriver(run_jobs_lst):
    """ Set vars for the ktpdriver run using the run_jobs_lst
    """
    write_messrate = False
    run_messrate = False
    run_fits = False
    if 'kinetics' in run_jobs_lst:
        write_messrate = True
        run_messrate = True
        run_fits = True
    else:
        if 'write_messrate' in run_jobs_lst:
            write_messrate = True
        if 'run_messrate' in run_jobs_lst:
            run_messrate = True
        if 'run_fits' in run_jobs_lst:
            run_fits = True

    return write_messrate, run_messrate, run_fits


def check_run_jobs_section(job_str, keyword_lst):
    """ Check the run jobs section
    """
    if not job_str:
        print('*ERROR: No "jobs" section given in run.dat')
        sys.exit()

    if not keyword_lst:
        print('*ERROR: No keyword given in jobs')
        sys.exit()

    chk = all(key in RUN_SUPPORTED_KEYWORDS for key in keyword_lst)
    if not chk:
        print('*ERROR: Not all keywords specified',
              'in run.dat "jobs" section are supported')
        print('*Allowed keys:')
        for key in RUN_SUPPORTED_KEYWORDS:
            print(key)
        sys.exit()


# PARSE THE ES_TSKS SECTION OF THE FILE #
def read_es_tsks(job_path):
    """ Build a dictionary for all the theory keywords
    """

    # Read the electronic structure tasks section
    run_str = ptt.read_inp_str(job_path, RUN_INP)
    es_tsks_str = es_tsks_block(run_str)

    # Check if section is there
    if es_tsks_str is None:
        print('*ERROR: No "es_tsks" section defined in run.dat')
        sys.exit()

    return es_tsks_str


def es_tsks_block(inp_str):
    """ Read the string that has the global model information
    """
    return remove_whitespace(
        apf.first_capture(ptt.end_section('es_tsks'), inp_str))


def build_run_es_tsks_lst(es_tsk_str, rxn_model_dct, thy_dct, saddle=False):
    """ Build the list of ES tasks, potentially w/ models
    """
    es_tsk_lst = tsks.es_tsk_lst(
        es_tsk_str, rxn_model_dct, thy_dct, saddle=saddle)

    return es_tsk_lst
