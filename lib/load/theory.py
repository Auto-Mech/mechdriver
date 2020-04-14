""" library of reader functions for the theory file
"""

import sys
import autoparse.find as apf
from lib.load import ptt
from lib.load.keywords import THY_REQUIRED_KEYWORDS, THY_SUPPORTED_KEYWORDS


THY_INP = 'inp/theory.dat'


def build_thy_dct(job_path):
    """ species input
    """
    # Obtain the species string
    thy_str = ptt.read_inp_str(job_path, THY_INP)
    
    # Obtain any thy_sections
    thy_sections = apf.all_captures(
        ptt.end_section_wname2('level'), thy_str)
    check_thy_sections_nonempty(thy_sections)

    # Build dictionary of theory methods
    thy_methods = {}
    for section in thy_sections:
        # Check formatting...
        name = section[0]
        keyword_dct = build_thy_keyword_dct(name, section[1])
        thy_methods[name] = keyword_dct

    return thy_methods


def build_thy_keyword_dct(name, thy_str):
    """ Build a dictionary for all the theory keywords
    """
    keyword_dct = ptt.build_keyword_dct(thy_str)
    check_thy_dct(name, keyword_dct)

    return keyword_dct


def check_thy_sections_nonempty(thy_sections):
    """ Make sure the theory dictionary keywords are all correct
    """
    if thy_sections is None:
        print('*ERROR: No level sections defined in theory.day')
        sys.exit()


def check_thy_dct(name, dct):
    """ Make sure the theory dictionary keywords are all correct
    """
    req_keys_def = all(key in dct.keys() for key in THY_REQUIRED_KEYWORDS)
    if not req_keys_def:
        print('*ERROR: Required keywords missing from thy section')
        print('level with issue: ', name)
        print('Required keys:')
        for key in THY_REQUIRED_KEYWORDS:
            print(key)
        sys.exit()
    sup_keys_def = all(key in THY_SUPPORTED_KEYWORDS for key in dct.keys())
    if not req_keys_def:
        print('*ERROR: Non-supported Required keywords missing from thy section')
        print('level with issue: ', name)
        print('Supported keys:')
        for key in THY_SUPPORTED_KEYWORDS:
            print(key)
        sys.exit()
