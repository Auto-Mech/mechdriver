""" Check the input info and set defaults
"""

# check dcts
def _new_check(fname, sname, dct, req_dct):
    """ New section check function

        :param sname: section name
        :param dct: dict to be assessed
        :param req_dct: dict of required keywords
    """

    print('Checking {} file...'.format(fname))
    print('Checking {} section...'.format(sname))
    
    # Check if nothing specified in section
    if not dct:
        print('*ERROR: No input given in {} section specified in {}'.format(sname, fname))
        sys.exit()

    # Check if required keywords provided
    kdef = all(key in dct.keys() for key in RUN_INP_REQUIRED_KEYWORDS)
    if not kdef:
        print('*ERROR: Not all required keywords specified',
              'in run.dat "input" section')
        print('*Required keys:')
        for key in RUN_INP_REQUIRED_KEYWORDS:
            print(key)
        sys.exit()



# run
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

    ksup = all(key in RUN_INP_SUPPORTED_KEYWORDS for key in dct.keys())
    if not ksup:
        print('*ERROR: Using unsupported ketwords',
              'in run.dat "input" section')
        print('*Supported keys:')
        for key in RUN_INP_SUPPORTED_KEYWORDS:
            print(key)
        sys.exit()

    if dct['mech'] not in RUN_INP_KEY_DCT['mech']:
        print('*ERROR: Unallowed value for mech keyword')
        sys.exit()
    if dct['spc'] not in RUN_INP_KEY_DCT['spc']:
        print('*ERROR: Unallowed value for spc keyword')
        sys.exit()


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

# thy
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
    if not sup_keys_def:
        print('*ERROR: unsupported required keywords missing from thy section')
        print('level with issue: ', name)
        print('Supported keys:')
        for key in THY_SUPPORTED_KEYWORDS:
            print(key)
        sys.exit()

