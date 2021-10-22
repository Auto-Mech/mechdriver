""" Libraries of supported keywords and their corresponding values

    Also includes functionalities for constructing dictionaries
    of default values as well as assessing the validity of user input.
"""

import sys


# Might just be easier to reinstitute required lists again.
# Right now it depends on if each val dct is used multiple times or once
# then a required might be added (think this fails for tasks dicts)

# Functions needed to build custom values
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


# Default Dictionary Builders
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


def defaults_with_dcts(dct):
    """ Build defaults for something that could have dcts

        now just dcts in MODPF, could also have lsts with MODKIN with
        additional lines
    """

    default_dct = {}
    for keywrd in dct.keys():
        val = dct[keywrd]
        if isinstance(val, dict):
            keywrds2 = tuple(val.keys())
            newv = dict(zip(keywrds2, (val[kwrd][2] for kwrd in keywrds2)))
        elif isinstance(val, tuple):
            newv = val[2]
        default_dct[keywrd] = newv

    return default_dct


# Dictionary Checkers
def check_dct1(inp_dct, val_dct, req_lst, section):
    """ Check all of the facets of a dictionary

        Works if section is just a val dct
    """
    _check_required_keys(inp_dct, req_lst, section)
    _check_supported_keys(inp_dct, val_dct, section)
    _check_supported_vals(inp_dct, val_dct, req_lst, section)


# def check_dct_with_dcts(dct, val_dct):
#     """ Check dictionaries (BROKEN)
#     """
#     for keywrd in dct.keys():
#         val = dct[keywrd]
#         if isinstance(val, dict):
#             keywrds2 = tuple(val.keys())
#             vals2 =
#             newv = dict(zip(keywrds2, (val[kwrd][2] for kwrd in keywrds2)))
#         else:
#             newv = val[2]
#     _check_required_keys(inp_dct, req_lst, section)
#     _check_supported_keys(inp_dct, val_dct, section)
#     _check_supported_vals(inp_dct, val_dct, req_lst, section)


def _check_supported_keys(inp_dct, val_dct, section):
    """ Check if all keywordws supplied in an input dictionary

        :param inp_dct: input dictionary to assess
        :type inp_dct: dict[]
        :param val_dct: internal dictionary containing supported keywords
        :type val_dct: dict[]
        :param section: Label for what section of the input is being checked
        :type section: str
    """

    inp_keys = set(inp_dct.keys())
    chk_keys = set(val_dct.keys())
    unsupported_keys = inp_keys - chk_keys

    if unsupported_keys:
        print(f'User defined unsupported keywords in {section}')
        print(f'Unsupported keywords: {(",".join(unsupported_keys))}')
        print(f'Accepted keywords: {(",".join(chk_keys))}')
        sys.exit()


def _check_supported_vals(inp_dct, val_dct, req_lst, section):
    """ check dct function 2

        if val in inp_dct is None? what to do?
        maybe check if None if it is required lst
    """

    # Assess if the keywords have the appropriate value
    for key, val in inp_dct.items():
        allowed_typs, allowed_vals, _ = val_dct[key]

        if val is not None:
            # Check val if one is given
            if type(val) not in allowed_typs:
                print(f'bad {section}'.format(section))
                print(f'{key}'.format(key))
                print(f'val {val} must be type {allowed_typs}')
                sys.exit()
            if allowed_vals:
                if val not in allowed_vals:
                    print(f'bad {section}'.format(section))
                    print(f'{key}'.format(key))
                    print(f'val is {val}, must be {allowed_vals}')
                    sys.exit()
        else:
            # If val is None, check if it is required
            if key in req_lst:
                print(f'bad {section}'.format(section))
                print(f'key {key} has no value defined even '
                      'though it is required')
                sys.exit()


def _check_required_keys(inp_dct, req_lst, section):
    """ Check if required keys are in the input dict
    """

    inp_keys = set(inp_dct.keys())
    req_keys = set(req_lst)
    undefined_required_keys = req_keys - inp_keys

    if undefined_required_keys:
        print(f'Required keywords have not been defined in {section}')
        for key in undefined_required_keys:
            print(key)
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
                print(f'User has not defined val in {section}')
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
