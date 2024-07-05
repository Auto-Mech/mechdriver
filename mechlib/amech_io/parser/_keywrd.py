""" Libraries of supported keywords and their corresponding values

    Also includes functionalities for constructing dictionaries
    of default values as well as assessing the validity of user input.
"""

import sys
from copy import deepcopy

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
    """Building the default dictionary for run inp dictionary

    prob works for spc as well
    """
    supp_keywrds = tuple(dct.keys())
    default_dct = dict(zip(supp_keywrds, (dct[key][2] for key in supp_keywrds)))

    return default_dct


def defaults_from_key_val_dcts(key, key_dct, val_dct):
    """Way of building the default dcts for various things, this is
    for tasks blocks

    only works for tsks since it's info is in two dcts
    """

    # Set all of the keywords that are allowed for a task
    keywrds = key_dct[key][1]

    # Now build a dct where all the keywords are defaulted to internal value
    default_dct = dict(zip(keywrds, (val_dct[kwrd][2] for kwrd in keywrds)))

    return default_dct


def defaults_with_dcts(dct):
    """Build defaults for something that could have dcts

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
    """Check all of the facets of a dictionary

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
    """Check if all keywordws supplied in an input dictionary

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
        print(f"User defined unsupported keywords in {section}")
        print(f'Unsupported keywords: {(",".join(unsupported_keys))}')
        print(f'Accepted keywords: {(",".join(chk_keys))}')
        sys.exit()


def _check_supported_vals(inp_dct, val_dct, req_lst, section):
    """check dct function 2

    if val in inp_dct is None? what to do?
    maybe check if None if it is required lst
    """

    # Assess if the keywords have the appropriate value
    for key, val in inp_dct.items():
        allowed_typs, allowed_vals, _ = val_dct[key]

        if val is not None:
            # Check val if one is given
            if type(val) not in allowed_typs:
                print(f"bad {section}".format(section))
                print(f"{key}".format(key))
                print(f"val {val} must be type {allowed_typs}")
                sys.exit()
            if allowed_vals:
                if val not in allowed_vals:
                    print(f"bad {section}".format(section))
                    print(f"{key}".format(key))
                    print(f"val is {val}, must be {allowed_vals}")
                    sys.exit()
        else:
            # If val is None, check if it is required
            if key in req_lst:
                print(f"bad {section}".format(section))
                print(f"key {key} has no value defined even " "though it is required")
                sys.exit()


def _check_required_keys(inp_dct, req_lst, section):
    """Check if required keys are in the input dict"""

    inp_keys = set(inp_dct.keys())
    req_keys = set(req_lst)
    undefined_required_keys = req_keys - inp_keys

    if undefined_required_keys:
        print(f"Required keywords have not been defined in {section}")
        for key in undefined_required_keys:
            print(key)
        sys.exit()


def check_thy_lvls(key_dct, method_dct, section=""):
    """For specific tasks, we need to a second level of cheeck to ensure
    that values of keywords that correspond to blocks defined in either
    the thy or model dat files.

    :param key_dct:
    :type key_dct: dict[]
    :param method_dct: thy or mod dct
    """

    thy_defined_methods = set(method_dct.keys())

    for key in ("runlvl", "inplvl", "var_splvl1", "var_splvl2", "var_scnlvl"):
        val = key_dct.get(key)
        if val is not None:
            if val not in thy_defined_methods:
                print(f"User has not defined val in {section}")
                print(key, val)
                sys.exit()


def check_model_combinations(pf_dct):
    """Check if a model combination is not implemented for PF routines"""
    if pf_dct["vib"] == "vpt2" and pf_dct["tors"] == "1dhr":
        print("*ERROR: VPT2 and 1DHR combination is not yet implemented")
        sys.exit()
    elif pf_dct["vib"] == "vpt2" and "tau" in pf_dct["tors"]:
        print("*ERROR: VPT2 and TAU combination is not yet implemented")
        sys.exit()


# helpers
def empty_if_none(obj):
    """Returns an empty dictionary if input object is None,
    otherwise return the object.

    :param obj: generic object
    :type obj: any
    """
    return {} if obj is None else obj


def without_nones(dct: dict) -> dict:
    """Remove `None` values from a nested dictionary

    :param dct: A dictionary
    :return: The dictionary, with `None` values removed
    """
    dct = dct.copy()
    for key, val in dct.items():
        if isinstance(val, dict):
            dct[key] = without_nones(val)
        elif val is not None:
            dct[key] = val
    return dct


def right_update(
    dct1: dict, dct2: dict, nested: bool = True, drop_none: bool = False
) -> dict:
    """Updates the entries of `dct1` with those of `dct2`.

    Defaults to recursive update for nested dictionaries

    :param dct1: dictionary1 that will be updated
    :param dct2: dictionary2 whose entries will override dct1
    :param drop_none: Drop `None` values from the updated dictionary?
    :return: The updated dictionary
    """

    dct = {}
    dct1 = empty_if_none(dct1)
    dct2 = empty_if_none(dct2)

    dct = dct1.copy()
    for key, val in dct2.items():
        if (
            nested
            and key in dct
            and isinstance(dct[key], dict)
            and isinstance(val, dict)
        ):
            dct[key] = right_update(dct[key], val)
        else:
            dct[key] = val

    if drop_none:
        dct = without_nones(dct)

    return dct


def separate_subdct(dct: dict, key="global"):
    """Pulls out a sub-dictionary indexed by the given key and returns it
    and the original dictioanry with the requested sub-dictionary removed
    """

    # Grab the sub-dictonary and
    sub_dct = dct.get(key, {})

    # Build a new copy of the input dictionary with the sub-dict removed
    dct2 = deepcopy(dct)
    dct2.pop(key, None)

    return dct2, sub_dct


def merge_subdct(dct: dict, key="global", keep_subdct=False):
    """Obtain a sub-dictionary indexed by a given key and merge its contents
    with all of the other sub-dictionaries in the main dictionary.

    :param dct: dictionary containing several sub-dictionaries
    :type dct: dict[str: dict]
    :param key: key for the sub-dictionary to be merged
    :type key: str, int, float, tuple
    :param keep_subdct: keep/remove the sub-dictionary following the merge
    :type keep_subdct: bool
    """

    if list(dct.keys()) == [key]:
        new_dct = {key: dct[key]}
    else:
        new_dct, sub_dct = separate_subdct(dct, key=key)
        for new_key in new_dct:
            new_dct[new_key] = right_update(
                sub_dct, new_dct[new_key], nested=False, drop_none=False
            )
        if keep_subdct:
            new_dct[key] = sub_dct

    return new_dct
