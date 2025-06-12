""" implements strings for option specification

thought: we could have default values specified through the keys like so:
    MAXITER_ = create('scf_maxiter', ['num=50'])
"""

from string import Formatter as _Formatter
import re
import yaml


FORMAT = '<{}>[{}]'
PATTERN = r'<(.*)>(\[.*\])'


def create(name_, keys_=()):
    """ Create an option specifier.

        :param name:  name of the option to create
        :type name: str
        :param keys: provide keys associated with the option
        :type keys: tuple(obj)
    """

    keys_str = ','.join(f'{{{key}}}' for key in keys_)

    return FORMAT.format(name_, keys_str)


def specify(osp_, *args, **kwargs):
    """ set values for an option specifier
        Function allow positional arguments, even though .format does not

        :param osp_: option specifier
    """
    # allow positional arguments, even though .format doesn't
    _keys = keys(osp_)
    nkeys = len(_keys)
    assert len(args) + len(kwargs) == nkeys
    _kwargs = dict(zip(_keys, args[:nkeys]))
    assert not set(_kwargs) & set(kwargs)
    kwargs.update(_kwargs)

    osp = osp_.format(**kwargs)
    return osp


def is_valid(osp_):
    """ is this an option specifier? (either)
    """
    match = re.fullmatch(PATTERN, osp_)
    return match is not None


def name(osp_):
    """ get the name from an option specifier (either)
    """
    assert is_valid(osp_)
    match = re.fullmatch(PATTERN, osp_)
    name_, _ = match.groups()
    return name_


def is_valueless(osp_):
    """ is this string option specifier free of values?
    """
    assert is_valid(osp_)
    val_str = _value_string(osp_)
    return not yaml.load(val_str, Loader=yaml.FullLoader)


def is_template(osp_):
    """ is this a string option specifier?
    """
    assert is_valid(osp_)
    val_str = _value_string(osp_)
    return bool(_value_string_keys(val_str))


def keys(osp_):
    """ get the keys from a option specifier (template)
    """
    assert is_template(osp_) or is_valueless(osp_)
    val_str = _value_string(osp_)
    return _value_string_keys(val_str)


def values(osp):
    """ get the values from an option specifier (string)
    """
    assert (not is_template(osp)) or is_valueless(osp)
    val_str = _value_string(osp)
    return tuple(yaml.load(val_str, Loader=yaml.FullLoader))


def _value_string(osp_):
    """ get the value string from an option specifier (either)
    """
    assert is_valid(osp_)
    match = re.fullmatch(PATTERN, osp_)
    _, val_str = match.groups()
    return val_str


def _value_string_keys(val_str):
    return tuple(fpar[1] for fpar in _Formatter().parse(val_str)
                 if fpar[1] is not None)
