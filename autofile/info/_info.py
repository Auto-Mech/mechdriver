""" implements an class for YAML-style information
"""
import numbers
from types import SimpleNamespace
try:
    from collections.abc import Collection as _Collection
except ImportError:
    from collections import Collection as _Collection
import yaml
from autofile.info._inspect import function_keys as _function_keys


def object_(inf_dct):
    """ create an information object from a dictionary
    """

    def _cast(obj):
        if isinstance(obj, dict):
            ret = {key: _cast(val) for key, val in obj.items()}
            ret = Info(**ret)
        elif _is_nonstring_sequence(obj):
            ret = _normalized_nonstring_sequence(map(_cast, obj))
        else:
            ret = obj
        return ret

    inf_obj = _cast(inf_dct)
    return inf_obj


def dict_(inf_obj):
    """ convert an information object back to a dictionary
    """

    def _cast(obj):
        if isinstance(obj, Info):
            keys = obj.keys_()
            ret = {key: _cast(val)
                   for key, val in vars(obj).items() if key in keys}
        elif _is_nonstring_sequence(obj):
            ret = _normalized_nonstring_sequence(map(_cast, obj))
        else:
            ret = obj
        return ret
    inf_dct = _cast(inf_obj)
    return inf_dct


def string(inf_obj):
    """ write an information object to a YAML string
    """
    inf_dct = dict(inf_obj)
    inf_str = yaml.dump(inf_dct, default_flow_style=False)
    return inf_str


def from_string(inf_str):
    """ read an information object from a YAML string
    """
    inf_dct = yaml.load(inf_str, Loader=yaml.FullLoader)
    inf_obj = object_(inf_dct)
    return inf_obj


def matches_function_signature(inf_obj, function):
    """ does the information object match this function signature?
    """
    assert isinstance(inf_obj, Info)
    return inf_obj.keys_() == _function_keys(function)


class Info(SimpleNamespace):
    """ information container class, implemented as a frozen namespace

    (values can change, but you can't add keys after initialization)
    """
    _frozen = False

    def _freeze(self):
        self._frozen = True

    def __init__(self, **kwargs):
        kwargs = {key: (_normalized_nonstring_sequence(val) if
                        _is_nonstring_sequence(val) else val)
                  for key, val in kwargs.items()}
        super(Info, self).__init__(**kwargs)
        self._freeze()

    def keys_(self):
        """ keys for this instance """
        keys = frozenset(key for key in vars(self).keys()
                         if not key == '_frozen')
        return keys

    def __iter__(self):
        """ used by the dict() function for conversion to dictionary """
        for key, val in dict_(self).items():
            yield key, val

    def __eq__(self, other):
        return self.__dict__ == other.__dict__

    def __repr__(self):
        dct = dict_(self)
        items = ("{}={!r}".format(k, dct[k]) for k in sorted(dct.keys()))
        return "Info({})".format(", ".join(items))

    def __setattr__(self, key, value):
        """ prevent adding new keys after the object is frozen """
        if self._frozen and not hasattr(self, key):
            raise TypeError("'{}' object does not support item assignment"
                            .format(self.__class__.__name__))
        object.__setattr__(self, key, value)


def _normalized_nonstring_sequence(seq):
    return [
        int(val) if isinstance(val, numbers.Integral) else
        float(val) if isinstance(val, numbers.Real) else val for val in seq]


def _is_nonstring_sequence(obj):
    return (isinstance(obj, _Collection)
            and not isinstance(obj, (dict, str, bytes, bytearray)))
