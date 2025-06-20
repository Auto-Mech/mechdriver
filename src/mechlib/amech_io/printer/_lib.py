"""
  Library of strings
"""

from ioformat import addchar
from mechlib.amech_io.printer._print import message


LINE_PLUS = (
    '++++++++++++++++++++++++++++++++++++++++++++++++' +
    '++++++++++++++++++++++++++++++++++++++'
)

LINE_DASH = (
    '------------------------------------------------' +
    '--------------------------------------'
)

DOUBLE_NEW_LINE = '\n\n'

LIB_DCT = {
    'line_plus': LINE_PLUS,
    'line_dash': LINE_DASH,
    'vspace': DOUBLE_NEW_LINE
}


def objs(keys, prechars=()):
    """ Print a set of standard objects
    """

    obj_str = ''
    for key, prechar in zip(keys, prechars):
        obj_str += obj(key, prechar=prechar)

    message(obj_str)


def obj(key, prechar=None):
    """ Print some standard object
    """

    assert key in LIB_DCT, (
        f'Object {key} not in library'
    )

    obj_str = LIB_DCT[key]
    if prechar is not None:
        obj_str = addchar(obj_str, prechar, side='pre')

    message(obj_str)
