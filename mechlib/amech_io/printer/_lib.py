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


def obj(key, prechar=None):
    """ Print some standard object
    """

    assert key in LIB_DCT, (
        'Object {} not in library'.format(key)
    )

    obj_str = LIB_DCT[key]
    if prechar is not None:
        obj_str = addchar(obj_str, prechar, side='pre')

    message(obj_str)

    return obj_str
