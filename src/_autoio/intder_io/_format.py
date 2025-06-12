""" Functions to format the data into strings appropriate
"""

import os
from itertools import chain
import automol.geom
import automol.zmat
import automol.util
from ioformat import build_mako_str
from ioformat import indent


# OBTAIN THE PATH TO THE DIRECTORY CONTAINING THE TEMPLATES #
SRC_PATH = os.path.dirname(os.path.realpath(__file__))
TEMPLATE_PATH = os.path.join(SRC_PATH, 'templates')


def header_format(geo):
    """ format the header
    """

    natom = automol.geom.count(geo)
    if not automol.geom.is_linear(geo):
        nintl = 3 * natom - 6
    else:
        nintl = 3 * natom - 5

    header_keys = {
        'natom': f'{natom:>5d}',
        'nintl': f'{nintl:>5d}',
        'comment': 'TED Calculation'
    }

    return build_mako_str(
        template_file_name='header.mako',
        template_src_path=TEMPLATE_PATH,
        template_keys=header_keys,
        remove_whitespace=False)


def internals_format(zma):
    """ Format the strings with the internal coordinates
    """

    key_mat = automol.zmat.key_matrix(zma)

    # Write the stretch, bend, and torsion coordinates
    istr = ''
    for i, row in enumerate(key_mat[1:]):
        istr += f'STRE{i+1+1:>6d}{row[0]+1:>5d}\n'
    for i, row in enumerate(key_mat[2:]):
        istr += f'BEND{i+2+1:>6d}{row[0]+1:>5d}{row[1]+1:>5d}\n'
    for i, row in enumerate(key_mat[3:]):
        istr += f'TORS{i+3+1:>6d}{row[0]+1:>5d}{row[1]+1:>5d}{row[2]+1:>5d}\n'

    # Remove final newline character
    istr = istr.rstrip()

    return istr


def geometry_format(geo):
    """ Format the geometry string
    """

    # Build geom str
    geo_str = ''
    for (_, xyz) in geo:
        geo_str += f'{xyz[0]:>14.5f}{xyz[1]:>14.5f}{xyz[2]:>14.5f}\n'

    # Indent the lines and remove final newline character
    geo_str = indent(geo_str, 4)
    geo_str = geo_str.rstrip()

    return geo_str


def symbols_format(geo):
    """ Format the symbols
    """

    symbs = automol.geom.symbols(geo)
    symb_str = automol.util.vector.string(
        symbs, num_per_row=6, val_format='{0:>12s}')

    return symb_str


def hessian_format(hess):
    """ Format a mass-weighted Hessian into a string for the
        auxiliary input file for INTDER.

        :param hess: mass-weighted Hessian
        :type hess: numpy.ndarray
        :rtype: str
    """

    # Flatten the Hessian out for ease of writing the string
    hess = tuple(chain.from_iterable(hess))
    hess_str = automol.util.vector.string(
        hess, num_per_row=3, val_format='{0:>20.10f}')

    return hess_str
