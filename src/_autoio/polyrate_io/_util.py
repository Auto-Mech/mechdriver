"""
utility functions
"""

from itertools import chain
import automol
from ioformat import indent


def pt_format(header, hess, vlabel, vval,
              slabel=None, sval=None, geo=None, grad=None):
    """ write a point string
    """

    # Intialize with header
    pt_str = f'*{header}'
    pt_str += '\n\n'

    # Write the energy and coordinate along reaction coordinate
    pt_str += energy_format(vlabel, vval)
    pt_str += '\n'
    if sval is not None:
        pt_str += energy_format(slabel, sval)
        pt_str += '\n'
    pt_str += '\n'

    # Write the structurea information
    if geo is not None:
        pt_str += geometry_format(geo)
        pt_str += '\n\n'
    if grad is not None:
        pt_str += gradient_format(grad)
        pt_str += '\n\n'
    pt_str += hessian_format(hess)
    pt_str += '\n'

    pt_str = indent(pt_str, 2)

    return pt_str


def energy_format(label, ene):
    """ write an energy
    """

    assert label in ('smep', 'vmep', 'svalue', 'vvalue'), (
        f'Label {label} != smep, vmep, svalue, or vvalue'
    )

    ene_str = f'{label:<8s}{ene:<14.12f}'

    return ene_str


def geometry_format(geo):
    """ Write geom
    """

    xyzs = tuple(xyz for _, xyz in geo)
    geo_str = automol.util.matrix.string(xyzs, val_format='{:>12.8f}')

    return _end_format('geom', 'end', geo_str)


def gradient_format(grad):
    """ format hessian
    """

    grad_str = automol.util.matrix.string(grad, val_format='{:>12.8f}')

    return _end_format('grads', 'end', grad_str)


def hessian_format(hess):
    """ format hessian
    """

    hess = list(chain.from_iterable(hess))
    hess_str = automol.util.vector.string(
        hess, num_per_row=6, val_format='{0:>12.8f}')

    return _end_format('hessian', 'end', hess_str)


def list_format(lst):
    """ Unpack a list of values and write them to a string
    """
    return '\n'.join(lst)


def _end_format(header, ender, dat_str):
    """ Write a block with an end
    """
    return (
        header + '\n' +
        dat_str + '\n' +
        ender
    )
