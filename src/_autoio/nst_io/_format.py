""" Format data for NST
"""

import automol.geom


def format_geometry(geo):
    """ Format all of the info related to the TS geometry.

        :param geo: guess geometry for a NST TS
        :type geo: automol.geom object
    """

    # Get number of atoms
    natoms = automol.geom.count(geo)

    # Format geometry into a string
    geo_str = automol.geom.string(geo)

    # Put all of the masses into the string
    mass_str = ''
    for mass in automol.geom.masses(geo):
        mass_str += f'{mass}\n'
    mass_str = mass_str.rstrip()

    return natoms, geo_str, mass_str


def format_hessian(hess):
    """ Format the hessian into a file to be
        read by the NST code.
    """

    hess_str = ''

    ndim = len(hess)
    for i in range(ndim):
        for j in range(ndim):
            hess_str += f'{hess[i][j]}\n'

    return hess_str
