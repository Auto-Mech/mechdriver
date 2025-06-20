"""
  Additional functions for formatting information for MESS strings
"""

import numpy
from phydat import phycon, ptab
from ioformat import remove_trail_whitespace


def write_data_str(geos, grads, hessians):
    """ Writes a string containing the geometry, gradient, and Hessian
        for either a single species or points along a reaction path
        that is formatted appropriately for the ProjRot input file.

        :param geos: geometries
        :type geos: list
        :param grads: gradients
        :type grads: list
        :param hessians: Hessians
        :type hessians: list
        :rtype: str
    """

    nsteps = len(geos)

    data_str = ''
    for i, (geo, grad, hess) in enumerate(zip(geos, grads, hessians)):
        data_str += f'Step    {str(i+1)}\n'
        data_str += 'geometry\n'
        data_str += _format_geo_str(geo)
        data_str += 'gradient\n'
        data_str += _format_grad_str(geo, grad)
        data_str += 'Hessian\n'
        data_str += _format_hessian_str(hess)
        if i < nsteps-1:
            data_str += '\n'

    return remove_trail_whitespace(data_str)


def _format_geo_str(geo):
    """ Formats a geometry into a string used for the ProjRot input file.

        :param geo: geometries (Angstrom)
        :type geo: list
        :rtype: str
    """

    # Format the strings for the xyz coordinates
    geo_str = ''
    for i, (sym, xyzs) in enumerate(geo):
        anum = int(ptab.to_number(sym))
        xyzs = [xyz * phycon.BOHR2ANG for xyz in xyzs]
        xyzs_str = f'{xyzs[0]:>14.8f}{xyzs[1]:>14.8f}{xyzs[2]:>14.8f}'
        geo_str += f'{i+1:2d}{anum:4d}{0:4d}{xyzs_str}\n'

    return remove_trail_whitespace(geo_str)


def _format_grad_str(geo, grad):
    """ Formats a gradient into a string used for the ProjRot input file.

        :param geom: geometries (Angstrom)
        :type geom: list
        :param grads: gradients (Eh/Bohr)
        :type grads: list
        :rtype: str
    """

    atom_list = []
    for i, (sym, _) in enumerate(geo):
        atom_list.append(int(ptab.to_number(sym)))

    # Format the strings for the xyz gradients
    full_grads_str = ''
    for i, grads in enumerate(grad):
        grads_str = f'{grads[0]:>14.8f}{grads[1]:>14.8f}{grads[2]:>14.8f}'
        full_grads_str += f'{i+1:2d}{atom_list[i]:4d}{grads_str}\n'

    return remove_trail_whitespace(full_grads_str)


def _format_hessian_str(hess):
    """ Formats a Hessian into a string used for the ProjRot input file.

        :param hess: hessians (Eh/Bohr^2)
        :type hess: list
        :rtype: str
    """

    # Format the Hessian
    hess = numpy.array(hess)
    nrows = numpy.shape(hess)[0]
    ncols = numpy.shape(hess)[1]

    if nrows % 5 == 0:
        nchunks = nrows // 5
    else:
        nchunks = (nrows // 5) + 1

    hess_str = '   ' + '    '.join([str(val) for val in range(1, 6)]) + '\n'
    cnt = 0
    while cnt+1 <= nchunks:
        for i in range(nrows):
            col_tracker = 1
            if i >= 5*cnt:
                hess_str += f'{str(i+1)}'
                for j in range(5*cnt, ncols):
                    if i >= j:
                        if col_tracker <= 5:
                            hess_str += f'  {hess[i][j]:5.8f}'
                            col_tracker += 1
                            if col_tracker == 6:
                                hess_str += '\n'
                        else:
                            continue
                    elif i < j and col_tracker != 6:
                        hess_str += '\n'
                        break
                    else:
                        break
            if i+1 == nrows and cnt+1 < nchunks-1:
                val_str = '     '.join(
                    [str(val)
                     for val in range(5*(cnt+1) + 1, 5*(cnt+1) + 6)])
                hess_str += '    ' + val_str + '\n'
            if i+1 == nrows and cnt+1 == nchunks-1:
                val_str = '     '.join(
                    [str(val)
                     for val in range(5*(cnt+1) + 1, nrows+1)])
                hess_str += '    ' + val_str + '\n'
        cnt += 1

    return remove_trail_whitespace(hess_str)
