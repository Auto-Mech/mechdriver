"""
 Utility functions for formatting
"""

from qcelemental import constants as qcc
from qcelemental import periodictable as ptab
from automol import geom

# Conversion factors
BOHR2ANG = qcc.conversion_factor('bohr', 'angstrom')


def determine_struct_type(geo):
    """ determines the linear string
    """

    # Use automol to determine the type of structure
    if geom.is_atom(geo):
        struct_type = 'Monoatomic'
    else:
        if geom.is_linear(geo):
            struct_type = 'Linear'
        else:
            struct_type = 'Nonlinear'

    return struct_type


def format_coords(geo):
    """ format the coords section
    """

    # Get the number of atoms
    natoms = len(geo)

    # Get the geometry information
    symbols = geom.symbols(geo)
    coordinates = geom.coordinates(geo)
    masses = [int(ptab.to_mass(symbol)) for symbol in symbols]

    # Build a string with the formatted coordinates string
    if geom.is_atom(geo):
        geo_str = '{0:<4s}{1:<6d}'.format(symbols[0], masses[0])
    else:
        geo_str = ''
        for symbol, mass, coords in zip(symbols, masses, coordinates):
            coords = [coord * BOHR2ANG for coord in coords]
            coords_str = '{0:>14.8f}{1:>14.8f}{2:>14.8f}'.format(
                coords[0], coords[1], coords[2])
            geo_str += '{0:<4s}{1:<6d}{2}\n'.format(
                symbol, mass, coords_str)
        # Remove final newline character from the string
        geo_str = geo_str.rstrip()

    return natoms, geo_str


def format_values_string(coord, values, conv_factor=1.0):
    """ format the values string for the divsur.inp file
    """
    if values:
        values = ', '.join('{0:.3f}'.format(val * conv_factor)
                           for val in values)
        values_string = '{0} = ({1})'.format(coord, values)
    else:
        values_string = ''

    return values_string


def format_pivot_xyz_string(idx, npivot, xyzP):
    """ format the pivot point xyz
    """

    assert npivot in (1, 2)

    atom_idx = idx
    if idx == 1:
        d_idx = 1
    else:
        d_idx = 2

    if npivot == 1:
        x_val = 'x{0} = {1}'.format(atom_idx, xyzP[0])
        y_val = '  y{0} = {1}'.format(atom_idx, xyzP[1])
        z_val = '  z{0} = {1}'.format(atom_idx, xyzP[2])
        pivot_xyz_string = (x_val + y_val + z_val)
    else:
        x_val1 = 'x{0} = {1} + d{2}*sin(p{0})*cos(t{0})'.format(
            atom_idx, xyzP[0], d_idx)
        y_val1 = '  y{0} = {1} + d{2}*sin(p{0})*sin(t{0})'.format(
            atom_idx, xyzP[1], d_idx)
        z_val1 = '  z{0} = {1} + d{2}*cos(p{0})'.format(
            atom_idx, xyzP[2], d_idx)
        x_val2 = 'x{0} = {1} - d{2}*sin(p{0})*cos(t{0})'.format(
            atom_idx+1, xyzP[0], d_idx)
        y_val2 = '  y{0} = {1} - d{2}*sin(p{0})*sin(t{0})'.format(
            atom_idx+1, xyzP[1], d_idx)
        z_val2 = '  z{0} = {1} + d{2}*cos(p{0})'.format(
            atom_idx+1, xyzP[2], d_idx)
        pivot_xyz_string = (x_val1 + y_val1 + z_val1 + '\n' +
                            x_val2 + y_val2 + z_val2)

    return pivot_xyz_string
