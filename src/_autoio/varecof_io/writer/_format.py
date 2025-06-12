"""
 Utility functions for formatting
"""

import os
import subprocess
from automol import geom
from phydat import phycon


# Obtain the path to a convert struct executable
SRC_PATH = os.path.dirname(os.path.realpath(__file__))


# Utility functions for the structure.inp writer
def determine_struct_type(geo):
    """ determines the linear string
    """

    # Remove dummy atoms
    geo = [coords for coords in geo
           if coords[0] != 'X']

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
    symbs = geom.symbols(geo)
    coords = geom.coordinates(geo)
    masses = tuple(round(mass) for mass in geom.masses(geo))

    # Build a string with the formatted coordinates string
    if geom.is_atom(geo):
        geo_str = f'{symbs[0]:<4s}{masses[0]:<6d}'
    else:
        geo_str = f'{str(natoms)} \n'
        for symb, mass, xyzs in zip(symbs, masses, coords):
            _xyzs = [x * phycon.BOHR2ANG for x in xyzs]
            xyzs_str = f'{_xyzs[0]:>14.8f}{_xyzs[1]:>14.8f}{_xyzs[2]:>14.8f}'
            geo_str += f'{symb:<4s}{mass:<6.0f}{xyzs_str}\n'
        # Remove final newline character from the string
        geo_str = geo_str.rstrip()

    return natoms, geo_str


# Utility functions for the tst.inp writer
def format_grids_string(grid, name, units):
    """ format the string using the grids for
        energy and angular momentum for tst.inp file
    """
    grid_str = (
        f'{name}_grid{grid[0]:>8d}{grid[1]:>9d}'
        f'{grid[2]:>11.2f}{grid[3]:>7d}'
    )
    grid_str += f'     {units:<8s}# {name} grid'

    return grid_str


def format_faces_string(faces):
    """ format faces keywords
    """
    faces_str = ' '.join((str(val) for val in faces))

    return faces_str


# Utility functions for the divsur.inp writer
def format_values_string(coord, values, conv_factor=1.0):
    """ format the values string for the divsur.inp file
    """
    if values:
        values = ', '.join(f'{val*conv_factor:.3f}'
                           for val in values)
        values_string = f'{coord} = ({values})'
    else:
        values_string = ''

    return values_string


def format_pivot_xyz_string(idx, npivot, xyzp, phi_dependence=False):
    """ format the pivot point xyz
    """

    assert npivot in (1, 2)

    atom_idx = idx
    if idx == 1:
        d_idx = 1
        t_idx = 1
    else:
        d_idx = 2
        t_idx = 2

    if npivot == 1:
        x_val = f'x{atom_idx} = {xyzp[0]:.3f}'
        y_val = f'  y{atom_idx} = {xyzp[1]:.3f}'
        z_val = f'  z{atom_idx} = {xyzp[2]:.3f}'
        pivot_xyz_string = (x_val + y_val + z_val)
    elif npivot > 1 and not phi_dependence:
        x_val1 = f'x{atom_idx} = {xyzp[0]:.3f} + d{d_idx}*cos(t{t_idx})'
        y_val1 = f'  y{atom_idx} = {xyzp[1]:.3f} + d{d_idx}*sin(t{t_idx})'
        z_val1 = f'  z{atom_idx} = 0.000'
        x_val2 = f'x{atom_idx+1} = {xyzp[0]:.3f} - d{d_idx}*cos(t{t_idx})'
        y_val2 = f'  y{atom_idx+1} = {xyzp[1]:.3f} - d{d_idx}*sin(t{t_idx})'
        z_val2 = f'  z{atom_idx+1} = 0.000'
        pivot_xyz_string = (x_val1 + y_val1 + z_val1 + '\n' +
                            x_val2 + y_val2 + z_val2)
    else:
        raise NotImplementedError
        # # Not sure if this implementation is any good
        # x_val1 = 'x{0} = {1:.0f} + d{2}*sin(p{0})*cos(t{0})'.format(
        #     atom_idx, xyzp[0], d_idx)
        # y_val1 = '  y{0} = {1:.0f} + d{2}*sin(p{0})*sin(t{0})'.format(
        #     atom_idx, xyzp[1], d_idx)
        # z_val1 = '  z{0} = {1:.0f} + d{2}*cos(p{0})'.format(
        #     atom_idx, xyzp[2], d_idx)
        # x_val2 = 'x{0} = {1:.0f} - d{2}*sin(p{0})*cos(t{0})'.format(
        #     atom_idx+1, xyzp[0], d_idx)
        # y_val2 = '  y{0} = {1:.0f} - d{2}*sin(p{0})*sin(t{0})'.format(
        #     atom_idx+1, xyzp[1], d_idx)
        # z_val2 = '  z{0} = {1:.0f} + d{2}*cos(p{0})'.format(
        #     atom_idx+1, xyzp[2], d_idx)
        # pivot_xyz_string = (x_val1 + y_val1 + z_val1 + '\n' +
        #                     x_val2 + y_val2 + z_val2)

    return pivot_xyz_string


# Utility functions for the species_corr.f correction potential writer
def format_corrpot_dist_string(aidx, bidx, asym, bsym):
    """ set distance string for two atoms for the file
    """

    lasym, lbsym = asym.lower(), bsym.lower()

    return (
        f"      n{lasym} = {aidx}\n" +
        f"      n{lbsym} = {bidx}\n" +
        f"      r{asym}{bsym} = dsqrt( (x(1,n{lbsym})-x(1,n{lasym}))**2 +\n" +
        f"     x             (x(2,n{lbsym})-x(2,n{lasym}))**2 +\n" +
        f"     x             (x(3,n{lbsym})-x(3,n{lasym}))**2)\n"  # +
        # f"      r{asym}{bsym} = r{asym}{bsym}*0.52917"
    )


def format_delmlt_string(asym, bsym):
    """ set distance string for two atoms for the file
    """

    return (
        "      delmlt = 1.0d0\n" +
        f"      if(r{asym}{bsym}.le.r{asym}{bsym}min) r{asym}{bsym} = " +
        f"r{asym}{bsym}min\n" +
        f"      if(r{asym}{bsym}.ge.r{asym}{bsym}max) then\n" +
        f"        delmlt = exp(-2.d0*(r{asym}{bsym}-r{asym}{bsym}max))\n" +
        f"        r{asym}{bsym}=r{asym}{bsym}max\n" +
        "      endif"
    )


def format_restrict_dist_string(sym1, sym2, name):
    """ build string that has the distance comparison
    """

    return (
        f"      if (r{sym1}{sym2}.lt.rAB) then\n" +
        f"        {name}_corr = 100.0\n" +
        "        return\n" +
        "      endif"
    )


def format_spline_strings(npot, sym1, sym2, species_name):
    """ spline fitting strings
    """

    corr_name = species_name+'_corr'

    spline_str = ''
    for i in range(npot):
        if i == 0:
            ifstr = 'if'
        else:
            ifstr = 'else if'
        spline_str += f'      {ifstr} (ipot.eq.{str(i+1)}) then\n'
        spline_str += (
            f'        call spline(rinp,dv{str(i+1)},nrin,dvp1,dvpn,dv20)\n' +
            f'        call splint(rinp,dv{str(i+1)},dv20,nrin,r{sym1}{sym2}' +
            f',{corr_name})\n'
        )
    spline_str += '      endif'

    return spline_str


def divsur_frame_geom_script():
    """ Run the VaReCoF utility script to calculate the fragment
        geometries contained in the divsur.out file
        (only requires the divsur.inp file)
    """

    conv_cmd = [
        '/lcrc/project/CMRP/amech/VaReCoF/build/convert_struct',
        'divsur.inp'
    ]
    with open(os.devnull, mode='w', encoding='utf-8') as nullfile:
        subprocess.check_call(
            conv_cmd, stdout=nullfile, stderr=nullfile)
