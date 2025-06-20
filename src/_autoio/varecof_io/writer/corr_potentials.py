"""
Writes the fortran files needed for the correction potential
"""

import os
import subprocess
from ioformat import build_mako_str
from varecof_io.writer import _format as vrcformat


# OBTAIN THE PATH TO THE DIRECTORY CONTAINING THE TEMPLATES #
SRC_PATH = os.path.dirname(os.path.realpath(__file__))
TEMPLATE_PATH = os.path.join(SRC_PATH, 'templates')


def species(rvalues, potentials, aidx, bidx,
            dist_restrict_idxs=(), pot_labels=(), species_name='mol'):
    """ Writes the string for a Fortran source file containing information
        used to build the correction potential for a species.

        :param rvalues: Intermolecular distance to scan over
        :type rvalues: list(float)
        :param potentials: list of potentials calculates along rvalues
        :type potentials: list(list(float))
        :param aidx: index for first radical atom that form bond
        :type aidx: int
        :param bidx: index for second radical atom that form bond
        :type bidx: int
        :param pot_labels: names of each of the potentials in .f file
        :type pot_labels: list(str)
        :param species_name: name given to mol_corr.f file
        :type species_name: str
        :rtype: str
    """

    npot = len(potentials)
    npot_terms = len(potentials[0])
    asym, bsym = 'A', 'B'

    assert npot > 0
    assert all(len(potential) == npot_terms for potential in potentials)

    # Put some comment lines giving a description of the correction potentials
    pot_labels_str = ''
    if pot_labels:
        for i, label in enumerate(pot_labels):
            pot_labels_str += f'c     dv{str(i+1)} = {label} correction\n'
    pot_labels_str = pot_labels_str.rstrip()

    # Strings to initialize the potential variables in the Fortran subroutine
    npot = len(potentials)
    npot_terms = len(potentials[0])
    dv_defs = ''
    for i in range(npot):
        dv_defs += f'dv{str(i+1)}({npot_terms}),'
    dv_defs = dv_defs[:-1]

    # Definitions of all of all the correction potential distances
    rvals = ''
    for i, rval in enumerate(rvalues):
        rvals += f'      data rinp({str(i+1)}) / {rval:.3f} /\n'
    rvals = rvals.rstrip()
    rmin = min(rvalues)
    rmax = max(rvalues)

    # Definitions of all of all the correction potential energies
    dv_vals = ''
    for i, potential in enumerate(potentials):
        for j, term in enumerate(potential):
            dv_vals += f'      data dv{str(i+1)}({str(j+1)}) / {term:.3f} /\n'
    dv_vals = dv_vals.rstrip()

    # Build principal distance string
    bond_distance_string = vrcformat.format_corrpot_dist_string(
        aidx, bidx, asym, bsym)

    # Build distance restriction strings
    restrict_distance_strings = ''
    for i, idxs in enumerate(dist_restrict_idxs):
        [idx1, idx2] = idxs
        sym1, sym2 = chr(67+2*i), chr(68+2*i)
        restrict_distance_strings += vrcformat.format_corrpot_dist_string(
            idx1, idx2, sym1, sym2)
        restrict_distance_strings += '\n'
        restrict_distance_strings += vrcformat.format_restrict_dist_string(
            sym1, sym2, species_name)
        restrict_distance_strings += '\n'

    # Build the delmlt string
    delmlt_string = vrcformat.format_delmlt_string(asym, bsym)

    # Build the spline fitting strings
    spline_strings = vrcformat.format_spline_strings(
        npot, asym, bsym, species_name)

    # Create dictionary to fill template
    corr_keys = {
        'species_name': species_name,
        'pot_labels': pot_labels_str,
        'npot': npot,
        'npot_terms': npot_terms,
        'dv_defs': dv_defs,
        'rvals': rvals,
        'dv_vals': dv_vals,
        'rmin': rmin,
        'rmax': rmax,
        'bond_distance_string': bond_distance_string,
        'restrict_distance_strings': restrict_distance_strings,
        'delmlt_string': delmlt_string,
        'spline_strings': spline_strings
    }

    return build_mako_str(
        template_file_name='species_corr.mako',
        template_src_path=TEMPLATE_PATH,
        template_keys=corr_keys)


def dummy():
    """ Writes string for the dummy correction potential Fortran file.

        :rtype: string
    """

    return build_mako_str(
        template_file_name='dummy_corr.mako',
        template_src_path=TEMPLATE_PATH,
        template_keys={})


def auxiliary():
    """ Writes string for the potential auxiliary functions Fortran file.

        :rtype: string
    """

    return build_mako_str(
        template_file_name='pot_aux.mako',
        template_src_path=TEMPLATE_PATH,
        template_keys={})


def makefile(fortran_compiler, pot_file_names=()):
    """ Writes string for a makefile to compile correction potentials.

        :param fortran_compiler: name of compiler to build potentials
        :type fortran_compiler: str
        :param pot_file_names: names of files with various potentials
        :type: pot_file_names: list(str)
        :return: string for the makefile
        :rtype: string
    """

    # Set species name
    corr_potential_names = ''
    if pot_file_names:
        for potential in pot_file_names:
            corr_potential_names += f'{potential}_corr.f '

    # Create dictionary to fill template
    make_keys = {
        'fc': fortran_compiler,
        'corr_potential_names': corr_potential_names
    }

    return build_mako_str(
        template_file_name='makefile.mako',
        template_src_path=TEMPLATE_PATH,
        template_keys=make_keys)


def compile_corr_pot(make_path):
    """ Compiles the correction potential using make.

        :param make_path: path to the makefile and correction potential src
        :type make_path: str
    """

    subprocess.check_call(
        ['make'], cwd=make_path,
        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
