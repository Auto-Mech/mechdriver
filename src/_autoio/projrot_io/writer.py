"""
  Functions to write ProjRot input files
"""

import os
from automol.util import tensor
from ioformat import build_mako_str
from ioformat import remove_trail_whitespace
from phydat import phycon
from projrot_io import util


# OBTAIN THE PATH TO THE DIRECTORY CONTAINING THE TEMPLATES #
SRC_PATH = os.path.dirname(os.path.realpath(__file__))
TEMPLATE_PATH = os.path.join(SRC_PATH, 'templates')


def rpht_input(geos, grads, hessians,
               saddle_idx=1,
               rotors_str='',
               coord_proj='cartesian',
               proj_rxn_coord=False):
    """ Writes a string for the input file for ProjRot.

        :param geos: geometry for single species or along a reaction path
        :type geos: list(list(float))
        :param grads: gradient for single species or along a reaction path
        :type grads: list(list(float))
        :param hessians: Hessian for single species or along a reaction path
        :type hessians: list(list(float))
        :param saddle_idx: idx denoting the saddle point along a reaction path
        :type saddle_idx: int
        :param rotors_str: ProjRot-format string with all the rotor definitions
        :type rotors_str: str
        :param coord_proj: choice of coordinate system to perform projections
        :type coord_proj: str
        :param proj_rxn_coord: whether to project out reaction coordinate
        :type proj_rxn_coord: bool
        :rtype: str
    """

    # Format the molecule info
    data_str = util.write_data_str(geos, grads, hessians)
    natoms = len(geos[0])
    # nsteps = len(geos)
    nrotors = rotors_str.count('pivotA')

    nsteps = len(geos)
    assert nsteps == len(hessians)
    if len(grads) != 0:
        assert len(grads) == nsteps
    assert coord_proj in ('cartesian', 'internal')

    # Create a fill value dictionary
    rpht_keys = {
        'natoms': natoms,
        'nsteps': nsteps,
        'saddle_idx': saddle_idx,
        'coord_proj': coord_proj,
        'proj_rxn_coord': proj_rxn_coord,
        'nrotors': nrotors,
        'rotors_str': rotors_str,
        'data_str': data_str
    }

    return build_mako_str(
        template_file_name='rpht_input.mako',
        template_src_path=TEMPLATE_PATH,
        template_keys=rpht_keys)


def rpht_path_coord_en(coords, energies, bnd1=(), bnd2=()):
    """ Writes a string for the auxiliary input file used for
        ProjRot SCT calculations, which contains information
        along the reaction path.

        :params coords: Values of reaction coordinate along reaction path
        :type coords: list(float)
        :params energies: Energies along reaction path
        :type energies: list(float)
        :params bnd1: values of bond in reactant side
        :type bnd1: list(float)
        :params bnd2: values of bond in product side
        :type bnd2: list(float)
        :rtype: str
    """

    nsteps = len(coords)

    # Check bnd1 and bnd2 lists and build the corresponding string lists
    assert len(bnd1) == len(bnd2)
    assert (bool(bnd1 and bnd2) or bool(not bnd1 and not bnd2))
    bnd_strs = []
    if bnd1 and bnd2:
        for bd1, bd2 in zip(bnd1, bnd2):
            bnd_strs.append(f'{bd1:>10.5f}{bd2:>10.5f}')
    else:
        for _ in range(len(coords)):
            bnd_strs.append(f'{1.0:>10.5f}{1.0:>10.5f}')

    # Check that all the lists are not empty and have the same length
    assert all(lst for lst in (coords, energies, bnd_strs))
    assert all(len(lst) == nsteps for lst in (coords, energies, bnd_strs))

    path_str = (
        f'{"Point":>7s}'
        f'{"Coordinate":>12s}'
        f'{"Energy":>10s}'
        f'{"Bond1":>10s}'
        f'{"Bond2":>10s}\n'
    )
    for i, (crd, ene, bnd_str) in enumerate(zip(coords, energies, bnd_strs)):
        path_str += f'{i+1:>7d}{crd:>12.5f}{ene:>10.5f}{bnd_str:>20s}'
        if i+1 != nsteps:
            path_str += '\n'

    return remove_trail_whitespace(path_str)


def rotors(axis, group):
    """ Write the sections that defines the rotors section

        :param group: idxs for the atoms of one of the rotational groups
        :type group: list(int)
        :param axis: idxs for the atoms that make up the rotational axis
        :type axis: list(int)
        :rtype: str
    """

    # Set up the keywords
    pivota, pivotb = axis
    atomsintopa = len(group)

    pivota += 1
    pivotb += 1
    topaatoms = '  '.join([str(val+1) for val in group])

    # Build the rotors_str
    rotors_str = f'\n{"pivotA":<32s}{pivota:<4d}\n'
    rotors_str += f'{"pivotB":<32s}{pivotb:<4d}\n'
    rotors_str += f'{"atomsintopA":<32s}{atomsintopa:<4d}\n'
    rotors_str += f'{"topAatoms":<32s}{topaatoms}\n'

    return util.remove_trail_whitespace(rotors_str)


def projection_distance_aux(dist_cutoff_dct=None):
    """ Write the auxiliary file for defining atom-atom distance cutoffs used
        to constuct parts of the molecule for rotors.

        :param dist_cutoff_dct: atom-atom cutoff distances (in Bohr)
        :type dist_cutoff_dct: dict[(str, str): float]
        :rtype: str
    """

    # Set dictionary of default cutoffs and update with the input dct
    # These distances are already in Angstrom since that is unit for file
    dist_dct = {
        ('C', 'H'): 1.50,
        ('H', 'C'): 1.25,
        ('C', 'C'): 2.40,
        ('C', 'N'): 1.80,
        ('C', 'O'): 1.60,
        ('O', 'H'): 1.60,
        ('H', 'O'): 1.50,
        ('O', 'O'): 1.60,
        ('S', 'H'): 1.50,
        ('N', 'H'): 1.30,
        ('H', 'H'): 1.20
    }
    if dist_cutoff_dct is not None:
        for key, val in dist_cutoff_dct.items():
            dist_cutoff_dct[key] = val * phycon.BOHR2ANG
        dist_dct.update(dist_cutoff_dct)

    # Set the header
    dist_str = (
        'projection distances\n'
        '\n'
        'distances to be read as Atom1 Atom2 distance\n'
        '\n'
    )

    # Write the distance strings
    for (atm1, atm2), dist in dist_dct.items():
        dist_str += f'{atm1:s} {atm2:s} {dist:.3f}\n'

    return dist_str


def bmatrix(bmat):  # pragma: no cover
    """ Write the B-Matrix to a ProjRot style input
    """

    nd1, nd2, nd3 = bmat.shape
    bmat_str = f'{nd1:>12d}{nd2*nd3:>12d}\n'
    bmat_str += tensor.string_submat_3d(bmat)

    return bmat_str


def cmatrix(cmat):  # pragma: no cover
    """ Write the C-Matrix to a ProjRot style input
    """

    _, nd2, nd3, nd4 = cmat.shape
    cmat_str = f'{nd2:>12d}{nd3*nd4:>12d}\n'
    cmat_str += tensor.string_submat_4d(cmat)

    return cmat_str
