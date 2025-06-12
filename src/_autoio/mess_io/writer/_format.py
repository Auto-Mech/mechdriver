"""
  Additional functions for formatting information for MESS strings
"""

import numpy
from ioformat import indent


# Format MESS labels
def mess_label_format(spc_label, aux_id_label=None, calc_dens=False):
    """ Format the full MESS label line that proceeds a reaction channel
        section keyword such as `Species`, `Well`, `Fragment`, etc.
        Function will add an additional id label (e.g., SMILES) and a
        density of state keyword, if requested.
    """

    lbl = spc_label
    if calc_dens:
        lbl = f'{lbl} density'
    if aux_id_label is not None:
        lbl = f'{lbl}   ! {aux_id_label}'

    return lbl


# Format various pieces of data into strings for MESS input files
def zero_energy_format(zero_ene):
    """ Formats the zero point energy into a string that
        is appropriate for a MESS input file.

        :param zero_ene: zero point energy value
        :type zero_ene: float
        :return zero_ene_str: MESS-format string containing energy
        :rtype string
    """
    return f'ZeroEnergy[kcal/mol]      {zero_ene:<8.2f}'


def elec_levels_format(elec_levels):
    """ Formats the list of electronic energy levels into a string that
        is appropriate for a MESS input file.

        :param elec_levels: levels, given as [[energy, degeneracy], ...]
        :type elec_levels: list(list(float))
        :return elec_levels_str: MESS-format string containing levels
        :rtype string
    """

    # Get the number of elec levles
    nlevels = len(elec_levels)

    # Build elec levels string
    elec_levels_str = ''
    for i, level in enumerate(elec_levels):
        elec_levels_str += '  '.join(map(str, level))
        if (i+1) != len(elec_levels):
            elec_levels_str += '\n'

    # Indent the lines
    elec_levels_str = indent(elec_levels_str, 4)

    return nlevels, elec_levels_str


def geometry_format(geo, indent_lines=True):
    """ Formats the geometry of a species into a string that
        is appropriate for a MESS input file.

        :param geo: geometry of a species
        :param indent_lines: indent the lines of the geometry
        :return natoms: number of atoms in the geometry
        :rtype int
        :return geo_str: MESS-format string containing geometry
        :rtype string
    """

    # Get the number of atoms
    natoms = len(geo)

    # Build geom string; converting the coordinates to angstrom
    gstr = ''
    for (symb, xyz) in geo:
        xyzc = tuple(val*0.529177 for val in xyz)
        gstr += f'{symb:<4s}{xyzc[0]:>14.5f}{xyzc[1]:>14.5f}{xyzc[2]:>14.5f}\n'

    # Remove final newline character and indent the lines
    if indent_lines:
        gstr = indent(gstr.rstrip(), 4)
    else:
        gstr = gstr.rstrip()

    return natoms, gstr


def mc_geometry_format(geo):
    """ Formats the geometry of a species into a string that
        is appropriate for a MESS MonteCarlo auxiliary data file.

        :param geo: geometry of a species
        :return natoms: number of atoms in the geometry
        :rtype int
        :return geo_str: MESS-format string containing geometry
        :rtype string
    """
    natoms, gstr = geometry_format(geo, indent_lines=False)
    return f'{natoms}\n{gstr}'


def freqs_format(freqs):
    """ Formats the vibrational frequencies of a species into a string that
        is appropriate for a MESS input file.

        :param freqs: vibrational frequencies of species
        :type freqs: list(float)
        :return nfreqs: number of frequences for the species
        :rtype int
        :return freq_str: MESS-format string containing frequencies
        :rtype string
    """

    # Get the number of freqs
    nfreqs = len(freqs)

    # Build freqs string
    freq_str = ''
    for i, freq in enumerate(freqs):
        if ((i+1) % 6) == 0 and (i+1) != len(freqs):
            freq_str += f'{int(freq):<8.0f}\n'
        else:
            freq_str += f'{int(freq):<8.0f}'

    # Indent the lines
    freq_str = indent(freq_str, 4)

    return nfreqs, freq_str


def intensities_format(intens):
    """ Formats the vibrational intenuencies of a species into a string that
        is appropriate for a MESS input file.

        :param intens: harmonic infrared intensities of species
        :type intens: list(float)
        :return nintens: number of harmonic infrared intensities for species
        :rtype int
        :return inten_str: MESS-format string containing infrared intensities
        :rtype string
    """

    # Get the number of intens
    nintens = len(intens)

    # Build intens string
    inten_str = ''
    for i, inten in enumerate(intens):
        if ((i+1) % 6) == 0 and (i+1) != len(intens):
            inten_str += f'{int(inten):<8.1f}\n'
        else:
            inten_str += f'{int(inten):<8.1f}'

    # Indent the lines
    inten_str = indent(inten_str, 4)

    return nintens, inten_str


def format_rotor_key_defs(rotor_keyword_vals):
    """ Formats strings that contain the 'Group', 'Axis', and 'Symmetry'
        keywords and values that are used to define hindered rotors and
        internal rotors in MESS input files.

        :param rotor_keyword_vals: values for the for some rotor keyword
        :type: rotor_keyword_vals: list(int)
        :return rotor_keyword_str: MESS-format string containing values
        :rtype str
    """

    # Build string containing the values of each keyword
    rotor_keyword_str = ''
    for vals in rotor_keyword_vals:
        rotor_keyword_str += f'{vals+1:<4d}'

    return rotor_keyword_str


def format_rotor_potential(potential):
    """ Formats the potential energy surface along a rotor into a string
        used to define hindered rotors and internal rotors in MESS input files.

        :param potential: value of the potential along torsion (kcal.mol-1)
        :type potential: list(float)
        :return npotential: number of values in the potential
        :rtype int
        :return potential_str: values of potential in a MESS-format string
        :rtype str
    """
    assert not any(map(numpy.isnan, potential.values())), (
        f"potential has nans: {potential}"
    )

    # Get the number of the terms in the potential
    npot = len(potential)

    # Build potentials string
    coord_str, ene_str = '', ''
    for i, (coord, energy) in enumerate(potential.items()):
        if ((i+1) % 6) == 0 and (i+1) != npot:
            coord_str += f'{coord[0]:<8.2f}\n'
            ene_str += f'{energy:<8.4f}\n'
        else:
            coord_str += f'{coord[0]:<8.2f}'
            ene_str += f'{energy:<8.4f}'

    # Indent the lines
    coord_str = indent(coord_str, 4)
    ene_str = indent(ene_str, 4)

    return npot, coord_str, ene_str


def format_rovib_coups(rovib_coups):
    """ Formats the matrix of rovibrational coupling terms for a species
        into a string appropriate for a MESS input file.

        :param rovib_coups: rovibrational coupling matrix
        :type rovib_coups: numpy.ndarray
        :return rovib_coups_str: values of potential in a MESS-format string
        :rtype str
    """

    # Join the values into a string
    rovib_coups_lst = []
    for line in rovib_coups:
        rovib_coups_lst.append(
            indent('\t'.join(str(val) for val in line), 4))

    rovib_coups_str = '\n'.join(rovib_coups_lst)

    return rovib_coups_str


def format_rot_dist_consts(rot_dists):
    """ Formats the list of rotational distortion constants
        into a string appropriate for a MESS input file.

        :param rot_dists: rotational distortion constants: [['aaa'], [val]]
        :type rot_dists: list(list(str), list(float))
        :return rot_dists_str: values of potential in a MESS-format string
        :rtype string
    """

    # Build rotational dists string
    rot_dists_str = ''
    for i, const in enumerate(rot_dists):
        rot_dists_str += '  '.join(map(str, const))
        if (i+1) != len(rot_dists):
            rot_dists_str += '\n'

    # Indent the lines
    rot_dists_str = indent(rot_dists_str, 4)

    return rot_dists_str


def format_xmat(xmat):
    """ Formats the anharmonicity (X) matrix for a species
        into a string appropriate for a MESS input file.

        :param xmat: anharmonicity matrix (cm-1)
        :type xmat: list(list(float))
        :return xmat_str: anharmonicity matrix in a MESS-format string
        :rtype string
    """

    xmat = numpy.array(xmat)

    # Loop over the rows of the anharm numpy array
    xmat_str = ''
    for i in range(xmat.shape[0]):
        xmat_str += ' '.join(
            [f'{val:>12.5f}' for val in list(xmat[i, :i+1])]
            # if val != 0.0]
        )
        if (i+1) != xmat.shape[0]:
            xmat_str += '\n'

    # Indent the lines
    xmat_str = indent(xmat_str, 2)

    return xmat_str


def molec_spec_format(geo):
    """ Parses out the atom labels of a Cartesian geometry and
        formats them into a string appropriate for definining
        molecular species for Monte Carlo calculations in MESS.

        :param geo: geometry
        :type geo: list
        :return atom_lst_str
        :rtype: string
    """

    # Build geom string; converting the coordinates to angstrom
    atom_lst_str = ''
    for (symb, _) in geo:
        atom_lst_str += f'{symb:s} '

    # Remove final newline character
    atom_lst_str = atom_lst_str.rstrip()

    # Indent the lines
    atom_lst_str = indent(atom_lst_str, 6)

    return atom_lst_str


def format_flux_mode_indices(atom_idxs):
    """ Formates the atom indices into a string that is used
        to define the fluxional (torsional) modes of a
        molecular species for Monte Carlo calculations in MESS.

        :param atom_idxs: idxs of atoms involved in fluxional mode
        :type atom_idxs: list(int)
        :return flux_mode_idx_str: formatted string of indices
        :rtype: string
    """

    # Build string containing the values of each keyword
    flux_mode_idx_str = ''
    for vals in atom_idxs:
        flux_mode_idx_str += f'{vals+1:<4d}'

    return flux_mode_idx_str


def format_ped_species(ped_spc_lst):
    """ Format the names of species to provided to the global
        PEDSpecies keyword of the input file.

        :param: ped_spc_lst: species to obtain PEDs for
        :type: tuple(str)
        :rtype: str
    """
    return '   '.join(ped_spc_lst)


def format_hot_enes(hot_enes_dct):
    """ Format the list hot energies for one or more species
        into a string appropriate for the HotEnergies keyword
        used at the top of the MESS keyword section

        :param hot_enes_dct: hot energies for each species in kcal/mol
        :type hot_enes_dct: dict[str: tuple(float)]
        :rtype: (int, str)
    """

    ene_str = ''
    n_enes = 0
    for spc, ene_lst in hot_enes_dct.items():
        _str = ':'.join((f'{ene:.1f}' for ene in ene_lst))
        ene_str += f'{spc:5s} {_str}\n'
        n_enes += 1

    return n_enes, ene_str.rstrip()


# Helpful checker to set MESS string writing
def is_atom_in_str(spc_str):
    """ Checks a MESS-formatted species data string to see
        if the species is, or contains, an Atom species definition.

        :param: spc_str: MESS species string
        :type spc_str: str
        rtype: bool
    """
    return bool('Atom' in spc_str)
