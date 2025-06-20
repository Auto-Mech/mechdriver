"""
Writes a section for an atom or molecule for the MESS input file.
Generally placed within a MESS "species", "well",
"bimolecular", or "barrier" section.

Takes in data from other MESS writer functions
"""

import os
from ioformat import build_mako_str
from mess_io.writer import _format as messformat


# OBTAIN THE PATH TO THE DIRECTORY CONTAINING THE TEMPLATES #
SRC_PATH = os.path.dirname(os.path.realpath(__file__))
TEMPLATE_PATH = os.path.join(SRC_PATH, 'templates')
SPECIES_PATH = os.path.join(TEMPLATE_PATH, 'species')


def atom(mass, elec_levels):
    """ Writes the string that defines the `Species` section
        for an atom for a MESS input file by
        formatting input information into strings and filling Mako template.

        :param mass: mass of the atom (of desired isotope)
        :type mass: int
        :param elec_levels: energy and degeneracy of atom's electronic states
        :type elec_levels: list(float)
        :rtype: str
    """

    # Build a formatted elec levels string
    nlevels, levels = messformat.elec_levels_format(elec_levels)

    # Create dictionary to fill template
    atom_keys = {
        'mass': mass,
        'nlevels': nlevels,
        'levels': levels
    }

    return build_mako_str(
        template_file_name='atom.mako',
        template_src_path=SPECIES_PATH,
        template_keys=atom_keys)


def molecule(core, elec_levels,
             freqs=(), hind_rot='', xmat=(),
             rovib_coups=(), rot_dists=(), inf_intens=(),
             freq_scale_factor=None,
             use_harmfreqs_key=False):
    """ Writes the string that defines the `Species` section
        for a molecule for a MESS input file by
        formatting input information into strings and filling Mako template.

        :param core: `Core` section string in MESS format
        :type core: str
        :param elec_levels: energy and degeneracy of atom's electronic states
        :type elec_levels: list(float)
        :param freqs: Harmonic/Anharmonic vibrational frequencies
        :type freqs: list(float)
        :param elec_levels: energy and degeneracy of atom's electronic states
        :param hind_rot: string of MESS-format `Rotor` sections for all rotors
        :type hind_rot: str
        :param xmat: anharmonicity matrix (cm-1)
        :type xmat: list(list(float))
        :param inf_intens: Harmonic Infrared Intensities for vibrational modes
        :type inf_intens: list(float)
        :param rovib_coups: rovibrational coupling matrix
        :type rovib_coups: numpy.ndarray
        :param rot_dists: rotational distortion constants: [['aaa'], [val]]
        :type rot_dists: list(list(str), list(float))
        :rtype: str
    """

    # Add in infrared intensities at some point

    # Build formatted frequencies, infrared intensities and elec levels string
    nfreqs, freqs = messformat.freqs_format(freqs)
    nintens, intens = messformat.intensities_format(inf_intens)
    nlevels, levels = messformat.elec_levels_format(elec_levels)
    interp_emax = None
    # Format the rovib couplings and rotational distortions if needed
    if rovib_coups:
        rovib_coups = messformat.format_rovib_coups(rovib_coups)
        interp_emax = 500
    else:
        rovib_coups = ''
    if rot_dists:
        rot_dists = messformat.format_rot_dist_consts(rot_dists)
    else:
        rot_dists = ''
    if xmat:
        anharm = messformat.format_xmat(xmat)
        interp_emax = 500
    else:
        anharm = ''

    # Indent various strings string if needed
    if hind_rot != '':
        hind_rot = messformat.indent(hind_rot, 2)

    # Create dictionary to fill template
    molec_keys = {
        'core': core,
        'nlevels': nlevels,
        'levels': levels,
        'nfreqs': nfreqs,
        'freqs': freqs,
        'hind_rot': hind_rot,
        'anharm': anharm,
        'nintens': nintens,
        'intens': intens,
        'rovib_coups': rovib_coups,
        'rot_dists': rot_dists,
        'interp_emax': interp_emax,
        'freq_scale_factor': freq_scale_factor,
        'use_harmfreqs_key': use_harmfreqs_key
    }

    return build_mako_str(
        template_file_name='molecule.mako',
        template_src_path=SPECIES_PATH,
        template_keys=molec_keys)
