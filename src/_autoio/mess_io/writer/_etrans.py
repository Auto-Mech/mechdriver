"""
Writes the energy transfer section of a MESS input file
"""

import os
from phydat import phycon
from ioformat import build_mako_str


# OBTAIN THE PATH TO THE DIRECTORY CONTAINING THE TEMPLATES #
SRC_PATH = os.path.dirname(os.path.realpath(__file__))
TEMPLATE_PATH = os.path.join(SRC_PATH, 'templates')
SECTION_PATH = os.path.join(TEMPLATE_PATH, 'sections')
ENE_TRANS_PATH = os.path.join(SECTION_PATH, 'energy_transfer')


def energy_down(exp_factor, exp_power, exp_cutoff):
    """ Writes the energy transfer section of the MESS input file by
        formatting input information into strings a filling Mako template.

        :param exp_factor: 300 K energy-down value for collision model (cm-1)
        :type exp_factor: float
        :param exp_power: n, for [energy-down * (T/300K)^n] for collision model
        :type exp_power: float
        :param exp_cutoff: cutoff for assuming transition probability is zero
        :type exp_cutoff: float
        :rtype: string
    """

    # Put the values into a string
    exp_factor_str = f'{exp_factor:<10.3f}'
    exp_power_str = f'{exp_power:<10.3f}'
    exp_cutoff_str = f'{exp_cutoff:<10.3f}'

    # Create dictionary to fill template
    etrans_keys = {
        'exp_factor': exp_factor_str,
        'exp_power': exp_power_str,
        'exp_cutoff': exp_cutoff_str
    }

    return build_mako_str(
        template_file_name='edown.mako',
        template_src_path=ENE_TRANS_PATH,
        template_keys=etrans_keys)


def collision_frequency(eps1, eps2, sig1, sig2, mass1, mass2):
    """ Writes the energy transfer section of the MESS input file by
        formatting input information into strings a filling Mako template.

        :param eps1: A+A Lennard-Jones epsilon parameter of spc 1 (cm-1)
        :type eps1: float
        :param eps2: A+A Lennard-Jones epsilon parameter of spc 2 (cm-1)
        :type eps2: float
        :param sig1: A+A Lennard-Jones sigma parameter of spc 1 (Angstrom)
        :type sig1: float
        :param sig2: A+A Lennard-Jones sigma parameter of spc 2 (Angstrom)
        :type sig2: float
        :param mass1: mass of Species 1 (amu)
        :type mass1: float
        :param mass2: mass of Species 2 (amu)
        :type mass2: float
        :rtype: string
    """

    # Convert LJ units and format them into string
    eps_str = f'{eps1*phycon.EH2WAVEN:<10.3f} {eps2*phycon.EH2WAVEN:<10.3f}'
    sig_str = f'{sig1*phycon.BOHR2ANG:<10.3f} {sig2*phycon.BOHR2ANG:<10.3f}'
    mass_str = f'{mass1:<10.3f} {mass2:<10.3f}'

    # Create dictionary to fill template
    etrans_keys = {
        'epsilons': eps_str,
        'sigmas': sig_str,
        'masses': mass_str
    }

    return build_mako_str(
        template_file_name='collid_freq.mako',
        template_src_path=ENE_TRANS_PATH,
        template_keys=etrans_keys)
