""" Write THERMP input file information
"""

import os
from ioformat import build_mako_str
from phydat import phycon
from thermp_io import util


# OBTAIN THE PATH TO THE DIRECTORY CONTAINING THE TEMPLATES #
SRC_PATH = os.path.dirname(os.path.realpath(__file__))
TEMPLATE_PATH = os.path.join(SRC_PATH, 'templates')


def input_file(ntemps, formula, delta_h, enthalpy_temp=0.0, break_temp=1000.0):
    """ Writes a string for the input file for ThermP.

        :param ntemps: number of temperatures
        :type ntemps: int
        :param formula: chemical formula for species
        :type formula: str
        :param delta_h: enthalpy of formation [Eh]
        :type delta_h: float
        :param enthalpy_temp: temperature corresponding to enthalpy
        :type enthalpy_temp: float
        :param break_temp: temperature delineating low-T and high-T for fits
        :type break_temp: float
        :rtype: str
    """

    # Get the stoichiometry of all elements to build composition string
    atom_dict = util.get_atom_counts_dict(formula)
    composition_str = ''
    for key, val in atom_dict.items():
        composition_str += f'{key}  {val}\n'
    composition_str = composition_str.rstrip()

    # Get the delta H in kcal
    delta_h = f'{delta_h * phycon.EH2KCAL:.3f}'

    # Create a fill value dictionary
    thermp_keys = {
        'ntemps': ntemps,
        'formula': formula,
        'deltaH': delta_h,
        'enthalpyT': enthalpy_temp,
        'breakT': break_temp,
        'composition_str': composition_str
    }

    return build_mako_str(
        template_file_name='thermp.mako',
        template_src_path=TEMPLATE_PATH,
        template_keys=thermp_keys)
