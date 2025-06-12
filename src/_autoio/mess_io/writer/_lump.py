""" Writes the section that details which wells are lumpted togeher
    on the PES during the Master Equation Simulation.
"""

import os
from ioformat import build_mako_str


# OBTAIN THE PATH TO THE DIRECTORY CONTAINING THE TEMPLATES #
SRC_PATH = os.path.dirname(os.path.realpath(__file__))
TEMPLATE_PATH = os.path.join(SRC_PATH, 'templates')
SECTION_PATH = os.path.join(TEMPLATE_PATH, 'sections')


def well_lump_scheme(well_merge_lst, separator='+'):
    """ Write the string containing all of the wells to be merged
        in MESS well lumping scheme.

        :param well_merge_lst: list of well groups
        :type well_merge_lst: tuple(tuple(str))
        :rtype: str
    """

    # Build string for the wells being merged
    well_str = '  '.join((separator.join(lst) for lst in well_merge_lst))

    # Create dictionary to fill template
    well_lump_keys = {
        'well_str': well_str,
        'separator': separator
    }

    return build_mako_str(
        template_file_name='well_lump.mako',
        template_src_path=SECTION_PATH,
        template_keys=well_lump_keys)
