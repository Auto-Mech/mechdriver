"""
  Determine the stability of radicals
"""

import chemkin_io


CKIN_FILES = [
    './abs_mig/l1_abs_mig.ckin',
#    './abs_mig/l2_abs_mig.ckin',
#    './abs_mig/l3_abs_mig.ckin',
#    './abs_mig/l3_cbs_abs_mig.ckin',
]

PRESSURES = [0.1, 1.0, 10.0, 100.0]


def read_radical_stabilities(mech1_str, pressures):
    """ Read the reaction sections of a CHEMKIN files and
        calculate the rate constants
    """

    # Build dict for the rate constants
    mech1_reaction_block = chemkin_io.parser.util.clean_up_whitespace(
        chemkin_io.parser.mechanism.reaction_block(
            mech1_str, remove_comments=False))
    # mech1_units = chemkin_io.parser.mechanism.reaction_units(mech1_str)

    mech1_rxnstr_dct = chemkin_io.parser.reaction.data_dct(
        mech1_reaction_block, data_entry='strings')

    # Add fancy sorting for radical searching with a csv file later

    # Print sorting by the additions
    stab_dct = {}
    radicals = []
    for rxn, dstr in mech1_rxnstr_dct.items():
        [reacs, prods] = rxn
        if len(reacs) == 2 and len(prods) == 1:
            temps = _get_max_temps(dstr, pressures)
            if temps:
                tstr = '   '.join((str(int(val)) for val in temps))
                radicals.append([prods[0], temps[0], temps])
            else:
                tstr = 'use isom'
                radicals.append([prods[0], 1000000])
            stab_dct[prods[0]] = '{}   {}'.format(prods[0], tstr)

    # Get list of radicals starting out at least stable
    radicals.sort(key=lambda x: x[1])

    # Prind the stabilities
    for radical in radicals:
        print(stab_dct[radical[0]])


def _get_max_temps(dstr, pressures):
    """ grab the temperature ranges for various pressures
    """
    tranges = {}
    dlines = dstr.splitlines()
    for line in dlines:
        if 'TempRange' in line and 'High' not in line:
            pressure = line.strip().split()[1]
            trange = line.strip().split()[5]
            maxtemp = trange.split('-')[1]
            tranges[float(pressure)] = maxtemp

    maxtemps = []
    if tranges:
        for pressure in pressures:
            maxtemps.append(float(tranges[pressure]))

    return maxtemps


if __name__ == '__main__':
    for ckin_file in CKIN_FILES:
        with open(ckin_file, 'r') as mech_file:
            mech_str = mech_file.read()
        print(ckin_file)
        read_radical_stabilities(mech_str, PRESSURES)
