"""
Writes strings containing NASA-7 polynomials in the Chemkin format
"""


def thermo_entry(spc_name, nasa7_params):
    """ Writes a string for a single thermo entry
    """
    notes = nasa7_params[0]
    composition = nasa7_params[1]
    phase = nasa7_params[2]
    tlimits = nasa7_params[3]
    try:
        high_coeffs = nasa7_params[4][0]
    except TypeError:
        print(f'TypeError for the species {spc_name}')
    low_coeffs = nasa7_params[4][1]

    line1 = (
        '{0:<16s}{1:<2s}{2:<6s}'.format(spc_name, '', notes) +
        '{0:>20s}{1:<1s}{2:>10.2f}'.format(composition, phase, tlimits[0]) +
        '{0:>10.2f}{1:>8.2f}'.format(tlimits[1], tlimits[2]) +
        '{0:>5s}{1:<2s}\n'.format('', ' 1')
    )
    # note: need to fix line 1 to incorporate
    # more detailed composition information
    line2 = (
        '{0:>15.8E}{1:>15.8E}{2:>15.8E}{3:>15.8E}{4:>15.8E}{5:>5s}\n'
        ).format(
            high_coeffs[0], high_coeffs[1], high_coeffs[2], high_coeffs[3],
            high_coeffs[4], '2')
    line3 = (
        '{0:>15.8E}{1:>15.8E}{2:>15.8E}{3:>15.8E}{4:>15.8E}{5:>5s}\n'
        ).format(
            high_coeffs[5], high_coeffs[6], low_coeffs[0],
            low_coeffs[1], low_coeffs[2], '3')
    line4 = (
        '{0:>15.8E}{1:>15.8E}{2:>15.8E}{3:>15.8E}{4:>20s}\n'
        ).format(
            low_coeffs[3], low_coeffs[4],
            low_coeffs[5], low_coeffs[6], '4')

    thermo_str = line1 + line2 + line3 + line4

    return thermo_str
