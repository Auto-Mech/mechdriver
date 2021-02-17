"""
Helper function
"""


def pull_highp_from_dct(param_dct):
    """ seperate the high pressure rates from the param
    """

    # Build Pressure dependent ktp dct
    pdep_dct = {}
    pressures = [pressure for pressure in param_dct.keys()
                 if pressure != 'high']
    for pressure in pressures:
        pdep_dct[pressure] = param_dct[pressure]

    # Get the high pressure parameters
    if 'high' in param_dct:
        highp_params = param_dct['high']
    else:
        highp_params = tuple()

    return highp_params, pdep_dct, pressures
