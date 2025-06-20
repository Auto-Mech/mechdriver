"""
  Read the the ThermP output
"""

import autoparse.pattern as app
import autoparse.find as apf
from phydat import phycon


def hf298k(output_str):
    """ Read the Heat of Formation at 298 K.

        :param output_str: string for output file of ThermP
        :type output_str: str
        :rtype: tuple(float)
    """

    ptt = (
        'h298' +
        app.SPACES +
        'final' +
        app.SPACES +
        app.capturing(app.NUMBER)
    )

    caps = apf.all_captures(ptt, output_str)

    if caps:
        hfs = tuple(float(val) for val in caps)
        hfs = hfs[0] * phycon.KCAL2EH
    else:
        hfs = None

    return hfs


def properties_temp_dct(output_str):
    """ Read the heat capacities, entropies, and enthalpy changes
        for the whole set of temperatures in a thermp.out and
        put them as a dictionary value (as a tuple) with temperature
        as their corresponding key.

        :param output_str: string for output file of ThermP
        :type output_str: str
        :rtype: dictionary(tuple)
    """
    ptt = (
        app.LINE_START + app.SPACES +
        app.capturing(
            app.FLOAT + app.SPACES +
            app.FLOAT + app.SPACES +
            app.FLOAT + app.SPACES +
            app.FLOAT + app.SPACES +
            app.EXPONENTIAL_FLOAT + app.SPACES +
            app.FLOAT)
    )
    caps = apf.all_captures(ptt, output_str)
    hsc_t_dct = {}
    if caps:
        for cap in caps:
            temp, heat_cap, entropy, _, _, enthalpy = cap.split()
            hsc_t_dct[float(temp)] = (
                float(heat_cap), float(entropy), float(enthalpy))
    return hsc_t_dct
