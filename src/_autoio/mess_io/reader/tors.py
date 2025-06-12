"""
Reads the outoput of a MESSPF calculation for the
frequencies and ZPEs for torsional modes

Right now, we assume only a single species is in the output
"""

import autoparse.pattern as app
import autoparse.find as apf
from phydat import phycon


def min_energy_value(output_str):
    """ Reads the minumum energy on the spline fit
    """

    pattern = (app.escape('potential spline minimal value[kcal/mol] =') +
               app.one_or_more(app.SPACE) +
               app.capturing(app.NUMBER))
    tors_mins = tuple(float(val)
                      for val in apf.all_captures(pattern, output_str))
    return tors_mins


def analytic_frequencies(output_str):
    """ Reads the analytic frequencies for each of the
        hindered rotors from MESS output file string.

        Frequency corresponds to the minimum from a Fourier fit to
        the user supplied hindered rotor potential in the input.

        :param output_str: string of lines of MESS output file
        :type output_str: str
        :return freqs: frequency for each of the rotors
        :rtype: list(float)
    """

    # Pattern for the frequency of a rotor
    pattern = (app.escape('analytic  frequency at minimum[1/cm] =') +
               app.one_or_more(app.SPACE) +
               app.capturing(app.NUMBER))

    # Obtain each frequency from the output string
    tors_freqs = tuple(float(val)
                       for val in apf.all_captures(pattern, output_str))

    return tors_freqs


def grid_minimum_frequencies(output_str):
    """ Reads the analytic frequencies for each of the
        hindered rotors from MESS output file string.

        Frequency corresponds to the minimum from the minimum on the grid
        of the user supplied hindered rotor potential in the input.

        :param output_str: string of lines of MESS output file
        :type output_str: str
        :return freqs: frequency for each of the rotors
        :rtype: list(float)
    """

    # Pattern for the frequency of a rotor
    pattern = (app.escape('first point frequency estimate =') +
               app.one_or_more(app.SPACE) +
               app.capturing(app.NUMBER) +
               app.one_or_more(app.SPACE) +
               app.escape('1/cm'))

    # Obtain each frequency from the output string
    tors_freqs = tuple(float(val)
                       for val in apf.all_captures(pattern, output_str))

    return tors_freqs


def first_point_harmonic_frequencies(output_str):
    """ Reads the analytic frequencies for each of the
        hindered rotors from MESS output file string.

        Frequency corresponds to the minimum from the minimum on the grid
        of the user supplied hindered rotor potential in the input.

        :param output_str: string of lines of MESS output file
        :type output_str: str
        :return freqs: frequency for each of the rotors
        :rtype: list(float)
    """

    # Pattern for the frequency of a rotor
    pattern = (app.escape('harmonic  frequency at first point =') +
               app.one_or_more(app.SPACE) +
               app.capturing(app.NUMBER) +
               app.one_or_more(app.SPACE) +
               app.escape('1/cm'))

    # Obtain each frequency from the output string
    tors_freqs = tuple(float(val)
                       for val in apf.all_captures(pattern, output_str))

    return tors_freqs


def zero_point_vibrational_energies(output_str):
    """ Reads the zero-point energies for each of the
        hindered rotors from MESS output file string.

        :param output_str: string of lines of MESS output file
        :type output_str: str
        :return zpves: zero-point energy for each of the rotors
        :rtype: list(float)
    """

    # Patterns for the ZPVE of a rotor
    num_patterns = (app.EXPONENTIAL_FLOAT, app.FLOAT, app.INTEGER)
    pattern1 = (app.escape('minimum energy[kcal/mol]') +
                app.one_or_more(app.SPACE) +
                '=' +
                app.one_or_more(app.SPACE) +
                app.capturing(app.one_of_these(num_patterns)))

    pattern2 = (app.escape('ground energy [kcal/mol]') +
                app.one_or_more(app.SPACE) +
                '=' +
                app.one_or_more(app.SPACE) +
                app.capturing(app.one_of_these(num_patterns)))

    # Obtain each ZPVE from the output string
    tmp1 = tuple(-float(val)
                 for val in apf.all_captures(pattern1, output_str))
    tmp2 = tuple(float(val)
                 for val in apf.all_captures(pattern2, output_str))
    tors_zpves = [sum(tmp) for tmp in zip(tmp1, tmp2)]

    # Convert to Hartrees
    for i, _ in enumerate(tors_zpves):
        tors_zpves[i] *= phycon.KCAL2EH

    return tuple(tors_zpves)
