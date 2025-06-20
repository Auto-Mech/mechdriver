""" Read output files of the NST program
"""

import autoread as ar
import autoparse.pattern as app
import autoparse.find as apf
import automol


def optimized_msx_geometry(output_str):
    """ Read the optimized geometry that exists on the
        minimum at the crossing seam for the two potential
        energy surfaces.

        :param output_str: string of the program's output file
        :type output_str: str
        :rtype: automol molecular geometry data structure
    """

    symbs, xyzs = ar.geom.read(
        output_str,
        start_ptt=app.escape('Optimized geometry (A)')+app.NEWLINE)

    if all(x is not None for x in (symbs, xyzs)):
        geo = automol.geom.from_data(symbs, xyzs, angstrom=True)
    else:
        geo = None

    return geo


def rotated_geometry(output_str):
    """ Read the geometry that has been rotated into the
        appropriate frame to do Hessian calculations
        that required for MSX vibrational frequencies.

        :param output_str: string of the program's output file
        :type output_str: str
        :rtype: automol molecular geometry data structure
    """

    symbs, xyzs = ar.geom.read(
        output_str,
        start_ptt=app.escape('Rotated geometry (A)')+app.NEWLINE)

    if all(x is not None for x in (symbs, xyzs)):
        geo = automol.geom.from_data(symbs, xyzs, angstrom=True)
    else:
        geo = None

    return geo


def msx_vibrational_frequencies(output_str):
    """ Read the vibrational frequencies for the geometry that exists on the
        minimum at the crossing seam for the two potential
        energy surfaces. The frequencies are computed internally
        using an effective Hessian combining the Hessians computed
        at the two PES spin states.

        :param output_str: string of the program's output file
        :type output_str: str
        :rtype: automol molecular geometry data structure
    """

    block_ptt = app.block_pattern(
        begin_pattern='MSX frequencies',
        end_pattern='Finished writing the nej.dat and ne.dat')
    block = apf.first_capture(block_ptt, output_str)

    if block is not None:
        vals = apf.all_captures(app.capturing(app.FLOAT), block)
        freqs = tuple(float(val) for val in vals)
    else:
        freqs = None

    return freqs
