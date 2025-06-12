""" Read useful information from the output file
"""

import ioformat
import autoparse.find as apf
import autoparse.pattern as app


def internal_coordinates(output_str):
    """ Read the internal coordinates defined in output.

        Coordinates returned as a list where each entry is
        (type, (idxs,)) where type is STRE, BEND, TORS, etc and
        idxs correspond to the atom indices of a Cartesian
        geometry associated with the internal coordinates.

        :param output_str: string the INTDER output file
        :type output_str: str
        :rtype: tuple(tuple(str, tuple(int)))
    """

    begin_ptt = 'DEFINITION OF INTERNAL COORDINATES'
    end_ptt = 'VALUES OF SIMPLE INTERNAL COORDINATES'
    block_ptt = app.block_pattern(begin_ptt, end_ptt)
    block = apf.last_capture(block_ptt, output_str)

    if block is not None:
        block = ioformat.remove_whitespace_from_string(block)

        intl_coords = ()
        for line in block.splitlines()[1:]:
            tmp = line.strip().split('=')[1]
            tmp2 = tmp.split()
            # Get the coord typ and idxs
            # only get non-zero idxs as zeros correspond to no atoms
            typ = tmp2[0]
            idxs = tuple(int(val) for val in tmp2[1:])
            idxs = tuple(val-1 for val in idxs if val > 0)

            intl_coords += ((typ, idxs),)
    else:
        intl_coords = None

    return intl_coords


def ted_assignments(output_str):
    """ Parse the TED assignments to each vibrational mode which include
        the vibrational frequency as well as the index of the internal
        coordinates that make up the mode. For each internal coordinate,
        the percentage of contribution to the vibrational mode is also
        returned.

        The output is a tuple for each mode given as (freq, {idx: percentage})
        where freq is in cm-1 and the percentage is given in N%.

        :param output_str: string the INTDER output file
        :type output_str: str
        :rtype: tuple(tuple(float, dict[int: float]))
    """

    # Block patterns
    begin_ptt = 'DOMINANT COMPONENTS OF TED'
    end_ptt = 'NORMAL MODE ANALYSIS IN CARTESIAN COORDINATES'
    block_ptt = app.block_pattern(begin_ptt, end_ptt)
    block = apf.last_capture(block_ptt, output_str)

    # Other patterns
    pct_ptt = app.block_pattern(app.escape('('), app.escape(')'))
    coord_ptt = app.INTEGER + app.SPACE + app.escape('(')

    if block is not None:
        block = ioformat.remove_whitespace_from_string(block)

        # Parse line for mode info
        freqs, idxs_lst, pcts_lst = (), (), ()
        for line in block.splitlines():
            # Get the mod frequency
            freq = float(line.split()[1])
            freqs += (freq,)
            idxs = apf.all_captures(coord_ptt, line)

            # Get the idx numbers of coordinates comprising mode (in 0-index)
            # When subtracting 1 for 0-idx, be careful of negative values
            idxs = tuple(int(val.replace(' (', '')) for val in idxs)
            idxs0 = ()
            for idx in idxs:
                if idx > 0:
                    idxs0 += (idx-1,)
                else:
                    idxs0 += (-1*(abs(idx)-1),)
            idxs_lst += (idxs0,)

            # Get the percentages of coordinates comprising mode
            # If TED %s are undefined (***) we replace with 1e10, ignoring it
            init_pcts = apf.all_captures(pct_ptt, line)
            conv_pcts = ()
            for val in init_pcts:
                if '***' not in val:
                    conv_pcts += (float(val.strip()),)
                else:
                    conv_pcts += (1e10,)
            pcts_lst += (conv_pcts,)

        # Build the TED list
        ted = ()
        for freq, idxs, pcts in zip(freqs, idxs_lst, pcts_lst):
            ted += ((freq, dict(zip(idxs, pcts))),)
    else:
        ted = None

    return ted
