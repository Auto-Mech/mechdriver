""" geometry parsers
"""

from autoparse import cast as _cast
import autoparse.find as apf
import autoparse.pattern as app
from autoread import par


def read(string,
         symb_ptt=par.Pattern.ATOM_SYMBOL,
         val_ptt=par.Pattern.NUMERIC_VALUE,
         start_ptt=None,
         line_sep_ptt=None,
         line_start_ptt=None,
         last=True,
         case=False):
    """ Read atom symbols and xyz coordinates of a Cartesian molecular geometry
        from a string by capturing symbols (from column 1) and
        coordinates (from columns 2-4).

        :param symb_ptt: matches atom symbol in the first column
        :type symb_ptt: str
        :param val_ptt: matches coordinate values in columns 2-4
        :type val_ptt: str
        :param start_ptt: pattern before the start of the geometry block
        :type start_ptt: str
        :param line_sep_ptt: matches separator between
            column 1 and columns 2-4 in each line
        :type line_sep_ptt: str
        :param line_start_ptt: matches at the start of each geometry line
        :type line_start_ptt: str
        :param last: capture the last match, instead of the first?
        :type last: bool
        :param case: make the match case-sensitive?
        :type case: bool
        :rtype: (tuple(str), tuple(tuple(float)))
    """

    line_ptt_ = line_pattern(
        symb_ptt=app.capturing(symb_ptt), val_ptt=app.capturing(val_ptt),
        sep_ptt=line_sep_ptt, start_ptt=line_start_ptt)
    block_ptt_ = app.capturing(block_pattern(
        symb_ptt=symb_ptt, val_ptt=val_ptt, line_sep_ptt=line_sep_ptt,
        line_start_ptt=line_start_ptt))

    block_ptt_ = block_ptt_ if start_ptt is None else start_ptt + block_ptt_

    block_str = (apf.last_capture(block_ptt_, string, case=case) if last else
                 apf.first_capture(block_ptt_, string, case=case))

    caps = apf.all_captures(line_ptt_, block_str)
    if caps is not None:
        symbs, xcomps, ycomps, zcomps = zip(*_cast(caps))
        xyzs = tuple(zip(xcomps, ycomps, zcomps))
    else:
        symbs, xyzs = None, None

    return symbs, xyzs


def read_xyz(string, symb_ptt=par.Pattern.ATOM_SYMBOL,
             val_ptt=par.Pattern.NUMERIC_VALUE):
    """ Reads a Cartesian molecular geometry from an .xyz-formatted string.

        :param string: xyz-string containing geometry
        :type string: str
        :param symb_ptt: matches atom symbol in the first column of xyz lines
        :type symb_ptt: str
        :param val_ptt: matches coordinate values in columns 2-4 of xyz lines
        :type val_ptt: str
        :rtype: (tuple(str), tuple(tuple(float)))
    """

    lines = string.splitlines()
    try:
        natms = int(lines[0])
    except ValueError as valerr:
        raise ValueError('invalid xyz string') from valerr

    geo_str = '\n'.join(lines[2:natms+2])
    symbs, xyzs = read(geo_str, symb_ptt=symb_ptt, val_ptt=val_ptt)

    return symbs, xyzs


def block_pattern(symb_ptt=par.Pattern.ATOM_SYMBOL,
                  val_ptt=par.Pattern.NUMERIC_VALUE,
                  line_sep_ptt=None,
                  line_start_ptt=None):
    """ Builds a pattern that can match a Cartesian molecular geometry from
        block of lines.

        :param symb_ptt: matches atom symbol in the first column of xyz lines
        :type symb_ptt: str
        :param val_ptt: matches coordinate values in columns 2-4 of xyz lines
        :type val_ptt: str
        :param line_sep_ptt: pattern that delimits columns of xyz line
        :type line_sep_ptt: str
        :param line_start_ptt: pattern preceding atom symbols of geometry block
        :type line_start_ptt: str
        :rtype: str
    """

    line_ptt = line_pattern(
        symb_ptt=symb_ptt, val_ptt=val_ptt, sep_ptt=line_sep_ptt,
        start_ptt=line_start_ptt)
    block_ptt = app.series(line_ptt, app.padded(app.NEWLINE))

    return block_ptt


def line_pattern(symb_ptt=par.Pattern.ATOM_SYMBOL,
                 val_ptt=par.Pattern.NUMERIC_VALUE,
                 sep_ptt=None,
                 start_ptt=None):
    """ Builds a pattern that can match an atom line within a block of text
        for a Cartesian molecular geometry.

        :param symb_ptt: matches atom symbol in the first column of xyz line
        :type symb_ptt: str
        :param val_ptt: matches coordinate values in columns 2-4 of xyz line
        :type val_ptt: str
        :param line_sep_ptt: pattern that delimits columns of xyz line
        :type line_sep_ptt: str
        :param line_start_ptt: pattern preceding atom symbols of geometry block
        :type line_start_ptt: str
        :rtype: str
    """

    parts = (
        ([] if start_ptt is None else [start_ptt]) +
        [symb_ptt] +
        ([] if sep_ptt is None else [sep_ptt]) +
        3 * [val_ptt]
    )
    ptt = app.LINE_START + app.padded(app.LINESPACES.join(parts))

    return ptt
