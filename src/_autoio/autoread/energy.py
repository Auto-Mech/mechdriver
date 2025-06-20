""" energy parsers
"""

from autoparse import cast as _cast
import autoparse.find as apf
import autoparse.pattern as app


VALUE_PATTERN = app.one_of_these([app.EXPONENTIAL_FLOAT, 
                                  app.EXPONENTIAL_FLOAT_D, 
                                  app.FLOAT])
SEP_PATTERN = app.LINESPACES


def read(string,
         start_ptt,
         val_ptt=VALUE_PATTERN,
         last=True,
         case=False):
    """ read energy from a string

        :param start_ptt: matches before the numeric value
        :type start_ptt: str
        :param val_ptt: matches the numeric value
        :type val_ptt: str
        :param last: capture the last match, instead of the first?
        :type last: bool
        :param case: make the match case-sensitive?
        :type case: bool
        :return: ene
        :rtype: float
    """

    ptt_ = pattern(start_ptt=start_ptt, val_ptt=app.capturing(val_ptt))
    ene_str = (apf.last_capture(ptt_, string, case=case) if last else
               apf.first_capture(ptt_, string, case=case))
    if ene_str is not None:
        ene = _cast(ene_str.replace('D', 'E'))
    else:
        ene = None

    return ene


def pattern(start_ptt,
            val_ptt=VALUE_PATTERN):
    """ Build a regex string for the energy pattern.

        :param start_ptt: matches before the numeric value
        :type start_ptt: str
        :param val_ptt: matches the numeric value
        :type val_ptt: str
        :rtype: str
    """
    return start_ptt + app.lpadded(val_ptt)
