"""
  Reads the output of a MESSPF calculation for the
  partition functions and their derivatives
  corresponding to one species.
"""


def partition_function(output_str):
    """ Parses the MESSPF output file string for the parition function
        and related information for a single species.

        :param output_str: string of lines for MESSPF output file
        :type output_str: str
        :return [temps: logq, dq_dt, dq2_dt2]:
            List of temperatures
            loq(Q) where Q is partition function
            1st derivative of Q w/r to temperature
            2nd derivative of Q w/r to temperature
        :rtype: dict[float: tuple(float)]
    """

    temps, logq, dq_dt, dq2_dt2 = tuple(), tuple(), tuple(), tuple()
    for i, line in enumerate(output_str.splitlines()):
        if i not in (0, 1, 2):
            tmp = line.strip().split()
            temps += (float(tmp[0]),)
            logq += (float(tmp[1]),)
            dq_dt += (float(tmp[2]),)
            dq2_dt2 += (float(tmp[3]),)

    # pf_dct = dict(zip(temps, zip(logq, dq_dt, dq2_dt2)))
    # return pf_dct
    return temps, logq, dq_dt, dq2_dt2


def entropy(output_str):
    """ Parses the MESSPF output file string for the entropy
        for a single species.

        :param output_str: string of lines for MESSPF output file
        :type output_str: str
        :return [temps: s_t]:
            List of temperatures
            Entropy
        :rtype: dict[float: tuple(float)]
    """

    temps, s_t = tuple(), tuple()
    for i, line in enumerate(output_str.splitlines()):
        if i not in (0, 1, 2):
            tmp = line.strip().split()
            temps += (float(tmp[0]),)
            s_t += (float(tmp[4]),)

    s_dct = dict(zip(temps, s_t))

    return s_dct


def heat_capacity(output_str):
    """ Parses the MESSPF output file string for the heat capacity
        for a single species.

        :param output_str: string of lines for MESSPF output file
        :type output_str: str
        :return [temps: cp_t]:
            List of temperatures
            Entropy
        :rtype: dict[float: tuple(float)]
    """

    temps, cp_t = tuple(), tuple()
    for i, line in enumerate(output_str.splitlines()):
        if i not in (0, 1, 2):
            tmp = line.strip().split()
            temps += (float(tmp[0]),)
            cp_t += (float(tmp[5]),)

    cp_dct = dict(zip(temps, cp_t))

    return cp_dct
