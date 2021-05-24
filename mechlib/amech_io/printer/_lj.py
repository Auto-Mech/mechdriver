""" lj stuff
"""

def print_lj_parms(sigmas, epsilons):
    """ Print the lj parameters out
    """
    if sigmas and epsilons:
        ioprinter.info_message(
            '{0:<14s}{1:<16s}'.format('\nSigma (Ang)', 'Epsilon (cm-1)'))
        for sig, eps in zip(sigmas, epsilons):
            ioprinter.info_message(
                '{0:<14.4f}{1:<16.4f}'.format(sig, eps))

