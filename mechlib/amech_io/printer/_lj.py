""" lj stuff
"""

from mechlib.amech_io.printer._print import info_message


def lennard_jones_params(sigmas, epsilons):
    """ Print Lennard-Jones parameters on individual lines.

        Need to convert?

        :param sigmas: list of sigma values
        :type sigmas: tuple(float)
        :param epsilons: list of epsilon values
        :type epsilons: tuple(float)
    """
    if sigmas and epsilons:
        info_message(
            '{0:<14s}{1:<16s}'.format('\nSigma (Ang)', 'Epsilon (cm-1)'))
        for sig, eps in zip(sigmas, epsilons):
            info_message(
                '{0:<14.4f}{1:<16.4f}'.format(sig, eps))
