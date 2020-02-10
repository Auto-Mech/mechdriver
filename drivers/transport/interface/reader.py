"""
reads info
"""

import autoparse.pattern as app
import autoparse.find as apf


def lennard_jones(output_string):
    """ reads the lennard jones params from the output
    """

    sigma_ptt = (app.SPACES + app.INTEGER + app.SPACES +
                 app.capturing(app.FLOAT) + app.SPACES +
                 app.FLOAT)
    epsilon_ptt = (app.SPACES + app.INTEGER + app.SPACES +
                   app.FLOAT + app.SPACES +
                   app.capturing(app.FLOAT))

    sigmas = apf.all_captures(sigma_ptt, output_string)
    epsilons = apf.all_captures(epsilon_ptt, output_string)
    if sigmas is not None:
        sigmas = [float(val) for val in sigmas]
    if epsilons is not None:
        epsilons = [float(val) for val in epsilons]

    return sigmas, epsilons


def combine_params(param1, param2, rule='default'):
    """ perform a combining rule for two parameters
    """

    if rule == 'default':
        combined_param = (param1 + param2) / 2.0
    else:
        raise NotImplementedError

    return combined_param
