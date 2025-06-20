""" Parses information from output files of OneDMin
"""

import autoparse.pattern as app
import autoparse.find as apf
from phydat import phycon


def lennard_jones(output_string):
    """ reads the lennard jones params from the output

        this function works the lj.out file
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
        sigmas = tuple(float(val) * phycon.ANG2BOHR for val in sigmas)
    if epsilons is not None:
        epsilons = tuple(float(val) for val in epsilons)

    return sigmas, epsilons


def program_version(output_string):
    """ Read the version
    """

    pattern = (
        'OneDMin' + app.SPACES +
        'version' + app.SPACES +
        app.capturing(app.FLOAT)
    )
    capture = apf.last_capture(pattern, output_string)

    return capture


def random_seed_value(output_string):
    """ Read the ranseed from the standard log file.
    """

    pattern = (
        'RANSEED' + app.SPACES +
        '=' + app.SPACES +
        app.capturing(app.INTEGER)
    )
    capture = apf.first_capture(pattern, output_string)
    if capture:
        ranseed = int(capture)
    else:
        ranseed = None

    return ranseed
