""" Auxiliary CFOUR energy readers
"""

import autoparse.pattern as app
import autoparse.find as apf


def relativistic_energy(output_str):
    """ Read a relativistic energy
    """

    # 'relativistic' for MVD1
    methods = app.one_of_these(['DPT', 'MVD2', 'relativistic'])

    ptt = (
        'Total energy with' +
        app.SPACE +
        methods +
        app.SPACE +
        'correction:' +
        app.SPACES +
        app.capturing(app.FLOAT) +
        app.SPACES +
        'Hartree'
    )

    cap = apf.last_capture(ptt, output_str)
    val = float(cap) if cap is not None else None

    return val


def diagonal_born_oppenheimer_correction(output_str):
    """ DBOC
    """

    ptt = (
        'The total diagonal Born-Oppenheimer correction (DBOC) is:' +
        app.SPACES +
        app.capturing(app.FLOAT) +
        app.SPACES +
        'a.u.'
    )

    cap = apf.last_capture(ptt, output_str)
    val = float(cap) if cap is not None else None

    return val
