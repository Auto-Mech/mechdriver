""" Handle bad files for mechs and spc (prob disperse)
"""

import mechanalyzer.parser.mech as mparser


def test__mech_err():
    """ test mechanalyzer.parser.mech.
    """

    try:
        mparser.parse_mechanism('', 'bad-mech-type')
    except NotImplementedError:
        pass


if __name__ == '__main__':
    test__mech_err()
