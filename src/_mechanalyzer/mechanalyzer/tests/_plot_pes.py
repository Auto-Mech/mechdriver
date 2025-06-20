"""
Test the rate plotting functionality for comparing two mechanisms
"""

import mechanalyzer


ENE_DCT = {
    'R1': 0.00,
    'B1': 15.00,
    'B2': 20.00,
    'P1': -5.00,
    'P2': -2.00
}

CONN_LST = (
    ('R1', 'B1'),
    ('R1', 'B2'),
    ('B1', 'P1'),
    ('B2', 'P2')
)

FORMULA = 'CH2O'


def test__plot_pes():
    """ test mechanalyzer.plotter.pes.build
    """
    mechanalyzer.plotter.pes.build(ENE_DCT, CONN_LST)  # , FORMULA)


if __name__ == '__main__':
    test__plot_pes()
