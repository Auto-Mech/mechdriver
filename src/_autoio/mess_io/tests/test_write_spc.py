""" Tests the writing of the energy transfer section
"""

import os
from ioformat import pathtools
import mess_io.writer


PATH = os.path.dirname(os.path.realpath(__file__))
INP_PATH = os.path.join(PATH, 'data', 'inp')


# Atom/Molecule Data
ddPATH = os.path.dirname(os.path.realpath(__file__))
CORE_STR = """  Geometry[angstrom]        8
    C         -0.75583       0.00710      -0.01604
    C          0.75582      -0.00710       0.01604
    H         -1.16275      -0.10176       0.99371
    H         -1.12255       0.94866      -0.43557
    H         -1.13499      -0.81475      -0.63070
    H          1.13499       0.81475       0.63070
    H          1.16275       0.10176      -0.99371
    H          1.12255      -0.94866       0.43557
  Core RigidRotor
    SymmetryFactor          3.0
  End"""
#GEO = (
#    ('C', (-0.75583,  0.00710, -0.01604)),
#    ('C', ( 0.75582, -0.00710,  0.01604)),
#    ('H', (-1.16275, -0.10176,  0.99371)),
#    ('H', (-1.12255,  0.94866, -0.43557)),
#    ('H', (-1.13499, -0.81475, -0.63070)),
#    ('H', ( 1.13499,  0.81475,  0.63070)),
#    ('H', ( 1.16275,  0.10176, -0.99371)),
#    ('H', ( 1.12255, -0.94866,  0.43557)))
GEO = (
    ('C', (-1.4283116974047902, 0.013417055490750581, -0.03031120705234356)), 
    ('C', (1.4282928001435358, -0.013417055490750581, 0.03031120705234356)), 
    ('H', (-2.1972790523760897, -0.1922985305265886, 1.8778397481286984)), 
    ('H', (-2.121312062132685, 1.7927075861768231, -0.8231080084656662)), 
    ('H', (-2.1448202551333804, -1.5396543607167654, -1.1918502673262523)), 
    ('H', (2.1448202551333804, 1.5396543607167654, 1.1918502673262523)), 
    ('H', (2.1972790523760897, 0.1922985305265886, -1.8778397481286984)), 
    ('H', (2.121312062132685, -1.7927075861768231, 0.8231080084656662)))
SYM = 3.
FREQS = (10.0, 100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0, 800.0, 900.0,
         1000.0, 1100.0, 1200.0, 1300.0, 1400.0, 1500.0, 1600.0, 1700.0)
FREQ_SCALE_FACTOR = 1.833
INF_INTENS = (1.0, 11.0, 22.0, 33.0, 44.0, 55.0, 66.0, 77.0, 88.0, 99.0,
              110.0, 120.0, 130.0, 140.0, 150.0, 160.0, 170.0, 180.0)
MASS = 16.0
ELEC_LEVELS1 = ((1, 0.00), (3, 150.0), (9, 450.0))
ELEC_LEVELS2 = ((1, 0.00),)
HR_STR = """Rotor  Hindered
  Group                11  10  9   8   7   6
  Axis                 3   2   1
  Symmetry             1
  Potential[kcal/mol]  12
    0.00    2.91    9.06    12.63   9.97    3.51
    0.03    3.49    9.96    12.63   9.08    2.93
End
"""
XMAT = ((10000, 2000, 4000),
        (2000, 3000, 5000),
        (4000, 5000, 6000))
ROVIB_COUPS = ((100, 200, 300),)
ROT_DISTS = (('aaaa', 1000), ('bbaa', 2000), ('bbbb', 3000))


def test__atom_writer():
    """ Writes a string containing all info for an atom in MESS style
    """

    atom_str = mess_io.writer.atom(MASS, ELEC_LEVELS1)
    assert atom_str == pathtools.read_file(INP_PATH, 'atom_data.inp')


def test__molecule_writer():
    """ Writes a string containing all info for a molecule in MESS style
    """

    # mol1_str = mess_io.writer.molecule(
    #     CORE_STR, ELEC_LEVELS2)
    # assert mol1_str == pathtools.read_file(INP_PATH, 'mol1_data.inp')

    # mol2_str = mess_io.writer.molecule(
    #     CORE_STR, ELEC_LEVELS2,
    #     freqs=FREQS,
    #     freq_scale_factor=FREQ_SCALE_FACTOR,
    #     use_harmfreqs_key=True)
    # assert mol2_str == pathtools.read_file(INP_PATH, 'mol2_data.inp')
    mol3_core_str = mess_io.writer.core_rigidrotor(
        geo=GEO,
        sym_factor=SYM,
        interp_emax=None,
        freqs=FREQS,
        xmat=XMAT,
        rovib_coups=ROVIB_COUPS,
        rot_dists=ROT_DISTS,
    )

    mol3_str = mess_io.writer.molecule(
        mol3_core_str, ELEC_LEVELS2,
        freqs=FREQS, hind_rot=HR_STR,
        xmat=XMAT,
        rovib_coups=ROVIB_COUPS,
        rot_dists=ROT_DISTS,
        inf_intens=INF_INTENS)
    assert mol3_str == pathtools.read_file(INP_PATH, 'mol3_data.inp')


if __name__ == '__main__':
    test__atom_writer()
    test__molecule_writer()
