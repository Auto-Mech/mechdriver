""" Tests the writing of the species section
"""

import os
from ioformat import pathtools
import mess_io.writer


PATH = os.path.dirname(os.path.realpath(__file__))
INP_PATH = os.path.join(PATH, 'data', 'inp')


# Various label strings
SPC1_LABEL = 'Mol1'
SPC2_LABEL = 'Mol2'
WELL_LABEL = 'W1'
BIMOL_LABEL = 'P1'
TS_LABEL = 'B1'
AUX1_ID_LABEL = 'Aux1'
AUX2_ID_LABEL = 'Aux2'
AUXBIMOL_ID_LABEL = 'AuxBimol'

# Data Strings
MOL_MESS_STR = """RRHO
  Core RigidRotor
    SymmetryFactor          1.0
  End
  Geometry[angstrom]        3
    O   1.911401284  0.16134481659  -0.05448080419
    N   4.435924209  0.16134481659  -0.05448080419
    N   6.537299661  0.16134481659  -0.05448080419
  Frequencies[1/cm]         9
    100.00  200.00  300.00  400.00  500.00
    600.00  700.00  800.00  900.00
  ElectronicLevels[1/cm]    2
    1  0.0
    3  50.0
"""

ATOM_MESS_STR = """Atom
  Name O
  ElectronicLevels[1/cm]    3
    1  0.0
    3  150.0
    9  450.0
"""

ENE = -35.0
PATH_ENES = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
TUNNEL_STR = """Tunneling  Eckart
  ImaginaryFrequency[1/cm]  2000
  WellDepth[kcal/mol]       10
  WellDepth[kcal/mol]       20"""
EDOWN_STR = """EnergyRelaxation
  Exponential
     Factor[1/cm]                     150.000
     Power                            50.000
     ExponentCutoff                   80.000
End"""
COLLID_STR = """CollisionFrequency
  LennardJones
     Epsilons[1/cm]                   100.000    200.000
     Sigmas[angstrom]                 10.000     20.000
     Masses[amu]                      15.000     25.000
End"""


def test__species_writer():
    """ Writes the MESS input for a Well
    """

    spc1_str = mess_io.writer.species(
        SPC1_LABEL, MOL_MESS_STR,
        aux_id_label=None, zero_ene=None)
    assert spc1_str == pathtools.read_file(INP_PATH, 'spc1.inp').rstrip()

    spc2_str = mess_io.writer.species(
        SPC1_LABEL, MOL_MESS_STR,
        aux_id_label=None, zero_ene=ENE)
    assert spc2_str == pathtools.read_file(INP_PATH, 'spc2.inp').rstrip()

    spc3_str = mess_io.writer.species(
        SPC1_LABEL, MOL_MESS_STR,
        aux_id_label=AUX1_ID_LABEL, zero_ene=None)
    assert spc3_str == pathtools.read_file(INP_PATH, 'spc3.inp').rstrip()


def test__well_writer():
    """ Writes the MESS input for a Well
    """

    well1_str = mess_io.writer.well(
        WELL_LABEL, MOL_MESS_STR)
    assert well1_str == pathtools.read_file(
        INP_PATH, 'well1.inp').rstrip()

    well2_str = mess_io.writer.well(
        WELL_LABEL, MOL_MESS_STR,
        zero_ene=ENE,
        edown_str=EDOWN_STR,
        collid_freq_str=COLLID_STR)
    assert well2_str == pathtools.read_file(
        INP_PATH, 'well2.inp').rstrip()

    well3_str = mess_io.writer.well(
        WELL_LABEL, MOL_MESS_STR,
        aux_id_label=AUX1_ID_LABEL)
    assert well3_str == pathtools.read_file(
        INP_PATH, 'well3.inp').rstrip()


def test__bimolecular_writer():
    """ Writes the MESS input for a bimolecular set
    """

    bimol1_str = mess_io.writer.bimolecular(
        BIMOL_LABEL,
        SPC1_LABEL, ATOM_MESS_STR,
        SPC2_LABEL, MOL_MESS_STR,
        ENE)
    assert bimol1_str == pathtools.read_file(
        INP_PATH, 'bimol1.inp').rstrip()

    bimol2_str = mess_io.writer.bimolecular(
        BIMOL_LABEL,
        SPC1_LABEL, ATOM_MESS_STR,
        SPC2_LABEL, MOL_MESS_STR,
        ENE,
        auxbimol_id_label=AUXBIMOL_ID_LABEL,
        aux1_id_label=AUX1_ID_LABEL,
        aux2_id_label=AUX2_ID_LABEL,
        calc_spc1_density=True,
        calc_spc2_density=True)
    assert bimol2_str == pathtools.read_file(
        INP_PATH, 'bimol2.inp').rstrip()


def test__ts_sadpt_writer():
    """ ts sadpt writer
    """

    ts_sadpt1_str = mess_io.writer.ts_sadpt(
        TS_LABEL, WELL_LABEL, BIMOL_LABEL, MOL_MESS_STR)
    assert ts_sadpt1_str == pathtools.read_file(
        INP_PATH, 'ts_sadpt1.inp').rstrip()

    ts_sadpt2_str = mess_io.writer.ts_sadpt(
        TS_LABEL, WELL_LABEL, BIMOL_LABEL, MOL_MESS_STR,
        zero_ene=ENE, tunnel=TUNNEL_STR)
    assert ts_sadpt2_str == pathtools.read_file(
        INP_PATH, 'ts_sadpt2.inp').rstrip()

    ts_sadpt3_str = mess_io.writer.ts_sadpt(
        TS_LABEL, WELL_LABEL, BIMOL_LABEL, MOL_MESS_STR,
        aux_id_label=AUX1_ID_LABEL)
    assert ts_sadpt3_str == pathtools.read_file(
        INP_PATH, 'ts_sadpt3.inp').rstrip()


def test__ts_variational_writer():
    """ ts mess_io.writer.ts_variational
    """
    # fix
    var_pt_strings = []
    for i in range(10):
        var_pt_string = '+++++++++++++++++++++++++++++++++++'
        var_pt_string += f'! Path Point {str(i+1)}\n'
        var_pt_string += MOL_MESS_STR
        var_pt_strings.append(var_pt_string)

    ts_var1_str = mess_io.writer.ts_variational(
        TS_LABEL, WELL_LABEL, BIMOL_LABEL,
        var_pt_strings,
        zero_enes=PATH_ENES, tunnel='')
    assert ts_var1_str == pathtools.read_file(INP_PATH, 'ts_var1.inp')

    ts_var2_str = mess_io.writer.ts_variational(
        TS_LABEL, WELL_LABEL, BIMOL_LABEL,
        var_pt_strings,
        zero_enes=PATH_ENES, tunnel=TUNNEL_STR)
    assert ts_var2_str == pathtools.read_file(INP_PATH, 'ts_var2.inp')


def test__dummy_writer():
    """ tests mess_io.writer.dummy
    """

    dummy1_str = mess_io.writer.dummy(
        BIMOL_LABEL)
    assert dummy1_str == pathtools.read_file(INP_PATH, 'dummy1.inp')

    dummy2_str = mess_io.writer.dummy(
        BIMOL_LABEL, zero_ene=ENE)
    assert dummy2_str == pathtools.read_file(INP_PATH, 'dummy2.inp')

    dummy3_str = mess_io.writer.dummy(
        BIMOL_LABEL, aux_id_label=AUX1_ID_LABEL)
    assert dummy3_str == pathtools.read_file(INP_PATH, 'dummy3.inp')


def test__configs_union_writer():
    """ tests mess_io.writer.configs_union
    """

    mol_strs = [MOL_MESS_STR for _ in range(3)]
    tunnel_strs = [TUNNEL_STR for _ in range(3)]
    zero_enes = (1.00, 2.00, 3.00)

    union1_str = mess_io.writer.configs_union(
        mol_strs, zero_enes)
    union2_str = mess_io.writer.configs_union(
        mol_strs, zero_enes, tunnel_strs=tunnel_strs)

    assert union1_str == pathtools.read_file(INP_PATH, 'union1.inp')
    assert union2_str == pathtools.read_file(INP_PATH, 'union2.inp')

if __name__ == '__main__':
    test__bimolecular_writer()
    test__configs_union_writer()
    test__dummy_writer()
    test__species_writer()
    test__ts_sadpt_writer()
    test__ts_variational_writer()
    test__well_writer()
