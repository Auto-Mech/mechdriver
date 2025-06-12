"""
  Tests the varecof_io.writer functions
"""


import os
from ioformat import pathtools
import varecof_io.writer


PATH = os.path.dirname(os.path.realpath(__file__))
DAT_PATH = os.path.join(PATH, 'data')


# TST parameters
NSAMP_MAX = 2000
NSAMP_MIN = 500
FLUX_ERR = 5
PES_SIZE = 1
FACES = (0,)
FACES_SYMM = 1
ENER_GRID = [0, 20, 1.05, 179]
AMOM_GRID = [0, 3, 1.10, 50]

# Divsur
RDISTS = [8., 6., 4., 3.8, 3.6, 3.4, 3.2, 3., 2.8, 2.6, 2.4, 2.2]
NPIVOT1 = 1
NPIVOT2 = 2
XYZ_PIVOT1 = (1., 2., 3.)
XYZ_PIVOT2 = (1., 2., 3.)
FRAME1 = (0, 1, 2, 3)
FRAME2 = (0, 1, 2, 3)
R2DISTS = [8., 6., 4., 3.8, 3.6, 3.4, 3.2, 3., 2.8, 2.6, 2.4, 2.2]
D1DISTS = [0.01, 0.5, 1.]
D2DISTS = [0.01, 0.5, 1.]
T1ANGS = []
T2ANGS = []
P1ANGS = []
P2ANGS = []
CONDITIONS_DCT = {'r1 > r2': 5}

#  Structure
C2H5_GEO = (
    ('C', (-1.4283035320563338, 0.013425343735546437, -0.030302158896694683)),
    ('C', (1.4283027358735494, -0.013425597530894248, 0.0303022919384165)),
    ('H', (-2.1972722614281355, -0.19229727219177065, 1.8778380427620682)),
    ('H', (-2.121310184939721, 1.792702413487708, -0.8231106338374065)),
    ('H', (-2.1448124562913287, -1.5396513482615042, -1.191852168914227)),
    ('H', (2.1972712765396953, 0.1922944277301287, -1.8778395029874426)),
    ('H', (2.121312248031497, -1.7927029137609576, 0.8231123911174519)))
OH_GEO = (
    ('O', (1.4283027358735494, -0.013425597530894248, 0.0303022919384165)),
    ('H', (-2.1972722614281355, -0.19229727219177065, 1.8778380427620682)))

# Elstruct
EXE_PATH = '/path/to/exe'
LIB_PATH = '/path/to/lib'
BASE_NAME = 'mol'
NPOT = 3
DUMMY_NAME = 'x_corr_'
LIB_NAME = 'corrpot.so'
GEO_PTT = 'geometry_here'
ENE_PTT = 'energy_here'

# ramdp, parameters
FORTRAN_COMPILER = 'gfortran'
SPC_NAME = 'mol'
MEMORY = 4.0
R1DISTS = [8., 6., 4., 3.8, 3.6, 3.4, 3.2, 3., 2.8, 2.6, 2.4, 2.2]
R2DISTS_SR = [4., 3.8, 3.6, 3.4, 3.2, 3., 2.8, 2.6, 2.4, 2.2]
CONDITIONS = {}
FLUX_ERR = 10
PES_SIZE = 2
EXE_PATH = '/blues/gpfs/home/sjklipp/bin/molpro'


def test__tst_writer():
    """ tests varecof_io.writer.input_file.tst
    """

    tst_inp1_str = varecof_io.writer.input_file.tst(
        NSAMP_MAX, NSAMP_MIN, FLUX_ERR, PES_SIZE)

    tst_inp2_str = varecof_io.writer.input_file.tst(
        NSAMP_MAX, NSAMP_MIN, FLUX_ERR, PES_SIZE,
        faces=FACES,
        faces_symm=FACES_SYMM,
        ener_grid=ENER_GRID,
        amom_grid=AMOM_GRID)

    assert tst_inp1_str == pathtools.read_file(DAT_PATH, 'tst1.inp')
    assert tst_inp2_str == pathtools.read_file(DAT_PATH, 'tst2.inp')


def test__divsur_writer():
    """ tests varecof_io.writer.input_file.divsur
    """

    divsur_inp1_str = varecof_io.writer.input_file.divsur(
        RDISTS, NPIVOT1, NPIVOT2, XYZ_PIVOT1, XYZ_PIVOT2)

    divsur_inp2_str = varecof_io.writer.input_file.divsur(
        RDISTS, NPIVOT1, NPIVOT2, XYZ_PIVOT1, XYZ_PIVOT2,
        frame1=FRAME1, frame2=FRAME2,
        r2dists=R2DISTS,
        d1dists=D1DISTS, d2dists=D2DISTS,
        t1angs=T1ANGS, t2angs=T2ANGS,
        p1angs=P1ANGS, p2angs=P2ANGS,
        phi_dependence=False,
        **CONDITIONS_DCT)

    assert divsur_inp1_str == pathtools.read_file(DAT_PATH, 'divsur1.inp')
    assert divsur_inp2_str == pathtools.read_file(DAT_PATH, 'divsur2.inp')


def test__els_writer():
    """ tests varecof_io.writer.input_file.elec_struct
    """

    # Write the els input string
    els_inp1_str = varecof_io.writer.input_file.elec_struct(
        LIB_PATH, BASE_NAME, NPOT)

    els_inp2_str = varecof_io.writer.input_file.elec_struct(
        LIB_PATH, BASE_NAME, NPOT,
        dummy_name=DUMMY_NAME, lib_name=LIB_NAME,
        geo_ptt=GEO_PTT, ene_ptt=ENE_PTT)

    assert els_inp1_str == pathtools.read_file(DAT_PATH, 'els1.inp')
    assert els_inp2_str == pathtools.read_file(DAT_PATH, 'els2.inp')


def test__structure_writer():
    """ tests varecof_io.writer.input_file.structure
    """

    struct_inp_str = varecof_io.writer.input_file.structure(
        C2H5_GEO, OH_GEO)

    assert struct_inp_str == pathtools.read_file(DAT_PATH, 'structure.inp')


def test__mcflux_writer():
    """ tests varecof_io.writer.input_file.mc_flux
    """

    mc_flux_inp_str = varecof_io.writer.input_file.mc_flux()

    assert mc_flux_inp_str == pathtools.read_file(DAT_PATH, 'mc_flux.inp')


def test__convert_writer():
    """ tests varecof_io.writer.input_file.convert
    """

    conv_inp_str = varecof_io.writer.input_file.convert()

    assert conv_inp_str == pathtools.read_file(DAT_PATH, 'conv.inp')


if __name__ == '__main__':
    test__tst_writer()
    test__divsur_writer()
    test__els_writer()
    test__structure_writer()
    test__mcflux_writer()
    test__convert_writer()
