""" test intder_io.writer
"""

import os
import automol
from ioformat import pathtools
import intder_io.writer


PATH = os.path.dirname(os.path.realpath(__file__))
DAT_PATH = os.path.join(PATH, 'data')

# ethanol species
C2H5OH_GEO = automol.geom.from_string(
    pathtools.read_file(DAT_PATH, 'c2h5oh.xyz'))
C2H5OH_ZMA = automol.zmat.from_string(
    pathtools.read_file(DAT_PATH, 'c2h5oh.zmat'))
C2H5OH_HESS = pathtools.read_numpy_file(DAT_PATH, 'c2h5oh.hess')

# transition state
# CH4_H_GEO = automol.geom.from_string(
#     pathtools.read_file(DAT_PATH, 'ch4_h.xyz'))
# CH4_H_ZMA = automol.zmat.from_string(
#     pathtools.read_file(DAT_PATH, 'ch4_h.zmat'))
# CH4_H_HESS = pathtools.read_numpy_file(DAT_PATH, 'ch4_h.hess')


def test__input():
    """ test intder_io.writer.input_file
    """

    inp_str = intder_io.writer.input_file(C2H5OH_GEO, C2H5OH_ZMA)

    assert inp_str == pathtools.read_file(DAT_PATH, 'intder.inp').rstrip()


def test__hess():
    """ test intder_io.writer.cart_hess_file
    """

    hess_str = intder_io.writer.cartesian_hessian_file(C2H5OH_HESS)

    assert hess_str == pathtools.read_file(DAT_PATH, 'file15')
