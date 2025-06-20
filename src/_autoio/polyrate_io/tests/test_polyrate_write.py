"""
  test an intder writer
"""

import os
import numpy
from ioformat import pathtools
import polyrate_io.writer


PATH = os.path.dirname(os.path.realpath(__file__))
DAT_PATH = os.path.join(PATH, 'data')


# Base data to build the info objects
GEO = (('C', (-4.0048955763, -0.3439866053, -0.0021431734)),
       ('O', (-1.3627056155, -0.3412713280, 0.0239463418)),
       ('H', (-4.7435343957, 1.4733340928, 0.7491098889)),
       ('H', (-4.7435373042, -1.9674678465, 1.1075144307)),
       ('H', (-4.6638955748, -0.5501793084, -1.9816675556)),
       ('H', (-0.8648060003, -0.1539639444, 1.8221471090)))
GRAD = ((-4.0048955763, -0.3439866053, -0.0021431734),
        (-1.3627056155, -0.3412713280, 0.0239463418),
        (-4.7435343957, 1.4733340928, 0.7491098889),
        (-4.7435373042, -1.9674678465, 1.1075144307),
        (-4.6638955748, -0.5501793084, -1.9816675556),
        (-0.8648060003, -0.1539639444, 1.8221471090))
HESS = ((0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
        (0.0, 0.959, 0.0, 0.0, -0.452, 0.519, 0.0, -0.477, -0.23),
        (0.0, 0.0, 0.371, 0.0, 0.222, -0.555, 0.0, -0.279, -0.128),
        (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
        (0.0, -0.479, 0.279, 0.0, 0.455, -0.256, 0.0, -0.017, 0.051),
        (0.0, 0.251, -0.185, 0.0, -0.247, 0.607, 0.0, -0.012, 0.09),
        (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
        (0.0, -0.479, -0.279, 0.0, -0.003, -0.263, 0.0, 0.494, 0.279),
        (0.0, -0.251, -0.185, 0.0, 0.025, 0.947, 0.0, 0.292, 0.137))
GEOS = [GEO for i in range(21)]
GRADS = [GRAD for i in range(21)]
HESSES = [HESS for i in range(21)]

SVALS = list(x for x in numpy.arange(1.0, 22.0, 1.0))
VVALS = list(x for x in numpy.arange(1.0, 22.0, 1.0))

# Build info objects
RCT_INFO = (HESS, -10.0, 0.0)
PRD_INFO = (HESS, -20.0, 22.0)
SADPT_INFO = (HESS, 10.0)
MEP_INFOS = ()
for i in range(21):
    MEP_INFOS += ((HESSES[i], VVALS[i], SVALS[i], GEOS[i], GRADS[i]),)


def test__pot40():
    """ test polyrate_io.writer.potential_file
    """

    inp_str = polyrate_io.writer.potential_file(
        RCT_INFO, SADPT_INFO, MEP_INFOS)
    assert inp_str == pathtools.read_file(DAT_PATH, 'pot.fu40')
