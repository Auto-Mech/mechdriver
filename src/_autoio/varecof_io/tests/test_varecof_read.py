"""
  Tests the varecof_io.writer functions
"""

import os
from ioformat import pathtools
import automol.geom
import varecof_io


PATH = os.path.dirname(os.path.realpath(__file__))
DAT_PATH = os.path.join(PATH, 'data')

OUT_STR = pathtools.read_file(DAT_PATH, 'divsur.out')


def test__divsur_frag_geoms_reader():
    """ tests varecof_io.reader.divsur.frag_geoms_divsur_frame
    """
    ref_geos = (
        (('C', (0.0, -0.22676713505493937, 0.0)),
         ('H', (3.533787854606139, 1.3039110265659013, 0.0)),
         ('H', (-3.533787854606139, 1.3039110265659013, 0.0))),
        (('H', (0.0, 0.0, 0.0)),)
    )

    geos = varecof_io.reader.divsur.frame_geometries(OUT_STR)

    for rgeo, geo in zip(ref_geos, geos):
        assert automol.geom.almost_equal_dist_matrix(rgeo, geo)


if __name__ == '__main__':
    test__divsur_frag_geoms_reader()
