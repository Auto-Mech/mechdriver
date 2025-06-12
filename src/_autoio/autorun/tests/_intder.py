""" test autorun.intder
"""

import os
import tempfile
import automol
from ioformat import pathtools
import autorun


PATH = os.path.dirname(os.path.realpath(__file__))
DAT_PATH = os.path.join(PATH, 'data')

GEO = automol.geom.from_string(
    pathtools.read_file(DAT_PATH, 'ts2.xyz'))
ZMA = automol.zmat.from_string(
    pathtools.read_file(DAT_PATH, 'ts2.zmat'))
HESS = pathtools.read_numpy_file(DAT_PATH, 'ts2.hess')


def test__ted_zmatrix_coordinates():
    """ test autorun.intder.frequencies
    """

    # ref_ted_zmat_names = ('R10', 'A4', 'A7')

    script_str = autorun.SCRIPT_DCT['intder']
    # with tempfile.TemporaryDirectory(dir=PATH) as run_dir:
    run_dir = tempfile.mkdtemp()
    print(run_dir)
    ted_zmat_names = autorun.intder.ted_zmatrix_coordinates(
        script_str, run_dir,
        GEO, ZMA, HESS, 0)
    #     assert ref_ted_zmat_names == ted_zmat_names
    print(ted_zmat_names)


if __name__ == '__main__':
    test__ted_zmatrix_coordinates()
