""" Read various OneDMin output files
"""

import os
import numpy
from ioformat import pathtools
import onedmin_io.reader


PATH = os.path.dirname(os.path.realpath(__file__))
DAT_PATH = os.path.join(PATH, 'data')

LJ_OUT_STR = pathtools.read_file(DAT_PATH, 'lj.out')
ONEDMIN_OUT_STR = pathtools.read_file(DAT_PATH, 'onedmin.out')


def test__lennard_jones():
    """ test onedmin_io.reader.lennard_jones
    """

    ref_sigmas = (6.587188430859643, 6.908800920151311, 6.609902938887646)
    ref_epsilons = (17.88616, 15.23453, 17.731)

    sigmas, epsilons = onedmin_io.reader.lennard_jones(LJ_OUT_STR)

    assert numpy.allclose(ref_sigmas, sigmas)
    assert numpy.allclose(ref_epsilons, epsilons)


def test__program_version():
    """ test onedmin_io.reader.program_version
    """

    assert onedmin_io.reader.program_version(ONEDMIN_OUT_STR) == '1.0'


def test__ranseed():
    """ test onedmin_io.reader.random_seed_value
    """

    assert onedmin_io.reader.random_seed_value(ONEDMIN_OUT_STR) == 153214316
    assert onedmin_io.reader.random_seed_value('') is None
