""" test mess_io.reader._lump
"""

import os
import numpy
from ioformat import pathtools
import mess_io.reader


PATH = os.path.dirname(os.path.realpath(__file__))
AUX_PATH = os.path.join(PATH, 'data', 'out')

AUX_STR = pathtools.read_file(AUX_PATH, 'c3h3_rate.auxf')
LOG_STR = pathtools.read_file(AUX_PATH, 'c3h3_rate.logf')

TEMP1 = 600
TEMP2 = 1200
TEMP3 = 2000
PRESSURE = 1.0


def test__merge():
    """ test mess_io.reader._wells.merged_wells
    """

    assert mess_io.reader.merged_wells(AUX_STR, PRESSURE, TEMP1) == (
        ('W1', 'W5'), ('W2', 'W6'))
    assert mess_io.reader.merged_wells(AUX_STR, PRESSURE, TEMP2) == (
        ('W1', 'W3', 'W5', 'W6'),)
    assert mess_io.reader.merged_wells(AUX_STR, PRESSURE, TEMP3) == (
        ('W1', 'W2', 'W4', 'W5', 'W3'),)


def test__energies():
    """ test mess_io.reader._wells.well_average_energy
    """
    assert numpy.isclose(
        mess_io.reader.well_thermal_energy(LOG_STR, 'W1', TEMP1),
        0.009402248486357753)
    assert numpy.isclose(
        mess_io.reader.well_thermal_energy(LOG_STR, 'W1', TEMP2),
        0.027091224452217254)
    assert numpy.isclose(
        mess_io.reader.well_thermal_energy(LOG_STR, 'W1', TEMP3),
        0.055776050342800226)
    assert numpy.isclose(
        mess_io.reader.well_thermal_energy(LOG_STR, 'W3', TEMP1),
        0.07330566616482315)
    assert numpy.isclose(
        mess_io.reader.well_thermal_energy(LOG_STR, 'W3', TEMP3),
        0.11904202744591934)

if __name__ == '__main__':
    test__energies()
    test__merge()