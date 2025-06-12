""" Test reading files for ThermP
"""

import os
import numpy
from ioformat import pathtools
import thermp_io


PATH = os.path.dirname(os.path.realpath(__file__))
DAT_PATH = os.path.join(PATH, 'data')

OUT_STR1 = pathtools.read_file(DAT_PATH, 'thermp.out', PATH)
OUT_STR2 = ''


def test__hf298k():
    """ test thermp_io.readar.hf298k
    """

    # All values from file, but reader only returning first
    ref_hf298k = (0.1146914955191809, 0.08766401512449826,
                  0.07924181529799032, 0.0810591889970729)

    hf298k_1 = thermp_io.reader.hf298k(OUT_STR1)
    assert numpy.allclose(hf298k_1, ref_hf298k[0])

    hf298k_2 = thermp_io.reader.hf298k(OUT_STR2)
    assert hf298k_2 is None


def test_property_dct():
    """ test thermp_io.readar.hf298k
    """
    ref_out = {
        200.0: (8.7172, 42.94, 1.6312), 300.0: (9.4123, 46.607, 2.5382),
        400.0: (10.115, 49.411, 3.5146), 500.0: (10.818, 51.744, 4.5616),
        600.0: (11.495, 53.776, 5.6776), 700.0: (12.145, 55.597, 6.86),
        800.0: (12.771, 57.26, 8.1063), 900.0: (13.369, 58.799, 9.4138),
        1000.0: (13.933, 60.237, 10.779), 1100.0: (14.458, 61.59, 12.2),
        1200.0: (14.939, 62.869, 13.67), 1300.0: (15.377, 64.083, 15.187),
        1400.0: (15.772, 65.237, 16.745), 1500.0: (16.127, 66.337, 18.34),
        1600.0: (16.446, 67.388, 19.969), 1700.0: (16.731, 68.394, 21.629),
        1800.0: (16.986, 69.358, 23.315), 1900.0: (17.215, 70.283, 25.026),
        2000.0: (17.42, 71.171, 26.758), 2100.0: (17.604, 72.025, 28.51),
        2200.0: (17.77, 72.848, 30.279), 2300.0: (17.919, 73.641, 32.064),
        2400.0: (18.054, 74.407, 33.863), 2500.0: (18.177, 75.147, 35.675),
        2600.0: (18.288, 75.862, 37.498), 2700.0: (18.389, 76.554, 39.333),
        2800.0: (18.481, 77.224, 41.177), 2900.0: (18.565, 77.874, 43.029),
        3000.0: (18.642, 78.505, 44.89), 298.15: (9.3995, 46.549, 2.5208)}
    # All values from file, but reader only returning first
    out1 = thermp_io.reader.properties_temp_dct(OUT_STR1)
    assert out1 == ref_out

    out2 = thermp_io.reader.properties_temp_dct(OUT_STR2)
    assert not out2
