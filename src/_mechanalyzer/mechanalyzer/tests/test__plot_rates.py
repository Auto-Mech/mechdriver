""" Test the mechanalyzer.plotter.rates functions
"""


import tempfile
import numpy as np
from mechanalyzer.plotter import rates
from mechanalyzer.plotter._util import build_pdf


TMP_DIR = tempfile.mkdtemp()
FILENAME = 'rate_comparison.pdf'
print('Temp Run Dir:', TMP_DIR)

TEMPS = np.array([500, 1000, 1500])
KTS = np.array([1e10, 1e11, 1e12])


ALGN_RXN_KTP_DCT = {
    (('H2', 'O'), ('OH', 'H'), (None,)): [
        {'high': (TEMPS, KTS), 1: (TEMPS, KTS), 10: (TEMPS, KTS)},
        {'high': (TEMPS, 2*KTS), 10: (TEMPS, 2*KTS)}],

    (('H', 'O2'), ('OH', 'O'), (None,)): [
        {'high': (TEMPS, KTS), 1: (TEMPS, KTS), 10: (TEMPS, KTS)},
        None]
}


def test_build_plots():
    """ Test the build_plots function
    """
    figs = rates.build_plots(ALGN_RXN_KTP_DCT)
    build_pdf(figs, FILENAME, TMP_DIR)


def test_build_plots_ratio_sort():
    """ Test the build_plots function with sorting by ratio
    """
    figs = rates.build_plots(ALGN_RXN_KTP_DCT, ratio_sort=True)
    build_pdf(figs, FILENAME, TMP_DIR)


if __name__ == '__main__':
    test_build_plots()
    test_build_plots_ratio_sort()
