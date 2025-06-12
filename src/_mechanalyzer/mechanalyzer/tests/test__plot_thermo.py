""" Test the mechanalyzer.plotter.rates functions
"""

import tempfile
import numpy as np
from mechanalyzer.plotter import thermo
from mechanalyzer.plotter._util import build_pdf

TMP_DIR = tempfile.mkdtemp()
FILENAME = 'thermo_comparison.pdf'
print('Temp Run Dir:', TMP_DIR)

TEMPS = np.array([500, 1000, 1500])
H_T = np.array([-5000, -4000, -3000])
CP_T = np.array([-5000, -4000, -3000])
S_T = np.array([-5000, -4000, -3000])
G_T = np.array([-5000, -4000, -3000])
lnq_t = np.array([-5000, -4000, -3000])

ALGN_SPC_THERM_DCT = {
    'H': [
        [TEMPS, H_T, CP_T, S_T, G_T, lnq_t],
        [TEMPS, 1.5 * H_T, 1.5 * CP_T, 1.5 * S_T, 1.5 * G_T, 1.5 * lnq_t],
    ]
}


def test_build_plots():
    """ Test the build_plots function
    """
    figs, _ = thermo.build_plots(ALGN_SPC_THERM_DCT)
    build_pdf(figs, FILENAME, TMP_DIR)


if __name__ == '__main__':
    test_build_plots()
