""" Tests the err.py module
"""

import numpy as np
from autoreact.params import RxnParams
from ratefit.fit import err


KTP_DCT1 = {
   0.01:  (np.array([1000., 1500.]),
           np.array([1e-5, 1.6e-1])),
   'high': (np.array([1000., 1500.]),
            np.array([2.8e-2, 9.9e2]))
}


LIND_DCT = {'highp_arr': [[1.26e12, 0, 62620]],
            'lowp_arr':  [[1.04e15, 0, 59810]]}
PARAMS1 = RxnParams(lind_dct=LIND_DCT)


def test_get_err_dct():
    """ Tests the get_err_dct function
    """

    err_dct = err.get_err_dct(KTP_DCT1, PARAMS1)
    for _, (_, errs) in err_dct.items():
        assert max(abs(errs)) < 10  # check that all errors are less than 10%


if __name__ == '__main__':
    test_get_err_dct()
